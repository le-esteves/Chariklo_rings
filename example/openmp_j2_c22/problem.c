#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <omp.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"

//Tempo maximo:
const double tmax = 10000;
int test_particle_flag = 1;
#ifdef LOG
int collisions_count = 0;
#endif

int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c){
    //if (r->particles[c.p1].last_collision==r->t || r->particles[c.p2].last_collision==r->t) return 0;

    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    unsigned int i = c.p1;
    unsigned int j = c.p2;   //want j to be removed particle
    if ((r->particles[i].m)<(r->particles[j].m)){
        swap = 1;
        i = c.p2;
        j = c.p1;
    }
    //printf("here1");

    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
                
    double invmass = 1.0/(pi->m + pj->m);
    
    //Scale out energy from collision - initial energy
    double Ei=0, Ef=0;
    if(r->track_energy_offset){
        {
            double vx = pi->vx;
            double vy = pi->vy;
            double vz = pi->vz;
            // Calculate energy difference in inertial frame
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            Ei += 0.5*pi->m*(vx*vx + vy*vy + vz*vz);
        }
        {
            double vx = pj->vx;
            double vy = pj->vy;
            double vz = pj->vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            Ei += 0.5*pj->m*(vx*vx + vy*vy + vz*vz);
        }
        // problem is not here!
        const unsigned int N_active = ((r->N_active==-1)?r->N-r->N_var: (unsigned int)r->N_active);
        // No potential energy between test particles
        if (i<N_active || j<N_active){
            double x = pi->x - pj->x;
            double y = pi->y - pj->y;
            double z = pi->z - pj->z;
            double _r = sqrt(x*x + y*y + z*z);

            Ei += - r->G*pi->m*pj->m/_r;
        }
    }
    //printf("here2");
    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = cbrt(pi->r*pi->r*pi->r + pj->r*pj->r*pj->r);
    pi->last_collision = r->t;
    

    // Keeping track of energy offst
    if(r->track_energy_offset){
        {
            double vx = pi->vx;
            double vy = pi->vy;
            double vz = pi->vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            Ef += 0.5*pi->m*(vx*vx + vy*vy + vz*vz);
        }
        r->energy_offset += Ei - Ef;
    }
    int swap2 = swap?1:2;
    //printf("\nhere3 | %d", swap2);
    return swap2; // Remove particle p2 from simulation
}

int reb_collision_resolve_partialmerge(struct reb_simulation* const r, struct reb_collision c){
    struct reb_particle* const particles = r->particles;
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2;

#ifdef MPI
    int isloc = reb_communication_mpi_rootbox_is_local(r, c.ri);
    if (isloc==1){
#endif // MPI
        p2 = particles[c.p2];
#ifdef MPI
    }else{
        int N_root_per_node = r->N_root/r->mpi_num;
        int proc_id = c.ri/N_root_per_node;
        p2 = r->particles_recv[proc_id][c.p2];
    }
#endif // MPI
    
    const double G = 1;
    const double mEarth = 1; // mudar

    struct reb_vec6d gb = c.gb;
    double x21  = p1.x + gb.x  - p2.x; 
    double y21  = p1.y + gb.y  - p2.y; 
    double z21  = p1.z + gb.z  - p2.z; 
    double rp   = p1.r+p2.r;
    double oldvyouter;
    if (x21>0){
        oldvyouter = p1.vy;
    }else{
        oldvyouter = p2.vy;
    }
    if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;
    double vx21 = p1.vx + gb.vx - p2.vx; 
    double vy21 = p1.vy + gb.vy - p2.vy; 
    double vz21 = p1.vz + gb.vz - p2.vz; 
    if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0; // not approaching

    double r1[3] = {p1.x, p1.y, p1.z};
    double r2[3] = {p2.x, p2.y, p2.z};
    double relPos[3];
    for (int k = 0; k < 3; ++k) {
         relPos[k] = r1[k] - r2[k];
    }
    double relPosmod = sqrt(relPos[0]*relPos[0] + relPos[1]*relPos[1] + relPos[2]*relPos[2]);
    
    double comPos[3];
    for (int k = 0; k < 3; ++k) {
        comPos[k] = p1.m / (p1.m + p2.m) * relPos[k] + r2[k];
    }
    double comPosmod = sqrt(comPos[0]*comPos[0] + comPos[1]*comPos[1] + comPos[2]*comPos[2]);

    // // Jacobi energy of p1 and p2
    double v1[3] = {p1.vx, p1.vy, p1.vz};
    double v2[3] = {p2.vx, p2.vy, p2.vz};
    double relVel[3];
    for (int k = 0; k < 3; ++k) {
        relVel[k] = v1[k] - v2[k];
    }
    double relVelmod = sqrt(relVel[0]*relVel[0] + relVel[1]*relVel[1] + relVel[2]*relVel[2]);

    // if ((relVel[0] * relPos[0] + relVel[1] * relPos[1] + relVel[2] * relPos[2]) > 0) {
    //     return 0;
    // }
    double keplerVel = sqrt(G * mEarth * pow(comPosmod, -3.));
    double rhill = pow((p1.m + p2.m) / (3 * mEarth), 1.0/3.0) * comPosmod;

    double Ej = 0.5 * relVelmod*relVelmod - 1.5 * relPos[0]*relPos[0] * keplerVel*keplerVel +
                0.5 * relPos[2]*relPos[2] * keplerVel*keplerVel -
                G * (p1.m + p2.m) / relPosmod + 4.5 * rhill*rhill * keplerVel*keplerVel;
    //printf("\nColl %d and %d | Ej %f | %f < %f", p1.hash, p1.hash, Ej, (p1.r + p2.r), 0.7*rhill);

    if ((p1.m>=1.0)||(p2.m>=1.0)){ 
#ifdef LOG
        collisions_count++;
#endif
        int aux = reb_collision_resolve_merge(r, c);
        //if (aux>0) { printf("\nMerge %f and %f | %d %d particles", p1.m, p2.m, aux,r->N); }
        return aux;
    }
    //Merge
    if (Ej < 0 && (p1.r + p2.r) < 0.7 * rhill) {
#ifdef LOG
        collisions_count++;
#endif
        //printf("\nColl %d and %d | Ej %f | %f < %f", p1.hash, p2.hash, Ej, (p1.r + p2.r), 0.7*rhill);
        int aux = reb_collision_resolve_merge(r, c);
        
        //printf("\nMerge %d and %d | %d", p1.hash, p2.hash, aux);
        return aux;

    } else if (relPosmod < (p1.r + p2.r)) {
#ifdef LOG
        collisions_count++;
#endif
        //printf("Bounce %d and %d", p1.hash, p1.hash);
        // Hardsphere
        //if (p1.last_collision==t || p2.last_collision==t) return;
        // Bring the to balls in the xy plane.
        // NOTE: this could probabely be an atan (which is faster than atan2)
        double theta = atan2(z21,y21);
        double stheta = sin(theta);
        double ctheta = cos(theta);
        double vy21n = ctheta * vy21 + stheta * vz21;    
        double y21n = ctheta * y21 + stheta * z21;    
        
        // Bring the two balls onto the positive x axis.
        double phi = atan2(y21n,x21);
        double cphi = cos(phi);
        double sphi = sin(phi);
        double vx21nn = cphi * vx21  + sphi * vy21n;

        // Coefficient of restitution
        double eps= 1; // perfect bouncing by default 
        if (r->coefficient_of_restitution){
            eps = r->coefficient_of_restitution(r, vx21nn);
        }
        double dvx2 = -(1.0+eps)*vx21nn;
        double minr = (p1.r>p2.r)?p2.r:p1.r;
        double maxr = (p1.r<p2.r)?p2.r:p1.r;
        double mindv= minr*r->minimum_collision_velocity;
        double _r = sqrt(x21*x21 + y21*y21 + z21*z21);
        mindv *= 1.-(_r - maxr)/minr;
        if (mindv>maxr*r->minimum_collision_velocity)mindv = maxr*r->minimum_collision_velocity;
        if (dvx2<mindv) dvx2 = mindv;
        // Now we are rotating backwards
        double dvx2n = cphi * dvx2;        
        double dvy2n = sphi * dvx2;        
        double dvy2nn = ctheta * dvy2n;    
        double dvz2nn = stheta * dvy2n;


        // Applying the changes to the particles.
    #ifdef MPI
        if (isloc==1){
    #endif // MPI
        const double p2pf = p1.m/(p1.m+p2.m);
        particles[c.p2].vx -=    p2pf*dvx2n;
        particles[c.p2].vy -=    p2pf*dvy2nn;
        particles[c.p2].vz -=    p2pf*dvz2nn;
        particles[c.p2].last_collision = r->t;
    #ifdef MPI
        }
    #endif // MPI
        const double p1pf = p2.m/(p1.m+p2.m);
        particles[c.p1].vx +=    p1pf*dvx2n; 
        particles[c.p1].vy +=    p1pf*dvy2nn; 
        particles[c.p1].vz +=    p1pf*dvz2nn; 
        particles[c.p1].last_collision = r->t;
            
        // Return y-momentum change
        if (x21>0){
            r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
            r->collisions_log_n ++;
        }else{
            r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
            r->collisions_log_n ++;
        }
        return 0;
    } else {
        printf("Warning! No option found in collision.");
        return 0;
    }
}

int reb_collision_resolve_hardsphere_merge_central(struct reb_simulation* const r, struct reb_collision c){
    struct reb_particle* const particles = r->particles;
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2;

#ifdef MPI
    int isloc = reb_communication_mpi_rootbox_is_local(r, c.ri);
    if (isloc==1){
#endif // MPI
        p2 = particles[c.p2];
#ifdef MPI
    }else{
        int N_root_per_node = r->N_root/r->mpi_num;
        int proc_id = c.ri/N_root_per_node;
        p2 = r->particles_recv[proc_id][c.p2];
    }
#endif // MPI

    struct reb_vec6d gb = c.gb;
    double x21  = p1.x + gb.x  - p2.x; 
    double y21  = p1.y + gb.y  - p2.y; 
    double z21  = p1.z + gb.z  - p2.z; 
    double rp   = p1.r+p2.r;
    double oldvyouter;
    if (x21>0){
        oldvyouter = p1.vy;
    }else{
        oldvyouter = p2.vy;
    }
    if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;
    double vx21 = p1.vx + gb.vx - p2.vx; 
    double vy21 = p1.vy + gb.vy - p2.vy; 
    double vz21 = p1.vz + gb.vz - p2.vz; 
    if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0; // not approaching

    double r1[3] = {p1.x, p1.y, p1.z};
    double r2[3] = {p2.x, p2.y, p2.z};
    double relPos[3];
    for (int k = 0; k < 3; ++k) {
         relPos[k] = r1[k] - r2[k];
    }
    double relPosmod = sqrt(relPos[0]*relPos[0] + relPos[1]*relPos[1] + relPos[2]*relPos[2]);

    if ((p1.m>=1.0)||(p2.m>=1.0)){ 
#ifdef LOG
        collisions_count++;
#endif
        int aux = reb_collision_resolve_merge(r, c);
        //if (aux>0) { printf("\nMerge %f and %f | %d %d particles", p1.m, p2.m, aux,r->N); }
        return aux;
    } else if (relPosmod < (p1.r + p2.r)) {
#ifdef LOG
        collisions_count++;
#endif
        //printf("Bounce %d and %d", p1.hash, p1.hash);
        // Hardsphere
        //if (p1.last_collision==t || p2.last_collision==t) return;
        // Bring the to balls in the xy plane.
        // NOTE: this could probabely be an atan (which is faster than atan2)
        double theta = atan2(z21,y21);
        double stheta = sin(theta);
        double ctheta = cos(theta);
        double vy21n = ctheta * vy21 + stheta * vz21;    
        double y21n = ctheta * y21 + stheta * z21;    
        
        // Bring the two balls onto the positive x axis.
        double phi = atan2(y21n,x21);
        double cphi = cos(phi);
        double sphi = sin(phi);
        double vx21nn = cphi * vx21  + sphi * vy21n;

        // Coefficient of restitution
        double eps= 1; // perfect bouncing by default 
        if (r->coefficient_of_restitution){
            eps = r->coefficient_of_restitution(r, vx21nn);
            //printf("Here! cr = %.1f | mdv_coll = %e", eps, r->minimum_collision_velocity);
        }
        double dvx2 = -(1.0+eps)*vx21nn;
        double minr = (p1.r>p2.r)?p2.r:p1.r;
        double maxr = (p1.r<p2.r)?p2.r:p1.r;
        double mindv= minr*r->minimum_collision_velocity;
        double _r = sqrt(x21*x21 + y21*y21 + z21*z21);
        mindv *= 1.-(_r - maxr)/minr;
        if (mindv>maxr*r->minimum_collision_velocity)mindv = maxr*r->minimum_collision_velocity;
        if (dvx2<mindv) dvx2 = mindv;
        // Now we are rotating backwards
        double dvx2n = cphi * dvx2;        
        double dvy2n = sphi * dvx2;        
        double dvy2nn = ctheta * dvy2n;    
        double dvz2nn = stheta * dvy2n;


        // Applying the changes to the particles.
    #ifdef MPI
        if (isloc==1){
    #endif // MPI
        const double p2pf = p1.m/(p1.m+p2.m);
        particles[c.p2].vx -=    p2pf*dvx2n;
        particles[c.p2].vy -=    p2pf*dvy2nn;
        particles[c.p2].vz -=    p2pf*dvz2nn;
        particles[c.p2].last_collision = r->t;
    #ifdef MPI
        }
    #endif // MPI
        const double p1pf = p2.m/(p1.m+p2.m);
        particles[c.p1].vx +=    p1pf*dvx2n; 
        particles[c.p1].vy +=    p1pf*dvy2nn; 
        particles[c.p1].vz +=    p1pf*dvz2nn; 
        particles[c.p1].last_collision = r->t;
            
        // Return y-momentum change
        if (x21>0){
            r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
            r->collisions_log_n ++;
        }else{
            r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
            r->collisions_log_n ++;
        }
        return 0;
    } else {
        printf("Warning! No option found in collision.");
        return 0;
    }
}

double coefficient_of_restitution_constant(const struct reb_simulation* const r, double v){
double cr = 0.1;
    return cr;
}

void test_particles_remove(struct reb_simulation* const r){
    const int N = r->N;
    for (int i = 1; i < N; i++) { // starting from i = 1 to skip CB
        struct reb_particle *p = &r->particles[i];
        double dist = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
        if (dist <= 1.2 && p->m == 0.0) { // test particle and close enough to the central body
            struct reb_orbit o = reb_orbit_from_particle(r->G, *p, r->particles[0]);
            double q = o.a*(1.0 - o.e); // pericenter distance
            if (q < 1.0) {
                reb_simulation_remove_particle(r, i, 0); // do not sort the particles yet (third arg)
            }
        }
    }
}

//const double C22 = 0.007480263569693367; // modified by Python script
//const double C20 = -0.12436657282151249; // modified by Python script
const double C22 = 0.007802419354838727; // calculado via Giuliatti Winter et al 2023
const double C20 = -0.1256234391259105; // calculado via Giuliatti Winter et al 2023
//const double w = 0.5328671868528622; // modified by Python script
// w = 2.4919036 rad/s, usando T = 7.004 h
const double w = 0.4212; // modified by Python script
const double T = 2.*M_PI/w; // Rotational period 

const double M_QUAOAR = 1.0; // mass of Quaoar in (kg) 
const double R_QUAOAR = 1.0; // equi radius of Quaoar (m)
const double T_QUAOAR = 2.*M_PI/w; // Rotational period (s)

void force_harmonics(struct reb_simulation* r){
    const int N = r->N;
    if (C20==0 && C22==0) return;
    // Central body (Quaoar)
    int i=0;
    for (int k=0; k<N; k++) {
        if (r->particles[k].m >= 1.0) {
            i=k;
            break;
        }
    }
    const struct reb_particle planet = r->particles[i]; // cache
    const double cte  = r->G*planet.m*R_QUAOAR*R_QUAOAR;
    const double theta = r->t*w;
#pragma omp parallel for
    for (int j=0;j<N;j++){
        if (j == i){
            continue;
        }
        const struct reb_particle p = r->particles[j]; // cache
        const double sprx = p.x-planet.x;
        const double spry = p.y-planet.y;
        const double sprz = p.z-planet.z;
        const double x2 = sprx*sprx;
        const double y2 = spry*spry;
        const double z2 = sprz*sprz;
        const double z3 = z2*sprz;
        const double pr  = sqrt(x2 + y2 + z2);
        const double pr2  = pr*pr;      // distance^2 relative to planet
        const double pr3  = pr2*pr;
        const double pr5  = pr3*pr2;
        const double GMR2_r5  = cte/pr5;
        const double Csmenos = C22*cos(2.0*theta);
        const double Csplus  = C22*sin(2.0*theta);
        
        const double term1x = C20*sprx*(1.5 - 7.5*z2/pr2);
        const double term1y = C20*spry*(1.5 - 7.5*z2/pr2);
        const double term1z = C20*sprz*(4.5 - 7.5*z3/pr2);

        const double term2x = 6.0*(Csmenos*sprx + Csplus*spry);
        const double term2y = 6.0*(Csplus*sprx - Csmenos*spry);

        const double term3x = 15.0*sprx*(Csmenos*(y2 - x2) - 2.0*sprx*spry*Csplus)/pr2;
        const double term3y = 15.0*spry*(Csmenos*(y2 - x2) - 2.0*sprx*spry*Csplus)/pr2;
        const double term3z = 15.0*sprz*(Csmenos*(y2 - x2) - 2.0*sprx*spry*Csplus)/pr2;
 

        const double pax  = GMR2_r5*(term1x + term2x + term3x);
        const double pay  = GMR2_r5*(term1y + term2y + term3y);
        const double paz  = GMR2_r5*(term1z + term3z);
        
        r->particles[j].ax += pax;
        r->particles[j].ay += pay;
        r->particles[j].az += paz;
    }
}

// Função para encontrar o arquivo de restart mais recente
char* find_latest_restart(const char* directory) {
    DIR *dir;
    struct dirent *ent;
    char* latest_restart = NULL;
    int max_number = -1;

    if ((dir = opendir(directory)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (strstr(ent->d_name, "restart_") != NULL && strstr(ent->d_name, ".bin") != NULL) {
                int number;
                if (sscanf(ent->d_name, "restart_%d.bin", &number) == 1) {
                    if (number > max_number) {
                        max_number = number;
                        if (latest_restart) free(latest_restart);
                        latest_restart = malloc(strlen(directory) + strlen(ent->d_name) + 2);
                        sprintf(latest_restart, "%s/%s", directory, ent->d_name);
                    }
                }
            }
        }
        closedir(dir);
    } else {
        perror("Could not open directory");
        return NULL;
    }

    return latest_restart;
}

// necessário para partículas teste, ou auto-gravidade off
void reb_simulation_central_to_zero(struct reb_simulation* const r){
    struct reb_particle* restrict const particles = r->particles;
    int central_index = -1;
    for (int i = 0; i < r->N; i++) {
        if (particles[i].m >= 1.0) {
            central_index = i;
            break;
        }
    }
    if (central_index > 0) {
        struct reb_particle temp = particles[0];
        particles[0] = particles[central_index];
        particles[central_index] = temp;
    }
}

void heartbeat(struct reb_simulation* const r);

int run_sim(){
    const char* output_directory = "Outputs";
    char* latest_restart = find_latest_restart(output_directory);

#ifdef RESTART
    if (latest_restart == NULL) {
        printf("No restart files found. Starting from scratch.\n");
        latest_restart = "Outputs/restart_0.bin"; // Assume que o primeiro restart é sem numero
    } else {
        printf("Creating simulation from binary file: %s\n", latest_restart);
    }

    struct reb_simulationarchive* archive = reb_simulationarchive_create_from_file(latest_restart);
    struct reb_simulation* const r = reb_simulation_create_from_simulationarchive(archive, -1);
#else
    if (latest_restart != NULL) {
        printf("A binary file was found in Outputs. To start a new simulation delete %s. Otherwise, continue from the latest restart by using ~$make restart\n", latest_restart);
        return EXIT_FAILURE;
    }
    printf("Starting new simulation \n");
    struct reb_simulation* const r = reb_simulation_create();
#endif
    
    //reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator       = REB_INTEGRATOR_LEAPFROG;
    r->gravity          = REB_GRAVITY_BASIC;
    r->boundary         = REB_BOUNDARY_OPEN;
    r->opening_angle2   = 0.5;          // referente a arvore da GRAVIDADE
    r->G                = 1.0;            // Gravitational constant
    r->dt               = 0.005;         // Timestep
    r->softening        = 0.02;         // Gravitational softening length
    r->collision        = REB_COLLISION_TREE;
    r->coefficient_of_restitution = coefficient_of_restitution_constant;
    r->collision_resolve = reb_collision_resolve_partialmerge;
    r->additional_forces = force_harmonics;
    r->central_to_zero  = 1;             // setar o corpo central na posição 0 do vetor particles

    const double boxsize = 25.0; // Distancia de ejecao
    reb_simulation_configure_box(r,boxsize,1,1,1);

#ifndef RESTART
    // Initial conditions
    struct reb_particle star = {0};
    star.m  = 1.0;
    star.r = 1.0 ;
    star.hash = 0;
    reb_simulation_add(r, star);

    FILE *entrada;
    double mass, density, a_i, a_f, radius, a, e, inc, w, Omg, lmd;
    int i=1;
    entrada = fopen("../../Inputs/disc_chariklo_tp.in", "r");
    if (entrada == NULL) {
        perror("Erro ao abrir o arquivo");
        return EXIT_FAILURE;
    }
    while (fscanf(entrada, "%le %le %le %le %le %le %le %le %le %le %le\n", &mass, &a_i, &a_f, &radius, &density, &a, &e, &inc, &w, &Omg, &lmd) == 11) {
        
        struct reb_particle pt = reb_particle_from_orbit(r->G, star, mass, a, e, inc, Omg, w, lmd);
        pt.hash = i;
        pt.r = radius;
        reb_simulation_add(r,pt);
        i=i+1;
    }
    fclose(entrada);
#endif

    r->N_active = 1;
    r->heartbeat = heartbeat;

    double output_interval = 1;

    int next_restart_number = 0;
    char next_restart_filename[256];
#ifdef RESTART
    // Determina o próximo número de restart
    if (latest_restart != NULL) {
        sscanf(latest_restart, "Outputs/restart_%d.bin", &next_restart_number);
        next_restart_number++;
    }
#endif
    sprintf(next_restart_filename, "Outputs/restart_%d.bin", next_restart_number);

    reb_simulation_save_to_file_interval(r, next_restart_filename, output_interval);
    reb_simulation_integrate(r, tmax);
    // Cleanup
    reb_simulation_free(r); 
    return 0;
}

void heartbeat(struct reb_simulation* const r){
    reb_simulation_update_tree(r);
    reb_simulation_move_to_com(r);
    
    if (test_particle_flag) {
        test_particles_remove(r);
    }    

#ifdef LOG
    static double last_time = 0.0;
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double current_time = tim.tv_sec + (tim.tv_usec / 1000000.0);
    double timestep_runtime = current_time - last_time;
    last_time = current_time;
    printf("Time: %.3f | N_part: %d | colls.: %d | part.m[0]: %.6f | dt runtime: %.6f s\n", r->t, r->N, collisions_count, r->particles[0].m, 
        timestep_runtime);
    collisions_count = 0;
#endif
}

int main(int argc, char* argv[]){
    // Set the number of OpenMP threads
    int np = 16;
    // Set the number of OpenMP threads to be the number of processors    
    omp_set_num_threads(np);
    run_sim();
}
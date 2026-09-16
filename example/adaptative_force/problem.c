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
	const double tmax = 80000; // altered by create_sim 
	int test_particle_flag = 0; // altered by create_sim 
#ifdef LOG
int collisions_count = 0;
int collisions_soft_count = 0;
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

int reb_collision_resolve_softsphere(struct reb_simulation* const r, struct reb_collision c){
    struct reb_particle* const particles = r->particles;
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2 = particles[c.p2];

    struct reb_vec6d gb = c.gb;
    double x21  = p1.x + gb.x  - p2.x; 
    double y21  = p1.y + gb.y  - p2.y; 
    double z21  = p1.z + gb.z  - p2.z; 
    double rp   = p1.r+p2.r;

    double d2 = x21*x21 + y21*y21 + z21*z21;
    if (rp*rp < d2) return 0; // not overlapping

    double vx21 = p1.vx + gb.vx - p2.vx; 
    double vy21 = p1.vy + gb.vy - p2.vy; 
    double vz21 = p1.vz + gb.vz - p2.vz; 

    if ((p1.m>=1.0)||(p2.m>=1.0)){ 
        // merge
        if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0; // not approaching
#ifdef LOG
        collisions_count++;
#endif
        int aux = reb_collision_resolve_merge(r, c);
        //if (aux>0) { printf("Merge %f and %f | %d-> to remove | N=%d\n", p1.m, p2.m, aux,r->N); }
        return aux;
    } else { // softsphere
        // ==== partículas em sobreposição ====
        // Vetor separação com correção de shear (já calculado: x21,y21,z21)
        #ifdef LOG
            collisions_soft_count++;
        #endif
        
        double d  = sqrt(d2);
    
        // Penetração alpha e unitário e_hat (da p2 -> p1)
        double alpha = rp - d;               // >0 garantido pelo if anterior

        double ex = x21/d, ey = y21/d, ez = z21/d;

        // Projeção normal da velocidade relativa e alpha_dot
        double v_n = vx21*ex + vy21*ey + vz21*ez; // d_dot
        double alpha_dot = -v_n;

        // Parâmetros do modelo
        double Tdur = r->ri_lfada.Tdur;   // duração "artificial" da colisão
        double eps = 1; // perfect bouncing by default 
        if (r->coefficient_of_restitution){
            eps = r->coefficient_of_restitution(r, v_n);
        }

        // Massas e m_eff
        double meff = (p1.m*p2.m)/(p1.m+p2.m);

        // k1, k2 conforme eq. (8) do paper
        // omega0 = pi/Tdur
        double omega0 = M_PI / Tdur;

        // Trata eps <= 0 de forma robusta (evita ln(0))
        if (eps < 1e-6) eps = 1e-6;
        double L = M_PI / log(eps); // pi / ln(eps)
        double root = sqrt(L*L + 1.0);

        double k1 = meff * (omega0*omega0);
        double k2 = meff * (2.0 * omega0 * root);

        // Força normal (direção +e_hat na partícula 1; oposta na 2)
        double fmag = (k1*alpha + k2*alpha_dot);
        // Acelerações
        double Fax = (fmag * ex);
        double Fay = (fmag * ey);
        double Faz = (fmag * ez);

        // Escreve no buffer de acelerações (kick) — +F/m em p1, -F/m em p2
        r->ri_lfada.par_soft_acc[c.p1].x +=  Fax / p1.m;
        r->ri_lfada.par_soft_acc[c.p1].y +=  Fay / p1.m;
        r->ri_lfada.par_soft_acc[c.p1].z +=  Faz / p1.m;

        r->ri_lfada.par_soft_acc[c.p2].x -=  Fax / p2.m;
        r->ri_lfada.par_soft_acc[c.p2].y -=  Fay / p2.m;
        r->ri_lfada.par_soft_acc[c.p2].z -=  Faz / p2.m;

        return 0; // não remove nenhuma partícula
    }
}

double coefficient_of_restitution(const struct reb_simulation* const r, double v){
	double cr = 0.1; // altered by create_sim 
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
    int i=-1;
    if (!r->central_to_zero) {
        for (int k=0; k<N; k++) {
            if (r->particles[k].m >= 1.0) {
                i=k;
                break;
            }
        }
    } else {
        i = 0; // assume the first particle is the central body if central_to_zero is true
    }
    if (i==-1) {
        printf("Warning: no central body found for harmonics!\n");
        return;
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
        
        if (r->integrator == REB_INTEGRATOR_LEAPFROG_ADA){
            r->ri_lfada.par_soft_acc[j].x += pax;
            r->ri_lfada.par_soft_acc[j].y += pay;
            r->ri_lfada.par_soft_acc[j].z += paz;
        } else {
            r->particles[j].ax += pax;
            r->particles[j].ay += pay;
            r->particles[j].az += paz;
        }
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
    
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator       = REB_INTEGRATOR_LEAPFROG_ADA;
	r->gravity          = REB_GRAVITY_TREE; // altered by create_sim 
    r->boundary         = REB_BOUNDARY_OPEN;
    r->opening_angle2   = 0.5;            // referente a arvore da GRAVIDADE e COLISÃO
    r->G                = 1.0;            // Gravitational constant
	r->ri_lfada.dt_global = 0.005; // altered by create_sim 
	r->ri_lfada.dt_factor = 20; // altered by create_sim 
	r->ri_lfada.mode    = 1; // altered by create_sim 
    r->dt               = r->ri_lfada.dt_global/r->ri_lfada.dt_factor; // Timestep
	r->softening        = 0.00263312388; // altered by create_sim 
    r->collision        = REB_COLLISION_TREE;
	r->collision_resolve = reb_collision_resolve_softsphere; // altered by create_sim 
    r->ri_lfada.Tdur    = r->dt * 50; // coeficiente da duração da colisão
    r->coefficient_of_restitution = coefficient_of_restitution;
    r->additional_forces = force_harmonics;
	r->central_to_zero  = 0; // altered by create_sim 

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
	entrada = fopen("../../Inputs/disc_chariklo_1e-3.in", "r"); // altered by create_sim 
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

	r->N_active = r->N; // altered by create_sim 
    r->heartbeat = heartbeat;

	double output_interval = 6; // altered by create_sim 

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
    if (r->ri_lfada.step_counter == r->ri_lfada.dt_factor-1) {
        reb_simulation_update_tree(r);
        reb_simulation_move_to_com(r);
    }
    
    // if (test_particle_flag) {
    //     test_particles_remove(r);
    // }    

#ifdef LOG
    static double last_time = 0.0;
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double current_time = tim.tv_sec + (tim.tv_usec / 1000000.0);
    double timestep_runtime = current_time - last_time;
    last_time = current_time;

    if (r->ri_lfada.step_counter == 0) {
        printf("Global || Time: %.3f | dt: %.3f | N_part: %d | Ncolls.: %d | colls.: %d | colls. soft: %d | dt runtime: %.6f s\n", r->t, r->dt, 
        r->N, r->collisions_N, collisions_count, collisions_soft_count, timestep_runtime);
    }
    // } else {
    //     printf("Time: %.3f | dt: %.3f | N_part: %d | Ncolls.: %d | colls.: %d | colls. soft: %d | dt runtime: %.6f s\n", r->t, r->dt, 
    //     r->N, r->collisions_N, collisions_count, collisions_soft_count, timestep_runtime);
    // }
    collisions_count = 0;
    collisions_soft_count = 0;
#endif
}

int main(int argc, char* argv[]){
    // Set the number of OpenMP threads
	int np = 16; // altered by create_sim 
    // Set the number of OpenMP threads to be the number of processors    
    omp_set_num_threads(np);
    run_sim();
}
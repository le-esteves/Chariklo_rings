/**
 * @file 	collisions.c
 * @brief 	Collision search using a line sweep algorithm, O(N log(N)).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The routines in this file implement a collision detection
 * method called line sweep. It is very fast if all dimensions except one 
 * are small. The algorithm is similar to the original algorithm proposed
 * by Bentley & Ottmann (1979) but does not maintain a binary search tree.
 * This is much faster as long as the number of particle trajectories
 * currently intersecting the plane is small.
 *
 * The sweeping direction in this implementation is x.
 * 
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "collision.h"
#include "collisions_sweep.h"
#include "rebound.h"
#include "tree.h"
#include "boundary.h"
#ifdef OPENMP
#include <omp.h>
#endif


double 	collisions_max_r	= 0;
double 	collisions_max2_r	= 0;
int	sweeps_proc		= 1;	/**< Number of processors used for seeping algorithm. */
int 	sweeps_init_done 	= 0;	/**< Used for initialisation of data structures. */

//static inline double min(double a, double b){ return (a>b)?b:a;}
//static inline double max(double a, double b){ return (b>a)?b:a;}
static inline double sgn(const double a){ return (a>=0 ? 1. : -1); }

/** 
 * This function checks if two particles colliding during one drift step.
 * @param pt1 reb_particle 1. 
 * @param pt2 reb_particle 2. 
 * @param proci Processor id (OpenMP) for this collision.
 * @param crossing Flag that is one if one of the particles crosses a boundary in this timestep.
 * @param ghostbox Ghostbox used in this collision.
 */
void detect_collision_of_pair(struct reb_simulation* const r, int pt1, int pt2, int proci, int crossing, struct reb_vec6d gb);

/**
 * Structure that stores a start or end point of a particle trajectory.
 */
struct xvalue {
	double 	x;		// position along sweep axis
	int 	inout;		// start or endpoint
	int	    nx;		
	int 	crossing;	// crosses boundary
	int 	pt;		// particle
};

/**
 * Structure that contains a list of xvalues.
 */
struct xvaluelist {
	struct xvalue* xvalues;
	int 	N;		/**< Current array size. */
	int 	Nmax;		/**< Maximum array size before realloc() is needed. */
};
struct  xvaluelist* sweepx;	/**< Pointers to the SWEEPX list of each processor. */

/**
 * Structure that contains a list of collisions.
 */
// struct reb_collisionlist {
// 	struct reb_collision* collisions;
// 	int 	N;		/**< Current array size. */
// 	int 	Nmax;		/**< Maximum array size before realloc() is needed. */
// };
//struct 	collisionlist* clist;	/**< Pointers to the collisions list of each processor. */

/**
 * Adds a line to the SWEEPX array of processor proci.
 */
void add_line_to_xvsublist(double x1, double x2, int pt, int n, int proci, int crossing){
	int N = sweepx[proci].N;
	
	if (N+2>sweepx[proci].Nmax){
		sweepx[proci].Nmax 	+= 1024;
		sweepx[proci].xvalues	= (struct xvalue*)realloc(sweepx[proci].xvalues,sweepx[proci].Nmax*sizeof(struct xvalue));
	}

	sweepx[proci].xvalues[N].x 		= x1;
	sweepx[proci].xvalues[N].pt 		= pt;
	sweepx[proci].xvalues[N].nx 		= n;
	sweepx[proci].xvalues[N].inout 		= 0;
	sweepx[proci].xvalues[N].crossing 	= crossing;
	sweepx[proci].xvalues[N+1].x 		= x2;
	sweepx[proci].xvalues[N+1].pt 		= pt;
	sweepx[proci].xvalues[N+1].nx 		= n;
	sweepx[proci].xvalues[N+1].inout	= 1;
	sweepx[proci].xvalues[N+1].crossing	= crossing;

	sweepx[proci].N += 2;
}

/**
 * Adds a line to the SWEEPX array and checks for crossings of processor boundaries.
 */
void add_line_to_xvlist(struct reb_simulation* const r, double x1, double x2, int pt, int n, int crossing){
	int procix1 = (int)(floor( (x1/r->boxsize.x+0.5) *(double)sweeps_proc));// %sweeps.xvlists;
	int procix2 = (int)(floor( (x2/r->boxsize.x+0.5) *(double)sweeps_proc));// %sweeps.xvlists;
	if (procix2>=sweeps_proc){
		procix2 = sweeps_proc-1;
	}
	if (procix1<0){
		procix1 = 0;
	}

	if (procix1!=procix2){
		double b = -r->boxsize.x/2.+r->boxsize.x/(double)sweeps_proc*(double)procix2; 
		add_line_to_xvsublist(x1,b,pt,n,procix1,1);
		add_line_to_xvsublist(b,x2,pt,n,procix2,1);
	}else{
		add_line_to_xvsublist(x1,x2,pt,n,procix1,crossing);
	}
}

/**
 * Adds a line to the SWEEPX array and checks for crossings of simulation boundaries.
 */
void add_to_xvlist(struct reb_simulation* const r, double x1, double x2, int pt){
	double xmin, xmax;
	if (x1 < x2){
		xmin = x1;
		xmax = x2;
	}else{
		xmin = x2;
		xmax = x1;
	}
	double radius = r->particles[pt].r*1.0001; //Safety factor to avoid floating point issues.
	xmin -= radius;
	xmax += radius;

	if (xmin < -r->boxsize.x/2.){
		add_line_to_xvlist(r, xmin+r->boxsize.x,r->boxsize.x/2.,pt,1,1);
		add_line_to_xvlist(r, -r->boxsize.x/2.,xmax,pt,0,1);
		return;
	}
	if (xmax > r->boxsize.x/2.){
		add_line_to_xvlist(r, -r->boxsize.x/2.,xmax-r->boxsize.x,pt,-1,1);
		add_line_to_xvlist(r, xmin,r->boxsize.x/2.,pt,0,1);
		return;
	}
	add_line_to_xvlist(r, xmin,xmax,pt,0,0);
}

/**
 * Compares the x position of two xvalues.
 */
int compare_xvalue (const void * a, const void * b){
	const double diff = ((struct xvalue*)a)->x - ((struct xvalue*)b)->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

/**
 * Compares the x position of two particles.
 */
int compare_particle (const void * a, const void * b){
	const double diff = ((struct reb_particle*)a)->x - ((struct reb_particle*)b)->x;
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

/**
 * Sorts the array xvl with insertion sort.
 */
void collisions_sweep_insertionsort_xvaluelist(struct xvaluelist* xvl){
	struct xvalue* xv = xvl->xvalues;
	int _N = xvl->N;
	for(int j=1;j<_N;j++){
		struct xvalue key = xv[j];
		int i = j - 1;
		while(i >= 0 && xv[i].x > key.x){
		    xv[i+1] = xv[i];
		    i--;
		}
		xv[i+1] = key;
	}
}

/**
 * Sorts the particle array with insertion sort.
 */
void collisions_sweep_insertionsort_particles(struct reb_simulation* const r){
	for(int j=1;j<r->N;j++){
		const struct reb_particle key = r->particles[j];
		int i = j - 1;
		while(i >= 0 && r->particles[i].x > key.x){
		    r->particles[i+1] = r->particles[i];
		    i--;
		}
		r->particles[i+1] = key;
	}
}



void reb_collision_search_sweep(struct reb_simulation* const r, const double dt_global){
	if (r->collision == REB_COLLISION_NONE) return;
	// To implement a REB_COLLISION_SWEEPX flag to this function.
	if (sweeps_init_done!=1){
		sweeps_init_done = 1;
#ifdef OPENMP
		sweeps_proc 		= omp_get_max_threads();
#endif // OPENMP
		sweepx		= (struct xvaluelist*)   calloc(sweeps_proc,sizeof(struct xvaluelist));
		
		//clist		= (struct reb_collisionlist*)calloc(sweeps_proc,sizeof(struct reb_collisionlist));
//#ifndef TREE
		if (!r->tree_root) {
			// Sort particles according to their x position to speed up sorting of lines.
			// Initially the particles are not pre-sorted, thus qsort is faster than insertionsort.
			// Note that this rearranges particles and will cause problems if the particle id is used elsewhere.
			qsort (r->particles, r->N, sizeof(struct reb_particle), compare_particle);
			//printf("not tree x-sort \n");
		}
	} else {
		if (!r->tree_root) {
			// Keep particles sorted according to their x position to speed up sorting of lines.
			collisions_sweep_insertionsort_particles(r);
		}
//#endif //TREE
	}
	for (int i=0;i<r->N;i++){
		double oldx = r->particles[i].x-0.5*dt_global*r->particles[i].vx;	
		double newx = r->particles[i].x+0.5*dt_global*r->particles[i].vx;	
		add_to_xvlist(r, oldx,newx,i);
	}
	
	// Precalculate most comonly used ghostboxes
	const struct reb_vec6d gb00  = reb_boundary_get_ghostbox(r, 0, 0, 0);
	const struct reb_vec6d gb0p1 = reb_boundary_get_ghostbox(r, 0, 1, 0);
	const struct reb_vec6d gb0m1 = reb_boundary_get_ghostbox(r, 0, -1, 0);

#pragma omp parallel for schedule (static,1)
	for (int proci=0;proci<sweeps_proc;proci++){
		struct xvaluelist* sweepxi = &(sweepx[proci]);
//#ifdef TREE
		if (r->tree_root) {
			// Use quicksort when there is a tree. reb_particles are not pre-sorted.
			qsort (sweepxi->xvalues, sweepxi->N, sizeof(struct xvalue), compare_xvalue);
			//printf("tree x-sort \n");
		} else {
//#else //TREE 
			// Use insertionsort when there is a tree. reb_particles are pre-sorted.
			collisions_sweep_insertionsort_xvaluelist(sweepxi);	
		}
//#endif //TREE
		//printf("proc i: %d | sweepl N: %d \n", proci, sweepxi->N);			

		// SWEEPL: List of lines intersecting the plane.
		struct xvaluelist sweepl = {NULL,0,0};

		for (int i=0;i<sweepxi->N;i++){
			struct xvalue xv = sweepxi->xvalues[i];
			if (xv.inout == 0){
				// Add event if start of line
				if (sweepl.N>=sweepl.Nmax){
					sweepl.Nmax += 32;
					sweepl.xvalues = realloc(sweepl.xvalues,sizeof(struct xvalue)*sweepl.Nmax); 
				}
				sweepl.xvalues[sweepl.N] = xv;
				// Check for collisions with other particles in SWEEPL
				for (int k=0;k<sweepl.N;k++){
					int p1 = xv.pt;
					int p2 = sweepl.xvalues[k].pt;
					if (p1==p2) continue;
					int gbnx = xv.nx;
					if (sweepl.xvalues[k].nx!=0){
						if (sweepl.xvalues[k].nx==xv.nx) continue;
						int tmp = p1;
						p1 = p2;
						p2 = tmp;
						gbnx = sweepl.xvalues[k].nx;
					}
					if (gbnx==0){
						// Use cached ghostboxes if possible
						detect_collision_of_pair(r,p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb00);
						detect_collision_of_pair(r,p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb0p1);
						detect_collision_of_pair(r,p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb0m1);
					}else{
						for (int gbny = -1; gbny<=1; gbny++){
							struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbnx,gbny,0);
							detect_collision_of_pair(r,p1,p2,proci,sweepl.xvalues[k].crossing||xv.crossing,gb);
						}
					}
				}
				sweepl.N++;
			}else{
				// Remove event if end of line
				for (int j=0;j<sweepl.N;j++){
					if (sweepl.xvalues[j].pt == xv.pt){
						sweepl.N--;
						sweepl.xvalues[j] = sweepl.xvalues[sweepl.N];
						j--;
						break;
					}
				}
			}
			//printf("sweepl N: %d \n", sweepl.N);
		}	
		//printf("N detected colls: %d \n", r->collisions_N);
		free(sweepl.xvalues);
	}
	//printf("end reb_coll_search_sweep \n");
}

void detect_collision_of_pair(struct reb_simulation* const r, int pt1, int pt2, int proci, int crossing, struct reb_vec6d gb){
	const struct reb_particle p1 = r->particles[pt1];
	const struct reb_particle p2 = r->particles[pt2];
	double x  = p1.x  + gb.x	- p2.x;
	double y  = p1.y  + gb.y	- p2.y;
	double z  = p1.z  + gb.z	- p2.z;
	double vx = p1.vx + gb.vx	- p2.vx;
	double vy = p1.vy + gb.vy   - p2.vy;
	double vz = p1.vz + gb.vz	- p2.vz;

	double a = vx*vx + vy*vy + vz*vz;
	double b = 2.*(vx*x + vy*y + vz*z);
	double rr = p1.r + p2.r;
	double c = -rr*rr + x*x + y*y + z*z;

	double root = b*b-4.*a*c;
	if (root>=0.){
		// Floating point optimized solution of a quadratic equation. Avoids cancelations.
		double q = -0.5*(b+sgn(b)*sqrt(root));
		double time1 = c/q;
		double time2 = q/a;
		if (time1>time2){
			double tmp = time2;
			time2=time1;
			time1=tmp;
		}
		if ( (time1>-r->dt/2. && time1<r->dt/2.) || (time1<-r->dt/2. && time2>r->dt/2.) ){
#pragma omp critical
            {
			if (r->N_allocated_collisions<=r->collisions_N){
				// Allocate memory if there is no space in array.
				// Init to 32 if no space has been allocated yet, otherwise double it.
				r->N_allocated_collisions = r->N_allocated_collisions ? r->N_allocated_collisions * 2 : 32;
				r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->N_allocated_collisions);
			}
			}
			//struct reb_collisionlist* clisti = &(clist[proci]);
			// if (clisti->N>=clisti->Nmax){
			// 	clisti->Nmax	 	+= 1024;
			// 	clisti->collisions	= (struct reb_collision*)realloc(clisti->collisions,clisti->Nmax*sizeof(struct reb_collision));
			// }
			//struct reb_collision* c = &(clisti->collisions[clisti->N]);
			struct reb_collision* c = &(r->collisions[r->collisions_N]);
			c->p1		= pt1;
			c->p2		= pt2;
			c->gb	 	= gb;
			if ( (time1>-r->dt/2. && time1<r->dt/2.)) {
				c->time 	= time1;
			}else{
				c->time 	= 0;
			}
			// Ignorar crossing até eu entender se tem como usar a arvore em conjunto com o plane sweep para colisões.
			//c->crossing 	= crossing;
			r->collisions_N++;
		}
	}
}

void collisions_resolve(struct reb_simulation* const r){
	//printf("start Resolve for N colls: %d \n", r->collisions_N);
	if (r->collision == REB_COLLISION_NONE) return;

	#ifdef OPENMP
	omp_lock_t boundarylock;
	omp_init_lock(&boundarylock);
#endif //OPENMP
	int (*resolve) (struct reb_simulation* const r, struct reb_collision c1) = r->collision_resolve;
	if (resolve==NULL){
        // Default is to throw an exception
        resolve = reb_collision_resolve_halt;
    }

	// Reset softsphere accelerations
	if (r->ri_lfada.mode) {
		for (int i=0; i<r->N; i++){
			r->ri_lfada.par_soft_acc[i].x = 0;
			r->ri_lfada.par_soft_acc[i].y = 0;
			r->ri_lfada.par_soft_acc[i].z = 0;
		}
	}

	// Randomize array.	
	for(int i=0; i<r->collisions_N; i++){
		int new = rand_r(&(r->rand_seed))%r->collisions_N;
		struct reb_collision ctemp = r->collisions[i];
		r->collisions[i]=r->collisions[new];
		r->collisions[new]=ctemp;
	}
		
//#pragma omp parallel for schedule (static,1)
	for(int i=0; i<r->collisions_N; i++){
		struct reb_collision* c1 = &(r->collisions[i]);

		// Skip collisions involving removed particles
		if (c1->p1 == -1 || c1->p2 == -1){
			continue;
		}

		int outcome = resolve(r, *c1);

		// Remove particles
		if (outcome & 1){
			// Remove p1
			int removedp1 = reb_simulation_remove_particle(r, c1->p1, 0); // 0 for not keeping sorted (in a Tree)
			if (removedp1){
				if (r->tree_root){ // In a tree, particles get removed later.
					for (int j=i+1;j<r->collisions_N;j++){ // Update other collisions
						struct reb_collision* cp = &(r->collisions[j]);
						// Skip collisions which involved the removed particle
						if (cp->p1==c1->p1 || cp->p2==c1->p1){
							cp->p1 = -1;
							cp->p2 = -1;
						}
					}
				}
				c1->p1 = -1; // Mark as removed
				c1->p2 = -1;
			}
		}
		if (outcome & 2){
			// Remove p2
			int removedp2 = reb_simulation_remove_particle(r, c1->p2, 0);
			if (removedp2){ // Update other collisions
				if (r->tree_root){ // In a tree, particles get removed later. 
					for (int j=i+1;j<r->collisions_N;j++){ // Update other collisions
						struct reb_collision* cp = &(r->collisions[j]);
						// Skip collisions which involved the removed particle
						if (cp->p1==c1->p2 || cp->p2==c1->p2){
							cp->p1 = -1;
							cp->p2 = -1;
						}
					}
				}
				c1->p1 = -1; // Mark as removed
				c1->p2 = -1;
			}
		}
	}
	// Caso o próximo dt seja o global, tal que dt+dt/2 % dt_global == 0, resetar as colisões
	if (r->ri_lfada.step_counter == r->ri_lfada.dt_factor-1) {
		//printf("entered\n");
		r->collisions_N = 0;
		for (int proci = 0; proci < sweeps_proc; proci++){
			sweepx[proci].N = 0;	// Reset the sweepx list for next run.
		}
	}
	//printf ("Coll_n: %d",r->collisions_N);
#ifdef OPENMP
	omp_destroy_lock(&boundarylock);
#endif //OPENMP
}


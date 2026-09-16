/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
 * 
 * @section 	LICENSE
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
#include "rebound.h"
// const int MAX_TIMESTEP_LEVELS = 6;

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void reb_integrator_leapfrog_ada_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
    const unsigned int N = r->N;

    struct reb_integrator_leapfrog_ada* const rlf = &(r->ri_lfada);
    if ((r->ri_lfada.dt_global/r->ri_lfada.dt_factor != r->dt) && (r->status == REB_STATUS_RUNNING)) {
        reb_exit("Error: dt_global and dt_factor do not match dt. Please set dt = dt_global/dt_factor to use adaptive leapfrog integrator.");
    }    
    const double dt = r->dt;
    if (rlf->par_soft_acc == NULL || rlf->par_jerk == NULL) {
        r->ri_lfada.par_jerk = realloc(r->ri_lfada.par_jerk, sizeof(struct reb_vec3d) * r->N_allocated);
        r->ri_lfada.par_soft_acc = realloc(r->ri_lfada.par_soft_acc, sizeof(struct reb_vec3d) * r->N_allocated);
    }
    // if (rlf->i_colliders == NULL) {
    //     rlf->i_colliders = malloc(sizeof(int) * N);
    //     rlf->N_allocated_colliders = N;
    // }

    struct reb_particle* restrict const particles = r->particles;
    
#pragma omp parallel for schedule(guided)
    for (unsigned int i=0;i<N;i++){
        particles[i].x  += 0.5* dt * particles[i].vx;
        particles[i].y  += 0.5* dt * particles[i].vy;
        particles[i].z  += 0.5* dt * particles[i].vz;
    }
    r->t += dt / 2.;
}

void reb_integrator_leapfrog_ada_part2(struct reb_simulation* r){
    const unsigned int N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    struct reb_integrator_leapfrog_ada* const rlf = &(r->ri_lfada);
    const double dt = r->dt;

    // t_last_acc é o último instante em que a força foi calculada via árvore
    if (r->ri_lfada.step_counter == 0) {
        rlf->t_last_acc = r->t;  // tempo em que é calculado a última acc: manter atualizado
    }
    
    const double t_rel = r->t - rlf->t_last_acc;
#pragma omp parallel for schedule(guided)
    for (unsigned int i = 0; i < N; i++){
        double ax_corr = particles[i].ax + rlf->par_jerk[i].x * t_rel;
        double ay_corr = particles[i].ay + rlf->par_jerk[i].y * t_rel;
        double az_corr = particles[i].az + rlf->par_jerk[i].z * t_rel;
        // if (particles[i].hash != rlf->hash_particles[i]) {
        //     printf("particle.hash = %u | jerk.particle.hash = %u\n", particles[i].hash, rlf->hash_particles[i]);
        // }

        particles[i].vx += dt * ax_corr;
        particles[i].vy += dt * ay_corr;
        particles[i].vz += dt * az_corr;

        if (rlf->mode) {
            particles[i].vx += dt * rlf->par_soft_acc[i].x;
            particles[i].vy += dt * rlf->par_soft_acc[i].y;
            particles[i].vz += dt * rlf->par_soft_acc[i].z;
        }

        particles[i].x  += 0.5 * dt * particles[i].vx;
        particles[i].y  += 0.5 * dt * particles[i].vy;
        particles[i].z  += 0.5 * dt * particles[i].vz;
    }
    r->t += dt / 2.;
    r->dt_last_done = r->dt;
}

// void reb_integrator_leapfrog_ada_part2(struct reb_simulation* r){
//     const unsigned int N = r->N;
//     struct reb_particle* restrict const particles = r->particles;
//     const double dt = r->dt;
// #pragma omp parallel for schedule(guided)
//     for (unsigned int i=0;i<N;i++){
//         particles[i].vx += dt * particles[i].ax;
//         particles[i].vy += dt * particles[i].ay;
//         particles[i].vz += dt * particles[i].az;
//         particles[i].x  += 0.5* dt * particles[i].vx;
//         particles[i].y  += 0.5* dt * particles[i].vy;
//         particles[i].z  += 0.5* dt * particles[i].vz;
//     }
//     r->t+=dt/2.;
//     r->dt_last_done = r->dt;
// }

void reb_integrator_leapfrog_ada_synchronize(struct reb_simulation* r){
    // Do nothing.
}


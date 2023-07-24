#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "defs.h"

#define Jval 1.0
#define PI2 6.283185307179586476925286766559

double IsingProp_CMCMOVE(System *s, Latt *l, double T, int RND_SITE, double RND_TRIAL) {
//CMCMOVE (CYCLIC).

        srand(time(NULL));

        double dE = 0.;
        double coup = 0.;
        int randsite = RND_SITE;
	double randy = RND_TRIAL;
	s->sim_sites[randsite].Npicks++;

        int s_val = l->coord_arr[randsite].spin_val;
        int sk_val = s_val*(-1);

        for (int i = 0; i < l->nb_cnt; ++i)
                coup += l->coord_arr[ l->coord_arr[randsite].nb_arr[i] ].spin_val;

        dE = (2.0*(-Jval)*sk_val)*coup;

        if (dE <= 0.) {

                l->coord_arr[randsite].spin_val = sk_val;
		s->sim_sites[randsite].Nflips_acpt++;
        }

        else {
                double bman = exp(-dE*T);

                if (randy < bman) {

                        l->coord_arr[randsite].spin_val = sk_val;
			s->sim_sites[randsite].Nflips_acpt++;
                }

                else {
                        l->coord_arr[randsite].spin_val = s_val;
			s->sim_sites[randsite].Nflips_rjct++;
                }
        }

        return dE;
}

double IsingProp_RMCMOVE(System *s, Latt *l, double T, int RND_SITE, double RND_TRIAL) {
//RMCMOVE_(RNDM).

	srand(time(NULL));

        double dE = 0.;
        double coup = 0.;
        int coup_cnt = 0;

        int randsite = RND_SITE;
	double randy = RND_TRIAL;
        s->sim_sites[randsite].Npicks++;

        int s_val = l->coord_arr[randsite].spin_val;
        int sk_val = s_val*(-1);

        //CENTRAL ATOM (NOT BOUNDARY).
        if (l->coord_arr[randsite].edge_val == 0) {

                for (int i = 0; i < l->nb_cnt; ++i) {
                        coup += l->coord_arr[ l->coord_arr[randsite].nb_arr[i] ].spin_val;
                        coup_cnt++;
                }

        //VERTEX ATOM
        } else if (l->coord_arr[randsite].edge_val == 5) {

                for (int i=0; i<l->nb_cnt; ++i) {

                        if (l->coord_arr[ l->coord_arr[randsite].nb_arr[i] ].edge_val == 5) {

                                int dummy = l->edge_states[rand() % l->num_edgeatms];
                                coup = coup + l->coord_arr[dummy].spin_val;
                                coup_cnt++;

                        } else {
                                coup = coup + l->coord_arr[l->coord_arr[randsite].nb_arr[i]].spin_val;
                                coup_cnt++;
                        }
                }

        // EDGE ATOM 
        } else {

                for (int i=0; i<l->nb_cnt; i++) {

                        //Check if neighbours have the same edgeval, i.e. do they lie in the same line.
                        if (l->coord_arr[randsite].edge_val == l->coord_arr[ l->coord_arr[randsite].nb_arr[i] ].edge_val) {

                                // exchange with sites in the same line.
                                coup = coup + l->coord_arr[l->coord_arr[randsite].nb_arr[i]].spin_val;
                                coup_cnt++;

                        } else {

                                // does the next site lie in the center of the lattice, or do we cross a BOUNDARY?
                                if (l->coord_arr[ l->coord_arr[randsite].nb_arr[i] ].edge_val == 0) {

                                        coup = coup + l->coord_arr[l->coord_arr[randsite].nb_arr[i]].spin_val;
                                        coup_cnt++;

                                // In this branch, the remaining site has to be across a lattice edge, so we spin the wheel.
                                } else {

                                        int dummy = l->edge_states[rand() % l->num_edgeatms];
                                        coup = coup + l->coord_arr[dummy].spin_val;
                                        coup_cnt++;
                                }
                        }
                }
        }

        dE = (2.0*(-Jval)*sk_val)*coup;

        //if E. is lowered, flip the spin.
        if (dE <= 0.) {

                l->coord_arr[randsite].spin_val = sk_val;
                s->sim_sites[randsite].Nflips_acpt++;
        }

        else {
                double bman = exp(-dE*T);

                if (randy < bman) {

                        l->coord_arr[randsite].spin_val = sk_val;
                        s->sim_sites[randsite].Nflips_acpt++;
                }

                else {
                        l->coord_arr[randsite].spin_val = s_val;
                        s->sim_sites[randsite].Nflips_rjct++;
                }
        }

        return dE;
}

double XYProp_MCMOVE(Latt *l, double T, int RND_SITE, double RND_TRIAL, int MPI_RANK) {

        srand(time(NULL));

        double dE = 0.;
        double coup1 = 0.;
        double coup2 = 0.;
        int randsite = RND_SITE;
	double randy = RND_TRIAL;

        double s_val = l->coord_arr[randsite].xy_spin_val;
        double sk_val = ((double)rand()*PI2)/(double)RAND_MAX;

        for (int j=0; j < l->nb_cnt; j++) {

                coup1 += cos(s_val - l->coord_arr[l->coord_arr[randsite].nb_arr[j]].xy_spin_val);
                coup2 += cos(sk_val - l->coord_arr[l->coord_arr[randsite].nb_arr[j]].xy_spin_val);
        }

        dE = (-Jval*coup2) + (Jval*coup1);
        //printf("dE = %.5f\n", dE);

        if (dE <= 0.) {

                l->coord_arr[randsite].xy_spin_val = sk_val;

        } else {

                double bman = exp(-dE*T);

                if (randy < bman) {

                        l->coord_arr[randsite].xy_spin_val = sk_val;

                } else {

                        l->coord_arr[randsite].xy_spin_val = s_val;
                }
        }

        return dE;
}

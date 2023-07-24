#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rstat.h>

#include "defs.h"
#include "utilities.h"

#define Jval 1.0
#define PI2 6.283185307179586476925286766559

void IsingSpinArr_init(Latt *l) {

        for (int i = 0; i < l->Nsites; ++i) {

                srand((int)time(NULL) + (i+1));

                if (rand()%2 == 0)
                        l->coord_arr[i].spin_val = 1;
                else
                        l->coord_arr[i].spin_val = -1;
        }

        return;
}

void XYSpinArr_init(Latt *l) {

        for (int i = 0; i < l->Nsites; i++) {
                srand((time(NULL)*i+1));
                l->coord_arr[i].xy_spin_val = ((double)rand()*PI2)/(double)RAND_MAX;
        }

        return;
}

void IsingBOUNDS_find_edges(Latt *l) {
/*      EDGE VALUES:
        1 - top, 2 - left, 3 - bottom, 4 - right, 5 - vertex */
        double x_max = 0.;
        double y_max = 0.;

// First loop finds all the positions on the left and bottom side of the square, and finds XMAX & YMAX.
// Second loop finds all the vertices of the square and other edges.

        for (int i=0; i<l->Nsites; i++) {

                if (l->coord_arr[i].pos[0] == 0.0 && l->coord_arr[i].pos[1] == 0.0) {

                        l->coord_arr[i].edge_val = 5;

                } else if (l->coord_arr[i].pos[0] == 0.0 && l->coord_arr[i].pos[1] > 0.) {

                        l->coord_arr[i].edge_val = 2;

                        if (l->coord_arr[i].pos[1] > y_max)
                                y_max = l->coord_arr[i].pos[1];

                } else if (l->coord_arr[i].pos[1] == 0.0 && l->coord_arr[i].pos[0] > 0.) {

                        l->coord_arr[i].edge_val = 3;

                        if (l->coord_arr[i].pos[0] > x_max)
                                x_max = l->coord_arr[i].pos[0];
                }
        }

        for (int i=0; i<l->Nsites; i++) {

                if (l->coord_arr[i].pos[0] == 0.0 && l->coord_arr[i].pos[1] == y_max) {

                        l->coord_arr[i].edge_val = 5;

                } else if (l->coord_arr[i].pos[0] == x_max && l->coord_arr[i].pos[1] == 0.0) {

                        l->coord_arr[i].edge_val = 5;

                } else if (l->coord_arr[i].pos[0] == x_max && l->coord_arr[i].pos[1] == y_max) {

                        l->coord_arr[i].edge_val = 5;

                } else if (l->coord_arr[i].pos[0] > 0.0 && l->coord_arr[i].pos[1] == y_max) {

                        l->coord_arr[i].edge_val = 1;

                } else if (l->coord_arr[i].pos[0] == x_max && l->coord_arr[i].pos[1] > 0.) {

                        l->coord_arr[i].edge_val = 4;
                }
        }

        return;
}

int IsingBOUNDS_update_edges(Latt *l) {

        l->edge_states = (int*)calloc(((2*l->length) + (2*(l->length-2))), sizeof(int));
        int cnt = 0;

        for (int i=0; i<l->Nsites; i++) {

                if (l->coord_arr[i].edge_val != 0) {
                        l->edge_states[cnt] = l->coord_arr[i].index;
                        cnt++;
                }
        }
        return cnt;
}

double IsingProp_En(Latt *l) {

        double energy = 0.;
        double e_site;
        double coup;

        for (int i = 0; i < l->Nsites; ++i) {

                e_site = 0.;
                coup = 0.;

                for (int j = 0; j < l->nb_cnt; ++j) {
                        coup += l->coord_arr[l->coord_arr[i].nb_arr[j]].spin_val;
                }
                e_site = (-Jval) * coup * l->coord_arr[i].spin_val;
                energy += e_site;
        }

        return energy/2.0;
}

double IsingProp_Mag(Latt *l, bool MODEL) {

        double mag = 0.;

        if (MODEL == 0) {
                for (int i = 0; i < l->Nsites; ++i)
                        mag += l->coord_arr[i].spin_val;

        } else {
                for (int i = 0; i < l->Nsites; ++i)
                        mag += l->coord_arr[i].xy_spin_val;
        }

        return mag;
}

void IsingProp_Corr(Latt *l, Task *Tsch, double dR, int dummy,  double am) {
/* Gc = <sisj> - <si><sj> */

        double pair_corr = 0.;
        double dist_cur = 0.;
        int cnt_r = 0;

        am = pow(am,2);

        double r_init = 0.;

        //<sisj> 
        for (int i=0; i <= dummy; i++) {

                r_init += dR;
                Tsch->dists[i] = r_init;
                pair_corr = 0.;

                for (int j=0; j<l->Nsites; j++) {

                        for (int k=1; k<l->Nsites-1; k++) {

                                if (j != k) {

                                        dist_cur = IsingLatt_dist(l, j, k);

                                        if (dist_cur <= r_init) {

                                                pair_corr += ((l->coord_arr[j].spin_val*l->coord_arr[k].spin_val) - am);
                                                cnt_r++;
                                        }
                                }
                        }
                }

                if (cnt_r == 0) {

                        Tsch->corr[i] = 0.;

                } else { Tsch->corr[i] = (1.0/cnt_r)*pair_corr; }
        }

        return;
}

double XYProp_En(Latt *l) {

        double energy = 0.;
        double e_site;

        for (int i = 0; i < l->Nsites; ++i) {

                e_site = 0.;

                for (int j = 0; j < l->nb_cnt; ++j)
                        e_site += cos(l->coord_arr[i].xy_spin_val - l->coord_arr[l->coord_arr[i].nb_arr[j]].xy_spin_val);

                e_site = e_site*(-Jval);
                energy += e_site;
        }

        return energy/2.0;
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "defs.h"
#include "utilities.h"

#define Jval 1.0
#define PI2 6.283185307179586476925286766559

// Memory allocation & File reading //
void IsingLatt_params(char *g_path, Latt *l) {
/* MEMORY LEAK HERE... */

        FILE *fp;
        char dummy[1024];
        l->fpth = g_path;
        l->PBC[0] = 1;
        l->PBC[1] = 1;
        l->PBC[2] = 0;

        fp = fopen(l->fpth, "r");
        fscanf(fp, "%s", dummy);
        fscanf(fp, "%d", &l->length);
        fscanf(fp, "%d", &l->Nsites);
        fscanf(fp, "%s", dummy);
        fscanf(fp, "%lf", &l->a);
        fscanf(fp, "%lf", &l->b);
        fscanf(fp, "%lf", &l->c);
        fscanf(fp, "%s", dummy);
        fclose(fp);

        l->num_edgeatms = (2*l->length) + (2*(l->length-2));
        //printf("%s\n", dummy);
        return;
}

void IsingSpinArr_alloc(int len, Latt *l) {

        l->coord_arr = (Atm*)malloc(l->Nsites*sizeof(Atm));

        return;
}

void IsingLatt_coords(Latt *l) {

        FILE *fp;
        char dummy[1024];
        int skip = 8;

        fp = fopen(l->fpth, "r");

        for (int i = 0; i < skip; ++i) 
                fscanf(fp, "%s", dummy);

        for(int j = 0; j < l->Nsites; ++j) {

                fscanf(fp, "%s", l->coord_arr[j].spec);
                fscanf(fp, "%lf", &l->coord_arr[j].pos[0]);
                fscanf(fp, "%lf", &l->coord_arr[j].pos[1]);
                fscanf(fp, "%lf", &l->coord_arr[j].pos[2]);
                l->coord_arr[j].index = j;
        }
        fclose(fp);

        return;
}

void IsingSpinArr_free(int len, Latt *l, bool flag) {

        if (flag == 1) {

                for (int i = 0; i < l->Nsites; ++i) {

                        free(l->coord_arr[i].nb_arr);
                }
        }

        free(l->coord_arr);
        free(l->edge_states);

        return;
}

// Nearest Neighbours //
int IsingLatt_nbchk(Latt *l, double cutrd, int atmIdx, bool flag) {
/* checks the surroundings of atmIdx for neighbours within cutrd. */

        double dist_cur;
        int nb_cnt = 0;

        for (int i = 0; i < l->Nsites; ++i) {
                if (l->coord_arr[i].index != l->coord_arr[atmIdx].index) {
                        dist_cur = IsingLatt_dist(l, l->coord_arr[atmIdx].index, l->coord_arr[i].index);               
                        if (dist_cur < cutrd + 0.1){
//                              printf("ATM (%d, %d), DIST = %.3f\n", l->coord_arr[atmIdx].index, l->coord_arr[i].index, dist_cur);
                                if (flag == 1)
                                        l->coord_arr[atmIdx].nb_arr[nb_cnt] = l->coord_arr[i].index;

                                nb_cnt++;
                                }
                        }
                }
        return nb_cnt;
}

void IsingSpinArr_alloc_nb(Latt *l) {

        for (int i = 0; i < l->Nsites; ++i)
                l->coord_arr[i].nb_arr = (int*)malloc(l->nb_cnt*sizeof(int));

        return;
}

// Geometrical Utilities //
double IsingLatt_dist(Latt *l, int atmIdx1, int atmIdx2) {

        double dist;
        double rx, ry, rz;

        rx = (l->coord_arr[atmIdx1].pos[0] - l->coord_arr[atmIdx2].pos[0]);
        ry = (l->coord_arr[atmIdx1].pos[1] - l->coord_arr[atmIdx2].pos[1]);
        rz = (l->coord_arr[atmIdx1].pos[2] - l->coord_arr[atmIdx2].pos[2]);

        if (l->PBC[0] && l->PBC[1] == 1) {

                rx = rx - (round(rx/l->a) * l->a);
                ry = ry - (round(ry/l->b) * l->b);
        }

        return dist = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));
}

double IsingLatt_shrtst(Latt *l) {
/* finds the shortest distance, from all the interatomic distances in lattice. */

        double *dist_arr = (double*)malloc(l->Nsites*sizeof(double));

        for (int i = 0; i < l->Nsites-1; ++i){

                dist_arr[i] = floorf((IsingLatt_dist(l, l->coord_arr[i].index, l->coord_arr[i+1].index))*1000)/1000;
                //printf("RNDED DOWN = %.4f\n", dist_arr[i]);
                //printf("last distance = (%d, %d)", l->coord_arr[i].index, l->coord_arr[i+1].index);
        }

        double mindist = min(l->Nsites-1, dist_arr);

        free(dist_arr);

        //printf("shortest dist = %.3f\n", mindist);
        return mindist;
}


void Ising_OutputReadings(System *s, char *OUT, double **mag_readings, double **ener_readings, double *T) {

	FILE *fp;

	fp = fopen(OUT, "w");

	for (int i=0; i < s->Tempsteps; i++) {

		fprintf(fp, "T = %.4f\n", 1.0/T[i]);
		fprintf(fp, "t\t        m(t)\t        e(t)\n");

		for (int j=0; j<s->Simtime; j++) {

			fprintf(fp, "%d          %.10f          %.10f\n", j, mag_readings[i][j], ener_readings[i][j]);

		}
	}

	fclose(fp);

	return;
}

void Ising_OutputOBS(char *OUT, double *E, double *err_E, double *M, double *err_M, double *X, double *err_X, double *C, double *err_C, int BETA, int Tsteps, double *T) {

	FILE *fp;

	fp = fopen(OUT, "w");
	fprintf(fp, "   T\t\t     <E>\t    err<E>\t    <M>\t\t    <|M|>\t  err<M>\t    <X>\t\t   err<X>\t     <C>\t    err<C>\n");

	for (int j=0; j < (Tsteps); j++) {

		if (BETA == 0)
			T[j] = 1.0/T[j];

		fprintf(fp, "%.7f       %.7f       %.7f       %.7f       %.7f      %.7f      %.7f      %.7f      %.7f      %.7f\n", T[j], E[j], err_E[j],  M[j], fabs(M[j]), err_M[j], X[j], err_X[j], C[j], err_C[j]);

	}

	fclose(fp);

	return;
}

// MISCELLANEOUS //
void Ising_SiteChk(Latt *l) {

	for (int i=0; i<l->Nsites; i++) {

			printf("site %d:\n", l->coord_arr[i].index);

		for (int j=0; j<l->nb_cnt; j++) {

			printf("NB[%d] = %d, EDGEVAL = %d ",j, l->coord_arr[l->coord_arr[i].nb_arr[j]].index, l->coord_arr[l->coord_arr[i].nb_arr[j]].edge_val);
		}

			printf("\n");
	}
	
	return;
}

void Ising_OutputCnfg(Latt *l, char *out) {

	int dim = l->length;

	int coords[l->Nsites][2];
	int spins[l->Nsites];

	for (int i=0; i<dim; i++) {

		for(int j=0; j<dim; j++) {

			coords[i*dim +j][0] = i;
			coords[i*dim +j][1] = j;
			spins[i*dim +j] = l->coord_arr[i*dim + j].spin_val;
		}
	}

	FILE *fp;

	fp = fopen(out, "w");
	fprintf(fp, "Site (x,y)\t Spin\n");

	for (int i=0; i<l->Nsites; i++)
		fprintf(fp, "%d,%d             %d\n", coords[i][0], coords[i][1], spins[i]);
	fclose(fp);

	return;
}

void IsingSpinArr_print(Latt *l, bool flag) {

        int array[l->Nsites];

        for (int i = 0; i < l->Nsites; ++i)
                array[i] = l->coord_arr[i].spin_val;

        reverse(array, l->Nsites);

        for (int i = 0; i < l->Nsites; ++i) {

                if (flag == 0) {

                        if (i%l->length == 0)
                                printf("\n");

                        printf("(%s, %d)\t", l->coord_arr[i].spin_val > 0 ? "+" : "-", l->coord_arr[i].index);
                }

                else {

                        if (i%l->length == 0)
                                printf("\n");

                        printf("%s\t", l->coord_arr[i].spin_val > 0 ? "+" : "-");
                }
        }

        printf("\n");

        return;
}

double min(int length, const double *arr) {

        double min = arr[0];

        for(int i = 0; i < length; ++i) {
                if (min > arr[i] && arr[i] != 0.)
                        min = arr[i];
        }

        return min;
}

void reverse(int *arr, int length) {

        int *s, c, d;

        s = (int*)malloc(length*sizeof(int));

        if (s == NULL)
                exit(1);

        for ( c = length - 1, d = 0; c >= 0; c--, d++ )
                *(s+d) = *(arr+c);

        for ( c = 0; c < length; c++ )
                *(arr+c) = *(s+c);

        free(s);

        return;
}

void Utils_Tset(double min, double max, int N, double *data, int BETA) {

	double spacing = 0.; 

	if ( max < min ) {

		spacing = (min - max)/(N-1);

		for (int i=0; i<=(N-1); i++) {

			if (i == 0) {

				data[i] = min;

			} else { data[i] = data[i-1] - spacing; }
		}

	} else {

		spacing = (max - min)/(N-1);

		for (int i=0; i<=(N-1); i++) {

			if (i == 0) {

				data[i] = min;

			} else { data[i] = data[i-1] + spacing; }
		}
	}

	//printf("Tmin = %.4f, Tmax = %.4f, spacing = %.4f\n", min, max, spacing);

	if (BETA == 0) {

		for (int i=0; i < (N); i++)
			data[i] = 1.0/data[i];

	}

	return;
}

void IsingSpinSum(Latt *l) {

        int sum_up = 0;
        int sum_dwn = 0;

        for (int i = 0; i < l->Nsites; ++i) {

                if (l->coord_arr[i].spin_val > 0)
                        sum_up++;
                else
                        sum_dwn++;
        }

        printf("SPINUP = %d, SPINDOWN = %d\n", sum_up, sum_dwn);

        return;
}

void swap(double *a, double *b) {

	double dummy = *a;

	*a = *b;

	*b = dummy;

	return;
}

void IsingSpinSnap(Latt *l, char *fname) {

        FILE *fp;
        fp = fopen(fname, "w");

        for (int i = 0; i < l->length; i++) {

                for (int j = 0; j < l->length; j++) {

                        fprintf(fp, ",%d", l->coord_arr[ i + j ].spin_val);
                }

        fprintf(fp, "\n");

        }

        fclose(fp);

        return;
}

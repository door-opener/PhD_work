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

void Statistics_Qalloc(System *s) {

		s->q_readings = (double**)malloc(s->Tempsteps*sizeof(double*));

		for (int i=0; i<s->Tempsteps; i++) {
//			printf("alloc indx = %d\n", i);
			s->q_readings[i] = (double*)malloc(s->Simtime*sizeof(double));
		}


	return;
}

void Statistics_Qfree(System *s, int MPI_RANK) {

		for (int i=0; i<s->Tempsteps; i++) {
//			printf("free indx = %d\n", i);
			free(s->q_readings[i]);
		}

		free(s->q_readings);

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

// MPI Routines //
void TskAllocateT_MPI(Task *Tsch, int numprocs, int Tsteps, int rank, double *Trange) {

        for (int i=0; i < numprocs; i++) {

                for (int j=0; j < Tsteps; j++) {

                        if (rank == Tsch->rank_id && i == rank) {

                                if (i < 1)
                                        Tsch->Trng[j] = Trange[j];
                                else
                                        Tsch->Trng[j] = Trange[(Tsteps*i)+j];
                        }
                }
        }

        return;
}

void TskMemAlloc_MPI(Task *Tsch, System *s, int Tsteps, int MPI_RANK, int Q) {

        Tsch->Trng = (double*)calloc(Tsteps, sizeof(double));
        Tsch->E = (double*)calloc(Tsteps, sizeof(double));
        Tsch->M = (double*)calloc(Tsteps, sizeof(double));
        Tsch->X = (double*)calloc(Tsteps, sizeof(double));
        Tsch->C = (double*)calloc(Tsteps, sizeof(double));

	if (Q == 1) {

		Tsch->rq_readings = malloc(Tsteps*sizeof(double*));

//		printf("simtime = %d\n", s->Simtime);

		for (int i=0; i < Tsteps; i++) {
				Tsch->rq_readings[i] = malloc(s->Simtime*sizeof(double));
		}

	}

        Tsch->rank_id = MPI_RANK;

        return;

}

void TskCollectRoot_MPI(Task *Tsch, int Tsteps, double *T, int MPI_RANK, MPI_Status *stat, int numprocs, int MODEL, int BETA) {

	if (MPI_RANK == 0) {

		double *E = (double*)calloc((numprocs*Tsteps), sizeof(double));
		double *M = (double*)calloc((numprocs*Tsteps), sizeof(double));
		double *X = (double*)calloc((numprocs*Tsteps), sizeof(double));
		double *C = (double*)calloc((numprocs*Tsteps), sizeof(double));

		for (int j=0; j<Tsteps; j++) {

			E[j] = Tsch->E[j];
			M[j] = Tsch->M[j];
			X[j] = Tsch->X[j];
			C[j] = Tsch->C[j];
		}

		for (int j=1; j < numprocs; j++) {

			MPI_Recv(&E[(j*Tsteps)-1], Tsteps, MPI_DOUBLE, j, (2020+j), MPI_COMM_WORLD, stat);
			MPI_Recv(&M[(j*Tsteps)-1], Tsteps, MPI_DOUBLE, j, (2120+j), MPI_COMM_WORLD, stat);
			MPI_Recv(&X[(j*Tsteps)-1], Tsteps, MPI_DOUBLE, j, (2220+j), MPI_COMM_WORLD, stat);
			MPI_Recv(&C[(j*Tsteps)-1], Tsteps, MPI_DOUBLE, j, (2320+j), MPI_COMM_WORLD, stat);

		}

		char *out = NULL;

		if (MODEL == 1)
			out = "XY_OBS.out";
		else
			out = "Ising_OBS.out";

		FILE *fp;

		fp = fopen(out, "w");
		fprintf(fp, "   T\t\t <E>\t\t<M>\t\t<|M|>\t\t<X>\t     <C>  \n");

		for (int j=0; j < (Tsteps*numprocs); j++) {

			if (BETA == 0)
				T[j] = 1.0/T[j];

			fprintf(fp, "%.5f       %.5f       %.5f       %.5f       %.5f      %.5f\n", T[j], E[j], M[j],fabs(M[j]), X[j], C[j]);

		}

		fclose(fp);

		free(E);
		free(M);
		free(X);
		free(C);
		
	} else {

		for (int j=1; j < numprocs; j++) {

			if (j == MPI_RANK) {

				MPI_Ssend(&Tsch->E[0], Tsteps, MPI_DOUBLE, 0, (2020+j), MPI_COMM_WORLD);
				MPI_Ssend(&Tsch->M[0], Tsteps, MPI_DOUBLE, 0, (2120+j), MPI_COMM_WORLD);
				MPI_Ssend(&Tsch->X[0], Tsteps, MPI_DOUBLE, 0, (2220+j), MPI_COMM_WORLD);
				MPI_Ssend(&Tsch->C[0], Tsteps, MPI_DOUBLE, 0, (2320+j), MPI_COMM_WORLD);
			}
		}
	}

	return;
}

void TskGatherStatRoot_MPI(System *s, Task *Tsch, int Tsteps, int MPI_RANK, MPI_Status *stat, int numprocs) {

	// Prepare Memory at Root.
	if (MPI_RANK == 0) {
		Statistics_Qalloc(s);
	}

	MPI_Datatype sub_arr, resized_type;

	// Dimensions of send type and recv type.
	int dim_full_arr[2] = {(Tsteps*numprocs), s->Simtime};
	int dim_sub_arr[2] = {Tsteps, s->Simtime};
	int start_coords[2] = {0, 0};

	if (MPI_RANK != 0)
		start_coords[0] = (MPI_RANK*Tsteps) - 1;

	printf("rank = %d Starting position in global = [%d, %d]\n", MPI_RANK, start_coords[0], start_coords[1]);

	MPI_Type_create_subarray(2, dim_full_arr, dim_sub_arr, start_coords, MPI_ORDER_C, MPI_DOUBLE, &sub_arr);
	MPI_Type_create_resized(sub_arr, 0, s->Simtime*sizeof(double), &resized_type);
	MPI_Type_commit(&resized_type);

	//rcounts => number of subarrays being sent.
	//displs => diplacement relative to (0,0) in global array in units of block extent (defined as s->Simtime*double)
	int displs[numprocs];
	int rcounts[numprocs];

	for (int i=0; i<numprocs; i++) {

		if (i == 0) {

			displs[i] = 0;
			rcounts[i] = 1;

		} else {

			displs[i] = (int)(i*Tsteps) - 1;
			rcounts[i] = 1;
		}

		if (MPI_RANK == 0)
			printf("rank[%d]: displacement = %d, rcount = %d\n", i, displs[i], rcounts[i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);	

	//MPI_Gatherv(&Tsch->rq_readings[0][0], (Tsteps*s->Simtime)*numprocs, MPI_DOUBLE, &s->q_readings[0][0], rcounts, displs, resized_type, 0, MPI_COMM_WORLD);

	return;
}

// MISCELLANEOUS //
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

void IsingTempset(double *array, double Tinit, double Tfin, int Tsteps, int MPI_RANK, bool BETA) {

        double dT = (Tinit - Tfin)/Tsteps;
        double Tcur = 0.;

        if (MPI_RANK == 0) {

                printf("Temp step = %.3f\n", dT);

        }

        for (int i = 0; i < Tsteps; ++i) {

                if (i == 0) {

                        if (BETA == 1) {
                                array[i] = Tinit;
                                Tcur = Tinit - dT;
                        }

                        else {
                                array[i] = 1.0/Tinit;
                                Tcur = Tinit - dT; }
                }

                else {

                        if (BETA == 1) {
                                array[i] = Tcur;
                        } else {

                        array[i] = 1.0/Tcur; }
                }

                Tcur = Tcur - dT;
        }
        return;
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rstat.h>

#include "defs.h"
#include "statistics.h"
#include "statmech.h"
#include "utilities.h"
#include "moveclass.h"

#define Jval 1.0
#define PI 3.14159265358979323846

/* NB. INDEX,COORDS & SPIN of AN ATOM = l->coord_arr[i].index, l->coord_arr[i].pos[0-2], l->coord_arr[i].spin_val,

   CORR. FOR NEIGHBOURS. = l->coord_arr[l->coord_arr[i].nb_arr[j]].index,l->coord_arr[l->coord_arr[i].nb_arr[j]].pos, l->coord_arr[l->coord_arr[i].nb_arr[j]].spin_val, */

/* NEED TO ADD: 

	- BINDER CUMULANT!
	- NB. WHEN USING XY, DON'T FORGET TO CHANGE BETA!.
*/

/* VERSION: 11/11/2019 */
/* MPI-VERSION: 12/10/2019 */

int main(int argc, char* argv[]) {

	int rank, numprocs;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	MPI_Status stat;

	char g_path[1024];
	strcpy(g_path, argv[1]);

	unsigned int seed = time(NULL);	

	Latt *l_test;
	l_test = (Latt*)malloc(sizeof(Latt));
	
	IsingLatt_params(g_path, l_test);

	//BOUNDARY CONDITIONS ARE SET HERE!
	l_test->boundaries = 1;

	/* if 1 - USE XY MODEL! */
	bool MODEL = 0; 

	if (rank == 0) {

		if(MODEL == 1)
			printf("Spin Model = XY\n");
		else
			printf("Spin Model = Ising\n");

		printf("Lattice Natms = %d\n", l_test->Nsites);
		printf("Lattice length = %d\n", l_test->length);
		printf("Lattice a,b,c = %.3lf, %.3lf, %.3lf\n", l_test->a, l_test->b, l_test->c);
		printf("PBC = [%d %d %d]\n", l_test->PBC[0], l_test->PBC[1], l_test->PBC[2]);
		printf("Number of atoms sitting on edges = %d\n", l_test->num_edgeatms);
		
		if (l_test->boundaries == 1) {
			printf("BOUNDS SETTING: Exchange is RANDOM at boundaries.\n");

		} else { printf("BOUNDS SETTING: Default\n"); }

	}

	IsingSpinArr_alloc(l_test->length, l_test);

	IsingLatt_coords(l_test);
	
	double nb_cutoff = IsingLatt_shrtst(l_test);

	// CHECK ONE SITE TO CALCULATE NEIGHBOURS
	l_test->nb_cnt = IsingLatt_nbchk(l_test, nb_cutoff, l_test->coord_arr[0].index, 0);

	if (rank == 0)	
		printf("Neighbour count = %d\n", l_test->nb_cnt);

	IsingSpinArr_alloc_nb(l_test);

	for (int i = 0; i < l_test->Nsites; ++i)
		IsingLatt_nbchk(l_test, nb_cutoff, i, 1);

	IsingBOUNDS_find_edges(l_test);
	IsingBOUNDS_update_edges(l_test);

	if (MODEL == 1) {
		XYSpinArr_init(l_test);
	} else { IsingSpinArr_init(l_test); }

	if (rank == 0) {

		printf("Initial lattice energy = %.3f\n", IsingProp_En(l_test));
		printf("Initial Magnetisation = %.4f\n", IsingProp_Mag(l_test, MODEL));
	}

	long long int mcsweeps = 1000;
	long long int eqmsweeps = 750;
	long long int simtime = 0;

	long long int mcsteps = l_test->Nsites*mcsweeps;
	long long int flush_thresh = mcsteps/10;
	long long int eqm_steps = l_test->Nsites*eqmsweeps;
	long long int ens = 3;

	/* if beta == 0, Tinit is in Kelvin, otherwise in units of Beta */
	bool beta = 0; 
	
	double Tinit = 5.0;
	double Tfin = 0.5;
	unsigned int Tsteps = 100;
	double *T = (double*)malloc(Tsteps*sizeof(double));

	/* Prepare Sim. Statistics */
	System *latt_stats;
        latt_stats = (System*)malloc(sizeof(System));
        Statistics_Init(latt_stats, l_test, mcsweeps, Tsteps);

	IsingTempset(T, Tinit, Tfin, Tsteps, rank, beta);

	/* MPI PART I - Partition & Allocation*/

	Tsteps = Tsteps/numprocs;

	Task *TaskSch;
	TaskSch = (Task*)malloc(sizeof(Task));

	/* If Q == 1, then alloc memory for statistics */
	TskMemAlloc_MPI(TaskSch, latt_stats, Tsteps, rank, 1);
	TskAllocateT_MPI(TaskSch, numprocs, Tsteps, rank, &T[0]);

	MPI_Barrier(MPI_COMM_WORLD);

	for (int k = 0; k < numprocs; k++) {

		for (int j = 0; j < Tsteps; j++) {

			if (rank == k)
				printf("PROC[%d] T[%d] = %.10f\n", TaskSch->rank_id, j, 1.0/TaskSch->Trng[j]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* MPI PART I END */

	if (rank == 0) {

		printf("Starting MC...\n"); 
		printf("EQM PERIOD = %lld\n", eqm_steps);
		printf("MC STEPS = %lld\n", mcsteps);
		printf("MPI_Tasks per proc = %d\n", Tsteps);

	}

	double E1, M1, E2, M2, AM;
	double Tsq = 0.;
	double mag, ener = 0.;
	double abs_mag = 0.;
	int i = 0;

	double norm1 = 1.0/(mcsteps*l_test->Nsites);
	double norm2 = 1.0/(pow(mcsteps,2)*l_test->Nsites);
	double norm3 = 1.0/ens;

	for (i = 0; i < Tsteps; ++i) {

		if (rank == 0)
			printf("Current MPI_Task %d of %d.\n", (i+1), Tsteps);

		for (int l = 0; l < ens; l++) {
			
			if (rank == 0)
				printf("Ensemble %d of %lld.\n", l+1, ens);

			E1 = 0.; M1 = 0.; E2 = 0.; M2 = 0.; Tsq = 0.; AM = 0.;
	
			for (long long int k = 0; k <= (eqm_steps + mcsteps); k++) {

				if (k== 0) { simtime = 0 - eqmsweeps; } else if (k%l_test->Nsites == 0) { simtime++; }

				if (k == flush_thresh)
					fflush(stdin);	

				IsingProp_RMCMOVE(latt_stats, l_test, TaskSch->Trng[i], (seed*k));
				//IsingProp_CMCMOVE(l_test, TaskSch->Trng[i], ((seed)*k*TaskSch->Trng[i]));
				//XYProp_MCMOVE(l_test, TaskSch->Trng[i], ((seed)*k*TaskSch->Trng[i]), rank);

				if (simtime >= 0) {

					mag = IsingProp_Mag(l_test, MODEL);
					abs_mag = fabs(mag);
					//ener = XYProp_En(l_test);
					ener = IsingProp_En(l_test);
	
					E1 += ener;
					M1 += mag;
					AM += abs_mag;
					M2 += ( pow(mag, 2) );
					E2 += ( pow(ener, 2) );

					//keeps track of sweep number
					if (k%l_test->Nsites == 0) {

						//DO NOT TOUCH THIS IF STATEMENT!
						if (simtime != mcsweeps) {
							TaskSch->rq_readings[i][simtime] = mag;
						}
					}

					}
				}
			
			Tsq = pow(TaskSch->Trng[i], 2);

			TaskSch->E[i] += (E1*norm1*norm3);
			TaskSch->M[i] += (M1*norm1*norm3);
			TaskSch->X[i] += (((M2*norm1) - (pow(AM,2)*norm2)) * TaskSch->Trng[i])*norm3;
			TaskSch->C[i] += (((E2*norm1) - (pow(E1,2)*norm2)) * Tsq)*norm3;
			//TaskSch->C[i] += (((E2*norm1) - (pow(E1,2)*norm2)))*norm3;

			}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* MPI Part II Agglomeration and Output */
	TskCollectRoot_MPI(TaskSch, Tsteps, &T[0], rank, &stat, numprocs, MODEL, beta);

	MPI_Barrier(MPI_COMM_WORLD);
	
	TskGatherStatRoot_MPI(latt_stats, TaskSch, Tsteps, rank, &stat, numprocs);

	if (rank == 0)
		printf("Communication finished!\n");

/*	if (rank == 0) {

		int dummy = (Tsteps*numprocs);

		for (int j=0; j<dummy; j++) {

			printf("T[%d] = %.2f\n",j, T[j]);
			printf("|  t  |\t |  q(t)  |\n");

			for (int k=0; k<latt_stats->Simtime; k++) {

				printf("|  %d  |\t |  %.3f  |\n", k+1, latt_stats->q_readings[j][k]);

			}

		}

	} */

	MPI_Barrier(MPI_COMM_WORLD);

	/* MPI_CLEAN */

	for (int j=0; j < Tsteps; j++)
		free(TaskSch->rq_readings[j]);

	free(TaskSch->rq_readings);

	free(TaskSch->Trng);
	free(TaskSch->E);
	free(TaskSch->M);
	free(TaskSch->X);
	free(TaskSch->C);

	MPI_Barrier(MPI_COMM_WORLD);

	/* Stat clean */
	//if (rank == 0)
	//	Statistics_Qfree(latt_stats, rank);

	free(latt_stats->sim_sites);
        free(latt_stats);

	IsingSpinArr_free(l_test->length, l_test, 1);

	free(T);
	free(l_test);

	if (rank == 0) {

		printf("DEALLOC SUCCESS!\n");
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

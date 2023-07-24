#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_rng.h>

#include "defs.h"
#include "statistics.h"
#include "statmech.h"
#include "utilities.h"
#include "moveclass.h"

/* NB. INDEX,COORDS & SPIN of AN ATOM = l->coord_arr[i].index, l->coord_arr[i].pos[0-2], l->coord_arr[i].spin_val,
 *
 *    CORR. FOR NEIGHBOURS. = l->coord_arr[l->coord_arr[i].nb_arr[j]].index,l->coord_arr[l->coord_arr[i].nb_arr[j]].pos, l->coord_arr[l->coord_arr[i].nb_arr[j]].spin_val, */

int main(int argc, char **argv) {

	int rank, numprocs;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	MPI_Status stat;

        char g_path[1024];
        strcpy(g_path, argv[1]);

        unsigned long int seed = time(NULL) + rank;

        //RNG FOR SITE RANDOMNESS
        gsl_rng * r_sites = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(r_sites, seed);

        Latt *l_test;
        l_test = (Latt*)malloc(sizeof(Latt));

        IsingLatt_params(g_path, l_test);

        bool MODEL = 0;

	if (rank == 0) {

		if (MODEL == 1)
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

	if (rank == 0) {

        	printf("Neighbour count = %d\n", l_test->nb_cnt);
	}

        IsingSpinArr_alloc_nb(l_test);

        for (int i = 0; i < l_test->Nsites; ++i)
                IsingLatt_nbchk(l_test, nb_cutoff, i, 1);

        IsingBOUNDS_find_edges(l_test);
        IsingBOUNDS_update_edges(l_test);

        //Ising_SiteChk(l_test);

        IsingSpinArr_init(l_test, r_sites);

	if (rank == 0) {

		Ising_OutputCnfg(l_test, "Initialstate.txt");
		fprintf(stdout, "Initial lattice energy = %.3f\n", IsingProp_En(l_test));
		fprintf(stdout, "Initial Magnetisation = %.4f\n", IsingProp_Mag(l_test, MODEL));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	long long int mcsweeps = 10000;
        long long int eqmsweeps = 150;
        long long int simtime = 0;

        long long int mcsteps = l_test->Nsites*mcsweeps;
        long long int flush_thresh = mcsteps/10;
        long long int eqm_steps = l_test->Nsites*eqmsweeps;
        long long int ens = 200;

        //RNG FOR METROPOLIS TRIALS.
        gsl_rng * r_trials = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(r_trials, seed);
        double *rnd_trials = NULL;
        rnd_trials = calloc((mcsteps + eqm_steps), sizeof(double));

        for (int i=0; i<(mcsteps + eqm_steps); i++)
                rnd_trials[i] = gsl_rng_uniform_pos(r_trials);

        gsl_rng_free(r_trials);

	MPI_Barrier(MPI_COMM_WORLD);

        /* if beta == 0, Tinit is in absT, otherwise in units of Beta */
        bool beta = 0;

        double Tinit = 0.2;
        double Tfin = 5.0;
        unsigned int Tsteps = 96;
        double *T = (double*)malloc(Tsteps*sizeof(double));
        Utils_Tset(Tinit, Tfin, Tsteps, &T[0], beta);

        /* Prepare Sim. Statistics */
        System *latt_stats;
        latt_stats = (System*)malloc(sizeof(System));
        Statistics_Init(latt_stats, l_test, mcsweeps, Tsteps);

	/* MPI PART 1 */
	Tsteps = Tsteps/numprocs;

	Task *TaskSch;
	TaskSch = (Task*)malloc(sizeof(Task));
	TskMemAlloc_MPI(TaskSch, latt_stats, Tsteps, rank, 1);
	TskAllocateT_MPI(TaskSch, numprocs, Tsteps, rank, &T[0]);

	MPI_Barrier(MPI_COMM_WORLD);

	TskShow_MPI(TaskSch, numprocs, Tsteps);

	if (rank == 0) {

		fprintf(stdout, "Starting MC...\n");
		fprintf(stdout, "EQM PERIOD = %lld\n", eqm_steps);
		fprintf(stdout, "MC STEPS = %lld\n", mcsteps);
	}

        double E1, M1, E2, M2, AM;
        double Tsq = 0.;
        double mag, ener = 0.;
        double abs_mag = 0.;
        int RND_SITE = 0;
        int i = 0;

        double nrm_site = 1.0/(l_test->Nsites);
        double norm1 = 1.0/(mcsteps*l_test->Nsites);
        double norm2 = 1.0/(pow(mcsteps,2)*l_test->Nsites);
        double norm3 = 1.0/ens;

	MPI_Barrier(MPI_COMM_WORLD);

        //GSL Running Statistics.
        gsl_rstat_workspace *wspace_mag_stat = gsl_rstat_alloc();
        gsl_rstat_workspace *wspace_ener_stat = gsl_rstat_alloc();

        for (i = 0; i < Tsteps; i++) {

		if (rank == 0) 
                	fprintf(stdout, "Current Tstep %d of %d:\n", (i+1), Tsteps);

                for (int l = 0; l < ens; l++) {

                        //printf("Ensemble %d of %lld.\n", l+1, ens);

                        E1 = 0.; M1 = 0.; E2 = 0.; M2 = 0.; Tsq = 0.; AM = 0.;

                        for (long long int k = 0; k <= (eqm_steps + mcsteps); k++) {

                                if (k == 0) { simtime = 0 - eqmsweeps; } else if (k%l_test->Nsites == 0) { simtime++; }

                                if (k == flush_thresh)
                                        fflush(stdin);

                                RND_SITE = gsl_rng_uniform_int(r_sites, l_test->Nsites);

                                //IsingProp_RMCMOVE(latt_stats, l_test, TaskSch->Trng[i], RND_SITE, rnd_trials[k]);
                                IsingProp_CMCMOVE(latt_stats ,l_test, TaskSch->Trng[i], RND_SITE, rnd_trials[k]);
                                //XYProp_MCMOVE(l_test, T[i], RND_SITE, rnd_trials[k], rank);

                                if (simtime >= 0) {

                                        mag = IsingProp_Mag(l_test, MODEL);
                                        abs_mag = fabs(mag);
                                        //ener = XYProp_En(l_test);
                                        ener = IsingProp_En(l_test);

                                        //GSL workspace.
                                        gsl_rstat_add((mag), wspace_mag_stat);
                                        gsl_rstat_add((ener), wspace_ener_stat);

                                        M1 += mag;
                                        AM += abs_mag;
                                        M2 += ( pow(mag, 2) );

                                        //keeps track of sweep number --> [rowindex*colindex + colindex]
                                        if (k%l_test->Nsites == 0) {

                                                TaskSch->q1_readings[i][simtime] += mag*norm3;
                                                TaskSch->q2_readings[i][simtime] += ener*norm3;

                                                }
                                        }
                                }

                        Tsq = pow(TaskSch->Trng[i], 2);

                        //GSL STATS
                        TaskSch->E[i] += gsl_rstat_mean(wspace_ener_stat)*norm3*nrm_site;
                        TaskSch->M[i] += gsl_rstat_mean(wspace_mag_stat)*norm3*nrm_site;
                        TaskSch->X[i] += (((M2*norm1) - (pow(AM,2)*norm2)) * TaskSch->Trng[i])*norm3;
                        TaskSch->C[i] += (gsl_rstat_variance(wspace_ener_stat) * Tsq)*norm3*nrm_site;

                        //TAU[i] += AutocorrelationIII(latt_stats, l_test, &latt_stats->q1_readings[i][0])*norm3;
                        //TaskSch->C[i] += (((E2*norm1) - (pow(E1,2)*norm2)))*norm3;

                        }

                        gsl_rstat_reset(wspace_mag_stat);
                        gsl_rstat_reset(wspace_ener_stat);
        }

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {

        	Ising_OutputCnfg(l_test, "Finalstate.txt");
	}

	// ERRORS ON EACH PROC //
	for (int i=0; i<Tsteps; i++) {

		TaskSch->err_E[i] = Statistics_ErrBSTRAP(latt_stats, l_test, &TaskSch->q1_readings[i][0], TaskSch->Trng[i], 0);
		TaskSch->err_X[i] = Statistics_ErrBSTRAP(latt_stats, l_test, &TaskSch->q2_readings[i][0], TaskSch->Trng[i], 2);
		TaskSch->err_C[i] = Statistics_ErrBSTRAP(latt_stats, l_test, &TaskSch->q1_readings[i][0], TaskSch->Trng[i], 1);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	// MPI PART II: COMM ROUTINES //
        /* Allocate memory for Observables and Statistics */

        double *E = NULL;
        double *err_E = NULL;
        double *M = NULL; 
        double *X = NULL;  
        double *err_X = NULL;  
        double *C = NULL;
        double *err_C = NULL;

	if (rank == 0) {

		E = (double*)calloc((numprocs*Tsteps), sizeof(double));
		err_E = (double*)calloc((numprocs*Tsteps), sizeof(double));
		M = (double*)calloc((numprocs*Tsteps), sizeof(double));
		X = (double*)calloc((numprocs*Tsteps), sizeof(double));
		err_X = (double*)calloc((numprocs*Tsteps), sizeof(double));
		C = (double*)calloc((numprocs*Tsteps), sizeof(double));
		err_C = (double*)calloc((numprocs*Tsteps), sizeof(double));
	
	}

	TskGatherOBS_MPI(TaskSch, Tsteps, &E[0], &err_E[0], &M[0], &X[0], &err_X[0], &C[0], &err_C[0]);

        // OUTPUT DATA + ERRORS //	
	if (rank == 0) {

	        FILE *fp;
		fp = fopen("IsingOBS.out", "w");
		fprintf(fp, "   T\t\t           <E>\t                err<E>                <M>                  <|M|>                  <X>                 err<X>                   <C>               err<C>\n");
		for (int i=0; i<(Tsteps*numprocs); i++) {

			fprintf(fp, "%12.8f          %12.8f          %12.12f          %12.8f            %12.8f         %12.8f          %12.12f         %12.8f          %12.12f\n", 1.0/T[i], E[i], err_E[i], M[i], fabs(M[i]), X[i], err_X[i], C[i], err_C[i]);

		}
		fclose(fp);
	}

        IsingSpinArr_free(l_test->length, l_test, 1);
        free(latt_stats->sim_sites);
        free(latt_stats);

	// LOCAL FREES //
	if (rank == 0) {

		free(E);
		free(err_E);
		free(M);
		free(X);
		free(err_X);
		free(C);
		free(err_C);
		free(T);
	}

	// MPI_FREES //	
	for (int i=0; i<Tsteps; i++) {

		free(TaskSch->q1_readings[i]);
		free(TaskSch->q2_readings[i]);
	}

	free(TaskSch->Trng);
	free(TaskSch->E);
	free(TaskSch->err_E);
	free(TaskSch->X);
	free(TaskSch->err_X);
	free(TaskSch->C);
	free(TaskSch->err_C);
	free(TaskSch->M);

        free(rnd_trials);
        gsl_rng_free(r_sites);
        gsl_rstat_free(wspace_mag_stat);
        gsl_rstat_free(wspace_ener_stat);

        free(l_test);

	MPI_Finalize();

        return 0;
}

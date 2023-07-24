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

int main(int argc, char **argv) {

	char g_path[1024];
	strcpy(g_path, argv[1]);
	
	unsigned long int seed = time(NULL);

	//RNG FOR SITE RANDOMNESS
	gsl_rng * r_sites = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r_sites, seed);

	Latt *l_test;
	l_test = (Latt*)malloc(sizeof(Latt));

	IsingLatt_params(g_path, l_test);

	bool MODEL = 0; 

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

	IsingSpinArr_alloc(l_test->length, l_test);

	IsingLatt_coords(l_test);
	
	double nb_cutoff = IsingLatt_shrtst(l_test);

	// CHECK ONE SITE TO CALCULATE NEIGHBOURS
	l_test->nb_cnt = IsingLatt_nbchk(l_test, nb_cutoff, l_test->coord_arr[0].index, 0);

	printf("Neighbour count = %d\n", l_test->nb_cnt);

	IsingSpinArr_alloc_nb(l_test);

	for (int i = 0; i < l_test->Nsites; ++i)
		IsingLatt_nbchk(l_test, nb_cutoff, i, 1);

	IsingBOUNDS_find_edges(l_test);
	IsingBOUNDS_update_edges(l_test);

	//Ising_SiteChk(l_test);

	if (MODEL == 1) {
		XYSpinArr_init(l_test);
	} else { IsingSpinArr_init(l_test, r_sites); }

	Ising_OutputCnfg(l_test, "Initialstate.txt");
	printf("Initial lattice energy = %.3f\n", IsingProp_En(l_test));
	printf("Initial Magnetisation = %.4f\n", IsingProp_Mag(l_test, MODEL));
	
	long long int mcsweeps = 1000000;
	long long int eqmsweeps = 150;
	long long int simtime = 0;

	long long int mcsteps = l_test->Nsites*mcsweeps;
	long long int flush_thresh = mcsteps/10;
	long long int eqm_steps = l_test->Nsites*eqmsweeps;
	long long int ens = 5;
	
	//RNG FOR METROPOLIS TRIALS.
	gsl_rng * r_trials = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r_trials, seed);

	double *rnd_trials = NULL;

	rnd_trials = calloc((mcsteps + eqm_steps), sizeof(double));

	for (int i=0; i<(mcsteps + eqm_steps); i++)
		rnd_trials[i] = gsl_rng_uniform_pos(r_trials);

	gsl_rng_free(r_trials); 

	/* if beta == 0, Tinit is in absT, otherwise in units of Beta */
	bool beta = 0; 
	
	double Tinit = 0.5;
	double Tfin = 5.0;
	unsigned int Tsteps = 25;
	double *T = (double*)malloc(Tsteps*sizeof(double));

	/* Prepare Sim. Statistics */
	System *latt_stats;
        latt_stats = (System*)malloc(sizeof(System));
        Statistics_Init(latt_stats, l_test, mcsweeps, Tsteps);

	Utils_Tset(Tinit, Tfin, Tsteps, &T[0], beta);

	/* Allocate memory for Observables and Statistics */
	double *E = calloc(Tsteps, sizeof(double));
	double *M = calloc(Tsteps, sizeof(double));
	double *X = calloc(Tsteps, sizeof(double));
	double *C = calloc(Tsteps, sizeof(double));

	latt_stats->q1_readings = calloc(Tsteps, sizeof(double*));
	latt_stats->q2_readings = calloc(Tsteps, sizeof(double*));

	for (int i=0; i<Tsteps; i++) {

		latt_stats->q1_readings[i] = calloc(mcsweeps, sizeof(double));
		latt_stats->q2_readings[i] = calloc(mcsweeps, sizeof(double));
	}
	
	printf("Starting MC...\n"); 
	printf("EQM PERIOD = %lld\n", eqm_steps);
	printf("MC STEPS = %lld\n", mcsteps);

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

	//GSL Running Statistics.
	gsl_rstat_workspace *wspace_mag_stat = gsl_rstat_alloc();
	gsl_rstat_workspace *wspace_ener_stat = gsl_rstat_alloc();

	for (i = 0; i < Tsteps; i++) {

		printf("Current Tstep %d of %d:\n", (i+1), Tsteps);

		for (int l = 0; l < ens; l++) {
			
			//printf("Ensemble %d of %lld.\n", l+1, ens);

			E1 = 0.; M1 = 0.; E2 = 0.; M2 = 0.; Tsq = 0.; AM = 0.;
	
			for (long long int k = 0; k <= (eqm_steps + mcsteps); k++) {

				if (k == 0) { simtime = 0 - eqmsweeps; } else if (k%l_test->Nsites == 0) { simtime++; }

				if (k == flush_thresh)
					fflush(stdin);

				RND_SITE = gsl_rng_uniform_int(r_sites, l_test->Nsites);

				//IsingProp_RMCMOVE(latt_stats, l_test, T[i], RND_SITE, rnd_trials[k]);
				IsingProp_CMCMOVE(latt_stats ,l_test, T[i], RND_SITE, rnd_trials[k]);
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
					if (k%l_test->Nsites == 0 && l == (ens-1)) {

						latt_stats->q1_readings[i][simtime] = mag;
						latt_stats->q2_readings[i][simtime] = ener;
					}

					}
				}
			
			Tsq = pow(T[i], 2);
		
			//GSL STATS
			E[i] += gsl_rstat_mean(wspace_ener_stat)*norm3*nrm_site;
			M[i] += gsl_rstat_mean(wspace_mag_stat)*norm3*nrm_site;
			X[i] += (((M2*norm1) - (pow(AM,2)*norm2)) * T[i])*norm3;
			C[i] += (gsl_rstat_variance(wspace_ener_stat) * Tsq)*norm3*nrm_site;
			//TaskSch->C[i] += (((E2*norm1) - (pow(E1,2)*norm2)))*norm3;

			}

			gsl_rstat_reset(wspace_mag_stat);
			gsl_rstat_reset(wspace_ener_stat);
	}

	Ising_OutputCnfg(l_test, "Finalstate.txt");

	/* Errors */
	double err_E[Tsteps];
	double err_X[Tsteps];
	double err_C[Tsteps];
	double err_M[Tsteps];

	for (int i=0; i < Tsteps; i++) {

		err_E[i] = Statistics_ErrBSTRAP(latt_stats, l_test, latt_stats->q2_readings[i], T[i], 0);
		err_C[i] = Statistics_ErrBSTRAP(latt_stats, l_test, latt_stats->q2_readings[i], T[i], 1);
		err_X[i] = Statistics_ErrBSTRAP(latt_stats, l_test, latt_stats->q1_readings[i], T[i], 2);
		err_M[i] = Statistics_ErrBSTRAP(latt_stats, l_test, latt_stats->q1_readings[i], T[i], 3);
	}

	/* OUTPUT DATA */
	FILE *fp;
	fp = fopen("IsingOBS.out", "w");
	fprintf(fp, "   T\t\t           <E>\t                err<E>                <M>                err<M>                 <|M|>                  <X>                 err<X>                   <C>               err<C>\n");
	for (int i=0; i<Tsteps; i++) {

		fprintf(fp, "%12.8f          %12.8f          %12.12f          %12.8f          %12.12f         %12.8f         %12.8f          %12.12f         %12.8f          %12.12f\n", 1.0/T[i], E[i], err_E[i], M[i], err_M[i], fabs(M[i]), X[i], err_X[i], C[i], err_C[i]);

	}
	fclose(fp);

	AlgorithmStats(latt_stats, l_test, mcsteps, eqm_steps, ens, 2);
	//Ising_OutputReadings(latt_stats, "simreadings.out", latt_stats->q1_readings, latt_stats->q2_readings, &T[0]);
	IsingSpinArr_free(l_test->length, l_test, 1);

	free(latt_stats->sim_sites);
	free(latt_stats);

	free(E);
	free(M);
	free(X);
	free(C);
	free(T);

	free(rnd_trials);
	gsl_rng_free(r_sites);
	gsl_rstat_free(wspace_mag_stat);
	gsl_rstat_free(wspace_ener_stat);

	free(l_test);

	return 0;
}

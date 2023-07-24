/* SIMULATON STATISTICS FUNCTIONS */

// 14/12/2019 //

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>

#include "defs.h"

#define Jval 1.0
#define PI2 6.283185307179586476925286766559

void Statistics_Init(System *s, Latt *l, long long int Simtime, int Tsteps) {

	s->Nsites = l->Nsites;
	s->sim_sites = (Site*)malloc(s->Nsites*sizeof(Site));
	s->Simtime = (int)Simtime;
	s->Tempsteps = Tsteps;
	s->a = l->a;
	s->b = l->b;
	s->c = l->c;

	for (int i=0; i<l->Nsites; i++) {

		s->sim_sites[i].index = l->coord_arr[i].index;
		s->sim_sites[i].edge_val = l->coord_arr[i].edge_val;
		s->sim_sites[i].Npicks = 0;
		s->sim_sites[i].Nflips_acpt = 0;
		s->sim_sites[i].Nflips_rjct = 0;
	}

	return;
}

//Simulation Stats.
double Statistics_Autocorr(System *s, double *q_arr, int OUT, char *out) {
//NEWMAN BARKEMA (3.21).

	long int tdmax = s->Simtime/2;

	double tau_est = 0.;
	double autocorr[tdmax];
	double M0MT, mt, mtt, norm_td;

	long int dt = 0.;

	for (long int i=0; i<tdmax; i++) {

		M0MT = 0.; mt = 0.; mtt = 0.; norm_td = 0;

		for (long long int j=0; j<s->Simtime; j++) {

			//EQN.2 (3.21)
			M0MT += (q_arr[j]*q_arr[j+dt]);
			mt += (q_arr[j]);
			mtt += (q_arr[j+dt]);

		}

		norm_td = 1.0/((double)(s->Simtime - dt));
		
		autocorr[dt] = (norm_td*M0MT) - ((norm_td*mt)*(norm_td*mtt));
		dt++;
	}

	for (int i=0; i<dt; i++){
		autocorr[i] = autocorr[i]/autocorr[0];
	}

	//compute the integrated correlation time.	
	for (int i=0; i<tdmax; i++)
		if (autocorr[i] > 0) 
			tau_est += autocorr[i];

	printf("Estimated TAU (E1) = %.8f\n", tau_est);

	if (OUT == 1) {

		FILE *fp;

		fp = fopen(out, "w");
		fprintf(fp, "Estimated TAU = %.10f\n", tau_est);
		fprintf(fp, "t\t\t X(t)\n");
		for (int i=0; i<dt; i++)
			fprintf(fp, "%d\t     %12.10f\n", i, autocorr[i]);
		fclose(fp);
	}

	return tau_est;
}

double Statistics_AutocorrII(System *s, double *q_arr, int OUT, char *out) {
//NEWMAN BARKEMA (3.17)

	long int tdmax = s->Simtime/2;

	double tau_est = 0.;
	double autocorr[tdmax];
	double M0MT, norm_td;
	double ams, varm, msq;

	varm = gsl_stats_variance(q_arr, 1, s->Simtime);

	long int dt = 0.;

	for (long int i=0; i<tdmax; i++) {

		M0MT = 0.; norm_td = 0; ams = 0.; msq = 0.;

		for (long long int j=0; j<s->Simtime; j++) {

			M0MT += q_arr[j]*q_arr[j+dt];
			ams += q_arr[j];
		}

		norm_td = 1.0/((double)(s->Simtime - dt));

		ams = pow(ams,2)/pow(((double)s->Simtime), 2);

		M0MT = M0MT*norm_td;
		
		autocorr[dt] = ((M0MT) - (ams))/varm;
		dt++;
	}

	for (int i=0; i<dt; i++){
		autocorr[i] = autocorr[i]/autocorr[0];
	}

	//compute the integrated correlation time.	
	for (int i=0; i<tdmax; i++)
		if (autocorr[i] > 0) 
			tau_est += autocorr[i];

	printf("Estimated TAU (E2) = %.8f\n", tau_est);

	if (OUT == 1) {

		FILE *fp;

		fp = fopen(out, "w");
		fprintf(fp, "Estimated TAU = %.10f\n", tau_est);
		fprintf(fp, "t\t\t X(t)\n");
		for (int i=0; i<dt; i++)
			fprintf(fp, "%d\t     %12.10f\n", i, autocorr[i]);
		fclose(fp);
	}

	return tau_est;
}

double Statistics_AutocorrIII(System *s, double *q_arr, int OUT, char *out) {
//BINNEY NEWMAN (4.47)

	long int tdmax = s->Simtime/2;

	double tau_est = 0.;
	double autocorr[tdmax];
	double norm_td, tdev;
	double am, varm;

	varm = gsl_stats_variance(q_arr, 1, s->Simtime);
	am = gsl_stats_mean(q_arr, 1, s->Simtime);

	long int dt = 0.;

	for (long int i=0; i<tdmax; i++) {

		norm_td = 0; tdev = 0.;

		for (long long int j=0; j<s->Simtime; j++) {

			tdev += (q_arr[j] - am) * (q_arr[j+dt] - am);
		}

		norm_td = (1.0/ ((double)(s->Simtime - dt)));

		autocorr[dt] = (norm_td*tdev)/varm;

		dt++;
	}

	for (int i=0; i<dt; i++){
		autocorr[i] = autocorr[i]/autocorr[0];
	}

	//compute the integrated correlation time.	
	for (int i=0; i<tdmax; i++)
		if (autocorr[i] > 0) 
			tau_est += autocorr[i];

	printf("Estimated TAU (E3) = %.8f\n", tau_est);

	if (OUT == 1) {

		FILE *fp;

		fp = fopen(out, "w");
		fprintf(fp, "Estimated TAU = %.10f\n", tau_est);
		fprintf(fp, "t\t\t X(t)\n");
		for (int i=0; i<dt; i++)
			fprintf(fp, "%d\t     %12.10f\n", i, autocorr[i]);
		fclose(fp);
	}

	return tau_est;
}

double Statistics_AutocorrIV(System *s, double *q_arr, int OUT, char *out) {
//BINNEY NEWMAN (4.47)
//Fit the time displacement with a tolerance on X(t).

	int tdmax = s->Simtime/2;
	double tau_est = 0.;
	double autocorr[tdmax];
	double norm_td, tdev;
	double am, varm;

	double dummy = 1;

	varm = gsl_stats_variance(q_arr, 1, s->Simtime);
	am = gsl_stats_mean(q_arr, 1, s->Simtime);

	long int dtmax = 0;
	long int dt = 0;	

	while (dummy >= 0.05) {

		norm_td = 0; tdev = 0.;

		for (long long int j=0; j<s->Simtime; j++) {

			tdev += (q_arr[j] - am) * (q_arr[j+dtmax] - am);
		}

		norm_td = (1.0/ ((double)(s->Simtime - dtmax)));

		dummy = (norm_td*tdev)/varm;
		dtmax++;
	}

	for (long long int i=s->Simtime; i==0; i = (i - dt)) {

		norm_td = 0.; tdev = 0.;

		for (long long int j=0; j<s->Simtime; j++) {

			tdev += (q_arr[j] - am) * (q_arr[j+dtmax] - am);
		}

		autocorr[dt] = (norm_td*tdev)/varm;

		if (dt < dtmax)
			dt++;
	}

	for (int i=0; i<dt; i++){
		autocorr[i] = autocorr[i]/autocorr[0];
	}

	//compute the integrated correlation time.	
	for (int i=0; i<tdmax; i++)
		if (autocorr[i] > 0) 
			tau_est += autocorr[i];

	printf("Estimated TAU (E3) = %.8f\n", tau_est);

	if (OUT == 1) {

		FILE *fp;

		fp = fopen(out, "w");
		fprintf(fp, "Estimated TAU = %.10f\n", tau_est);
		fprintf(fp, "t\t\t X(t)\n");
		for (int i=0; i<dt; i++)
			fprintf(fp, "%d\t     %.10f\n", i, autocorr[i]);
		fclose(fp);
	}

	return tau_est;
}

//Error Estimation
double Statistics_ErrBSTRAP(System *s, Latt *l, double *q_arr, double T, int FLAG) {

	// if FLAG is 1, compute errors for <C>
	//# of samples is s->Simtime.
	unsigned long int seed = time(NULL);
	int ncycles = 200;

	double Tsq = pow(T,2);

	double Q_ESTIMATES[ncycles];
	double Q_EST_RIGID[ncycles];

	double Q, Qsq;
	double mean = gsl_stats_mean(q_arr, 1, s->Simtime);

	double norm1 = 1.0/(s->Simtime*l->Nsites);
	double norm2 = 1.0/(pow(s->Simtime, 2)*l->Nsites);
	double nrm_site = 1.0/l->Nsites;

	// GSL RNG
	gsl_rng	* r = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r, seed);

	//GSL WSPACE
	gsl_rstat_workspace *error_space = gsl_rstat_alloc();

	for (int i=0; i<ncycles; i++) {

		Q = 0; Qsq = 0;

		for (int j=0; j<s->Simtime; j++) {

			int RND = gsl_rng_uniform_int(r, s->Simtime);

			if (FLAG == 3 || FLAG == 2) { 

				gsl_rstat_add(fabs(q_arr[RND]), error_space); 
				Q += fabs(q_arr[RND]);
				Qsq += pow(q_arr[RND], 2);

			} else { 

				gsl_rstat_add(q_arr[RND], error_space); 
				Q += q_arr[RND];
				Qsq += pow(q_arr[RND], 2);
			}
		}

		if (FLAG == 1) {

			Q_ESTIMATES[i] = (gsl_rstat_variance(error_space)*Tsq)*nrm_site;
			Q_EST_RIGID[i] = ( (Qsq*norm1) - ( pow(Q,2)*norm2 ) ) * Tsq;

		} else if (FLAG == 2) {
	
			Q_ESTIMATES[i] = (gsl_rstat_variance(error_space)*T)*nrm_site;
			Q_EST_RIGID[i] = ( (Qsq*norm1) - ( pow(Q,2)*norm2 ) ) * T;

		} else { 

			Q_ESTIMATES[i] = gsl_rstat_mean(error_space)*nrm_site;
			Q_EST_RIGID[i] = Q*norm1;
		}

		gsl_rstat_reset(error_space);
	}

	double error = gsl_stats_sd(Q_ESTIMATES, 1, ncycles);
	//double error = gsl_stats_sd_m(Q_ESTIMATES, 1, ncycles, mean);
	double error2 = gsl_stats_sd_m(Q_EST_RIGID, 1, ncycles, mean);

	gsl_rng_free(r);
	gsl_rstat_free(error_space);

	return error;
}

//Dynamics Algorithm Statistics.
void AlgorithmStats(System *s, Latt *l, long long int mcsteps, long long int eqm_steps, long long int ens, int FLAG) {

	int coords[l->Nsites][2];
	long long int flip_tot = s->Tempsteps*(ens*(mcsteps + eqm_steps));
	long long int sumNpicks = 0;

	//Normalised statistics.
	double Npicks[l->Nsites];
	double Nacpt[l->Nsites];
	double Nrjct[l->Nsites];

	for (int i=0; i<l->length; i++) {

		for (int j=0; j<l->length; j++) {

			coords[i*l->length +j][0] = i;
			coords[i*l->length +j][1] = j;

			sumNpicks += s->sim_sites[i*l->length + j].Npicks;
		}
	}

	//Absolute stats.
	if (FLAG == 1) {

		fprintf(stdout, "Site(x,y)\t\t Npicks\t\t Naccepted\t\t Nrejected\t\t Edgeval\n");

		for (int i=0; i<l->Nsites; i++) {


			printf("%d,%d\t\t         %ld\t\t   %ld\t\t   %ld\t\t     %d\n", coords[i][0], coords[i][1], s->sim_sites[i].Npicks, s->sim_sites[i].Nflips_acpt, s->sim_sites[i].Nflips_rjct, s->sim_sites[i].edge_val);
		}

		printf("%s = %lld\n", "Total Flips", flip_tot);
		printf("%s = %lld\n", "Sum Npicks", sumNpicks);

	//Normalised stats.
	} else if (FLAG == 2) {

		FILE *fp;
		fp = fopen("algstats_norm.out", "w");

		for (int i=0; i<l->Nsites; i++) {

			Npicks[i] = ((double)s->sim_sites[i].Npicks/(double)sumNpicks);
			Nacpt[i] = ((double)s->sim_sites[i].Nflips_acpt/(double)s->sim_sites[i].Npicks);
			Nrjct[i] = ((double)s->sim_sites[i].Nflips_rjct/(double)s->sim_sites[i].Npicks);

		}

		fprintf(fp, "Total Flips = %lld\n", sumNpicks);
		fprintf(fp, "Site(x,y)\t\t   Npicks\t\t   Naccepted\t\t   Nrejected\t\t Edgeval\n");

		for (int i=0; i<l->Nsites; i++)
			fprintf(fp, "%d,%d\t\t         %8.6f\t\t   %8.6f\t\t   %8.6f\t\t     %d\n", coords[i][0], coords[i][1], Npicks[i], Nacpt[i], Nrjct[i], s->sim_sites[i].edge_val);

		fclose(fp);
	}

	return;
}


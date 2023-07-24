/* SIMULATON STATISTICS FUNCTIONS */

// 28/11/2019 //

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_statistics_double.h>

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

/* MPI Task Scheduler */

typedef struct task_sch Task;

struct task_sch {

	int rank_id;

	double *Trng;
	double *E; 
	double *M;
	double *X;
	double *C;

	double **rq_readings;

	double *dists;
	double *corr;
};

/* Lattice and Geometry */

typedef struct atom Atm;

struct atom {

	char spec[2];
	char type;
	double pos[3];
	int index;
	int spin_val;
	double xy_spin_val;

	int *nb_arr;
	short int edge_val;

};

typedef struct lattice Latt;

struct lattice {

	char *fpth;

	double a,b,c;
	short int PBC[3];
	int length;
	int Nsites;
	int nb_cnt;
	
	Atm *coord_arr;
	int *edge_states;

	short int boundaries;
	int num_edgeatms;

};

/* Simulation Statistics */

typedef struct site Site;

struct site {

        int index;

        long int Npicks;
        long int Nflips_acpt;
        long int Nflips_rjct;
        short int edge_val;

};

typedef struct system System;

struct system {

        int Nsites;
        double a, b, c;
	long long int Simtime;
	int Tempsteps;
        Site *sim_sites;

        double **q_readings;
};


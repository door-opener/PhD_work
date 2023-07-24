/* HEADER FILE FOR MOVECLASS ROUTINES */
// 27/11/19 //

// Ising Model Routines //
double IsingProp_CMCMOVE(Latt *l, double T, unsigned int seed); 				/* MC IS SENSITIVE TO RNDMSEED! (UPDATED 2/9/19) */
double IsingProp_RMCMOVE(System *s, Latt *l, double T, unsigned int seed); 			/* (11/11/19) */

// XY Model Routines //
double XYProp_MCMOVE(Latt *l, double T, unsigned int seed, int MPI_RANK); 			/* 23/10/19 */

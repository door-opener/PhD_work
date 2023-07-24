/* HEADER FILE FOR MOVECLASS ROUTINES */
// 27/11/19 //

// Ising Model Routines //
double IsingProp_CMCMOVE(System *s, Latt *l, double T, int RND_SITE, double RND_TRIAL); 				/* MC IS SENSITIVE TO RNDMSEED! (UPDATED 2/9/19) */
double IsingProp_RMCMOVE(System *s, Latt *l, double T, int RND_SITE, double RND_TRIAL); 			/* (11/12/19) */

// XY Model Routines //
double XYProp_MCMOVE(Latt *l, double T, int RND_SITE, double RND_TRIAL, int MPI_RANK); 			/* 23/10/19 */

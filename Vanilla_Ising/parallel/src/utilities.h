// Magnetics Utilities //
// 27/11/2019 //

// Memory Allocation & Free Routines //
void IsingLatt_params(char *fpth, Latt *l);
void IsingSpinArr_alloc(int len, Latt *l);
void IsingLatt_coords(Latt *l);
void IsingSpinArr_free(int len, Latt *l, bool flag);
void Statistics_Qalloc(System *s);						/* 28/11/19 */
void Statistics_Qfree(System *s, int MPI_RANK);

// Nearest Neighbours //
int IsingLatt_nbchk(Latt *l, double cutrd, int atmIdx, bool flag); 					/* updated 27/8  IF FLAG == 1, the function will append nbs to the correct array. */
void IsingSpinArr_alloc_nb(Latt *l);

// Geometrical Utilities //
double IsingLatt_dist(Latt *l, int atmIdx1, int atmIdx2);
double IsingLatt_shrtst(Latt *l);

// MPI Routines //
void TskAllocateT_MPI(Task *Tsch, int numprocs, int Tsteps, int rank, double *Trange);
void TskMemAlloc_MPI(Task *Tsch, System *s, int Tsteps, int MPI_RANK, int Q);
void TskCollectRoot_MPI(Task *Tsch, int Tsteps, double *T, int MPI_RANK, MPI_Status *stat, int numprocs, int MODEL, int BETA);
void TskGatherStatRoot_MPI(System *s, Task *Tsch, int Tsteps, int MPI_RANK, MPI_Status *stat, int numprocs);

// MISCELLANEOUS //
void IsingSpinArr_print(Latt *l, bool flag);
double min(int length, const double *arr);
void IsingTempset(double *array, double Tinit, double Tfin, int Tsteps, int MPI_RANK, bool BETA); 	/* if Beta is True, Temperature is read in as Beta not Kelvin */
void reverse(int *arr, int length);
void IsingSpinSum(Latt *l);
void IsingSpinSnap(Latt *l, char *fname);

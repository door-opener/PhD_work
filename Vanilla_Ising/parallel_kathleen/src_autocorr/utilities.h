// Magnetics Utilities //
// 27/11/2019 //

// Memory Allocation & Free Routines //
void IsingLatt_params(char *fpth, Latt *l);
void IsingSpinArr_alloc(int len, Latt *l);
void IsingLatt_coords(Latt *l);
void IsingSpinArr_free(int len, Latt *l, bool flag);

// Nearest Neighbours //
int IsingLatt_nbchk(Latt *l, double cutrd, int atmIdx, bool flag); 					/* updated 27/8  IF FLAG == 1, the function will append nbs to the correct array. */
void IsingSpinArr_alloc_nb(Latt *l);
void IsingShowNBarr(Latt *l);

// Geometrical Utilities //
double IsingLatt_dist(Latt *l, int atmIdx1, int atmIdx2);
double IsingLatt_shrtst(Latt *l);

// MPI Routines //
void TskAllocateT_MPI(Task *Tsch, int numprocs, int Tsteps, int MPI_RANK, double *Trange);
void TskMemAlloc_MPI(Task *Tsch, System *s, int Tsteps, int MPI_RANK, int Q);
void TskCollectRoot_MPI(Task *Tsch, int Tsteps, int MPI_RANK, MPI_Status *stat, int numprocs, int MODEL, int BETA, double *E, double *M, double *X, double *C);
void TskGatherStatRoot_MPI(System *s, Task *Tsch, int Tsteps, int MPI_RANK, MPI_Status *stat, int numprocs);
void TskGatherOBS_MPI(Task *Tsch, int Tsteps, double *E, double *err_E, double *M, double *X, double *err_X, double *C, double *err_C);
void TskShow_MPI(Task *Tsch, int numprocs, int Tsteps);

// MISCELLANEOUS //
void Ising_SiteChk(Latt *l);
void Ising_OutputCnfg(Latt *l, char *out);
void Ising_OutputReadings(System *s, char *OUT, double **mag_readings, double **ener_readings, double *T);
void IsingSpinArr_print(Latt *l, bool flag);
void Ising_OutputOBS(char *OUT, double *E, double *err_E, double *M, double *err_M, double *X, double *err_X, double *C, double *err_C, int BETA, int Tsteps, double *T);
double min(int length, const double *arr);
void Utils_Tset(double min, double max, int N, double* data, int BETA); // if BETA == 0, T = KELVIN.
void reverse(int *arr, int length);
void IsingSpinSum(Latt *l);
void IsingSpinSnap(Latt *l, char *fname);
void swap(double *a, double *b);

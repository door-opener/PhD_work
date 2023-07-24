/* DATA STRUCTURES FOR CALCULATING SIMULATION STATISTICS & ERRORS */
// 27/11/19 //

// Memory Allocation //
void Statistics_Init(System *s, Latt *l, long long int Simtime, int Tsteps);
void Statistics_Qalloc(System *s);
void Statistics_QRealloc(System *s);
void Statistics_Qfree(System *s);

//Sim stats.
double Statistics_Autocorr(System *s, double *q_arr, int OUT, char *out);
double Statistics_AutocorrII(System *s, double *q_arr, int OUT, char *out);
double Statistics_AutocorrIII(System *s, double *q_arr, int OUT, char *out);
double Statistics_AutocorrIV(System *s, double *q_arr, int OUT, char *out);
double Statistics_ErrBSTRAP(System *s, Latt *l, double *q_arr, double T, int FLAG);

//Algorithm stats
void AlgorithmStats(System *s, Latt *l, long long int mcsteps, long long int eqm_steps, long long int ens, int FLAG);

//MISC.
void Stats_ShwReadings(System *s, double *T, int FLAG);

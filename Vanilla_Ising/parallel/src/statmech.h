/* HEADER FILE FOR SPIN LATTICE & STAT MECH FUNCTIONS */

// 27/11/2019 //

// Spin-Lattice Initialisation Routines //
void IsingSpinArr_init(Latt *l);
void XYSpinArr_init(Latt *l); 								/* ADDED 16/10/19 */

// Periodic boundary routines //
void IsingBOUNDS_find_edges(Latt *l);
int IsingBOUNDS_update_edges(Latt *l);

// Ising Model Routines //
double IsingProp_En(Latt *l); 								/* working, updated 28/8 */
double IsingProp_Mag(Latt *l, bool MODEL); 						/* working, checked 30/8 */
void IsingProp_Corr(Latt *l, Task *Tsch, double dR, int dummy,  double am); 		/* ADDED 14/11/19 */

// XY Model Routines //
double XYProp_En(Latt *l);								/* /16/10/19 */

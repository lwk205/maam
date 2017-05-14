//read.cc
void initializeprm();
int readprm(ifstream &file_prm);
void fixprm();
int read_restart(ifstream &file_in);

//write.cc
int showprm(const molecule &chn,const int min_cell4atom,const int max_cell4atom);
int write_volchg(ofstream &file_out);
int write_volchn(ofstream &file_out_profile,ofstream &file_out_radius);
int writemolecule();
	//Output potenziale e concentrazioni ioniche su file_dat
int writedata(double *phi, double *k, double *cl);
int writebackup(double *phi, double *k, double *cl);
int writerestart(double *phi, double *k, double *cl);
	//Output dati per eventuale ripresa calcolo
int outputDS();	//Output in cout di tutte le matrici A e vettori b
int outputresults(double *phi, double *k, double *cl);	
	//Output in cout di potenziale e concentrazioni ioniche
int outputvector(double *v);
int outputmatrix(double *A);

//discretize.cc
int getslice(double z, double zmin, double zmax, double dz);
int grid(int nr, int nt,int nz,double zmin, double zmax, double rmax,//Dimensioni griglia
		double rad_probe,//Raggio di probe utilizzato nella discretizzazione
		double zmem_min, double zmem_max,//Limiti regione di membrana
		double epsH2O, double epsMOL,//Costanti dielettriche relative
		double phimem, double phiterm, double *phibnd,
		double *Dk, double *Dcl,
		double rad_k, double rad_cl,
		double zminpaine, double zmaxpaine,
		molecule &chn,//Molecola da discretizzare
		int disc_chg,//Metodo discretizzazione
		int PNP,// = 0 Solve Poisson  = 1 Solve PNP
		double kext, double kint, double clext, double clint,//Condizioni al contorno
		double *chg_mol,//Carica discretizzata sulla gliglia
		int *inside_end, int *outside_begin,//Limiti molecola
		int *is_sol,//Elementi in soluzione
		double *sum_phi,
		double *Acs_k,double  *bcs_k,double *Acs_cl,double  *bcs_cl,
		double *A_phi,double *bcs_phi,
		int &min_cell4atom,//Minimo Numero di celle usate per discretizzare un atomo
		int &max_cell4atom,//Massimo Numero di celle usate per discretizzare un atomo
		double &distancemed,//Distanza media tra centro reale atomo e centro carica discretizzata
		double &distancemax//Distanza media tra centro reale atomo e centro carica discretizzata
		);
int setchg(int nr, int nt,int nz,double zmin, double zmax, double rmax,//Dimensioni griglia
		int disc_chg,//Metodo discretizzazione
		molecule &chn,//Molecola da discretizzare
		double *chg_mol,//Carica discretizzata sulla gliglia
		int &min_cell4atom,//Minimo Numero di celle usate per discretizzare un atomo
		int &max_cell4atom,//Massimo Numero di celle usate per discretizzare un atomo
		double &distancemed,//Distanza media tra centro reale atomo e centro carica discretizzata
		double &distancemax//Distanza media tra centro reale atomo e centro carica discretizzata
		);
int mkDSP(int nr, int nt,int nz, double dr, double dt, double dz,
		double phiterm,double phimem, double *phibnd,
		double *sum_phi,
		double *chg_mol,
		double *A_phi,
		double *bcs_phi
		);
int mkDSN(int nr, int nt,int nz, double dr, double dt, double dz,
		int *inside_end,int *outside_begin,
		double *D, double Cext, double Cint, double *Acs, double *bcs);
int updateDSN(int nr,int nt, int nz, double dr, double dt, double dz,
		int *inside_end,int *outside_begin,
		double *phi,
		double phimem,double phiterm,
		double *Acs_k, double *Acs_cl, double *bcs_k, double *bcs_cl,
		double *A_k, double *A_cl, double *b_k, double *b_cl
	     );
int updateDSN_exp(int nr, int nt,int nz,double dr, double dt, double dz,
		int *inside_end,int *outside_begin,
		double phiterm, double phimem,
		double *phiexp,
		double *Acs_k, double *Acs_cl, double *bcs_k, double *bcs_cl,
		double *A_k, double *A_cl, double *b_k, double *b_cl
		);
	//Costruisce le strutture dati per la risoluzione di Nernst nel campo esponenziale

//core.cc
int CG_P(double *x, double &tol, int &it, double *A, double *b, double *x0, int max_iter, double tol_aim);
int CG_N(double *x, double &tol, int &it, double *A, double *b, double *x0, int max_iter, double tol_aim);
int BiCGSTAB_N(double *x,double &tol,int &it,double *A,double *b,double *x0,int max_iter,double tol_aim);
int BiCGSTAB_P(double *x,double &tol,int &it,double *A,double *b,double *x0,int max_iter,double tol_aim);
int CGS_N(double *x,double &tol,int &it,double *A,double *b,double *x0,int max_iter,double tol_aim);

//analysis.cc
int fluxes(double* k, double *cl, double *phi);	//Calcolo & Output flussi
int get_flux(double *D, double *c, double cint, double cext, double *phi, int val,	
		double *Jbk,double *Jfw,double *Jdx,double *Jsx,double *Jdw,double *Jup);
	//Subroutine calcolo flussi
	
//auxiliary.cc
double volume_cell(int ir,double dr,double dt,double dz);
int printhelp(char *prm_name);

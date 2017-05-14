#ifndef GLOBAL
#define GLOBAL

#ifndef MAIN
#define PREF extern
#else
#define PREF 
#endif

#include <include/common.h>
#include <lib/etc.h>
#include <lib/molecule.h>

////////////////////////////////////////////////////////////////////////////////
//VARIABILI GLOBALI
//
////////////////////////////////////////////////////////////////////////////////
  
//INPUT FILES
  PREF char filepdb[301];	//File pdb della molecola
  PREF ifstream file_pdb;
  PREF char filestr[301];	//File contenente dati a inizio calcolo
  PREF ifstream file_str;

//OUTPUT FILES
  PREF char filegrd[301];	//Molecola discretizzata sulla griglia
  PREF ofstream file_grd;
  PREF char filedat[301];	//Potenziale e concentrazione
  PREF ofstream file_dat;
  PREF char fileflx[301];	//Dati di flusso
  PREF ofstream file_flx;
  PREF char fileend[301];	//File contenente dati alla fine del calcolo
  PREF ofstream file_end;	

//PARAMETRI GRIGLIA  
  PREF int nr,nt,nz;		//Dimensioni griglia
  PREF double rmax,zmin,zmax;	//Estensione del dominio di integrazione
  PREF double dr,dt,dz;		//Step della griglia
  PREF unsigned long int nrntnz;//Numero elementi di griglia
  PREF int ntnz,nrnt,nrnt_nr;		
  PREF double *cell_volume;	//Volume elemendo di griglia [0:nr-1]
  PREF int disc_chg;		//Metodo di discretizzazione della carica:
  				//	0	Uniforme sul volume dell'atomo
				//	1	Inversamente proporzionale alla distanza
				//	2	Concentrata nell'elemento di griglia contenente l'atomo
  
//PARAMETRI SISTEMA 
  PREF double zmem_min,zmem_max;	//Limiti strato lipidico 
  PREF int izmem_min,izmem_max;

//CONDIZIONI AL CONTORNO
  PREF double phimem,*phibnd;			//Condizioni al contorno del potenziale [0:nz-1]
  PREF double kext, kint, clext, clint;	//Condizioni al contorno concentrazioni

//PARAMETRI FISICI
  PREF double phiterm;	//Potenziale termico [mv]
  PREF double beta,beta2;//1/phiterm - 1/(2*phiterm)
  PREF double epsH2O, epsMOL;	//Costanti dielettriche soluzione e molecola
  PREF double *Dk, *Dcl;	//Coefficienti di diffusione	[0:nz-1]
  PREF double zminpaine,zmaxpaine;	//Limiti entro i cui applicare eventualmene Paine-Scherr
  PREF double rad_k,rad_cl,vol_k,vol_cl;
  PREF double rad_probe;	//Raggio molecola di probe

//STRUTTURE DATI
  PREF double *sum_phi, *A_phi, *bcs_phi, *b_phi;	
  					//sum_phi [0:nrntnz-1] Denominatore dei coefficienti
  					//A_phi   [0:6*nrntnz-1] Matrice dei coefficienti
					//bcs_phi [0:nrntnz-1] Parte costante del termine noto
					//b_phi [0:nrntnz-1] Termine noto
  PREF double *Acs_k, *Acs_cl, *bcs_k, *bcs_cl, *A_k, *A_cl, *b_k, *b_cl;
  				//Acs_k/cl [0:6*nrntnz-1] Parte costante della matrice dei coefficienti
				//bcs_k/cl [0:nrntnz-1] Parte costante del termine noto
  				//A_k/cl [0:6*nrntnz-1] Matrice dei coefficienti
				//b_k/cl [0:nrntnz-1] Termine noto
  PREF int *inside_end, *outside_begin;	
  	//inside_end [0:ntnz-1] numero dell'ultimo elemento interno di soluzione
	//outside_begin [0:ntnz-1] numero del'ultimo elemento esterno di molecola
  PREF int *is_sol;		//is_sol[0:nrntnz-1]
  				//1	Elemento di soluzione
				  
  PREF double *chg_mol;	//chg_mol[0:nrntnz-1] Cariche fisse molecola
  
//VARIABILI
  PREF double *phiodd, *phieven, *phiexp;	//Vettori potenziale [0:nrntnz-1]
  PREF double *kodd, *clodd, *keven, *cleven;	//Vettori concentrazioni [0:nrntnz-1]

//VARIABILI DI CONTROLLO
  PREF int verbose;		//Livello output
  				//0	Livello base
				//	
				//
				//1	+ Risultati a ogni iterazione in file.dat
				//	+ Stop a ogni iterazione
				//
				//2	writedata: Output potenziale e concentrazioni su tutte le slice
				//	analysis: Campo di flusso
				//	discretize: Carica proteina nella griglia
				//
				//3	+ Risultati a ogni iterazione in cout	
				//
				//4	+ Strutture date a ogni iterazione in cout	
				//
  PREF int restart;		//Flag per lettura/scrittura dati di restart
  				//0	Nessuna lettura
				//1 	Lettura
				//2 	Scrittura
				//3	Lettura/Scrittura

//CONTROLLO ITERAZIONE
  PREF int PNP;	//0	Soluzione Poisson
  		//1	Soluzione PNP
  PREF int lin_exp;	//0	Soluzione nel campo lineare
  			//1	Soluzione nel campo esponenziale
  PREF int max_it,max_itP, max_itN;	//Numero massimo iterazioni
  PREF int it_restart;
  PREF double tol,tol_P, tol_N;	//Tolleranze
  PREF double wP, wN, mwP, mwN;	//Fattori si smorzamento delle soluzioni
  PREF int smooth;			//Tipo di aggiornamento iterativo:
  					//	0	aggiornamento diretto
					//	1	xnew = (1-w) * xold + w * xnew 
					//		w da input
					//	2	xnew = (1-w) * xold + w * xnew 
					//		w da variazioni massime alla prima iterazione
					//		+ Aggiornamento incrementale a ogni iterazione
					//	3	xnew = (1-w) * xold + w * xnew 
					//		w da input 
					//		+ Aggiornamento incrementale a ogni iterazione
  PREF double MAX_DELTA_PHI,MAX_DELTA_ION;
  		//Variazioni massime consentite di potenziale e concentrazioni ioniche
  		//nel caso smooth sia 2
						  

#endif

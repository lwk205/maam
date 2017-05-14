////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
  
#define MAIN

#include <include/common.h>
#include <src/nepal/global.h>
#include <src/nepal/functions.h>
#include <lib/etc.h>
#include <lib/molecule.h>

bool got_SIGINT=false;	//Used to control SIGINT (kill -s 2 PID) 
			//in the iterative cycle
			//See: /usr/include/bits/signum.h For signal definitions
void postpone_SIGINT(int sig){
	got_SIGINT=true;
	return;
}

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],fileout_volchg[201];
  	ifstream file_in;//File parametri - Input
  	ofstream file_out;//File log - Output
  	ofstream file_out_volchg;//File log - Output - Charge for VMD
	bool out_volchg_flag;
	char fileout_volchn_profile[201];
	char fileout_volchn_radius[201];
  	ofstream file_out_volchn_profile;//File log - Output - Channel Profile for VMD
  	ofstream file_out_volchn_radius;//File log - Output - Channel Profile for VMD
	bool out_volchn_flag;

	molecule chn;//Molecola

	int min_cell4atom, max_cell4atom;//Controllo qualita' della discretizzazione
	double distancemed,distancemax;

	//Variabili temporanee e di debug
	unsigned long int im,iv;
	int iz,itmp;

	//Controllo convergenza algoritmo
	int n_it, reached,it_tmpP,it_tmpNk,it_tmpNcl;
	double tol_tmpP,tol_tmpNk,tol_tmpNcl;
	double rmsd_phi,rmsd_k,rmsd_cl;
	double max_delta_phi;				//Smoothing delle soluzioni iterative
	double max_delta_ion;			
	double dwP,dwN;


//---------------------------------------------------------
//INIZIALIZZAZIONE
//---------------------------------------------------------
#ifdef HAVE_FEENABLEEXCEPT
//Setting dei flag per eccezioni floating point
  feenableexcept(FE_DIVBYZERO);
  feenableexcept(FE_OVERFLOW);
  feenableexcept(FE_INVALID);
#endif  

//To control SIGIN, used in the iterative routine
  signal(SIGINT,postpone_SIGINT);

//To exclude the path from the program name
  if(strrchr(argv[0],'/'))prm_name=(strrchr(argv[0],'/')+1);
  else prm_name=argv[0];

//Timer
  time(&time_start);

//Input-Output file names
  strcpy(filein,prm_name);
  strcat(filein,".in.prm");
  strcpy(fileout,prm_name);
  strcat(fileout,".out.log");

//Initialize
  out_volchg_flag=out_volchn_flag=false;

//Parameters reading
  if(argc==1)
  {
          printhelp(prm_name);
          return 0;
  }
  for(ind_arg=1;ind_arg<(argc);ind_arg++)
  {
          if(argv[ind_arg][0] == '-')
          {
        	  if(strcmp(argv[ind_arg]+1,"h")==0){//Help
        		  printhelp(prm_name);
        		  return 0;
        	  }
        	  else if(strcmp(argv[ind_arg]+1,"i")==0){//Input File
        		  strcpy(filein,argv[ind_arg+1]);
  		          ind_arg++;
        	  }
        	  else if(strcmp(argv[ind_arg]+1,"o")==0){//Output File
        		  strcpy(fileout,argv[ind_arg+1]);
  		          ind_arg++;
        	  }
        	  else if(strcmp(argv[ind_arg]+1,"out_volchg")==0){//Output File - Charge for VMD
        		  strcpy(fileout_volchg,argv[ind_arg+1]);
			  out_volchg_flag=true;
  		          ind_arg++;
        	  }
        	  else if(strcmp(argv[ind_arg]+1,"out_volchn")==0){//Output File - Channel Profile for VMD
        		  strcpy(fileout_volchn_profile,argv[ind_arg+1]);
        		  strcpy(fileout_volchn_radius,argv[ind_arg+1]);
  			  strcat(fileout_volchn_profile,".profile.dx");
  			  strcat(fileout_volchn_radius,".radius.dx");
			  out_volchn_flag=true;
  		          ind_arg++;
        	  }
        	  else{//Wrong Option
        		  cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
  		  	  for(ind_arg=0;ind_arg<(argc);ind_arg++)cerr<<argv[ind_arg]<<" ";
        		  cerr<<endl;
        		  printhelp(prm_name);
        		  return 1;
        	  }
          }
  }

//Input-Output control
  file_in.open(filein);
  if(!file_in){
  	cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein<<endl;
  	return 1;
  }
  if(strcmp(fileout,"nepal.out.log")!=0){
	file_out.open(fileout);
  	if(!file_out)
  	{
  		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
  		return 1;
  	}
  	cout.rdbuf(file_out.rdbuf());
  }
  if(out_volchg_flag){
	file_out_volchg.open(fileout_volchg);
  	if(!file_out_volchg){
  		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout_volchg<<endl;
  		return 1;
  	}
  }
  if(out_volchn_flag){
	file_out_volchn_profile.open(fileout_volchn_profile);
	file_out_volchn_radius.open(fileout_volchn_radius);
  	if(!file_out_volchn_profile){
  		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout_volchn_profile<<endl;
  		return 1;
  	}
  	if(!file_out_volchn_radius){
  		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout_volchn_radius<<endl;
  		return 1;
  	}
  }

//I saluti iniziali  
  cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;

//Inizializzazione parametri
  readprm(file_in);

//Allocazione memoria
  sum_phi = new double [nrntnz];	//Strutture dati Poisson
  A_phi = new double [6*(nrntnz)];
  bcs_phi = new double [nrntnz];
  b_phi = new double [nrntnz];

  inside_end = new int [ntnz];		//Strutture dati molecola
  outside_begin = new int [ntnz];
  is_sol = new int [nrntnz];
  chg_mol = new double [nrntnz];

  Acs_k = new double [6*(nrntnz)];	//Strutture dati Nernst
  Acs_cl = new double [6*(nrntnz)];
  A_k = new double [6*(nrntnz)];
  A_cl = new double [6*(nrntnz)];
  bcs_k = new double [nrntnz];
  bcs_cl = new double [nrntnz];
  b_k = new double [nrntnz];
  b_cl = new double [nrntnz];

  phiodd = new double [nrntnz];		//Variabili Poisson
  phieven = new double [nrntnz];
  phiexp = new double [nrntnz];
  
  kodd = new double [nrntnz];		//Variabili Nernst
  clodd = new double [nrntnz];
  keven = new double [nrntnz];
  cleven = new double [nrntnz];

//Controllo allocazione memoria
  if ( (sum_phi==NULL)||(A_phi==NULL)||(bcs_phi==NULL)||(b_phi==NULL)
  		  ||(inside_end==NULL)||(outside_begin==NULL)||(is_sol==NULL)||(chg_mol==NULL)
		  ||(Acs_k==NULL)||(Acs_cl==NULL)||(A_k==NULL)||(A_cl==NULL)
		  ||(bcs_k==NULL)||(bcs_cl==NULL)||(b_k==NULL)||(b_cl==NULL)
		  ||(phiodd==NULL)||(phieven==NULL)||(phiexp==NULL)
		  ||(kodd==NULL)||(clodd==NULL)||(keven==NULL)||(cleven==NULL)) {
	  cerr<<"ERROR: it is not possible to allocate the requested memory"<<endl;
	  return 1;
  }

  //Inizializzazione di tutte le variabili allocate
  for(iv=0;iv<nrntnz;iv++){
	  bcs_phi[iv]=b_phi[iv]=bcs_k[iv]=bcs_cl[iv]=b_k[iv]=b_cl[iv]=0.0;
	  sum_phi[iv]=phiodd[iv]=phieven[iv]=phiexp[iv]=0.0;
	  kodd[iv]=clodd[iv]=keven[iv]=cleven[iv]=0.0;
	  chg_mol[iv]=0.0;
	  is_sol[iv]=0;
  }
  for(itmp=0;itmp<ntnz;itmp++){
	  inside_end[itmp]=outside_begin[itmp]=0;
  }
  for(im=0;im<6*nrntnz;im++){
          A_phi[im]=Acs_k[im]=Acs_cl[im]=A_k[im]=A_cl[im]=0.0;
  }

//Lettura molecola
  file_pdb.open(filepdb);
  if (!file_pdb)
  {
    cerr<<"ERROR "<<filepdb<<" DOES NOT EXIST\n";
    return 1;
  }
  chn.readpdb(file_pdb);
  file_pdb.close();

//Discretizzazione
  grid(nr,nt,nz,zmin,zmax,rmax,//Dimensioni griglia
		rad_probe,//Raggio di probe utilizzato nella discretizzazione
		zmem_min,zmem_max,//Limiti regione di membrana
		epsH2O,epsMOL,//Costanti dielettriche relative
		phimem,phiterm,phibnd,
		Dk,Dcl,
		rad_k,rad_cl,
		zminpaine,zmaxpaine,
		chn,//Molecola da discretizzare
		disc_chg,//Metodo discretizzazione
		PNP,// = 0 Solve Poisson = 1 Solve PNP
		kext,kint,clext,clint,//Condizioni al contorno
		chg_mol,//Carica discretizzata sulla gliglia
		inside_end,outside_begin,//Limiti molecola
		is_sol,//Elementi in soluzione
		sum_phi,
		Acs_k,bcs_k,Acs_cl,bcs_cl,
		A_phi,bcs_phi,
		min_cell4atom,//Minimo Numero di celle usate per discretizzare un atomo
		max_cell4atom,//Massimo Numero di celle usate per discretizzare un atomo
		distancemed,//Distanza media tra centro reale atomo e centro carica discretizzata
		distancemax//Distanza media tra centro reale atomo e centro carica discretizzata
		);

//Output parametri
  showprm(chn,min_cell4atom,max_cell4atom);

//Output Discretizzazione Molecola sulla griglia  
  if(strcmp(filegrd,"NEPAL.grd")!=0)
  {
  	file_grd.open(filegrd);
  	if (!file_grd)
  	{
  		cerr<<"ERROR IT IS NOT POSSIBLE TO OPEN "<<filegrd<<"\n";
  		exit(1);
  	}
  	writemolecule();
  	file_grd.close();
  }
  if(out_volchg_flag)write_volchg(file_out_volchg);
  if(out_volchn_flag)write_volchn(file_out_volchn_profile,file_out_volchn_radius);

//Inizializzazioni pre-iterazione (necessarie per azzerrare strutture dati precedentemente sporcate,
//alcune inserite solo per sicurezza)
  for(iv=0;iv<nrntnz;iv++)
  {
	  keven[iv]=0;
	  cleven[iv]=0;
	  kodd[iv]=0;
	  clodd[iv]=0;
	  phiodd[iv]=0;
	  phieven[iv]=0;
  }

//Lettura dati di restart
  if((restart==1)||(restart==3))
  {
  	file_str.open(filestr,ios::binary);
  	if (!file_str){
  	  cerr<<"ERROR "<<filestr<<" DOES NOT EXIST\n";
  	  return 1;
  	}
	read_restart(file_str);
	//Controllo modalita' di smorzamento
	if ((smooth!=1)&&(smooth!=3)){
		cerr<<"WARNING: Restart subroutine works better with smooth set to 1 or 3\n";
	}
	//Conversione del potenziale da mV a phiterm
	for(iv=0;iv<nrntnz;iv++)phieven[iv]=(phieven[iv]/phiterm);
	file_str.close();
  }

//Temporizzazione Inizializzazione
  time(&time_end);
  cout<<"INPUT READING ENDED AFTER "<<difftime(time_end,time_start)<<" SECONDS\n\n";

////////////////////////////////////////////////////////////////////////////////  
//ALGORITMO ITERATIVO  
  
//Inizializzazioni
  n_it=0;
  reached=0;
  dwP=1-wP;
  dwN=1-wN;
  rmsd_phi=rmsd_k=rmsd_cl=INF;
  time(&time_start);
  cout<<"ITERATIVE PROCEDURE IS STARTING ON "<<ctime(&time_start)<<"\n";

  if(PNP==0){//SOLUZIONE POISSON
	  //Per poter tenere conto di una distribuzione di carica ionica
	  //calcolata in precedenza
	for(iv=0;iv<nrntnz;iv++)
		b_phi[iv]=bcs_phi[iv]+((keven[iv]-cleven[iv])*(sum_phi[iv]));
  	BiCGSTAB_P(phieven,tol_tmpP,it_tmpP,A_phi,b_phi,phiodd,max_itP,tol_P);
  	cout<<setw(6)<<"POISSON SOLVED IN "<<it_tmpP<<" steps, TOL = "<<tol_tmpP<<endl;
  }

  else{//SOLUZIONE POISSON-NERNST-PLANCK
  cout<<setw(15)<<"n_it"<<setw(15)<<"Poisson"<<setw(15)<<"Nernst-Cation"<<setw(15)<<"Nernst-Anion"<<"\n";
  //Inizio iterazione
  while ((n_it<max_it)&&(!reached))
  {
	  //Aggiornamento b_phi f(keven,cleven)
	  for(iv=0;iv<nrntnz;iv++)
	 	b_phi[iv]=bcs_phi[iv]+((keven[iv]-cleven[iv])*(sum_phi[iv]));

	  //Soluzione Poisson -->phiodd
	  BiCGSTAB_P(phiodd,tol_tmpP,it_tmpP,A_phi,b_phi,phieven,max_itP,tol_P);

	  //Aggiornamento soluzione iterativa Poisson
	  switch(smooth)
	  {
		  case 3:
			  if(n_it==0)dwP=wP*1e-4;
			  if( ((wP+dwP)<1.0) ) {
			          wP=wP+dwP;
			          mwP=1.0-wP;
			  }
			  //phiodd = mwP * phieven + wP * phiodd
			  for(iv=0;iv<nrntnz;iv++)phiodd[iv]=((wP*phiodd[iv])+(mwP*phieven[iv]));
			  break;
		  case 2:
			  if(n_it==0){
	  		  	max_delta_phi=0.0;
				max_delta_ion=kext;
				if(kint>max_delta_ion)max_delta_ion=kint;
				if(clext>max_delta_ion)max_delta_ion=clext;
				if(clint>max_delta_ion)max_delta_ion=clint;
	  		  	for(iv=0;iv<nrntnz;iv++){
	  		  	        if( abs(phiodd[iv]) > max_delta_phi )max_delta_phi=abs(phiodd[iv]);
	  		  	}
			  	if(max_delta_phi>MAX_DELTA_PHI)
					wP=(MAX_DELTA_PHI/max_delta_phi)*(MAX_DELTA_ION/max_delta_ion);
			  	else wP=1.0;
				dwP=wP*1e-3;
				mwP=1.0-wP;
			  }
			  if( ((wP+dwP)<1.0) ) {
			          wP=wP+dwP;
			          mwP=1.0-wP;
			  }
			  //phiodd = mwP * phieven + wP * phiodd
			  for(iv=0;iv<nrntnz;iv++)phiodd[iv]=((wP*phiodd[iv])+(mwP*phieven[iv]));
			  break;
		  case 1:
			  //phiodd = mwP * phieven + wP * phiodd
			  for(iv=0;iv<nrntnz;iv++)phiodd[iv]=((wP*phiodd[iv])+(mwP*phieven[iv]));
			  break;
		  case 0:
			  break;
	  }

	  //Aggiornamento A_k/b_k A_cl/b_cl in funzione di phiodd
	  if(lin_exp==0){
		updateDSN(nr,nt,nz,dr,dt,dz,
			inside_end,outside_begin,
			phiodd,
			phimem,phiterm,
			Acs_k,Acs_cl,bcs_k,bcs_cl,
			A_k,A_cl,b_k,b_cl
	     	);
	  }
	  else{
		for(iv=0;iv<nrntnz;iv++)phiexp[iv]=exp(beta2 * phiodd[iv]);
		updateDSN_exp(nr,nt,nz,dr,dt,dz,
			inside_end,outside_begin,
			phiterm,phimem,
			phiexp,
			Acs_k,Acs_cl,bcs_k,bcs_cl,
			A_k,A_cl,b_k,b_cl
			);
		for(iv=0;iv<nrntnz;iv++){
			phiexp[iv]=exp(beta * phieven[iv]);
			keven[iv]=keven[iv]*phiexp[iv];
			cleven[iv]=cleven[iv]/phiexp[iv];
		}
	  }

	  //Soluzione Nernst -->kodd/clodd
	  BiCGSTAB_N(kodd,tol_tmpNk,it_tmpNk,A_k,b_k,keven,max_itN,tol_N);
	  BiCGSTAB_N(clodd,tol_tmpNcl,it_tmpNcl,A_cl,b_cl,cleven,max_itN,tol_N);
	  if(lin_exp==1)
	  {
		for(iv=0;iv<nrntnz;iv++)
		{
			keven[iv]=keven[iv]/phiexp[iv];
			cleven[iv]=cleven[iv]*phiexp[iv];
			phiexp[iv]=exp(beta * phiodd[iv]);
			kodd[iv]=kodd[iv]/phiexp[iv];
			clodd[iv]=clodd[iv]*phiexp[iv];
		}
	  }

	  //Aggiornamento soluzione iterativa Nernst
	  switch(smooth)
	  {
		  case 3:
			  if(n_it==0)dwN=wN*1e-2;
			  if( ((wN+dwN)<1.0) ){
			          wN=wN+dwN;
			          mwN=1.0-wN;
			  }
		  case 2:
		  case 1:
			  //kodd = mwN * keven + wN * kodd
			  //clodd = mwN * cleven + wN * clodd
			  for(iv=0;iv<nrntnz;iv++){
				  kodd[iv]=((wN*kodd[iv])+(mwN*keven[iv]));
				  clodd[iv]=((wN*clodd[iv])+(mwN*cleven[iv]));
			  }
			  break;
		  case 0:
			  break;
	  }

	  //Output
	  switch(verbose){
		  case 4:
			  outputDS();
		  case 3:
			  outputresults(phiodd,kodd,clodd);
		  case 1:
			  writebackup(phiodd,kodd,clodd);
			  cout<<"Press any key...\n";
			  cin.get();
		  case 0:
  	  		  cout<<setw(6)<<"ODD :"<<setw(9)<<n_it<<setw(15)<<it_tmpP
		  	      <<setw(15)<<it_tmpNk<<setw(15)<<it_tmpNcl<<setw(15)<<"# iterations"
		  	      <<setw(10)<<" || wP = "<<setw(5)<<wP<<"\n";
  	  		  cout<<setw(15)<<""<<setw(15)<<tol_tmpP<<setw(15)<<tol_tmpNk<<setw(15)
		              <<tol_tmpNcl<<setw(15)<<"Tollerance"<<setw(10)<<" || wN = "<<setw(5)<<wN<<"\n";
			  break;
	  }
  
	  //Aggiornamento b_phi f(kodd,clodd)
	  for(iv=0;iv<nrntnz;iv++)
	 	b_phi[iv]=bcs_phi[iv]+((kodd[iv]-clodd[iv])*(sum_phi[iv]));
	  
	  //Soluzione Poisson -->phieven
	  BiCGSTAB_P(phieven,tol_tmpP,it_tmpP,A_phi,b_phi,phiodd,max_itP,tol_P);

	  //Aggiornamento soluzione iterativa Poisson
	  switch(smooth){
		  case 3:
		  case 2:
		  case 1:
			  //phieven = mwP * phiodd + wP * phieven
			  for(iv=0;iv<nrntnz;iv++)phieven[iv]=((wP*phieven[iv])+(mwP*phiodd[iv]));
			  break;
		  case 0:
			  break;
	  }

	  //Aggiornamento A_k/b_k A_cl/b_cl in funzione di phieven
	  if(lin_exp==0){
		updateDSN(nr,nt,nz,dr,dt,dz,
			inside_end,outside_begin,
			phieven,
			phimem,phiterm,
			Acs_k,Acs_cl,bcs_k,bcs_cl,
			A_k,A_cl,b_k,b_cl
			);
	  }
	  else{
		for(iv=0;iv<nrntnz;iv++)phiexp[iv]=exp(beta2 * phieven[iv]);
		updateDSN_exp(nr,nt,nz,dr,dt,dz,
			inside_end,outside_begin,
			phiterm,phimem,
			phiexp,
			Acs_k,Acs_cl,bcs_k,bcs_cl,
			A_k,A_cl,b_k,b_cl
			);
		for(iv=0;iv<nrntnz;iv++){
			phiexp[iv]=exp(beta * phiodd[iv]);
			kodd[iv]=kodd[iv]*phiexp[iv];
			clodd[iv]=clodd[iv]/phiexp[iv];
		}
	  }

	  //Soluzione Nernst -->keven/cleven
	  BiCGSTAB_N(keven,tol_tmpNk,it_tmpNk,A_k,b_k,kodd,max_itN,tol_N);
	  BiCGSTAB_N(cleven,tol_tmpNcl,it_tmpNcl,A_cl,b_cl,clodd,max_itN,tol_N);
	  if(lin_exp==1){
		for(iv=0;iv<nrntnz;iv++){
			kodd[iv]=kodd[iv]/phiexp[iv];
			clodd[iv]=clodd[iv]*phiexp[iv];
			phiexp[iv]=exp(phieven[iv]*beta);
			keven[iv]=keven[iv]/phiexp[iv];
			cleven[iv]=cleven[iv]*phiexp[iv];
		}
	  }

	  //Aggiornamento soluzione iterativa Nernst
	  switch(smooth)
	  {
		  case 3:
		  case 2:
		  case 1:
			  //keven = mwN * kodd + wN * keven
			  //cleven = mwN * clodd + wN * cleven
			  for(iv=0;iv<nrntnz;iv++){
				  keven[iv]=((wN*keven[iv])+(mwN*kodd[iv]));
				  cleven[iv]=((wN*cleven[iv])+(mwN*clodd[iv]));
			  }
			  break;
		  case 0:
			  break;
	  }
	  
	  //Output
	  switch(verbose)
	  {
		  case 4:
			  outputDS();
		  case 3:
			  outputresults(phieven,keven,cleven);
		  case 1:
			  writebackup(phieven,keven,cleven);
			  cout<<"Press any key...\n";
			  cin.get();
		  case 0:
  	  		  cout<<setw(6)<<"EVEN :"<<setw(9)<<n_it<<setw(15)<<it_tmpP
		  	      <<setw(15)<<it_tmpNk<<setw(15)<<it_tmpNcl<<setw(15)<<"# iterations"<<setw(10)
		  	      <<" || wP = "<<setw(5)<<wP<<"\n";
  	  		  cout<<setw(15)<<""<<setw(15)<<tol_tmpP<<setw(15)<<tol_tmpNk<<setw(15)
		  	      <<tol_tmpNcl<<setw(15)<<"Tollerance"<<setw(10)<<" || wN = "<<setw(5)<<wN<<"\n";
			  break;
	  }

	  //Test convergenza
	  rmsd_phi=rmsd_k=rmsd_cl=0.0;
	  for(iv=0;iv<nrntnz;iv++){
		  rmsd_phi+=((phieven[iv]-phiodd[iv])*(phieven[iv]-phiodd[iv]));
		  rmsd_k+=((keven[iv]-kodd[iv])*(keven[iv]-kodd[iv]));
		  rmsd_cl+=((cleven[iv]-clodd[iv])*(cleven[iv]-clodd[iv]));
	  }
	  rmsd_phi=sqrt(rmsd_phi/nrntnz);
	  rmsd_k=sqrt(rmsd_k/nrntnz);
	  rmsd_cl=sqrt(rmsd_cl/nrntnz);
  	  cout<<setw(15)<<""<<setw(15)<<rmsd_phi<<setw(15)<<rmsd_k<<setw(15)<<rmsd_cl<<setw(15)<<"RMSD"<<"\n";
	  cout<<"--------------------------------------------------------------------------------\n";
	  if((rmsd_phi<tol)&&(rmsd_k<tol)&&(rmsd_cl<tol)){
		  if((wP==1)&&(wN==1)&&(smooth==0)){
			  reached=1;
			  break;
		  }else{
		 cout<<"Convergence reached using the smoothing factors wP = "<<wP<<" and wN = "<<wN<<"\n";
		 cout<<"Now cycles with wP = 1 and wN = 1 will test the solution\n";
	  	 cout<<"--------------------------------------------------------------------------------\n";
			  wP=1.0; mwP=0.0;
			  wN=1.0; mwN=0.0;
			  smooth=0;
		  }
	  }

	  if((got_SIGINT)
	      ||((it_restart!=0)&&((n_it%it_restart)==0))
	      ||((rmsd_phi<tol)&&(rmsd_k<tol)&&(rmsd_cl<tol))){
		  // Received interrupt --> Compute Fluxes and Write Output
		  // or periodic writing of the restart files 
		  // or reached convergence, but maybe the wP=1 wN=1 test is still to do
		writebackup(phieven,keven,cleven);
		got_SIGINT=false;
	  }

	  cout.flush();
	  n_it++;
  }
  //Controllo di modalita' termine iterazione
  if(reached){
          cout<<"CONVERGENCE REACHED IN "<<n_it<<" STEPS \n";
          cout<<"--------------------------------------------------------------------------------\n\n";
  }else{
          cout<<"CONVERGENCE WAS NOT REACHED AFTER "<<n_it<<" STES\n";
          cout<<"--------------------------------------------------------------------------------\n\n";
  }

  }
////////////////////////////////////////////////////////////////////////////////  

//---------------------------------------------------------
//FINE CALCOLO	  
//CHIUSURA PROGRAMMA
//---------------------------------------------------------
//Output finale
  time(&time_end);
  cout<<"\nITERATIVE PROCEDURE ENDED AFTER "<<difftime(time_end,time_start)<<" SECONDS\n\n";

  //Espressione di phi in [mV]
  for(iv=0;iv<nrntnz;iv++)phieven[iv]=(phieven[iv]*phiterm);
  for(iz=0;iz<nz;iz++)phibnd[iz]=(phibnd[iz]*phiterm);
  phimem=phimem*phiterm;
  //Output
  switch(verbose)
  {
	  case 3:
  		  outputresults(phieven,keven,cleven);
	  case 0:
		  if(PNP==1)fluxes(keven,cleven,phieven);
		  if(strcmp(filedat,"NEPAL.dat")!=0)writedata(phieven,keven,cleven);
		  break;
  }

//Output per restart
  if((restart==2)||(restart==3))
  {
	if(strcmp(fileend,"NEPAL.end")!=0)
	{
  		file_end.open(fileend,ios::binary);
  		if (!file_end){
  		  cerr<<"ERROR it is not possible to write "<<fileend<<"\n";
  		  return 1;
  		}
		writerestart(phieven,keven,cleven);
		file_end.close();
	}
  }

//Deallocazione memoria
  delete[] cell_volume;

  delete[] phibnd;
  delete[] Dk;
  delete[] Dcl;

  delete[] sum_phi; 
  delete[] A_phi;
  delete[] bcs_phi; 
  delete[] b_phi;

  delete[] phiodd;
  delete[] phieven;
  delete[] phiexp;
  delete[] kodd;
  delete[] clodd;
  delete[] keven;
  delete[] cleven;

  delete[] Acs_k; 
  delete[] Acs_cl;
  delete[] bcs_k; 
  delete[] bcs_cl;
  delete[] A_k;
  delete[] A_cl;
  delete[] b_k;
  delete[] b_cl;

  delete[] outside_begin;
  delete[] inside_end;
  delete[] is_sol;
  delete[] chg_mol;

  file_in.close();
  file_out.close();
  if(out_volchg_flag) file_out_volchg.close();
  if(out_volchn_flag){
	  file_out_volchn_profile.close();
	  file_out_volchn_radius.close();
  }

  //Temporizzazione codice
  time(&time_end);
  cout<<"\nNEPAL SAYS GOODBYE ON "<<ctime(&time_end)<<"\n";

  return 0;
}

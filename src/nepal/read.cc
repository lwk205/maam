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
  
#include <include/common.h>
#include <lib/etc.h>
#include <lib/molecule.h>
#include <src/nepal/global.h>
#include <src/nepal/functions.h>

#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

////////////////////////////////////////////////////////////////////////////////
//INIZIALIZZAZIONE AI VALORI DI DEFAULT
void initializeprm()
{
  //FILE
  strcpy(filepdb,"NEPAL.pdb");
  strcpy(filestr,"NEPAL.str");
  strcpy(filegrd,"NEPAL.grd");
  strcpy(filedat,"NEPAL.dat");
  strcpy(fileflx,"NEPAL.flx");
  strcpy(fileend,"NEPAL.end");

  //PUNTATORI
  cell_volume=NULL;
  phibnd=NULL;
  sum_phi=A_phi=bcs_phi=b_phi=NULL;
  Acs_k=Acs_cl=bcs_k=bcs_cl=A_k=A_cl=b_k=b_cl=NULL;
  inside_end=outside_begin=is_sol=NULL;
  chg_mol=NULL;
  phiodd=phieven=kodd=keven=clodd=cleven=NULL;

//PARAMETRI DI GRIGLIA
  nr=nt=nz=0;
  rmax=zmin=zmax=0.0;
  disc_chg=0;

//PARAMETRI DI MEMBRANA
  zmem_min=zmem_max=0.0;

//CONDIZIONI AL CONTORNO
  phimem=kext=kint=clext=clint=0.0;

//PARAMETRI FISICI
  phiterm=25;
  epsH2O=80;
  epsMOL=2;
  rad_probe=0.0;
  rad_k=1;
  rad_cl=1;
  Dk=Dcl=NULL;
  zminpaine=zmaxpaine=0.0;

//VARIABILI CONTROLLO OUTPUT
  verbose=0;

//VARIABILI CONTROLLO ITERAZIONE  
  PNP=1;
  restart=0;
  lin_exp=0;
  tol=tol_P=tol_N=1e-3;
  max_it=max_itP=max_itN=100;
  it_restart=0;
  smooth=0;
  wP=wN=1.0;
  MAX_DELTA_ION=100.0;	

  return;
}

////////////////////////////////////////////////////////////////////////////////
//READPRM
//
////////////////////////////////////////////////////////////////////////////////
int readprm(ifstream &file_prm)
{
  char line[301],stmp[301];
  double Dtmp,Dzmin,Dzmax;
  int Dizmin,Dizmax;
  int ind_start,ind_restart,len,iz;
  int flag_zread;

//Inizializzazioni
  initializeprm();
  flag_zread=0;

  while (!file_prm.eof())
  {
    file_prm.getline(line,300);
    if((line[0]!='#')&&(strchr(line,'='))){//La linea contiene il carattere =
            ind_start=ind_restart=len=0;
	    while(line[ind_start]==' ')ind_start++;//Per eliminare eventuali spazi iniziali
            while((line[len+ind_start]!=' ')&&(line[len+ind_start]!='='))len++;
	    while((line[ind_start+len+ind_restart]==' ')
		    ||(line[ind_start+len+ind_restart]=='='))ind_restart++;
	    if(strncmp(line+ind_start,"filepdb",len)==0){
	    	strcpy(filepdb,line+ind_start+len+ind_restart);
	    }
	    else if(strncmp(line+ind_start,"filegrd",len)==0){
	    	strcpy(filegrd,line+ind_start+len+ind_restart);
	    }
	    else if(strncmp(line+ind_start,"filestr",len)==0){
				    strcpy(filestr,line+ind_start+len+ind_restart);
	    }
	    else if(strncmp(line+ind_start,"fileend",len)==0){
				    strcpy(fileend,line+ind_start+len+ind_restart);
	    }
	    else if(strncmp(line+ind_start,"filedat",len)==0){
				    strcpy(filedat,line+ind_start+len+ind_restart);
	    }
	    else if(strncmp(line+ind_start,"fileflx",len)==0){
				    strcpy(fileflx,line+ind_start+len+ind_restart);
	    }
	    else if(strncmp(line+ind_start,"PNP",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    PNP=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"nr",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    nr=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"nt",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    nt=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"nz",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    nz=atoi(stmp);
				    if(Dk!=NULL)delete[] Dk;
				    if(Dcl!=NULL)delete[] Dcl;
				    Dk = new double [nz];
				    Dcl = new double [nz];
				    for(iz=0;iz<nz;iz++)Dk[iz]=Dcl[iz]=0.0;
				    flag_zread++;
	    }
	    else if(strncmp(line+ind_start,"rmax",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    rmax=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"zmin",len)==0){
				    if(flag_zread==0)
				    {
					    cerr<<"ERROR IT IS NOT POSSIBLE TO READ ZMIN BEFORE NZ\n";
					    exit(1);
				    }
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    zmin=atof(stmp);
				    flag_zread++;
				    if(flag_zread==3)
				    {
  					dz=(zmax-zmin)/nz;
				    }
	    }
	    else if(strncmp(line+ind_start,"zmax",len)==0){
				    if(flag_zread==0)
				    {
					    cerr<<"ERROR IT IS NOT POSSIBLE TO READ ZMAX BEFORE NZ\n";
					    exit(1);
				    }
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    zmax=atof(stmp);
				    flag_zread++;
				    if(flag_zread==3)
				    {
  					dz=(zmax-zmin)/nz;
				    }
	    }
	    else if(strncmp(line+ind_start,"zmem_min",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    zmem_min=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"zmem_max",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    zmem_max=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"phimem",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    phimem=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"kext",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    kext=atof(stmp);

	    }
	    else if(strncmp(line+ind_start,"kint",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    kint=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"clext",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    clext=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"clint",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    clint=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"epsH2O",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    epsH2O=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"epsMOL",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    epsMOL=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"Dk",len)==0){
				    if(flag_zread!=3)
				    {
				    	cerr<<"ERROR IT IS NOT POSSIBLE TO READ THE DIFFUSION "
				            <<"COEFFICIENTS BEFORE NZ, ZMIN AND ZMAX\n";
				    	exit(1);
				    }
				    if(
					(sscanf(line+ind_start+len+ind_restart,"%lg %s",&Dtmp,stmp)==2)
					&&
					(strcmp(stmp,"PaineScherr")==0)
				      )
				    {//Definizione dei coefficienti di diffusione con modalita' PaineSherr
			             //Il -Dtmp come valore di diffusione e' utilizzato come messaggio per
				     //le routine di discretizzazione che dovranno occuparsi
				     //della definizione dei coefficienti di diffusione sulla base
				     //del raggio della particolare slice
				    	   for(iz=0;iz<nz;iz++)Dk[iz]=-Dtmp;
				    }
				    else
				    {//Definizione manuale dei coefficienti di diffusione
				    	   if ((sscanf(line+ind_start+len+ind_restart,"%lg %lg %lg",&Dtmp,
						    &Dzmin,&Dzmax))!=3)
				    	   {
					    	cerr<<"ERROR IN DIFFUSION COEFFICIENTS READING\n";
						cerr<<line+ind_start+len+ind_restart<<"\n";
					    	exit(1);
				    	   }
				    	   Dizmin=getslice(Dzmin,zmin,zmax,dz);
				    	   Dizmax=getslice(Dzmax,zmin,zmax,dz);
				    	   for(iz=Dizmin;iz<Dizmax;iz++)Dk[iz]=Dtmp;
				    }
	    }
	    else if(strncmp(line+ind_start,"Dcl",len)==0){
				    if(flag_zread!=3)
				    {
					    cerr<<"ERROR IT IS NOT POSSIBLE TO READ THE DIFFUSION "
						    <<"COEFFICIENTS BEFORE NZ, ZMIN AND ZMAX\n";
					    exit(1);
				    }
				    if(
					(sscanf(line+ind_start+len+ind_restart,"%lg %s",&Dtmp,stmp)==2)
					&&
					(strcmp(stmp,"PaineScherr")==0)
				      )
				    {//Definizione dei coefficienti di diffusione con modalita' PaineSherr
			             //Il -Dtmp come valore di diffusione e' utilizzato come messaggio per
				     //le routine di discretizzazione che dovranno occuparsi
				     //della definizione dei coefficienti di diffusione sulla base
				     //del raggio della particolare slice
				    	   for(iz=0;iz<nz;iz++)Dcl[iz]=-Dtmp;
				    }
				    else
				    {//Definizione manuale dei coefficienti di diffusione
				    	   if ((sscanf(line+ind_start+len+ind_restart,"%lg %lg %lg",&Dtmp,
				    	       	    &Dzmin,&Dzmax))!=3)
				    	   {
				    	           cerr<<"ERROR IN DIFFUSION COEFFICIENTS READING\n";
				    	           exit(1);
				    	   }
				    	   Dizmin=getslice(Dzmin,zmin,zmax,dz);
				    	   Dizmax=getslice(Dzmax,zmin,zmax,dz);
				    	   for(iz=Dizmin;iz<Dizmax;iz++)Dcl[iz]=Dtmp;
				    }
	    }
	    else if(strncmp(line+ind_start,"tol",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    tol=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"max_it",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    max_it=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"tol_P",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    tol_P=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"max_itP",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    max_itP=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"tol_N",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    tol_N=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"max_itN",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    max_itN=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"it_restart",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    it_restart=atoi(stmp);
				    if(it_restart<0){
					    cerr<<"ERROR Negative value for it_restart\n";
					    exit(1);
				    }
	    }
	    else if(strncmp(line+ind_start,"phiterm",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    phiterm=atof(stmp);
				    if(phiterm==0)
				    {
					    cerr<<"ERROR Termic Potential can not be 0\n";
					    exit(1);
				    }
	    }
	    else if(strncmp(line+ind_start,"wP",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    wP=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"wN",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    wN=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"smooth",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    smooth=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"rad_probe",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    rad_probe=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"lin_exp",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    lin_exp=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"verbose",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    verbose=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"restart",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    restart=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"disc_chg",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    disc_chg=atoi(stmp);
	    }
	    else if(strncmp(line+ind_start,"rad_k",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    rad_k=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"rad_cl",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    rad_cl=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"zminpaine",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    zminpaine=atof(stmp);
	    }
	    else if(strncmp(line+ind_start,"zmaxpaine",len)==0){
				    strcpy(stmp,line+ind_start+len+ind_restart);
				    stmp[strlen(line)-(ind_start+len+ind_restart)]='\0';
				    zmaxpaine=atof(stmp);
	    }

    }
  }

  fixprm();

  return 0;
}

void fixprm(){
  int ir,iz,flagk,flagcl;

//PARAMETRI DI GRIGLIA
  dr=rmax/nr;
  dt=(2*PI)/nt;
  dz=(zmax-zmin)/nz;
  nrntnz=nr*nt*nz;
  nrnt=nr*nt;
  nrnt_nr=nrnt-nr;
  ntnz=nt*nz;
  if(cell_volume==NULL)cell_volume = new double [nr];	//Volume elementi di griglia
  for(ir=0;ir<nr;ir++)
  {
	  cell_volume[ir]=volume_cell(ir,dr,dt,dz);
  }

//PARAMETRI DI MEMBRANA
  izmem_min=(int)round((zmem_min-zmin)/dz);
  izmem_min++;	//Numero primo elemento di membrana
  izmem_max=(int)round((zmem_max-zmin)/dz); //Numero ultimo elemento di membrana
  if(izmem_min>izmem_max)
  {
	  cerr<<"ERROR WRONG MEMBRANE THICKNESS\n";
	  exit(1);
  }

//PARAMETRI CONTROLLO CONVERGENZA
  switch(smooth)
  {
	  case 0:
		  if(wP!=1)
		  {
			  wP=1;
			  cerr<<"WARNING: no smoothing and wP != 1 ---> wP was set to 1\n";
		  }
		  if(wN!=1)
		  {
			  wN=1;
			  cerr<<"WARNING: no smoothing and wN != 1 ---> wN was set to 1\n";
		  }
		  break;
	  case 1:
	  case 2:
	  case 3:
		  break;
	  default:
		  cerr<<"ERROR SMOOTH = "<<smooth<<"\n";
		  exit(1);
			  
  }
  MAX_DELTA_PHI=2.0;	//Espresso in potenziali termici
  
    
//CONDIZIONI AL CONTORNO
  phimem=phimem/phiterm;
  if(phibnd==NULL)phibnd = new double [nz];
  for(iz=0;iz<nz;iz++)
  {
	  phibnd[iz]=phimem-(phimem/(nz+1))*(iz+1);
  }

//PARAMETRI FISICI
  beta=1.0;
  beta2=1.0/2.0;
  vol_k=rad_k*rad_k*rad_k*(4.0/3.0)*PI;
  vol_cl=rad_cl*rad_cl*rad_cl*(4.0/3.0)*PI;
  if(Dk[0]<0)flagk=0;
  else flagk=1;
  if(Dcl[0]<0)flagcl=0;
  else flagcl=1;
  for(iz=0;iz<nz;iz++)	//Controllo definizione corretta coefficienti di diffusione
  {
	  if((Dk[iz]==0)||(Dcl[iz]==0))
	  {
		  cerr<<"ERROR IN DIFFUSION COEFFICIENTS DEFINITION\n";
		  cerr<<"Dk[iz="<<iz<<"] = "<<Dk[iz]<<"\t(z = "<<iz*dz+(dz/2.0)<<")\n";
		  cerr<<"Dcl[iz="<<iz<<"] = "<<Dk[iz]<<"\t(z = "<<iz*dz+(dz/2.0)<<")\n";
		  exit(1);
	  }
	  if(((flagk)&&(Dk[iz]<0))||((flagcl)&&(Dcl[iz]<0)))
	  {
		  cerr<<"ERROR IN DIFFUSION COEFFICIENTS DEFINITION\n";
		  cerr<<"Dk[iz="<<iz<<"] = "<<Dk[iz]<<"\t(z = "<<iz*dz+(dz/2.0)<<")\n";
		  cerr<<"Dcl[iz="<<iz<<"] = "<<Dk[iz]<<"\t(z = "<<iz*dz+(dz/2.0)<<")\n";
		  exit(1);
	  }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
//READRESTART
//
//Legge un file contenenti i dati per la ripresa dei calcoli
////////////////////////////////////////////////////////////////////////////////
int read_restart(ifstream &file_in)
{
	int nr_tmp,nt_tmp,nz_tmp;

	//Lettura dimensione griglia dal file binario e confronto con i valori attuali
	file_in.read((char *)&nr_tmp,sizeof(int));
	if((nr_tmp!=nr)&&(nr!=0)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)&nt_tmp,sizeof(int));
	if((nt_tmp!=nt)&&(nt!=0)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)&nz_tmp,sizeof(int));
	if((nz_tmp!=nz)&&(nz!=0)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	nr=nr_tmp;
	nt=nt_tmp;
	nz=nz_tmp;
	nrntnz=nr*nt*nz;

	//Allocate memory if needed
  	if(phieven==NULL) phieven = new double [nrntnz];
  	if(keven==NULL) keven = new double [nrntnz];
  	if(cleven==NULL) cleven = new double [nrntnz];
  	if(Dk==NULL) Dk = new double [nz];
  	if(Dcl==NULL) Dcl = new double [nz];
	
	//Lettura dati di potenzile e concentrazioni ioniche
	file_in.read((char *)phieven,nrntnz*sizeof(double));
	if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)keven,nrntnz*sizeof(double));
	if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)cleven,nrntnz*sizeof(double));
	if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	//Limiti molecola  
	file_in.read((char *)&zmin,sizeof(double));
	file_in.read((char *)&zmax,sizeof(double));
	file_in.read((char *)&rmax,sizeof(double));
	//Raggio ione
	file_in.read((char *)&rad_k,sizeof(double));
	file_in.read((char *)&rad_cl,sizeof(double));
	//Coefficienti di diffusione
	file_in.read((char *)Dk,nz*sizeof(double));
	if((unsigned int)file_in.gcount()!=nz*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)Dcl,nz*sizeof(double));
	if((unsigned int)file_in.gcount()!=nz*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}

	return 0;
}

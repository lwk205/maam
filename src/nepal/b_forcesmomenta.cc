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
#include <include/cpplapack.h>

#include <lib/etc.h>
#include <lib/molecule.h>
#include <src/nepal/functions.h>

//Dichiarazione funzioni

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],filepdb[201];
	ifstream file_in,file_pdb;
	ofstream file_out;

	double x_hinge,y_hinge,z_hinge;
	bool not_x_hinge,not_y_hinge,not_z_hinge;

	int nr,nt,nz;
	unsigned long int nrntnz,nrnt;
	double *phi,*k,*cl;
	double zmin,zmax,rmax;
	double dr,dt,dz;

//---------------------------------------------------------
//START-UP	BEGIN
//Floatint point exception
#ifdef HAVE_FEENABLEEXCEPT
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_OVERFLOW);
	feenableexcept(FE_INVALID);
#endif

//To exclude the path from the program name
	if(strrchr(argv[0],'/'))prm_name=(strrchr(argv[0],'/')+1);
	else prm_name=argv[0];

//Timer
	time(&time_start);

//Input-Output file names
	strcpy(filein,prm_name);
	strcat(filein,".in.end");
	strcpy(filepdb,prm_name);
	strcat(filepdb,".in.pdbrq");
	strcpy(fileout,prm_name);
	strcat(fileout,".out.dat");

//Initialization
	x_hinge=y_hinge=z_hinge=0.0;
	not_x_hinge=not_y_hinge=not_z_hinge=true;

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
	      	  if(strcmp(argv[ind_arg]+1,"h")==0)//Help
	      	  {
	      		  printhelp(prm_name);
	      		  return 0;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"i")==0)//Input File
	      	  {
	      		  strcpy(filein,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"str")==0)//Structure File - Where to compute the force
	      	  {
	      		  strcpy(filepdb,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"x_hinge")==0)//X - Hinge Point
	      	  {
			  not_x_hinge=false;
			  x_hinge=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"y_hinge")==0)//Y - Hinge Point
	      	  {
			  not_y_hinge=false;
			  y_hinge=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"z_hinge")==0)//Z - Hinge Point
	      	  {
			  not_z_hinge=false;
			  z_hinge=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"o")==0)//Output File
	      	  {
	      		  strcpy(fileout,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else//Wrong Option
	      	  {
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
  	file_pdb.open(filepdb);
	if(!file_pdb){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filepdb<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}
	if(not_x_hinge){
		cout<<"x coordinate of the hinge point = ";
		cin>>x_hinge;
		cout<<endl;
	}
	if(not_y_hinge){
		cout<<"y coordinate of the hinge point = ";
		cin>>y_hinge;
		cout<<endl;
	}
	if(not_z_hinge){
		cout<<"z coordinate of the hinge point = ";
		cin>>z_hinge;
		cout<<endl;
	}

	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//START-UP	END
//---------------------------------------------------------

//---------------------------------------------------------
//INPUT		BEGIN
  
//Dimensioni griglia  
  file_in.read((char *)&nr,sizeof(int));
  file_in.read((char *)&nt,sizeof(int));
  file_in.read((char *)&nz,sizeof(int));
  nrntnz=nr*nt*nz;
  nrnt=nr*nt;

//Allocazione memoria
  phi = new double [nrntnz];
  k = new double [nrntnz];
  cl = new double [nrntnz];

//Potenziale  
  file_in.read((char *)phi,nrntnz*sizeof(double));
  if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double)){
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }

//COncentrazioni
  file_in.read((char *)k,nrntnz*sizeof(double));
  if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double)){
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }
  file_in.read((char *)cl,nrntnz*sizeof(double));
  if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double)){
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }

//Limiti molecola  
  file_in.read((char *)&zmin,sizeof(double));
  file_in.read((char *)&zmax,sizeof(double));
  file_in.read((char *)&rmax,sizeof(double));

//Calcolo parametri derivati  
  dr=rmax/nr;
  dt=(2*PI)/nt;
  dz=(zmax-zmin)/nz;
//INPUT		END
//---------------------------------------------------------

  
//---------------------------------------------------------
//ANALYSIS	BEGIN
  
  molecule inpdb;

  double *chg_pdb;//Carica su cui calcolare la forza
  int ir,it,iz,ind_atm;
  unsigned long int iv,iv_dx,iv_sx;
  double r,t,z,x_cell,y_cell,z_cell,x_atm,y_atm,z_atm,rad2,dist2;
  int min_cell4atom,max_cell4atom;//Qualita' procedure discretizzazione carica
  double distancemed,distancemax;//Qualita' procedure discretizzazione carica

  //Allocazione memoria
  chg_pdb = new double [nrntnz];

  //-----------------------------------------
  //Individuazione regione in cui calcolare forza
  inpdb.readpdb(file_pdb);
  setchg(nr,nt,nz,zmin,zmax,rmax,1,inpdb,chg_pdb,min_cell4atom,max_cell4atom,distancemed,distancemax);
  for(iv=0;iv<nrntnz;iv++){
	  	chg_pdb[iv]=chg_pdb[iv]*(8.8544/1.6022); //Conversione delle cariche a [a10-6]
		//if(chg_pdb[iv]!=0)cout<<chg_pdb[iv]<<endl;
  }
  //-----------------------------------------
    
  //-----------------------------------------
  double mass;
  double cnt_mol[3];
  double I[9];
  //calcolo massa
  mass=inpdb.mass();
  //calcolo baricentro
  inpdb.centermass(cnt_mol);
  //calcolo momenti di inerzia
  inpdb.inertia(I,x_hinge,y_hinge,z_hinge);
  //Calcolo assi principali di inerzia
  int n=3;	//Dimensione matrice n*n
  char jobz='V'; //Calcola autovettori
  char uplo='U'; //Leggi dalla diagonale alta della matrice
  double lambda[n];	//Vettore di autovalori
  int lwork=1+6*n+2*n*n; //Vettori a uso interno
  int liwork=3+5*n;
  double work[lwork];
  int iwork[liwork];
  int info; // = 0 Se calcolo riuscito
  //file_out<<"#Inertia matrix Before Diagonalization = "<<endl;
  //file_out<<setw(15)<<I[0]<<"\t"<<setw(15)<<I[1]<<"\t"<<setw(15)<<I[2]<<endl;
  //file_out<<setw(15)<<I[3]<<"\t"<<setw(15)<<I[4]<<"\t"<<setw(15)<<I[5]<<endl;
  //file_out<<setw(15)<<I[6]<<"\t"<<setw(15)<<I[7]<<"\t"<<setw(15)<<I[8]<<endl;
  dsyevd_(&jobz,&uplo,&n,I,&n,lambda,work,&lwork,iwork,&liwork,&info);
  if(info!=0){
	  cerr<<"ERROR it is not possible to linearize the inertia matrix"<<endl;
	  exit(1);
  }
  //-----------------------------------------
  //-----------------------------------------
    
  //-----------------------------------------
  double dphidr_dx,dphidr_sx,dphidt_dx,dphidt_sx,dphidz_dx,dphidz_sx;
  double fr,ft,fz,fx,fy;
  double F[3],m[3];
  double lx,ly,lz;

  //Calcolo Forza - Momento
  m[0]=m[1]=m[2]=0.0;
  F[0]=F[1]=F[2]=0.0;
  //Nella prima slice non si puo' calcolare la derivata per cui la escludiamo
  //lo stesso faccio per l'ultima fetta
  for(iz=1;iz<nz;iz++){//cicli per prendere tutti gli elementi di griglia su cui calcolare la forza
	z_cell=zmin+(dz/2.0)+iz*dz;
  	for(it=0;it<nt;it++){
		t=(dt/2.0)+it*dt;
  		for(ir=1;ir<nr;ir++){	//Lo stesso motivo di prima si saltano il primo e' l'ultimo anello
			iv=iz*nrnt+it*nr+ir;
			r=(dr/2.0)+ir*dr;
			x_cell=r*cos(t);
			y_cell=r*sin(t);
			if(chg_pdb[iv]!=0){//Elementi di griglia occupati dalla regione in cui calcolare la forza
					   //E' la molecola contenuta nel file -str
				lx=x_cell-x_hinge;
				ly=y_cell-y_hinge;
				lz=z_cell-z_hinge;
				//Derivate direzione radiale
				dphidr_dx=(phi[iv+1]-phi[iv])/dr;//Derivata Dx
				dphidr_sx=(phi[iv]-phi[iv-1])/dr;//Derivata Sx
				//Derivate direzione angolare
				iv_dx=iv-nr;
				iv_sx=iv+nr;
				if(it==0)iv_dx=iv+nr*(nt-1);
				if(it==(nt-1))iv_sx=iv-nr*(nt-1);
				dphidt_dx=(phi[iv_sx]-phi[iv])/(dt*dr*(0.5+ir));//Derivata Dx
				dphidt_sx=(phi[iv]-phi[iv_dx])/(dt*dr*(0.5+ir));//Derivata Dx
				//Derizate direzione assiale
				dphidz_dx=(phi[iv+nrnt]-phi[iv])/dz;//Derivata Dx
				dphidz_sx=(phi[iv]-phi[iv-nrnt])/dz;//Derivata Sx
				//Calcolo Forze
				fr=-chg_pdb[iv]*((dphidr_dx+dphidr_sx)/2.0);
				ft=-chg_pdb[iv]*((dphidt_dx+dphidt_sx)/2.0);
				fz=-chg_pdb[iv]*((dphidz_dx+dphidz_sx)/2.0);
				fx=fr*cos(t)-ft*sin(t);
				fy=fr*sin(t)+ft*cos(t);
				//Calcolo momenti
				m[0]+=(ly*fz-lz*fy);
				m[1]+=(lz*fx-lx*fz);
				m[2]+=(lx*fy-ly*fx);
				//Forza Risultante
				F[0]+=fx;
				F[1]+=fy;
				F[2]+=fz;
				//file_out<<dphidt_dx<<'\t'<<dphidt_sx<<'\t'<<endl;
				//file_out<<iz<<'\t'<<it<<'\t'<<ir<<'\t'
				//	<<phi[iv]<<'\t'<<dphidr_sx<<'\t'<<dphidr_dx<<'\t'<<endl;
			}
		}
		//file_out<<endl;
	}
  }
  F[0]=F[0]*1.609e-3;//Conversione a fN (moltiplicata per la carica elementare in Coulomb)
  F[1]=F[1]*1.609e-3;
  F[2]=F[2]*1.609e-3;
  m[0]=m[0]*1.609e-6;//Conversione a pN*A (moltiplicata per la carica elementare in Coulomb)
  m[1]=m[1]*1.609e-6;
  m[2]=m[2]*1.609e-6;
  //Passaggio dei momenti nel sistema delle direzioni principali
  //l'algotimi di diagonalizzazione della matrice I restituisce gli autovettori come
  //righe all'interno di I stessa
  //Per il cambiamento di base basta quindi moltiplicare per questa matrice
  double m_eig[n];
  char trans='N';
  double done=1.0,dzero=0.0;
  int ione=1;
  dgemv_(&trans,&n,&n,&done,I,&n,m,&ione,&dzero,m_eig,&ione);
  double normF;
  normF=dnrm2_(&n,F,&ione);
  //-----------------------------------------

    
  //-----------------------------------------
  //		Output
  file_out<<"# Principal Axis"<<endl;
  file_out<<setw(15)<<x_hinge<<setw(15)<<y_hinge<<setw(15)<<z_hinge<<endl;
  file_out<<setw(15)<<x_hinge+I[0]<<setw(15)<<y_hinge+I[1]<<setw(15)<<z_hinge+I[2]<<endl<<endl<<endl;
  file_out<<setw(15)<<x_hinge<<setw(15)<<y_hinge<<setw(15)<<z_hinge<<endl;
  file_out<<setw(15)<<x_hinge+I[3]<<setw(15)<<y_hinge+I[4]<<setw(15)<<z_hinge+I[5]<<endl<<endl<<endl;
  file_out<<setw(15)<<x_hinge<<setw(15)<<y_hinge<<setw(15)<<z_hinge<<endl;
  file_out<<setw(15)<<x_hinge+I[6]<<setw(15)<<y_hinge+I[7]<<setw(15)<<z_hinge+I[8]<<endl<<endl<<endl;
  file_out<<"# Force direction"<<endl;
  file_out<<setw(15)<<x_hinge<<setw(15)<<y_hinge<<setw(15)<<z_hinge<<endl;
  file_out<<setw(15)<<x_hinge+(F[0]/normF)<<setw(15)
	  <<y_hinge+(F[1]/normF)<<setw(15)
	  <<z_hinge+(F[2]/normF)<<endl<<endl<<endl;
  file_out<<"# Principal Axis - Centered in zero"<<endl;
  file_out<<setw(15)<<0.0<<setw(15)<<0.0<<setw(15)<<0.0<<endl;
  file_out<<setw(15)<<0.0+I[0]<<setw(15)<<0.0+I[1]<<setw(15)<<0.0+I[2]<<endl<<endl<<endl;
  file_out<<setw(15)<<0.0<<setw(15)<<0.0<<setw(15)<<0.0<<endl;
  file_out<<setw(15)<<0.0+I[3]<<setw(15)<<0.0+I[4]<<setw(15)<<0.0+I[5]<<endl<<endl<<endl;
  file_out<<setw(15)<<0.0<<setw(15)<<0.0<<setw(15)<<0.0<<endl;
  file_out<<setw(15)<<0.0+I[6]<<setw(15)<<0.0+I[7]<<setw(15)<<0.0+I[8]<<endl<<endl<<endl;

  file_out<<"# GRID nr = "<<nr<<" nt = "<<nt<<" nz = "<<nz<<endl;
  file_out<<"# GRID rmax = "<<rmax<<" zmin = "<<zmin<<" zmax = "<<zmax<<endl<<endl;

  file_out<<"# Mass = "<<mass<<" [au]"<<endl;
  file_out<<"# Center of mass : "<<cnt_mol[0]<<" "<<cnt_mol[1]<<" "<<cnt_mol[2]<<" "<<endl;
  file_out<<"# Eigenvectors = "<<endl;
  file_out<<"#"<<setw(14)<<I[0]<<"\t"<<setw(15)<<I[3]<<"\t"<<setw(15)<<I[6]<<endl;
  file_out<<"#"<<setw(14)<<I[1]<<"\t"<<setw(15)<<I[4]<<"\t"<<setw(15)<<I[7]<<endl;
  file_out<<"#"<<setw(14)<<I[2]<<"\t"<<setw(15)<<I[5]<<"\t"<<setw(15)<<I[8]<<endl;
  file_out<<"# Inertia Matrix = [au*A]"<<endl;
  file_out<<"#"<<setw(14)<<lambda[0]<<"\t"<<setw(15)<<0.0<<"\t"<<setw(15)<<0.0<<endl;
  file_out<<"#"<<setw(14)<<0.0<<"\t"<<setw(15)<<lambda[1]<<"\t"<<setw(15)<<0.0<<endl;
  file_out<<"#"<<setw(14)<<0.0<<"\t"<<setw(15)<<0.0<<"\t"<<setw(15)<<lambda[2]<<endl<<endl;

  file_out<<"# Force on the center of mass [fN]  : "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
  file_out<<"# Force Norm = "<<normF<<" fN"<<endl;

/*
  //Momenti
  file_out<<"#Momenta = "<<endl;
  file_out<<"mx = "<<m[0]<<" pN*A"<<endl;
  file_out<<"my = "<<m[1]<<" pN*A"<<endl;
  file_out<<"mz = "<<m[2]<<" pN*A"<<endl;
  file_out<<"#Momenta along principal axes = "<<endl;
  file_out<<"m_eig_x = "<<m_eig[0]<<" pN*A"<<endl;
  file_out<<"m_eig_y = "<<m_eig[1]<<" pN*A"<<endl;
  file_out<<"m_eig_z = "<<m_eig[2]<<" pN*A"<<endl;
  file_out<<endl<<endl;

  file_out<<"#Grid point involved in force calculation"<<endl;
  file_out<<"#Usage: splot '"<<prm_name<<".out.dat' i 1 u 1:2:3 w p"<<endl;
  iv=0;
  for(iz=0;iz<nz;iz++){
	z_cell=zmin+(dz/2.0)+iz*dz;
  	for(it=0;it<nt;it++){
		t=(dt/2.0)+it*dt;
  		for(ir=0;ir<nr;ir++){
			r=(dr/2.0)+ir*dr;
			x_cell=r*cos(t);
			y_cell=r*sin(t);
			if(chg_pdb[iv]!=0){
				file_out<<x_cell<<'\t'<<y_cell<<'\t'<<z_cell<<endl;
			}
			iv++;
		}
	}
  }
  file_out<<endl<<endl;
*/
  //-----------------------------------------
    
//ANALYSIS	END
//---------------------------------------------------------

  time(&time_end);
  cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;
  
  delete[] phi;
  delete[] k;
  delete[] cl;
  
  delete[] chg_pdb;
  
  file_in.close();
  file_pdb.close();
  file_out.close();
  
  return 0;
}

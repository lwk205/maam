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
#include <lib/etc.h>
#include <lib/molecule.h>
#include <src/nepal/functions.h>

//Dichiarazione funzioni

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],filestr2[201],filepdb[201],filepdbfix[201];
	ifstream file_in,file_str2,file_pdb,file_pdbfix;
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
	strcpy(filepdbfix,prm_name);
	strcat(filepdbfix,".in.fix.pdbrq");
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
	      	  else if(strcmp(argv[ind_arg]+1,"str2")==0)//Structure File - Where we're going
	      	  {
	      		  strcpy(filestr2,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"fixchg")==0)//Structure File - Fixed charges
	      	  {
	      		  strcpy(filepdbfix,argv[ind_arg+1]);
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
	if(!file_in)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein<<endl;
		return 1;
	}
  	file_pdb.open(filepdb);
	if(!file_pdb)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filepdb<<endl;
		return 1;
	}
  	file_str2.open(filestr2);
	if(!file_str2)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filepdb<<endl;
		return 1;
	}
  	file_pdbfix.open(filepdbfix);
	if(!file_pdbfix)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filepdbfix<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out)
	{
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
  if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double))
  {
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }

//COncentrazioni
  file_in.read((char *)k,nrntnz*sizeof(double));
  if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double))
  {
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }
  file_in.read((char *)cl,nrntnz*sizeof(double));
  if((unsigned long int)file_in.gcount()!=nrntnz*sizeof(double))
  {
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

//Individua gli elementi di griglia appartenenti alla molecola nel filepdb

  molecule str2,inpdb,infix;

  double *chg_pdb;//Carica su cui calcolare la forza
  double *chg_fix;//Carica fissa da considerare nel calcolo della forza
  int ir,it,iz,ind_atm;
  unsigned long int iv,iv_dx,iv_sx;
  double r,t,z,x_cell,y_cell,z_cell,x_atm,y_atm,z_atm,rad2,dist2;
  int min_cell4atom,max_cell4atom;//Qualita' procedure discretizzazione carica
  double distancemed,distancemax;//Qualita' procedure discretizzazione carica

  //Allocazione memoria
  chg_pdb = new double [nrntnz];
  chg_fix = new double [nrntnz];

  //-----------------------------------------
  //Lettura Molecola con cariche fisse
  infix.readpdb(file_pdbfix);
  setchg(nr,nt,nz,zmin,zmax,rmax,1,infix,chg_fix,min_cell4atom,max_cell4atom,distancemed,distancemax);
  //-----------------------------------------

  //-----------------------------------------
  //Lettura Molecola direzione
  str2.readpdb(file_str2);
  //-----------------------------------------

  //-----------------------------------------
  //Individuazione regione in cui calcolare forza
  inpdb.readpdb(file_pdb);
  setchg(nr,nt,nz,zmin,zmax,rmax,1,inpdb,chg_pdb,min_cell4atom,max_cell4atom,distancemed,distancemax);
  //-----------------------------------------
    
  //-----------------------------------------
  double dphidr_dx,dphidr_sx,dphidt_dx,dphidt_sx,dphidz_dx,dphidz_sx;
  double fr,ft,fz,fx,fy;
  double mx,my,mz;
  double lx,ly,lz;
  double Fx,Fy,Fz;
  //Calcolo Forza

  mx=my=mz=0.0;
  Fx=Fy=Fz=0.0;
 	//Nella prima slice non si puo' calcolare la derivata per cui la escludiamo
  	//lo stesso faccio per l'ultima fetta
  for(iz=1;iz<nz;iz++){//cicli per prendere tutti gli elementi di griglia su cui calcolare la forza
  //for(iz=50;iz<60;iz++){//cicli per prendere tutti gli elementi di griglia su cui calcolare la forza
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
				fr=-chg_pdb[iv]*((dphidr_dx+dphidr_sx)/2.0);
				ft=-chg_pdb[iv]*((dphidt_dx+dphidt_sx)/2.0);
				fz=-chg_pdb[iv]*((dphidz_dx+dphidz_sx)/2.0);
				fx=fr*cos(t)-ft*sin(t);
				fy=fr*sin(t)+ft*cos(t);
				mx+=(ly*fz-lz*fy);
				my+=(lz*fx-lx*fz);
				mz+=(lx*fy-ly*fx);
				Fx+=fx;
				Fy+=fy;
				Fz+=fz;
				//file_out<<dphidt_dx<<'\t'<<dphidt_sx<<'\t'<<endl;
				//file_out<<iz<<'\t'<<it<<'\t'<<ir<<'\t'
				//	<<phi[iv]<<'\t'<<dphidr_sx<<'\t'<<dphidr_dx<<'\t'<<endl;
			}
		}
		//file_out<<endl;
	}
  }
  mx=mx*1.609;//Conversione a pN*A (moltiplicata per la carica elementare in Coulomb)
  my=my*1.609;
  mz=mz*1.609;
  //-----------------------------------------

  //-----------------------------------------
  //Output
  file_out<<"#GRID nr = "<<nr<<" nt = "<<nt<<" nz = "<<nz<<endl;
  file_out<<"#GRID rmax = "<<rmax<<" zmin = "<<zmin<<" zmax = "<<zmax<<endl<<endl;

  int pot_mx,pot_my,pot_mz;
  int pot_F,pot;
  double F,dirF;

  
  //Forze
  F=sqrt((Fx*Fx)+(Fy*Fy)+(Fz*Fz));
  dirF=(((str2.xcenter-inpdb.xcenter)*Fx)+
        ((str2.ycenter-inpdb.ycenter)*Fy)+
        ((str2.zcenter-inpdb.zcenter)*Fz));
  pot_F=0;
  while((F/10.0)>10.0){F=F/10.0;pot_F++;}
  for(pot=0;pot<pot_F;pot++){
	  Fx=Fx/10.0;
	  Fy=Fy/10.0;
	  Fz=Fz/10.0;
  }
  file_out<<"#F = "<<F<<"*10^"<<pot_F<<" pN*A"<<endl;
  file_out<<"#directionF = "<<dirF<<endl;
  file_out<<inpdb.xcenter<<'\t'<<inpdb.ycenter<<'\t'<<inpdb.zcenter<<endl;
  file_out<<inpdb.xcenter+Fx<<'\t'<<inpdb.ycenter+Fy<<'\t'<<inpdb.zcenter+Fz<<endl;
  file_out<<endl<<endl;

  //Momenti
  pot_mx=pot_my=pot_mz=0;
  //while((abs(mx)/10.0)>10.0){mx=mx/10.0;pot_mx++;}
  //while((abs(my)/10.0)>10.0){my=my/10.0;pot_my++;}
  //while((abs(mz)/10.0)>10.0){mz=mz/10.0;pot_mz++;}
  mx=mx/10000000;
  my=my/10000000;
  mz=mz/10000000;
  file_out<<"#mx = "<<mx<<"*10^"<<pot_mx<<" pN*A"<<endl;
  file_out<<x_hinge<<'\t'<<y_hinge<<'\t'<<z_hinge<<endl;
  file_out<<x_hinge+mx<<'\t'<<y_hinge<<'\t'<<z_hinge<<endl;
  file_out<<endl;
  file_out<<"#my = "<<my<<"10^"<<pot_my<<" pN*A"<<endl;
  file_out<<x_hinge<<'\t'<<y_hinge<<'\t'<<z_hinge<<endl;
  file_out<<x_hinge<<'\t'<<y_hinge+my<<'\t'<<z_hinge<<endl;
  file_out<<endl;
  file_out<<"#mz = "<<mz<<"10^"<<pot_mz<<" pN*A"<<endl;
  file_out<<x_hinge<<'\t'<<y_hinge<<'\t'<<z_hinge<<endl;
  file_out<<x_hinge<<'\t'<<y_hinge<<'\t'<<z_hinge+mz<<endl;
  file_out<<endl;
  file_out<<x_hinge<<'\t'<<y_hinge<<'\t'<<z_hinge<<endl;
  file_out<<x_hinge+mx<<'\t'<<y_hinge+my<<'\t'<<z_hinge+mz<<endl;
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
  file_out<<"#Grid point containing fixed charge"<<endl;
  file_out<<"#Usage: splot '"<<prm_name<<".out.dat' i 2 u 1:2:3 w p"<<endl;
  iv=0;
  for(iz=0;iz<nz;iz++){
	z_cell=zmin+(dz/2.0)+iz*dz;
  	for(it=0;it<nt;it++){
		t=(dt/2.0)+it*dt;
  		for(ir=0;ir<nr;ir++){
			r=(dr/2.0)+ir*dr;
			x_cell=r*cos(t);
			y_cell=r*sin(t);
			if(chg_fix[iv]!=0){
				file_out<<x_cell<<'\t'<<y_cell<<'\t'<<z_cell<<endl;
			}
			iv++;
		}
	}
  }
  //-----------------------------------------
    
//ANALYSIS	END
//---------------------------------------------------------

  time(&time_end);
  cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;
  
  delete[] phi;
  delete[] k;
  delete[] cl;
  
  delete[] chg_pdb;
  delete[] chg_fix;
  
  file_in.close();
  file_pdb.close();
  file_str2.close();
  file_pdbfix.close();
  file_out.close();
  
  return 0;
}

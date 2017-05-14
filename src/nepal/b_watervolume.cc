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

	char filein[201],fileout[201];
	ifstream file_in;
	ofstream file_out;

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
	strcpy(fileout,prm_name);
	strcat(fileout,".out.dat");

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
  	file_out.open(fileout);
	if(!file_out)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
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
  int ir,it,iz;
  unsigned long int iv;
  double volume_water;

  volume_water=0.0;
  iv=0;
  for(iz=0;iz<nz;iz++){
  	for(it=0;it<nt;it++){
  		for(ir=0;ir<nr;ir++){
			if((k[iv]!=0.0)||(cl[iv]!=0.0)){
				volume_water+=volume_cell(ir,dr,dt,dz);
			}
			iv++;
		}
	}
  }

//Output
  file_out<<"Volume_water "<<volume_water<<endl;
    
//ANALYSIS	END
//---------------------------------------------------------

  time(&time_end);
  cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;
  
  delete[] phi;
  delete[] k;
  delete[] cl;
  
  file_in.close();
  file_out.close();
  
  return 0;
}

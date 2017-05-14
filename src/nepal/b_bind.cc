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
#include <src/nepal/global.h>
#include <src/nepal/functions.h>

#define MAX_BINDING_SITES 10

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],fileinbs[201];
	ifstream file_in;
	ifstream file_inbs;
	ofstream file_out;

	bool matlab;
	char matname[201];
	char cmt;

	char line[201];
	int n_read;
	double ztrsl;
	char bs_name[MAX_BINDING_SITES][201];
	double z_bot[MAX_BINDING_SITES],z_top[MAX_BINDING_SITES];
	double bs_k[MAX_BINDING_SITES],bs_cl[MAX_BINDING_SITES];
	double bs_k_rad[MAX_BINDING_SITES],bs_cl_rad[MAX_BINDING_SITES];
	int iz_bot[MAX_BINDING_SITES],iz_top[MAX_BINDING_SITES];
	int i,bs_ind,bs_num,ind_rad_k;

	int ir,it,iz;
	unsigned long int iv;
	bool inside_bnd;

	double *phi,*k,*cl,*rad_int;


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
	strcpy(fileinbs,prm_name);
	strcat(fileinbs,".in.bs");
	strcpy(fileout,prm_name);
	strcat(fileout,".out.dat");

//Data Initialization
	ztrsl=0.0;
	matlab=false;
	strcpy(matname,"D");

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
	      	  else if(strcmp(argv[ind_arg]+1,"bs")==0)//Input binding site boundaries
	      	  {
	      		  strcpy(fileinbs,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"o")==0)//Output File
	      	  {
	      		  strcpy(fileout,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"ztrsl")==0)//Z axes Translation
	      	  {
			  ztrsl=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"matlab")==0)//Output for Matlab
	      	  {
			  matlab=true;
	      		  strcpy(matname,argv[ind_arg+1]);
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
  	file_inbs.open(fileinbs);
	if(!file_inbs)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileinbs<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}
	if(matlab)cmt='%';
	else cmt='#';

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
  Dk = new double [nz];
  Dcl = new double [nz];
  rad_int = new double [nz];

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

//Raggio ione
  file_in.read((char *)&rad_k,sizeof(double));
  file_in.read((char *)&rad_cl,sizeof(double));
  rad_k=1.33;
  rad_cl=1.81;

//Coefficienti di diffusione
  file_in.read((char *)Dk,nz*sizeof(double));
  if((unsigned int)file_in.gcount()!=nz*sizeof(double))
  {
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }
  file_in.read((char *)Dcl,nz*sizeof(double));
  if((unsigned int)file_in.gcount()!=nz*sizeof(double))
  {
  	cerr<<"ERROR IN RESTART READING: "<<filein<<" ENDS TOO SOON\n";
  	exit(1);
  }

//Calcolo parametri derivati  
  dr=rmax/nr;
  dt=(2*PI)/nt;
  dz=(zmax-zmin)/nz;
  vol_k=(4.0/3.0)*PI*rad_k*rad_k*rad_k;
  vol_cl=(4.0/3.0)*PI*rad_k*rad_k*rad_k;

	//Binding site boundaries File
	bs_ind=0;
	while(!file_inbs.eof())
	{
		file_inbs.getline(line,200);
		if(line[0]!='#')
		{
			n_read=sscanf(line,"%s %lf %lf\n"
					,bs_name[bs_ind],z_bot+bs_ind,z_top+bs_ind);
			if(n_read==3)
			{
				if(z_bot[bs_ind]>=z_top[bs_ind])
				{
			cerr<<"ERROR: BINDING SITE UPPER BOUNDARY < LOWER BOUNDARY"<<endl;
			exit(1);
				}
			bs_ind++;
			}
		}
	}
	bs_num=bs_ind;
//INPUT		END
//---------------------------------------------------------

//---------------------------------------------------------
//ANALYSIS	BEGIN
        file_out<<cmt<<setw(14)<<"NAME"
		<<setw(15)<<"z_bot"<<setw(15)<<"z_bot on grid"
		<<setw(15)<<"z_top"<<setw(15)<<"z_top on grid"<<endl;
	for(bs_ind=0;bs_ind<bs_num;bs_ind++)
	{
		bs_k[bs_ind]=bs_cl[bs_ind]=bs_k_rad[bs_ind]=bs_cl_rad[bs_ind]=0.0;
		iz_bot[bs_ind]=(int)((z_bot[bs_ind]-zmin-(dz/2.0))/dz);
		iz_top[bs_ind]=(int)((z_top[bs_ind]-zmin-(dz/2.0))/dz);
		file_out<<cmt<<setw(14)<<bs_name[bs_ind]
		    <<setw(15)<<z_bot[bs_ind]<<setw(15)<<zmin+(dz/2.0)+iz_bot[bs_ind]*dz
		    <<setw(15)<<z_top[bs_ind]<<setw(15)<<zmin+(dz/2.0)+iz_top[bs_ind]*dz
		    <<endl;
	}
	ind_rad_k=(int)floor(rad_k/dr);
	file_out<<cmt<<setw(14)<<"Cation Radius = "<<rad_k<<endl;
	file_out<<cmt<<setw(14)<<"* = Number of ions inside a radius equal to Cation Radius"<<endl;
	if(matlab)file_out<<cmt<<setw(14)<<"";
	else file_out<<cmt<<setw(14)<<"NAME";
	file_out<<setw(15)<<"<Z>"
	    <<setw(15)<<"CATION"<<setw(15)<<"CATION*"
	    <<setw(15)<<"ANION"<<setw(15)<<"ANION*"<<endl;
	if(matlab)file_out<<matname<<" = ["<<endl;
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
	       	{

	ir=0;
	inside_bnd=false;
	while((ir<nr)&&(!inside_bnd))
	{
		iv=(iz*nrnt)+(it*nr)+ir;
		//In this way cell outside the channel wall are not
		//taken into account
		if((k[iv]!=0)||(cl[iv]!=0))
		{
			for(bs_ind=0;bs_ind<bs_num;bs_ind++)
			{
			  if((iz>iz_bot[bs_ind])&&(iz<=iz_top[bs_ind]))
			  {
			  	bs_k[bs_ind]+=(k[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz));
			  	bs_cl[bs_ind]+=(cl[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz));
				if(ir<ind_rad_k)
				{
			  	  bs_k_rad[bs_ind]+=(k[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz));
			  	  bs_cl_rad[bs_ind]+=(cl[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz));
				}
			  }
			}
		}
		else inside_bnd=true;
		ir++;
	} 

		}
	}
	for(bs_ind=0;bs_ind<bs_num;bs_ind++)
	{
		if(matlab)file_out<<setw(15)<<"";
		else file_out<<setw(15)<<bs_name[bs_ind];
		file_out<<setw(15)<<((z_bot[bs_ind]+z_top[bs_ind])/2.0)+ztrsl
		    <<setw(15)<<bs_k[bs_ind]<<setw(15)<<bs_k_rad[bs_ind]
		    <<setw(15)<<bs_cl[bs_ind]<<setw(15)<<bs_cl_rad[bs_ind]<<endl;
	}

	if(matlab)file_out<<"];"<<endl;
//ANALYSIS	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_inbs.close();
	file_out.close();

  	return 0;
}

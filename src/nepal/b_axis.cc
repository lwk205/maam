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

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201];
	ifstream file_in;
	ofstream file_out;

	bool matlab;
	char matname[201];
	double ztop,zbot;
	char cmt;
	int ind_rad_k,iz,it,num_sol,ir,inside_bnd;
	double ztrsl,z,phi_tmp,k_tmp,cl_tmp,k_occ_tmp,cl_occ_tmp;
	double vol_k,vol_cl,vol_slice;
	unsigned long int iv;
	int N_int;
	double *rad_int,rad_intmed,rad_intmin,rad_botmin;
	double *phi,*k,*cl;

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

//Data Initialization
	ztrsl=0.0;
	matlab=false;
	strcpy(matname,"D");
	ztop=+INF;
	zbot=-INF;

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
	      	  else if(strcmp(argv[ind_arg]+1,"ztrsl")==0)//Z axes Translation
	      	  {
			  ztrsl=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"ztop")==0)//Z Top Analysis Radius
	      	  {
			  ztop=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"zbot")==0)//Z Bottom Analysis Radius
	      	  {
			  zbot=atof(argv[ind_arg+1]);
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
//INPUT		END
//---------------------------------------------------------

  
//---------------------------------------------------------
//ANALYSIS	BEGIN
	file_out<<cmt<<"Data along the channel axis"<<endl;
	file_out<<cmt<<"GRID nr = "<<nr<<" nt = "<<nt<<" nz = "<<nz<<endl;
	file_out<<cmt<<"GRID rmax = "<<rmax<<" zmin = "<<zmin<<" zmax = "<<zmax<<endl;
	file_out<<cmt<<"Cation radius = "<<rad_k<<"[A] (volume = "<<vol_k<<"[A^3])"<<endl;
	file_out<<cmt<<"Anion radius = "<<rad_cl<<"[A] (volume = "<<vol_cl<<"[A^3])"<<endl;
	file_out<<cmt<<"Colunm 1	z Axes"<<endl;
	file_out<<cmt<<"Colunm 2	Electrostatic Potential ir=0 it=0"<<endl;
	file_out<<cmt<<"Colunm 3	Cation concentration ir=0 it=0"<<endl;
	file_out<<cmt<<"Colunm 4	Cation occupancy ir=0 it=0"<<endl;
	file_out<<cmt<<"Colunm 5	Anion concentration ir=0 it=0"<<endl;
	file_out<<cmt<<"Colunm 6	Anion occupancy ir=0 it=0"<<endl;
	file_out<<cmt<<"Colunm 7	Electrostatic Potential ir=0 it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 8	Cation concentration ir=0 it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 9	Cation occupancy ir=0 it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 10	Anion concentration ir=0 it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 11	Anion occupancy ir=0 it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 12	Electrostatic Potential ir=0:inside_bnd it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 13	Cation concentration ir=0:inside_bnd it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 14	Cation occupancy ir=0:inside_bnd it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 15	Anion concentration ir=0:inside_bnd it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 16	Anion occupancy ir=0:inside_bnd it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 17	Electrostatic Potential ir=0:rad_k it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 18	Cation concentration ir=0:rad_k it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 19	Cation occupancy ir=0:rad_k it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 20	Anion concentration ir=0:rad_k it=0:nt-1"<<endl;
	file_out<<cmt<<"Colunm 21	Anion occupancy ir=0:rad_k it=0:nt-1"<<endl;
	file_out<<cmt<<""<<endl;
	if(matlab)file_out<<matname<<" = ["<<endl;
	ind_rad_k=(int)floor(rad_k/dr);
	if(((rad_k/dr)-ind_rad_k)>0.5)ind_rad_k++;
	for(iz=0;iz<nz;iz++)
	{
		//Value along the z-axis ir=it=0
		z=zmin+(dz/2.0)+iz*dz;
		z+=ztrsl;
		iv=(iz*nrnt);
		file_out<<setw(15)<<z<<setw(15)<<phi[iv]
	      		<<setw(15)<<k[iv]<<setw(15)<<k[iv]*mmol2ionA3*vol_k
	      		<<setw(15)<<cl[iv]<<setw(15)<<cl[iv]*mmol2ionA3*vol_cl;

		//Value along the z-axiz mediated along it (ir=0)
		phi_tmp=k_tmp=cl_tmp=k_occ_tmp=cl_occ_tmp=vol_slice=0.0;
		num_sol=0;
		for(it=0;it<nt;it++)
		{
			iv=(iz*nrnt)+(it*nr);
			phi_tmp+=phi[iv];
			k_tmp+=k[iv];
			k_occ_tmp+=k[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz);
			cl_tmp+=cl[iv];
			cl_occ_tmp+=cl[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz);
			if((k[iv]!=0)||(cl[iv]!=0))
			{
				num_sol++;
				vol_slice+=volume_cell(ir,dr,dt,dz);
			}
		}
		if(num_sol!=nt)
		{
			cerr<<"ERROR: "<<num_sol<<" GRID ELEMENTS AROUND THE AXIS\n";
			exit(1);
		}
		phi_tmp=phi_tmp/nt;
		k_tmp=k_tmp/num_sol;
		k_occ_tmp=k_occ_tmp*vol_k/vol_slice;
		cl_tmp=cl_tmp/num_sol;
		cl_occ_tmp=cl_occ_tmp*vol_cl/vol_slice;
		file_out<<setw(15)<<phi_tmp
			<<setw(15)<<k_tmp<<setw(15)<<k_occ_tmp
			<<setw(15)<<cl_tmp<<setw(15)<<cl_occ_tmp;


		//Value along the z-axiz mediated along it and ir (ir=0:inside_bnd)
		phi_tmp=k_tmp=cl_tmp=k_occ_tmp=cl_occ_tmp=vol_slice=0.0;
		num_sol=0;
		for(it=0;it<nt;it++)
		{
			ir=0;
			inside_bnd=0;
			while((ir<nr)&&(!inside_bnd))
			{
		  		iv=(iz*nrnt)+(it*nr)+ir;
				phi_tmp+=(phi[iv]*volume_cell(ir,dr,dt,dz));
				k_tmp+=k[iv];
				k_occ_tmp+=(k[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz));
				cl_tmp+=cl[iv];
				cl_occ_tmp+=(cl[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz));
				if((k[iv]!=0)||(cl[iv]!=0))
				{
					num_sol++;
					vol_slice+=volume_cell(ir,dr,dt,dz);
				}
				else
				{
					inside_bnd=1;
					if((ir*dr)<rad_int[iz])rad_int[iz]=ir*dr;
				}
				ir++;
			}
		}
		phi_tmp=phi_tmp;
		k_tmp=k_tmp/num_sol;
		k_occ_tmp=k_occ_tmp*vol_k/vol_slice;
		cl_tmp=cl_tmp/num_sol;
		cl_occ_tmp=cl_occ_tmp*vol_cl/vol_slice;
		file_out<<setw(15)<<phi_tmp
			<<setw(15)<<k_tmp<<setw(15)<<k_occ_tmp
			<<setw(15)<<cl_tmp<<setw(15)<<cl_occ_tmp;
		
		
		//Value along the z-axiz mediated along it and ir (ir=0:rad_k)
		phi_tmp=k_tmp=cl_tmp=k_occ_tmp=cl_occ_tmp=vol_slice=0.0;
		num_sol=0;
		for(it=0;it<nt;it++)
		{
			ir=0;
			inside_bnd=0;
			while((ir<ind_rad_k)&&(!inside_bnd))
			{
		  		iv=(iz*nrnt)+(it*nr)+ir;
				phi_tmp+=(phi[iv]*volume_cell(ir,dr,dt,dz));
				k_tmp+=k[iv];
				k_occ_tmp+=k[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz);
				cl_tmp+=cl[iv];
				cl_occ_tmp+=cl[iv]*mmol2ionA3*volume_cell(ir,dr,dt,dz);
				if((k[iv]!=0)||(cl[iv]!=0))
				{
					num_sol++;
					vol_slice+=volume_cell(ir,dr,dt,dz);
				}
				else
				{
					inside_bnd=1;
				}
				ir++;
			}
		}
		phi_tmp=phi_tmp;
		k_tmp=k_tmp/num_sol;
		k_occ_tmp=k_occ_tmp*vol_k/vol_slice;
		cl_tmp=cl_tmp/num_sol;
		cl_occ_tmp=cl_occ_tmp*vol_cl/vol_slice;
		file_out<<setw(15)<<phi_tmp
			<<setw(15)<<k_tmp<<setw(15)<<k_occ_tmp
			<<setw(15)<<cl_tmp<<setw(15)<<cl_occ_tmp<<endl;
	}
	if(matlab)file_out<<"];"<<endl;

	N_int=0;
	rad_intmin=INF;
	rad_botmin=INF;
	rad_intmed=0.0;
	for(iz=0;iz<nz;iz++){
		z=(iz*dz)+zmin+(dz/2.0);
		//cout<<"z = "<<z<<" rad = "<<rad_int[iz]<<endl;
		//cout<<z<<" "<<zbot<<" "<<ztop;
		if((z>=zbot)&&(z<=ztop)){
			//cout<<" OK!";
			if(rad_int[iz]<rad_intmin){
				rad_intmin=rad_int[iz];
				//cout<<" MIN = "<<rad_intmin;
			}
			N_int++;
			rad_intmed+=rad_int[iz];
		}
		if((z>=(zbot-dz))&&(z<=(zbot+dz))){
			if(rad_int[iz]<rad_botmin){
				rad_botmin=rad_int[iz];
				//cout<<" MIN = "<<rad_intmin;
			}
		}
		//cout<<endl;
	}
	rad_intmed=rad_intmed/N_int;
	file_out<<cmt<<"Internal Mean Radius in ["<<zbot<<":"<<ztop<<"] = "<<rad_intmed<<" [A]"<<endl;
	file_out<<cmt<<"Internal Min Radius in ["<<zbot<<":"<<ztop<<"] = "<<rad_intmin<<" [A]"<<endl;
	file_out<<cmt<<"Bottom Min Radius = "<<rad_botmin<<" [A]"<<endl;

//ANALYSIS	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();

  	return 0;
}

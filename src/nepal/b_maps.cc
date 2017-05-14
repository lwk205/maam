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

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],stmp[201];
	ifstream file_in;
	ofstream file_out;

	int bnd_kind;

	bool matlab,flag_phi;
	char matname[201];
	char cmt;
	int ind_rad_k,iz,it,num_sol,ir,inside_bnd,num_ele,sgn;
	double ztrsl,z,phi_tmp,k_tmp,cl_tmp,k_occ_tmp,cl_occ_tmp,r;
	double vol_k,vol_cl,vol_slice;
	unsigned long int iv;
	int N_int;
	double *rad_int,rad_intmed,rad_intmin,rad_botmin;
	double *phi,*k,*cl;
	double phimed,kmed,clmed;

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
	bnd_kind=1;

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
		  //cerr<<argv[ind_arg]<<endl;
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
	      	  else if(strcmp(argv[ind_arg]+1,"bnd")==0){//Boundary in the map
	      		  strcpy(stmp,argv[ind_arg+1]);
			  bnd_kind=atoi(stmp);
			  if((bnd_kind<0)||(bnd_kind>3)){
	      		  	cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
				exit(1);
			  }
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"phi")==0){//Write only potential in one slice
			  flag_phi = true;
	      	  }
	      	  else if(strncmp(argv[ind_arg]+1,"matlab",6)==0)//Output for Matlab
	      	  {
			  matlab=true;
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
	if(matlab){
		if((bnd_kind!=1)||(flag_phi)){
			cerr<<"=ERROR: MATLAB OUTPUT IMPLEMENTED JUST FOR BND = 1 AND flag_phi FALSE"<<endl;
			exit(1);
		}
		cmt='%';
	}
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
	file_out<<setiosflags(ios::fixed);
	file_out<<setprecision(5);
	if(flag_phi){
		file_out<<cmt<<setw(14)<<"z [A]"<<setw(15)<<"r [A]"<<setw(15)<<"phi(s) [mV]"<<endl;
	} else {
		file_out<<cmt<<setw(14)<<"z [A]"<<setw(15)<<"r [A]"
		    <<setw(15)<<"phi [mV]"<<setw(15)<<"K [mM]"<<setw(15)<<"Cl [mM]"
		    <<endl;
	}
	if(matlab){//Output for matlab
		  double *R,*PHI,*K,*CL;
		  R= new double [nz*nr];
		  PHI= new double [nz*nr];
		  K= new double [nz*nr];
		  CL= new double [nz*nr];
		  file_out<<setiosflags(ios::fixed);
		  file_out<<setprecision(5);
		  file_out<<"BND = ["<<endl;
		  for(iz=0;iz<nz;iz++){
		    	z=zmin+(dz/2.0)+iz*dz;
		        ir=0;
		  	inside_bnd=0;
		  	while((ir<nr)&&(!inside_bnd)){
		  	      	inside_bnd=0;
		  	  	for(it=0;it<nt;it++){
		  	  		iv=(iz*nrnt)+(it*nr)+ir;
		  	      		if((k[iv]==0)&&(cl[iv]==0))inside_bnd=1;
		  	      	}
				ir++;
			}
			r=ir*dr;
			file_out<<setw(15)<<r<<setw(15)<<z<<endl;
		  }
		  file_out<<"];"<<endl;
		  for(iz=0;iz<nz;iz++){
		    	z=zmin+(dz/2.0)+iz*dz;
		        ir=0;
		  	inside_bnd=0;
		        while(ir<nr){
				//Siamo oltre il bordo ?
		        	for(it=0;it<nt;it++){
		    			iv=(iz*nrnt)+(it*nr)+ir;
		          		if((k[iv]==0)&&(cl[iv]==0))inside_bnd=1;
		        	}
		        	if(inside_bnd){	//SI: continua a dare in output
						//l'ultimo valore di r, phi, k
					R[(iz*nr)+ir]=r;
					PHI[(iz*nr)+ir]=phimed;
					K[(iz*nr)+ir]=kmed;
					CL[(iz*nr)+ir]=clmed;
				}
				else{	//NO: dai in output il corretto valore
		          		r=(ir*dr);
		          		phimed=kmed=clmed=0.0;
		          		for(it=0;it<nt;it++){
		    	  	      		iv=(iz*nrnt)+(it*nr)+ir;
		          	      		phimed+=phi[iv];
		          	      		kmed+=k[iv];
		          	      		clmed+=cl[iv];
		          		}
		          		phimed=(phimed/nt);
		          		kmed=(kmed/nt);
		          		clmed=(clmed/nt);
					R[(iz*nr)+ir]=r;
					PHI[(iz*nr)+ir]=phimed;
					K[(iz*nr)+ir]=kmed;
					CL[(iz*nr)+ir]=clmed;
		  		}		
		        	ir++;
		  	}
		  }
		  file_out<<"R = ["<<endl;
		  for(iz=0;iz<nz;iz++){
		        for(ir=(nr-1);ir>=0;ir--)file_out<<setw(15)<<-R[(iz*nr)+ir];
		        for(ir=0;ir<nr;ir++)file_out<<setw(15)<<R[(iz*nr)+ir];
			file_out<<endl;
		  }
		  file_out<<"];"<<endl;
		  file_out<<"Z = ["<<endl;
		  for(iz=0;iz<nz;iz++){
		    	z=zmin+(dz/2.0)+iz*dz;
		        for(ir=0;ir<nr;ir++)file_out<<setw(15)<<z<<setw(15)<<z;
			file_out<<endl;
		  }
		  file_out<<"];"<<endl;
		  file_out<<"PHI = ["<<endl;
		  for(iz=0;iz<nz;iz++){
		        for(ir=(nr-1);ir>=0;ir--)file_out<<setw(15)<<PHI[(iz*nr)+ir];
		        for(ir=0;ir<nr;ir++)file_out<<setw(15)<<PHI[(iz*nr)+ir];
			file_out<<endl;
		  }
		  file_out<<"];"<<endl;
		  file_out<<"K = ["<<endl;
		  for(iz=0;iz<nz;iz++){
		        for(ir=(nr-1);ir>=0;ir--)file_out<<setw(15)<<K[(iz*nr)+ir];
		        for(ir=0;ir<nr;ir++)file_out<<setw(15)<<K[(iz*nr)+ir];
			file_out<<endl;
		  }
		  file_out<<"];"<<endl;
		  file_out<<"CL = ["<<endl;
		  for(iz=0;iz<nz;iz++){
		        for(ir=(nr-1);ir>=0;ir--)file_out<<setw(15)<<CL[(iz*nr)+ir];
		        for(ir=0;ir<nr;ir++)file_out<<setw(15)<<CL[(iz*nr)+ir];
			file_out<<endl;
		  }
		  file_out<<"];"<<endl;
		  delete [] R;
		  delete [] PHI;
		  delete [] K;
		  delete [] CL;
	}
	else{//Output for gnuplot
		if(flag_phi){//Potential on a slice
			for(iz=0;iz<nz;iz++){
				z=zmin+(dz/2.0)+iz*dz;
				for(ir=0;ir<nr;ir++){
					r=(ir*dr);
				      	file_out<<setw(15)<<z<<setw(15)<<r;
					for(it=0;it<nt;it++){
				  	      	iv=(iz*nrnt)+(it*nr)+ir;
				      		file_out<<setw(15)<<phi[iv];
					}
					file_out<<endl;
				}
				file_out<<endl;
			}
		} else {
			for(sgn=-1;sgn<=1;sgn+=2){
				for(iz=0;iz<nz;iz++){
				  	z=zmin+(dz/2.0)+iz*dz;
				        ir=0;
					inside_bnd=0;
				        while((ir<nr)&&(!inside_bnd)){
						if(bnd_kind==1){
			      				//Considera solo fino al primo elemento fuori
			      				//soluzione
				        		inside_bnd=0;
				      			for(it=0;it<nt;it++){
				  				iv=(iz*nrnt)+(it*nr)+ir;
				      		  		if((k[iv]==0)&&(cl[iv]==0))inside_bnd=1;
				      			}
			      				num_ele=nt;
						} else if(bnd_kind==2){
			      				//Considera fino a quando meta' degli elementi
			      				//sono in soluzione
				        		inside_bnd=1;
				      			num_ele=0;
				      			for(it=0;it<nt;it++){
				  				iv=(iz*nrnt)+(it*nr)+ir;
				      		  		if((k[iv]!=0)||(cl[iv]!=0))num_ele++;
				      				if(num_ele>(nt/2))inside_bnd=0;
				      			}
						} else if(bnd_kind==3){
			      				//Tutti gli elementi
			      				//inside_bnd=0;
			      				num_ele=nt;
						}
				      		if(!inside_bnd){
				      			r=(ir*dr);
				      			phimed=kmed=clmed=0.0;
				      			for(it=0;it<nt;it++){
				  		  	      	iv=(iz*nrnt)+(it*nr)+ir;
				      		      		phimed+=phi[iv];
				      		      		kmed+=k[iv];
				      		      		clmed+=cl[iv];
				      			}
				      			phimed=(phimed/num_ele);
				      			kmed=(kmed/num_ele);
				      			clmed=(clmed/num_ele);
				      			file_out<<setw(15)<<z<<setw(15)<<sgn*r
				      		    		<<setw(15)<<phimed
				      		    		<<setw(15)<<kmed
				      		    		<<setw(15)<<clmed
				      		    		<<endl;
						}		
				      		ir++;
			      		}
					file_out<<endl;
				}
				file_out<<endl;
			}
		}
	}
//ANALYSIS	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();

  	return 0;
}

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

	char filein[201],fileout[201],fileoutbar[201],stmp[201];
	ifstream file_in;
	ofstream file_out;
	ofstream file_outbar;

	int data_type;  //Data written to the output file: 
			// 0 Electrostatic Potential
			// 1 Potassium concentration
			// 2 Chloride concentration
			// 3 Protein-Charge distribution
			// 4 Channel profile
	
	double satp,satm;
	bool satp_flag,satm_flag;
	double min,max;
	bool outbar_flag;
	bool meanteta_flag;
	double meanteta;

	unsigned long int iv;
	int dummy_cnt;
	double xmax,ymax,dx,dy,x,y,dtmp;
	double r_grd,t_grd,r_cln,t_cln;
	int nx,ny,ix,iy,iz,it,ir,ir_cln,it_cln;

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
	strcpy(fileoutbar,prm_name);
	strcat(fileoutbar,".out.dat");

//Data Initialization
	data_type=0;
	satp=satm=0.0;
	satp_flag=satm_flag=outbar_flag=meanteta_flag=false;

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
	      	  else if(strcmp(argv[ind_arg]+1,"outbar")==0){ //Output colorbar
	      		  strcpy(fileoutbar,argv[ind_arg+1]);
			  outbar_flag=true;
			  ind_arg++;
	      	  }
	      	  else if(strncmp(argv[ind_arg]+1,"type",6)==0){ //Data written in the output file
			  data_type=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strncmp(argv[ind_arg]+1,"satp",6)==0){ //Positive Saturation
			  satp=atof(argv[ind_arg+1]);
			  satp_flag=true;
			  ind_arg++;
	      	  }
	      	  else if(strncmp(argv[ind_arg]+1,"satn",6)==0){ //Negative Saturation
			  satm=atof(argv[ind_arg+1]);
			  satm_flag=true;
			  ind_arg++;
	      	  }
	      	  else if(strncmp(argv[ind_arg]+1,"meanteta",6)==0){ //Average on angular coordiantes
			  meanteta_flag=true;
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
	if(!file_in) {
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}
	if(outbar_flag){
  		file_outbar.open(fileoutbar);
		if(!file_outbar){
			cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileoutbar<<endl;
			return 1;
		}
	}
	if((data_type<0)||(data_type>4)){
		cerr<<"ERROR: Data type "<<data_type<<" does not exist"<<endl;
		return 1;
	}

	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;

//START-UP	END
//---------------------------------------------------------

//---------------------------------------------------------
//INPUT		BEGIN
  
//Lettura Dati di griglia
  read_restart(file_in);

//Calcolo parametri derivati  
  nrntnz=nr*nt*nz;
  nrnt=nr*nt;
  dr=rmax/nr;
  dt=(2*PI)/nt;
  dz=(zmax-zmin)/nz;
  vol_k=(4.0/3.0)*PI*rad_k*rad_k*rad_k;
  vol_cl=(4.0/3.0)*PI*rad_k*rad_k*rad_k;
  
  xmax=ymax=rmax*((sqrt(2.0))/2.0);
  nx=ny=(2*((int)(nr*((sqrt(2.0))/2.0))) + 1);
  dx=dy=(2.0*xmax)/(nx-1);
//INPUT		END
//---------------------------------------------------------

  
//---------------------------------------------------------
//ANALYSIS	BEGIN

  	//Average on angular coordinate
	if(meanteta_flag){
		for(iz=0;iz<nz;iz++){
			for(ir=0;ir<nr;ir++){
				meanteta=0;
				for(it=0;it<nt;it++){
					iv=ir+it*nr+iz*nrnt;
					switch (data_type){
						case 0:
							meanteta+=phieven[iv];
							break;
						case 1:
							meanteta+=keven[iv];
							break;
						case 2:
							meanteta+=cleven[iv];
							break;
					}
				}
				meanteta=meanteta/nt;
				for(it=0;it<nt;it++){
					iv=ir+it*nr+iz*nrnt;
					switch (data_type){
						case 0:
							phieven[iv]=meanteta;
							break;
						case 1:
							keven[iv]=meanteta;
							break;
						case 2:
							cleven[iv]=meanteta;
							break;
					}
				}
			}
		}
	}

	//Produce Output
	dummy_cnt=1;
	max=-INF;
	min=INF;
	file_out<<setiosflags(ios::fixed);
	file_out<<setprecision(3);
	file_out<<"object 1 class gridpositions counts "<<nz<<" "<<nx<<" "<<ny<<endl;
	file_out<<"origin "<<-ymax<<" "<<-xmax<<" "<<zmin<<endl;
	file_out<<"delta 0.000000e+00 0.000000e+00 "<<dz<<endl;
	file_out<<"delta 0.000000e+00 "<<dy<<" 0.000000e+00"<<endl;
	file_out<<"delta "<<dx<<" 0.000000e+00 0.000000e+00"<<endl;
	file_out<<"object 2 class gridconnections counts "<<nz<<" "<<ny<<" "<<nx<<endl;
	file_out<<"object 3 class array type double rank 0 items "<<nx*ny*nz<<" data follows"<<endl;
	for(iz=0;iz<nz;iz++){
		for(iy=0;iy<ny;iy++){
			y=-ymax+(iy*dy);
			for(ix=0;ix<nx;ix++){
				x=-xmax+(ix*dx);
				crt2cln(x,y,r_grd,t_grd);
				if(t_grd < 0.0)t_grd+=2*PI;
				ir_cln=(int)(r_grd/dr);
				r_cln=(dr/2.0)+ir_cln*dr;
				it_cln=(int)(t_grd/dt);
				t_cln=(dt/2.0)+it_cln*dt;
				//cout<<x<<" "<<y<<" "<<r_cln*cos(t_cln)<<" "<<r_cln*sin(t_cln)
				//	<<" "<<r_grd<<" "<<t_grd<<" "<<ir_cln<<" "<<it_cln<<" "<<r_cln<<" "<<t_cln<<endl;
				iv=ir_cln+it_cln*nr+iz*nrnt;
				switch (data_type){
					case 0:
						if((satp_flag)&&(phieven[iv]>satp))dtmp=satp;
						else dtmp=phieven[iv];
						if((satm_flag)&&(dtmp<satm))dtmp=satm;
						break;
					case 1:
						if((satp_flag)&&(keven[iv]>satp))dtmp=satp;
						else dtmp=keven[iv];
						if((satm_flag)&&(dtmp<satm))dtmp=satm;
						break;
					case 2:
						if((satp_flag)&&(cleven[iv]>satp))dtmp=satp;
						else dtmp=cleven[iv];
						if((satm_flag)&&(dtmp<satm))dtmp=satm;
						break;
				}
				if(dtmp>max)max=dtmp;
				if(dtmp<min)min=dtmp;
				file_out<<setw(8)<<dtmp<<"\t";
				if((dummy_cnt%3)==0)file_out<<endl;
				dummy_cnt++;
			}
		}
	}
	file_out<<endl;
	file_out<<"attribute \"dep\" string \"positions\""<<endl;
	file_out<<"object \"regular positions regular connections\" class field"<<endl;
	file_out<<"component \"positions\" value 1"<<endl;
	file_out<<"component \"connections\" value 2"<<endl;
	file_out<<"component \"data\" value 3"<<endl;
	file_out<<endl;

	if(outbar_flag){
		dummy_cnt=1;
		if(satp_flag)max=satp;
		if(satm_flag)min=satm;
		nz=ny=2;
		file_outbar<<setiosflags(ios::fixed);
		file_outbar<<setprecision(3);
		file_outbar<<"object 1 class gridpositions counts "<<nz<<" "<<ny<<" "<<nx<<endl;
		file_outbar<<"origin "<<-ymax-4*dy<<" "<<-xmax-4*dx<<" "<<zmin<<endl;
		file_outbar<<"delta 0.000000e+00 0.000000e+00 "<<dz<<endl;
		file_outbar<<"delta 0.000000e+00 "<<ymax*2<<" 0.000000e+00"<<endl;
		file_outbar<<"delta "<<dx<<" 0.000000e+00 0.000000e+00"<<endl;
		file_outbar<<"object 2 class gridconnections counts "<<nz<<" "<<ny<<" "<<nx<<endl;
		file_outbar<<"object 3 class array type double rank 0 items "<<nx*ny*nz<<" data follows"<<endl;
		for(iz=0;iz<nz;iz++){
			for(iy=0;iy<ny;iy++){
				for(ix=0;ix<nx;ix++){
					file_outbar<<setw(8)<<min+((max-min)/(nx-1))*ix<<"\t";
					if((dummy_cnt%3)==0)file_outbar<<endl;
					dummy_cnt++;
				}
			}
		}
		file_outbar<<endl;
		file_outbar<<"attribute \"dep\" string \"positions\""<<endl;
		file_outbar<<"object \"regular positions regular connections\" class field"<<endl;
		file_outbar<<"component \"positions\" value 1"<<endl;
		file_outbar<<"component \"connections\" value 2"<<endl;
		file_outbar<<"component \"data\" value 3"<<endl;
	}

//ANALYSIS	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();
	file_outbar.close();

  	return 0;
}

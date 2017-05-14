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
#include <lib/molecule.h>

int printhelp(char *prm_name);

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],filein2[201];
	ifstream file_in;
	ifstream file_in2;
	ofstream file_out;

	molecule inpdb,in2pdb, tmpmol;
	double teta;
	int dir;
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
	strcat(filein,".in.pdb");
	strcpy(filein2,prm_name);
	strcat(filein2,".in2.pdb");
	strcpy(fileout,prm_name);
	strcat(fileout,".out.pdb");

//Data initialization
	teta=0.0;

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
	      	  else if(strcmp(argv[ind_arg]+1,"i2")==0){//Input File - For heterooligomers
	      		  strcpy(filein2,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"o")==0){//Output File
	      		  strcpy(fileout,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"teta")==0){//Theoretical Rotation
	      		  teta=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"cw")==0){//Rotate Clockwise
			  dir=1;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"acw")==0){//Rotate Anticlockwise
			  dir=-1;
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
  	file_in2.open(filein2);
	if(!file_in2){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein2<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}
	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//START-UP	END
//---------------------------------------------------------
  
//---------------------------------------------------------
//SEARCHING AXIS 	BEGIN

	inpdb.readpdb(file_in);
	in2pdb.readpdb(file_in2);

	double range,delta,rmsd,rmsd_min;
	double x,y,z,x_min,y_min,z_min;
	rmsd_min=INF;
	range=0.3;
	delta=0.01;

	x=-range;
	while(x<range+delta){
		y=-range;
		while(y<range+delta){
			tmpmol.copy(inpdb);
			z=sqrt(1.0-x*x-y*y);
			tmpmol.rotate_axis(x,y,z,dir*teta);
			rmsd=tmpmol.rmsd(in2pdb);
			if(rmsd<rmsd_min){
				rmsd_min=rmsd;
				x_min=x;
				y_min=y;
				z_min=z;
			}
			//cout<<"x "<<x<<" y "<<y<<" z "<<z<<" rmsd = "<<rmsd<<endl;
			y+=delta;
		}
		x+=delta;
	}
	tmpmol.copy(inpdb);
	tmpmol.rotate_axis(x_min,y_min,z_min,dir*teta);
	//tmpmol.rotate_axis(0.0,0.0,1.0,dir*teta);
	rmsd=tmpmol.rmsd(in2pdb);
	cout<<"x_min "<<x_min<<" y_min "<<y_min<<" z_min "<<z_min<<" rmsd_min = "<<rmsd<<endl;
	file_out<<"REMARK x_min "<<x_min<<" y_min "<<y_min<<" z_min "<<z_min<<" rmsd_min = "<<rmsd_min<<endl;
	file_out<<tmpmol<<endl;
//SEARCHING AXIS 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_in2.close();
	file_out.close();

  	return 0;
}

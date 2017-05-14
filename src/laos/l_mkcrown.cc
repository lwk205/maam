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

	molecule inpdb,in2pdb;
	bool hetero_flag;
	int np,dir;
	double r,tcenter,tcenter2;
	double x_axis,y_axis,z_axis,norm_axis;
	int i;
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
	hetero_flag=false;
	np=6;
	r=0.0;
	x_axis=y_axis=0.0;
	z_axis=1.0;

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
			  hetero_flag=true;
	      		  strcpy(filein2,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"o")==0){//Output File
	      		  strcpy(fileout,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"r")==0){//Crown Radius
	      		  r=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"np")==0){//Number of Promoters
	      		  np=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"cw")==0){//Rotate Clockwise
			  dir=1;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"acw")==0){//Rotate Anticlockwise
			  dir=-1;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"x_axis")==0){//X direction  Rotation axis
	      		  x_axis=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"y_axis")==0){//Y direction  Rotation axis
	      		  y_axis=atof(argv[ind_arg+1]);
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
	z_axis=sqrt(1-x_axis*x_axis-y_axis*y_axis);
	norm_axis=sqrt(x_axis*x_axis+y_axis*y_axis+z_axis*z_axis);
        if((norm_axis>1+1e-10)||(norm_axis<1-1e-10)){
                        cerr<<"WARNING: Rotation axis is not a versor norm = "<<norm_axis<<endl;
	}
  	file_in.open(filein);
	if(!file_in){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein<<endl;
		return 1;
	}
	if(hetero_flag){
  		file_in2.open(filein2);
		if(!file_in2){
			cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein2<<endl;
			return 1;
		}
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
//MAKING POLYMER 	BEGIN

	inpdb.readpdb(file_in);
	tcenter=atan2(inpdb.ycenter,inpdb.xcenter);
	cout<<"Wild-type radius = "<<inpdb.xcenter<<endl;
	if(r!=0){
		inpdb.rotate_axis(x_axis,y_axis,z_axis,tcenter);
		inpdb.translate(r-inpdb.xcenter,0.0,0.0);
		cout<<"Setting radius to = "<<r<<endl;
	}

	if(hetero_flag){
		in2pdb.readpdb(file_in2);
		tcenter2=atan2(in2pdb.ycenter,in2pdb.xcenter);
		inpdb.rotate_axis(x_axis,y_axis,z_axis,tcenter2);
		cout<<"Wild-type radius (Promoter B) = "<<in2pdb.xcenter<<endl;
		in2pdb.translate(r-in2pdb.xcenter,0.0,0.0);
	}

	inpdb.setchain('A');
	file_out<<inpdb<<endl; //Output first Promoter
	for(i=1;i<np;i++){
		if(hetero_flag){
			if((i%2)==1){ // Promoter B
				cout<<"Promoter B, Rotation = "<<4.0*PI/np<<endl;
				if(i==1) in2pdb.rotate_axis(x_axis,y_axis,z_axis,dir*2.0*PI/np);
				else in2pdb.rotate_axis(x_axis,y_axis,z_axis,dir*4.0*PI/np);
				in2pdb.setchain(i+65);
				file_out<<in2pdb<<endl; //Output Successive Promoters
			} else { // Promoter A
				cout<<"Promoter A, Rotation = "<<4.0*PI/np<<endl;
				inpdb.rotate_axis(x_axis,y_axis,z_axis,dir*4.0*PI/np);
				inpdb.setchain(i+65);
				file_out<<inpdb<<endl; //Output Successive Promoters
			}
		} else {
			inpdb.rotate_axis(x_axis,y_axis,z_axis,dir*2.0*PI/np);
			inpdb.setchain(i+65);
			file_out<<inpdb<<endl; //Output Successive Promoters
		}
	}

//MAKING POLYMER 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();
	if(hetero_flag)file_in2.close();

  	return 0;
}

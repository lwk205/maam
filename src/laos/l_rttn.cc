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

int printhelp(char *prm_name);

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201];
	ifstream file_in;
	ofstream file_out;

	molecule inpdb;

	double alpha,beta,gamma;
	double rad_alpha,rad_beta,rad_gamma;
	bool rad,deg;

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
	strcpy(fileout,prm_name);
	strcat(fileout,".out.pdb");

//Data initialization
	alpha=beta=gamma=0.0;
	rad_alpha=rad_beta=rad_gamma=0.0;
	rad=false;
	deg=true;

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
	      	  else if(strcmp(argv[ind_arg]+1,"rad")==0)//Radiant
	      	  {
			  rad=true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"deg")==0)//Degree
	      	  {
			  deg=true;
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
	      	  else if(strcmp(argv[ind_arg]+1,"x")==0)//x movement
	      	  {
	      		  alpha=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"y")==0)//y movement
	      	  {
	      		  beta=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"z")==0)//z movement
	      	  {
	      		  gamma=atof(argv[ind_arg+1]);
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
//ROTATION 	BEGIN
	inpdb.readpdb(file_in);
	if(rad&&deg)
	{
		cerr<<"ERROR: ANGLES CAN NOT BE BOTH IN RADIANTS AND DEGREE"<<endl;
		exit(1);
	}
	if(((alpha*beta)!=0.0)||((alpha*beta)!=0.0)||((alpha*beta)!=0.0))
	{
		cerr<<"WARNING: MORE THAN ONE ROTATION REQUESTED !"<<endl;
		cerr<<"         ROTATIONS WILL BE APPLIED IN THE ORDER x - y - z"<<endl; 
	}
	if(deg)
	{
		rad_alpha=deg2rad(alpha);
		rad_beta=deg2rad(beta);
		rad_gamma=deg2rad(gamma);
	}
	else
	{
		rad_alpha=alpha;
		rad_beta=beta;
		rad_gamma=gamma;
	}
	if(rad_alpha!=0.0)inpdb.rotate_x(rad_alpha);
	if(rad_beta!=0.0)inpdb.rotate_y(rad_beta);
	if(rad_gamma!=0.0)inpdb.rotate_z(rad_gamma);
	file_out<<inpdb;
//ROTATION 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();

  	return 0;
}

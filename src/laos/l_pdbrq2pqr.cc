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

	char filein[201],fileout[201];
	ifstream file_in;
	ofstream file_out;

	int ind_atm;
	double dtmp;

	molecule inpdb;

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
//CONVERSION 	BEGIN
	inpdb.readpdb(file_in);
	for(ind_atm=0;ind_atm<inpdb.num_atm;ind_atm++){
		dtmp=inpdb.mol[ind_atm].rad;
		inpdb.mol[ind_atm].rad=inpdb.mol[ind_atm].chg;
		inpdb.mol[ind_atm].chg=dtmp;
	}
	file_out<<inpdb;
//CONVERSION 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();

  	return 0;
}

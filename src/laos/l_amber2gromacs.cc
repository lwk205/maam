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
#include <lib/atom.h>

int printhelp(char *prm_name);

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201],fileout[201],fileconv[201];
	ifstream file_in;
	ifstream file_conv;
	ofstream file_out;

	molecule inpdb;

	int ind_atm;
	atom in_atom;

	char line[201];
	char amberatmname[5];
	char amberresname[5];
	char gromacsatmname[5];
	char gromacsresname[5];
	bool find;


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

//Temporary solution	
	strcpy(fileconv,"/home/simone/work/MySoftware/maam-1.0/dat/charmm2amber.dat");

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
  	file_conv.open(fileconv);
	if(!file_conv)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileconv<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}

	inpdb.readpdb(file_in);
	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//START-UP	END
//---------------------------------------------------------
  

//---------------------------------------------------------
//CONVERT ATOM NAME IN GROMACS FORMAT 	BEGIN

	for(ind_atm=0;ind_atm<inpdb.num_atm;ind_atm++)//For every atom...
	{
		//Find the equivalent nomenclature in gromacs
		in_atom=inpdb[ind_atm];
		file_conv.clear();
		file_conv.seekg(0,ios::beg);
		file_conv.getline(line,200);
		find=false;
		while((!file_conv.eof())&&(!find))
		{
			strncpy(amberatmname,line,4);
			amberatmname[4]='\0';
			strncpy(amberresname,line+5,4);
			amberresname[4]='\0';
			strncpy(gromacsatmname,line+10,4);
			gromacsatmname[4]='\0';
			strncpy(gromacsresname,line+15,4);
			gromacsresname[4]='\0';
			if(in_atom.cmpatomname(amberatmname)&&
			   in_atom.cmpresname(amberresname))find=true;
			file_conv.getline(line,200);
		}
		if(find){
			//cout<<"FOUND: "<<ind_atm<<endl;
			if(strcmp(gromacsatmname,"NONE")==0){
				cout<<"deleting: "<<in_atom<<endl;
			}
			else{
				strcpy(in_atom.name,gromacsatmname);
				strcpy(in_atom.res,gromacsresname);
				file_out<<in_atom<<endl;
			}
		}
		else{
			cerr<<"ERROR ATOM NOT FOUND"<<endl;
			cerr<<in_atom<<endl;
			exit(1);
		}
	}

//CONVERT ATOM NAME IN GROMACS FORMAT 	BEGIN
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_conv.close();
	file_out.close();

  	return 0;
}

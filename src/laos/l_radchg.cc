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

	char filein[201],fileout[201],fileprm[201];
	ifstream file_in;
	ifstream file_prm;
	ofstream file_out;

	molecule inpdb;
	molecule prmpdb;

	int res_str,res_last;
	char chn;
	bool flag_0chg,flag_nochg;
	int ind_atm;
	int term;
	atom in_atom, prm_atom;

//---------------------------------------------------------
//START-UP	BEGIN
//Floatint point exception
#ifdef FLOATINGPOINT
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

//Initialize
	res_str=-1;
	res_last=-1;
	chn='z';
	flag_nochg = false;

//Data file	
	strcpy(fileprm,DATADIR);
	strcat(fileprm,"/radchg_Nina_Charmm.pdb");

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
	      	  else if(strcmp(argv[ind_arg]+1,"c_0chg")==0)//Chain with 0 charge
	      	  {
			  chn=argv[ind_arg+1][0];
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"nochg")==0)//All charges set to zero
	      	  {
			  flag_nochg = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"res_str_0chg")==0)//First Residue with 0 charge
	      	  {
			  res_str=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"res_last_0chg")==0)//Last Residue with 0 charge
	      	  {
			  res_last=atoi(argv[ind_arg+1]);
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

//Controll parameters
	if((chn!='z')||(res_str!=-1)||(res_last!=-1)){
		if((chn=='z')||(res_str==-1)||(res_str==-1)){
			cerr<<"ERROR in parameters"<<endl;
			cerr<<"chn = "<<chn<<" res_str = "<<res_str<<" res_last = "<<res_last<<endl;
			exit(1);
		}
		else flag_0chg=true;
	}
	else flag_0chg=false;

//Input-Output control
  	file_in.open(filein);
	if(!file_in)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein<<endl;
		return 1;
	}
  	file_prm.open(fileprm);
	if(!file_prm)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileprm<<endl;
		return 1;
	}
  	file_out.open(fileout);
	if(!file_out)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}

//Molecule reading
	inpdb.readpdb(file_in);
	prmpdb.readpdb(file_prm);

	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//START-UP	END
//---------------------------------------------------------
  

//---------------------------------------------------------
//ADDING RADII AND CHARGE PARMETERS 	BEGIN
	
	//Per convertire il file di parametri
	//for(ind_atm=0;ind_atm<prmpdb.num_atm;ind_atm++){
	//	prmpdb.mol[ind_atm].x=prmpdb.mol[ind_atm].chg;
	//	prmpdb.mol[ind_atm].y=prmpdb.mol[ind_atm].chg;
	//	prmpdb.mol[ind_atm].z=0.0;
	//}
	//cout<<prmpdb;
	  
	for(ind_atm=0;ind_atm<inpdb.num_atm;ind_atm++)//For every atom...
	{
		//Find the parameters
		in_atom=inpdb[ind_atm];
		prm_atom=prmpdb.find(in_atom.res,in_atom.name);
		if(strcmp(prm_atom.name,"NULL")==0)
		{
			cerr<<"WARNING ATOM NOT FOUND IN THE PARAMETER FILE"<<endl;
			cerr<<in_atom<<endl;
		}
		else
		{
			term=inpdb.ister(ind_atm);
			switch (term){
				case -1:
					inpdb[ind_atm].chg=prm_atom.x;
					inpdb[ind_atm].rad=prm_atom.rad;
					break;
				case 0:
					inpdb[ind_atm].chg=prm_atom.chg;
					inpdb[ind_atm].rad=prm_atom.rad;
					break;
				case 1:
					inpdb[ind_atm].chg=prm_atom.y;
					inpdb[ind_atm].rad=prm_atom.rad;
					break;
				default:
					cerr<<"ERROR IN l_radchg: term = "<<term<<" ind_atm = "<<ind_atm<<endl;
					exit(1);

			}
			if(inpdb[ind_atm].rad==0)
			{
				cerr<<"WARNING Radius value equal to zero !\n";
				cerr<<in_atom<<endl;
			}
		}
		if((flag_0chg)&&(inpdb[ind_atm].chn==chn)&&
			(inpdb[ind_atm].res_num>=res_str)&&
			(inpdb[ind_atm].res_num<=res_last))inpdb[ind_atm].chg=0.0;
		if(flag_nochg)inpdb[ind_atm].chg=0.0;
	}

	file_out<<inpdb;
//ADDING RADII AND CHARGE PARMETERS 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_prm.close();
	file_out.close();

  	return 0;
}

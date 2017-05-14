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

	char filein[201],fileout[201],stmp[201];
	ifstream file_in;
	ofstream file_out;

	molecule inpdb;
	int ind_atm;
	double rbarrel;
	int nCA,N;
	double a,b;

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

//Parameter definition
	//cout<<"Distance between CA along the sheet [3.3 A] = ";
	//cin>>a;
	a=3.3;
	//cout<<"Distance among sheets [4.4 A] = ";
	//cin>>b;
	b=4.4;
	N=0;

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
	      	  else if(strcmp(argv[ind_arg]+1,"N")==0)//Number of b-sheets
	      	  {
	      		  strcpy(stmp,argv[ind_arg+1]);
			  N=atoi(stmp);
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
	if(N<=0){
		cerr<<"ERROR: l_flatbarrel Number of beta sheet <= 0"<<endl;
		exit(1);
	}

//Molecule reading
	inpdb.readpdb(file_in);

	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//START-UP	END
//---------------------------------------------------------

//---------------------------------------------------------
//FLATTING THE BETA-BARREL 	BEGIN
	rbarrel=0.0;
	nCA=0;
	//First cycle: define barrel radius
	//Radii of alpha-Carbons are used
	for(ind_atm=0;ind_atm<inpdb.num_atm;ind_atm++){
		if(inpdb[ind_atm].cmpatomname("CA")){
			if((inpdb[ind_atm].r)<=0.0){
				cerr<<"ERROR: l_flatbarrel Negative barrel radius"<<endl;
				exit(1);
			}
			rbarrel+=inpdb[ind_atm].r;
			nCA++;
			//DEBUG
			//file_out<<inpdb[ind_atm]<<" "<<inpdb[ind_atm].r<<" "<<inpdb[ind_atm].t<<endl;
		}
	}
	rbarrel=rbarrel/nCA;
	file_out<<"REMARK Beta-Barrel Radius = "<<rbarrel<<endl;
	int ind_chn,new_chn;
	double oldt,ttmp;
	ind_chn=0;
	new_chn=1;
	for(ind_atm=0;ind_atm<inpdb.num_atm;ind_atm++){
		if(ind_atm==inpdb.chn_end[ind_chn]){
			ind_chn++;
			new_chn=1;
			cerr<<inpdb[ind_atm]<<endl;
		}
		if(new_chn){
			cout<<endl<<endl;
			oldt=inpdb[ind_atm].t;
			new_chn=0;
		}else{
			if((oldt<-(PI/2))&&(inpdb[ind_atm].t>=0)){
				ttmp=(inpdb[ind_atm].t)-(2*PI);
				oldt=ttmp;
			}else if((oldt>(PI/2))&&(inpdb[ind_atm].t<=0)){
				ttmp=(inpdb[ind_atm].t)+(2*PI);
				oldt=ttmp;
			}else{
				ttmp=(inpdb[ind_atm].t);
				oldt=ttmp;
			}
			cout<<ttmp<<endl;
		}
		inpdb[ind_atm].y=rbarrel*ttmp;
		inpdb[ind_atm].x=(inpdb[ind_atm].r)-rbarrel;
		file_out<<inpdb[ind_atm]<<endl;
	}
	//Utilizzare catena 2 per aH epta
	//for(ind_atm=inpdb.chn_end[0];ind_atm<inpdb.chn_end[1];ind_atm++){
	//Utilizzare catena 6 per aH esa
	for(ind_atm=inpdb.chn_end[4];ind_atm<inpdb.num_atm;ind_atm++){
		if(ind_atm==0){
			oldt=inpdb[ind_atm].t;
		}else{
			if((oldt<-(3*PI/4))&&(inpdb[ind_atm].t>=0)){
				(inpdb[ind_atm].t)-=(2*PI);
				oldt=inpdb[ind_atm].t;
			}else if((oldt>(3*PI/4))&&(inpdb[ind_atm].t<=0)){
				(inpdb[ind_atm].t)+=(2*PI);
				oldt=inpdb[ind_atm].t;
			}else oldt=inpdb[ind_atm].t;
			//cout<<inpdb[ind_atm].t<<endl;
		}
		(inpdb[ind_atm].t)+=2*PI;
		inpdb[ind_atm].y=rbarrel*inpdb[ind_atm].t;
		inpdb[ind_atm].x=(inpdb[ind_atm].r)-rbarrel;
		file_out<<inpdb[ind_atm]<<" "<<inpdb[ind_atm].t<<endl;
	}
	cerr<<sqrt(N*N*( 4.0*rbarrel*rbarrel*sin(PI/N)*sin(PI/N) - b*b )/(a*a))<<endl;

//FLATTING THE BETA-BARREL 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();

  	return 0;
}

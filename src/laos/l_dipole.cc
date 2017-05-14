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

	char filein[201],fileout[201];
	ifstream file_in;
	ofstream file_out;

	molecule inpdb;

	int ind_atm;
	atom pos_atoms(1,1,"POT","POT");
	atom neg_atoms(2,1,"CLA","CLA");
	bool flag_vdw=false;
	double pos_norm=0.0;
	double neg_norm=0.0;
	double w;
	double dist, dipole;
	double phi,psi,distxy;

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
	      	  else if(strcmp(argv[ind_arg]+1,"vdw")==0)//Charges weighted by Van der Walls radii
	      	  {
			  flag_vdw=true;
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

//Molecule reading
	inpdb.readpdb(file_in);

	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//START-UP	END
//---------------------------------------------------------
  

//---------------------------------------------------------
//COMPUTE DIPOLE 	BEGIN

	for(ind_atm=0;ind_atm<inpdb.num_atm;ind_atm++)//For every atom...
	{
		if(inpdb[ind_atm].chg>0){
			pos_atoms.chg+=(inpdb[ind_atm].chg);
			if(flag_vdw) w=(inpdb[ind_atm].rad)*(flag_vdw);
			else w=1.0;
			pos_atoms.x+=((inpdb[ind_atm].x)*w);
			pos_atoms.y+=((inpdb[ind_atm].y)*w);
			pos_atoms.z+=((inpdb[ind_atm].z)*w);
			if(flag_vdw)pos_norm+=(inpdb[ind_atm].rad);
			else pos_norm+=1.0;
			//cerr<<"REMARK Positive = "<<pos_atoms.chg<<endl;
			//cerr<<inpdb[ind_atm]<<endl;

		}
		else{
			neg_atoms.chg+=(inpdb[ind_atm].chg);
			if(flag_vdw) w=(inpdb[ind_atm].rad)*(flag_vdw);
			else w=1.0;
			neg_atoms.x+=((inpdb[ind_atm].x)*w);
			neg_atoms.y+=((inpdb[ind_atm].y)*w);
			neg_atoms.z+=((inpdb[ind_atm].z)*w);
			if(flag_vdw)neg_norm+=(inpdb[ind_atm].rad);
			else neg_norm+=1.0;
			//cerr<<"REMARK Negative = "<<neg_atoms.chg<<endl;
			//cerr<<inpdb[ind_atm]<<endl;
		}
	}
	pos_atoms.x=((pos_atoms.x)/(pos_norm));
	pos_atoms.y=((pos_atoms.y)/(pos_norm));
	pos_atoms.z=((pos_atoms.z)/(pos_norm));
	neg_atoms.x=((neg_atoms.x)/(neg_norm));
	neg_atoms.y=((neg_atoms.y)/(neg_norm));
	neg_atoms.z=((neg_atoms.z)/(neg_norm));
	dist=pos_atoms.distance(neg_atoms);
	if( (pos_atoms.chg > (abs(neg_atoms.chg)+1.0)) || (pos_atoms.chg < (abs(neg_atoms.chg)-1.0)) ){
		file_out<<"REMARK WARNING: Positive charge = "<<pos_atoms.chg
			<<" Negative charge = "<<neg_atoms.chg<<endl;
		cerr<<"WARNING: Positive charge = "<<pos_atoms.chg
			<<" Negative charge = "<<neg_atoms.chg<<endl;
	}
	dipole=(pos_atoms.chg)*dist;

	if(flag_vdw)file_out<<"REMARK Atom positions weighted by Van der Walls radii"<<endl;
	file_out<<"REMARK DIPOLE = "<<dipole<<" eA"<<endl;
	file_out<<"REMARK DIPOLE = "<<dipole*eA2Deybe<<" Debye"<<endl;
	file_out<<"REMARK DIPOLE = "<<dipole*1.602<<"*10^-29 Cm"<<endl;

	//Output dipole strength, moving negative pole
	phi=acos((neg_atoms.z-pos_atoms.z)/dist);
	distxy=dist*sin(phi);
	psi=acos((neg_atoms.x-pos_atoms.x)/distxy);
	neg_atoms.x=pos_atoms.x+(dipole*eA2Deybe*sin(phi)*cos(psi));
	neg_atoms.y=pos_atoms.y+(dipole*eA2Deybe*sin(phi)*sin(psi));
	neg_atoms.z=pos_atoms.z+(dipole*eA2Deybe*cos(phi));
	  
	//Output dipole strength, moving positive pole
	//phi=acos((pos_atoms.z-neg_atoms.z)/dist);
	//distxy=dist*sin(phi);
	//psi=acos((pos_atoms.x-neg_atoms.x)/distxy);
	//pos_atoms.x=neg_atoms.x+(dipole*eA2Deybe*sin(phi)*cos(psi));
	//pos_atoms.y=neg_atoms.y+(dipole*eA2Deybe*sin(phi)*sin(psi));
	//pos_atoms.z=neg_atoms.z+(dipole*eA2Deybe*cos(phi));

	file_out<<pos_atoms<<endl;
	file_out<<neg_atoms<<endl;
//COMPUTE DIPOLE 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_in.close();
	file_out.close();

  	return 0;
}

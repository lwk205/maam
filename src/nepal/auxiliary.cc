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

////////////////////////////////////////////////////////////////////////////////
//VOLUME CELL
//Volume elemento di griglia
double volume_cell(int ir,double dr,double dt,double dz)
{
	return dr*dr*dz*dt*(ir+0.5);
}
////////////////////////////////////////////////////////////////////////////////
  
////////////////////////////////////////////////////////////////////////////////
//PRINTHELP
//Output messaggi di help
int printhelp(char *prm_name)
{
	cout<<"--------------------------------------------------"<<endl;
	if(strcmp(prm_name,"b_axis")==0)
	{
		cout<<"USAGE: b_axis [-h] [-i b_axis.in.end] [-o b_axis.out.dat]"<<endl<<endl;
		cout<<"Extract electrostatic potential, ion concentratios and occupancy indexes along the channel axis from a binary file"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option";
		cout<<setw(25)<<"Default Value";
		cout<<setw(35)<<endl;
		cout<<setw(15)<<"-ztrsl";
		cout<<setw(25)<<0.0;
		cout<<setw(35)<<"Translate the z axis"<<endl;
	}
	else if(strcmp(prm_name,"b_bind")==0)
	{
		cout<<"USAGE: b_axis [-h] [-i b_axis.in.end] [-bs b_bind.in.bs] [-o b_axis.out.dat]"<<endl<<endl;
		cout<<"Compute the number of ions in the binding sites stated in b_bind.in.bs"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option";
		cout<<setw(25)<<"Default Value";
		cout<<setw(35)<<endl;
		cout<<setw(15)<<"-ztrsl";
		cout<<setw(25)<<0.0;
		cout<<setw(35)<<"Translate the z axis"<<endl;
	}
	else if(strcmp(prm_name,"b_forces")==0)
	{
		cout<<"USAGE: b_forces [-h] [-i b_forces.in.end] [-str b_forces.in.pdbrq] [-o b_forces.out.dat]"<<endl<<endl;
		cout<<"Compute the electrical force on a the molecule -str"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option";
		cout<<setw(25)<<"Default Value";
		cout<<setw(35)<<endl;
	}
	else if(strcmp(prm_name,"b_forcesmomenta")==0)
	{
		cout<<"USAGE: b_forcesmomenta [-h] [-i b_forces.in.end] [-str b_forces.in.pdbrq] [-x_hinge 0.0]  [-y_hinge 0.0] [-z_hinge 0.0] [-o b_forces.out.dat]"<<endl<<endl;
		cout<<"Compute the electrical force on a the molecule -str"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option";
		cout<<setw(25)<<"Default Value";
		cout<<setw(35)<<endl;
	}
	else if(strcmp(prm_name,"nepal")==0)
	{
		cout<<"USAGE: nepal [-h] [-i nepal.prm] [-o nepal.log]"<<endl<<endl;
		cout<<"Solve the Poisson-Nerst-Planck equations"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option";
		cout<<setw(25)<<"Default Value";
		cout<<setw(35)<<endl;
	}
	else {
		cerr<<"ERROR IN SUBROUTINE: int printhelp(char *prm_name)"<<endl;
		cerr<<"PROGRAM "<<prm_name<<" does not exist"<<endl;
		exit(1);
	}
	cout<<"--------------------------------------------------"<<endl;

  	return 0;
}

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
//int printhelp(char *prm_name)
//Print the help for the program prm_name
int printhelp(char *prm_name)
{
	cout<<"--------------------------------------------------"<<endl;
	if(strcmp(prm_name,"l_trsl")==0)
	{
		cout<<"USAGE: l_trsl [-h] [-i l_trsl.in.pdb] [-o l_trsl.out.pdb] [-x xmove] [-y ymove] [-z zmove] [-zind]"
			<<endl<<endl;
		cout<<"Move the molecule read in l_trsl.in.pdb by {xmove ymove zmove}"<<endl;
		cout<<"If no translation is specified the molecule is moved to its geometric center"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_trsl.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_trsl.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-x"<<setw(25)<<"0.0"<<setw(35)<<"x translation"<<endl;
		cout<<setw(15)<<"-y"<<setw(25)<<"0.0"<<setw(35)<<"y translation"<<endl;
		cout<<setw(15)<<"-z"<<setw(25)<<"0.0"<<setw(35)<<"z translation"<<endl;
		cout<<setw(15)<<"-zind"<<setw(25)<<"false"<<setw(35)<<"center z of atom index"<<endl;
	}
	else if(strcmp(prm_name,"l_rttn")==0)
	{
		cout<<"USAGE: l_rttn [-h] [-i l_rttn.in.pdb] [-o l_rttn.out.pdb] [-x alpha] [-y beta] [-z gamma]"
			<<endl<<endl;
		cout<<"Rotate the molecule counter-clockwise around the selected axes"<<endl;
		cout<<"If more than one rotatiotion is given, these are applied in the order x - y - z"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_rttn.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_rttn.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-rad"<<setw(25)<<"false"<<setw(35)<<"Angle in radiants"<<endl;
		cout<<setw(15)<<"-deg"<<setw(25)<<"true"<<setw(35)<<"Angle in degrees"<<endl;
		cout<<setw(15)<<"-x"<<setw(25)<<"0.0"<<setw(35)<<"x rotation"<<endl;
		cout<<setw(15)<<"-y"<<setw(25)<<"0.0"<<setw(35)<<"y rotation"<<endl;
		cout<<setw(15)<<"-z"<<setw(25)<<"0.0"<<setw(35)<<"z rotation"<<endl;
	}
	else if(strcmp(prm_name,"l_radchg")==0)
	{
		cout<<"USAGE: l_radchg [-h] [-i l_rttn.in.pdb] [-o l_rttn.out.pdb]"
			<<endl<<endl;
		cout<<"Add radii and charges to a pdb file"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_rttn.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_rttn.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
	}
	else if(strcmp(prm_name,"l_fixnumber")==0)
	{
		cout<<"USAGE: l_fixnumber [-h] [-i l_fixnumber.in.pdb] [-o l_fixnumber.out.pdb]"
			<<endl<<endl;
		cout<<"Adjust the residue and atom number"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_fixnumber.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_fixnumber.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-bres"<<setw(25)<<"First Residue"<<setw(35)<<"1"<<endl;
		cout<<setw(15)<<"-jatom"<<setw(25)<<"Fix only atoms"<<setw(35)<<"False"<<endl;
	}
	else if(strcmp(prm_name,"l_amber2gromacs")==0)
	{
		cout<<"USAGE: l_amber2gromacs [-h] [-i l_amber2gromacs.in.pdb] [-o l_amber2gromacs.out.pdb]"
			<<endl<<endl;
		cout<<"Rotate the molecule counter-clockwise around the selected axes"<<endl;
		cout<<"If more than one rotatiotion is given, these are applied in the order x - y - z"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_amber2gromacs.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_amber2gromacs.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
	}
	else if(strcmp(prm_name,"l_dipole")==0)
	{
		cout<<"USAGE: l_dipole [-h] [-i l_dipole.in.pdb] [-o l_dipole.out.pdb] [-vdw]"
			<<endl<<endl;
		cout<<"Compute the electric dipole of the input molecule"<<endl;
		cout<<"Input molecule must have charges values in the occupancy colunm"<<endl;
		cout<<"Option -vdw = Charges weighted by Van der Walls radii"<<endl;
		cout<<"Produce in output a pdb file with the following format:"<<endl;
		cout<<"REMARK Molecule Dipole = X Deybe = Y Cm"<<endl;
		cout<<"ATOM 1: Position and value (in the occupancy colunm) for the baricenter of the positive charges"
			<<endl;
		cout<<"ATOM 1: Position and value (in the occupancy colunm) for the baricenter of the negative charges"
			<<endl;
		cout<<"--------------------------------------------------"<<endl;
	}
	else if(strcmp(prm_name,"l_betabarrel")==0)
	{
		cout<<"USAGE: l_betabarrel [-h] [-i l_betabarrel.in.pdb] [-o l_betabarrel.out.pdb] [-N 0] [-S 0] [-len_sheet 0]"
			<<endl<<endl;
		cout<<"Define a beta-barrel backbone structure"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_betabarrel.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_betabarrel.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-N"<<setw(25)<<"0"<<setw(35)<<"Number of beta-sheets"<<endl;
		cout<<setw(15)<<"-S"<<setw(25)<<"0"<<setw(35)<<"Sheer number"<<endl;
		cout<<setw(15)<<"-len-_sheet"<<setw(25)<<"0"<<setw(35)<<"Number of residues"<<endl;
	}
	else if(strcmp(prm_name,"l_flatbarrel")==0)
	{
		cout<<"USAGE: l_flatbarrel [-h] [-i l_flatbarrel.in.pdb] [-o l_flatbarrel.out.pdb] [-N 0]"
			<<endl<<endl;
		cout<<"Project a beta-barrel on a flat-surface"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_flatbarrel.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_flatbarrel.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-N"<<setw(25)<<"0"<<setw(35)<<"Number of beta-sheets"<<endl;
	}
	else if(strcmp(prm_name,"l_expandrad")==0)
	{
		cout<<"USAGE: l_expandrad-h] [-i l_expandrad.in.pdb] [-o l_expandrad.out.pdb] [-exp 1.0]"
			<<endl<<endl;
		cout<<"Expand a molecule in the radial direction"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_expandrad.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_expandrad.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-exp"<<setw(25)<<"1.0"<<setw(35)<<"Expansion"<<endl;
	}
	else if(strcmp(prm_name,"l_resetchg")==0)
	{
		cout<<"USAGE: l_resetchg [-h] [-i l_resetchg.in.pdb] [-o l_resetchg.out.pdb]"
			<<endl<<endl;
		cout<<"Set all the atomic charges to zero"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_resetchg.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_resetchg.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
	}
	else if(strcmp(prm_name,"l_mkcrown")==0)
	{
		cout<<"USAGE: l_mkcrown [-h] [-i l_mkcrown.in.pdb] [-o l_mkcrown.out.pdb]"
			<<endl<<endl;
		cout<<"A circular polymer (hexameric) is built using the monomer read from -i"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_mkcrown.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_mkcrown.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
		cout<<setw(15)<<"-r"<<setw(25)<<"0.0"<<setw(35)<<"Radius of the crown"<<endl;
	}
	else if(strcmp(prm_name,"l_pdbrq2pqr")==0) {
		cout<<"USAGE: l_pdbrq2pqr [-h] [-i l_pdbrq2pqr.pdb] [-o l_pdbrq2pqr.pdb]"
			<<endl<<endl;
		cout<<"Convert a pdbrq file to a pqr file [VMD usage]"<<endl;
		cout<<"--------------------------------------------------"<<endl;
		cout<<setw(15)<<"Option"<<setw(25)<<"Default Value"<<setw(35)<<endl;
		cout<<setw(15)<<"-h"<<setw(25)<<""<<setw(35)<<"Print help"<<endl;
		cout<<setw(15)<<"-i"<<setw(25)<<"l_pdbrq2pqr.in.pdb"<<setw(35)<<"Input Pdb File"<<endl;
		cout<<setw(15)<<"-o"<<setw(25)<<"l_pdbrq2pqr.out.pdb"<<setw(35)<<"Output Pdb File"<<endl;
	}
	else {
		cerr<<"ERROR IN SUBROUTINE: int printhelp(char *prm_name)"<<endl;
		cerr<<"PROGRAM "<<prm_name<<" does not exist"<<endl;
		exit(1);
	}
	cout<<"--------------------------------------------------"<<endl;

  	return 0;
}

#ifndef LIB_MOLECULE_H
#define LIB_MOLECULE_H

#include <include/common.h>
#include <lib/atom.h>

class molecule{
	public:
		atom  *mol;	//A molecule is a vector of atom (mol[0:num_atm-1])
		int *chn_end;	//Serial number of chain last atom (chn_end[0:num_chn])
		//Cicle on i-chain atoms:
		//for(ind_atm=chn_end[chn_i-1](0 if chn_i==0)
		//	;ind_atm<chn_end[chn_i];ind_atm++)
		int *res_end;	//Serial number of chain last atom (chn_end[0:num_chn])
		//Cicle on i-residue atoms:
		//for(ind_atm=res_end[res_i-1](0 if res_i==0)
		//	;ind_atm<res_end[res_i];ind_atm++)
		
		double *res_chg; //Residue charge res_chg[0:num_res-1]
		int *res_ter; 	//res_ter[0:num_res-1] = -1 C-terminal
		                //res_ter[0:num_res-1] = 0  Normal residue
		                //res_ter[0:num_res-1] = 1  N-terminal

		char name[101];	//Molecule Name
		int num_chn;	//Number of chains
		int num_res;	//Number of residues
		int num_atm;	//Number of atoms
		//[xyz]center used by l_trsl
		double xcenter,xmin,xmax;	//Geometrical characteristics
		double ycenter,ymin,ymax;
		double zcenter,zmin,zmax;
		double rmax;
		double chg;	//Charge

	//CONSTRUCTORS
		molecule(int numatm = 0, int numres = 0, int numchn = 0);
		molecule(ifstream &pdb_in);
		molecule(const molecule &obj);

	//INPUT FUNCTIONS	
		void clean();	//Reset a molecule
		void readpdb(ifstream &pdb_in);	//Read molecule from pdb
		void copy(const molecule &obj);	//Copy molecule

	//SELECTION FUNCTION
		atom find(char *res_name, char *atm_name);
			//Return the first atom (atm_name)
			//inside a residue (res_name) in the molecule
		int ister(int ind_atm);
			//Return -1 if in N-terminal residue
			//Return  1 if in C-terminal residue
			//Return  0 otherwise
			//Return  2 if ERROR
			  
	//TEXT EDITING
		void setchain(char chn);

	//GEOMETRICAL FUNCTIONS
		void fix();	//Center - Limits - Charge
		//translate used by l_trsl
		void translate(double *T);//Translation
		void translate(double xt, double yt, double zt);
		//rotate used by l_rtt
		void rotate(double *R);	//Rotation
		void rotate_x(double teta);
		void rotate_y(double teta);
		void rotate_z(double teta);
		void rotate_axis(double x1, double x2, double x3, double teta);
		double rmsd(const molecule &A); //RMSD between molecules
		int cmp(const molecule &A);	//cmp = 0 if the two molecule are identical
						//	1 Otherwise
	//MOLECULE PROPRIETIES
		double mass();//Molecule mass
		double centermass(double *bar);//Molecule center of mass - Return mass
		int inertia(double *I, double x0=0.0, double y0=0.0, double z0=0.0);//Inertia matrix
		
	//FRIEND OPERATORS
		atom& operator[](int ind_atm);

	//FRIEND OPERATORS
		friend ostream &operator<<(ostream &stream, molecule &obj);//Output Molecule in pdb

	//DESTRUCTOR
		~molecule();	//Distruttore
};

#endif //LIB_MOLECULE_H

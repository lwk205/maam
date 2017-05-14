#ifndef LIB_ATOM_H 
#define LIB_ATOM_H

class atom{
	public:
		int num;	//Serial number
		char name[5];	//Atom type
		char res[5];	//Residue name
		int res_num;	//Residue serial number
		char chn;	//Chain identifier
		double x,y,z;	//Cartesian coordinates
		double r,t;	//Cilindrical coordinates
		double rad;	//Radius                  
		double chg;	//Charge

		//Constructor
		atom();
		atom(int NUM,int RES_NUM,
			char *NAME, char *RES,
			char CHN='A',
			double X=0.0,double Y=0.0,double Z=0.0,
			double R=0.0,double T=0,
			double RAD=0.0,double CHG=0.0);

		//Function
		int clear();//Reset atom

		atom readpdb(char *line);//Read atom from pdb line

		int cmpatomname(char *atmname);//cmpname = 1 If (atmname == atm.name)
					       //        = 0 Otherwise		
					  
		int cmpresname(char *resname);//cmpname = 1 If (resname == atm.res)
					      //        = 0 Otherwise		
					  
		int isbackbone();//isbackbone = 1 If Backbone atom
				 //		0 otherwise
				 //In order to always have an elctrically neutral
				 //backbone:
				 //
				 //1)In Proline aminoacid H atoms bound 
				 //to the carbon atom in the ring are counted
				 //as backbone atoms
				 //
				 //2)In Glycine aminoacid all atoms are counted
				 //as backbone atoms

		int cmp(atom atm);//cmp = 1 If atoms are identical (name, residue, coordinates, charge)
				  //	= 0 Otherwise

		double distance(atom atm);//Compute the distance between two atoms
		
		double mass();

		//Friend functions
		friend double distance(atom atm1,atom atm2);//Atom distance

		//Friend operator
		friend ostream &operator<<(ostream &stream, atom atm);//Output pdb line

		//Destructor
		~atom();
};

#endif //LIB_ATOM_H

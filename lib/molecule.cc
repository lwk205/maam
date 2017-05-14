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
#include <lib/atom.h>
#include <lib/molecule.h>

molecule::molecule(int numatm, int numres, int numchn)
{
	clean();
	return;
}

molecule::molecule(ifstream &pdb_in)
{
	clean();
	readpdb(pdb_in);
	return;
}

molecule::molecule(const molecule &obj)
{
	int i;

	strcpy(name,obj.name);	
	num_chn=obj.num_chn;	
	num_res=obj.num_res;	
	num_atm=obj.num_atm;	
	xcenter=obj.xcenter; xmin=obj.xmin; xmax=obj.xmax;	
	ycenter=obj.ycenter; ymin=obj.ymin; ymax=obj.ymax;
	zcenter=obj.zcenter; zmin=obj.zmin; zmax=obj.zmax;
	rmax=obj.rmax;
	chg=obj.chg;		

	mol = new atom [num_atm];
	chn_end = new int [num_chn];
	res_end = new int [num_res];
	res_chg = new double [num_res];
	res_ter = new int [num_res];

	for(i=0;i<num_atm;i++)mol[i]=obj.mol[i];
	for(i=0;i<num_chn;i++)chn_end[i]=obj.chn_end[i];
	for(i=0;i<num_res;i++)res_end[i]=obj.res_end[i];
	for(i=0;i<num_res;i++)res_chg[i]=obj.res_chg[i];
	for(i=0;i<num_res;i++)res_ter[i]=obj.res_ter[i];

	return;
}

//INPUT FUNCTIONS
void molecule::clean()
{
	strcpy(name,"Unknown Molecule");
	mol=NULL;
	chn_end=NULL;
	res_end=NULL;
	res_chg=NULL;
	res_ter=NULL;
	num_chn=0;	
	num_res=0;
	num_atm=0;
	xcenter=xmin=xmax=0;
	ycenter=ymin=ymax=0;
	zcenter=zmin=zmax=0;
	rmax=0;
	chg=0;	

	return;
}

void molecule::readpdb(ifstream &pdb_in)
{
	char line[401];
	bool end_mol,end_chn; //Used to detect the end of molecule and chains
	atom atm,last_atm;

	if((mol)||(chn_end)||(res_end)||(num_atm)||(res_chg)||(res_ter))
	{
		cerr<<"WARNING - OVERWRITING MOLECULE: "<<name<<endl;
	}
	clean();

	//Molecule name
	pdb_in.getline(line,400);
	if(strncmp(line,"REMARK",6)==0)
	{
		strcpy(name,line+6);
		pdb_in.getline(line,400);
	}
	else strcpy(name,"Unknown Molecule");
	pdb_in.clear();
	pdb_in.seekg(0,ios::beg);	

	//Preliminary reading
	//It is used to detect the number of atoms, residues and chains
	end_mol=end_chn=false;
	num_atm=num_chn=num_res=0;
	while((strncmp(line,"ATOM",4)!=0)&&(!pdb_in.eof()))
	{
		pdb_in.getline(line,400);
	}
	if(pdb_in.eof())
	{
		cerr<<"WARNING - NO ATOM IN MOLECULE: "<<name<<endl;
		return;
	}
	atm.readpdb(line);
	last_atm=atm;
	num_atm++;
	pdb_in.getline(line,400);
	while( (!pdb_in.eof()) && (!end_mol) )
	{
		if(strncmp(line,"ATOM",4)==0)	//Atom line
		{
			atm.readpdb(line);	
			if(	//New chain
			    (atm.chn!=last_atm.chn)&&
			    (!end_chn)
			  )
			{
				num_chn++;
				end_chn=true;
			}
			if(	//New residue 
			    (strcmp(atm.res,last_atm.res)!=0)||
			    (atm.res_num!=last_atm.res_num)||
			    (end_chn)
			  )
			{
				num_res++;
				end_chn=false;
			}
			last_atm=atm;
			num_atm++;
		}
		if(strncmp(line,"TER",3)==0)	//New chain
		{
			num_chn++;
			end_chn=true;
		}
		if(strncmp(line,"END",3)==0)//Molecule end
		{
			num_res++;
			if(!end_chn)num_chn++;
			end_mol=true;		
		}
		pdb_in.getline(line,400);//Read new line
	}
	if((pdb_in.eof())&&(!end_mol))
	{
		if(!end_chn)num_res++;
		num_chn++;
	}

	//Allocate memory
	mol = new atom [num_atm];
	chn_end = new int [num_chn];
	res_end = new int [num_res];
	res_chg = new double [num_res];
	res_ter = new int [num_res];

	//Actuall pdb reading
	//Data structure are defined in this second reading of the file.pdb
	pdb_in.clear();
	pdb_in.seekg(0,ios::beg);	
	end_mol=end_chn=false;
	num_atm=num_chn=num_res=0;
	while((strncmp(line,"ATOM",4)!=0)&&(!pdb_in.eof()))
	{
		pdb_in.getline(line,400);
	}
	atm.readpdb(line);
	mol[num_atm]=atm;
	last_atm=atm;
	num_atm++;
	pdb_in.getline(line,400);
	while( (!pdb_in.eof()) && (!end_mol) )
	{
		if(strncmp(line,"ATOM",4)==0)	//Atom line
		{
			atm.readpdb(line);	
			if(	//New chain
			    (atm.chn!=last_atm.chn)&&
			    (!end_chn)
			  )
			{
				//cout<<atm.chn<<" "<<last_atm.chn<<endl;
				chn_end[num_chn]=num_atm;
				num_chn++;
				end_chn=true;
			}
			if(	//New residue 
			    (strcmp(atm.res,last_atm.res)!=0)||
			    (atm.res_num!=last_atm.res_num)||
			    (end_chn)
			  )
			{
				res_end[num_res]=num_atm;
				num_res++;
				end_chn=false;
			}
			last_atm=atm;
			mol[num_atm]=atm;
			num_atm++;
		}
		if(strncmp(line,"TER",3)==0)	//New chain
		{
			chn_end[num_chn]=num_atm;
			num_chn++;
			end_chn=true;
		}
		if(strncmp(line,"END",3)==0)//Molecule end
		{
			res_end[num_res]=num_atm;
			num_res++;
			if(!end_chn)
			{
				chn_end[num_chn]=num_atm;
				num_chn++;
			}
			end_mol=true;		
		}
		//cout<<line<<" "<<num_atm<<" "<<num_res<<" "<<num_chn<<endl;
		pdb_in.getline(line,400);//Read new line
	}
	if((pdb_in.eof())&&(!end_mol))
	{
		if(!end_chn)
		{
			chn_end[num_chn]=num_atm;
			num_chn++;
		}
		res_end[num_res]=num_atm;
		num_res++;
	}

	fix();

	return;
}

void molecule::copy(const molecule &obj)
{
	int i;

	clean();
	strcpy(name,obj.name);	
	num_chn=obj.num_chn;	
	num_res=obj.num_res;	
	num_atm=obj.num_atm;	
	xcenter=obj.xcenter; xmin=obj.xmin; xmax=obj.xmax;	
	ycenter=obj.ycenter; ymin=obj.ymin; ymax=obj.ymax;
	zcenter=obj.zcenter; zmin=obj.zmin; zmax=obj.zmax;
	rmax=obj.rmax;
	chg=obj.chg;		

	mol = new atom [num_atm];
	chn_end = new int [num_chn];
	res_end = new int [num_res];
	res_chg = new double [num_res];
	res_ter = new int [num_res];

	for(i=0;i<num_atm;i++)mol[i]=obj.mol[i];
	for(i=0;i<num_chn;i++)chn_end[i]=obj.chn_end[i];
	for(i=0;i<num_res;i++)res_end[i]=obj.res_end[i];
	for(i=0;i<num_res;i++)res_chg[i]=obj.res_chg[i];
	for(i=0;i<num_res;i++)res_ter[i]=obj.res_ter[i];

	fix();

	return;
}

//SELECTION FUNCTION
atom molecule::find(char *res_name, char *atm_name)
{
	int ind_atm;
	atom dummy;

	for(ind_atm=0;ind_atm<num_atm;ind_atm++)
	{
		if((mol[ind_atm].cmpresname(res_name))   
			&&(mol[ind_atm].cmpatomname(atm_name)))
			return mol[ind_atm];
	}
	return dummy;
}
int molecule::ister(int ind_atm)
{
	int ind_atm_tmp,ind_res;

	ind_res=0;
	for(ind_atm_tmp=0;ind_atm_tmp<num_atm;ind_atm_tmp++)
	{
		if(ind_atm_tmp==res_end[ind_res]){
			ind_res++;
		}
		if(ind_atm_tmp==ind_atm)return res_ter[ind_res];
	}
	return 2;
}

//TEXT EDITING
void molecule::setchain(char chn){
	int ind_atm;
	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		mol[ind_atm].chn=chn;
	}
	return;
}

//GEOMETRICAL FUNCTIONS
void molecule::fix(){
	int ind_atm,ind_res,ind_chn;

	xcenter=ycenter=zcenter=0.0;
	xmax=ymax=zmax=rmax=-INF;
	xmin=ymin=zmin=INF;
	chg=0;
	for(ind_atm=0;ind_atm<num_atm;ind_atm++)
	{
		xcenter+=mol[ind_atm].x;
		ycenter+=mol[ind_atm].y;
		zcenter+=mol[ind_atm].z;
		if(mol[ind_atm].x>xmax)xmax=mol[ind_atm].x;
		if(mol[ind_atm].y>ymax)ymax=mol[ind_atm].y;
		if(mol[ind_atm].z>zmax)zmax=mol[ind_atm].z;
		if(mol[ind_atm].x<xmin)xmin=mol[ind_atm].x;
		if(mol[ind_atm].y<ymin)ymin=mol[ind_atm].y;
		if(mol[ind_atm].z<zmin)zmin=mol[ind_atm].z;
		crt2cln(mol[ind_atm].x,mol[ind_atm].y,mol[ind_atm].r,mol[ind_atm].t);	
		if(mol[ind_atm].r>rmax)rmax=mol[ind_atm].r;
		chg+=mol[ind_atm].chg;
	}
	xcenter=xcenter/num_atm;
	ycenter=ycenter/num_atm;
	zcenter=zcenter/num_atm;
	for(ind_res=0;ind_res<num_res;ind_res++){
		res_chg[ind_res]=0;
		res_ter[ind_res]=0;
	}
	ind_res=0;
	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		if(ind_atm==res_end[ind_res])ind_res++;
		res_chg[ind_res]+=(mol[ind_atm].chg);
	}
	ind_res=ind_chn=0;
	res_ter[0]=-1;
	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		//cout<<ind_atm<<" "<<res_end[ind_res]<<" "<<chn_end[ind_chn]<<endl;
		if(ind_atm==chn_end[ind_chn]){
			res_ter[ind_res]=1;
			ind_chn++;
			if(ind_chn<num_chn){
				res_ter[ind_res+1]=-1;
			}
		}
		if(ind_atm==res_end[ind_res]){
			ind_res++;
		}
	}
	res_ter[num_res-1]=1;

	return;
}

void molecule::translate(double *T)
{
	int ind_atm;

	for(ind_atm=0;ind_atm<num_atm;ind_atm++)
	{
		mol[ind_atm].x+=T[0];
		mol[ind_atm].y+=T[1];
		mol[ind_atm].z+=T[2];
	}
	fix();

	return;
}
void molecule::translate(double xt, double yt, double zt)
{
	double T[3];

	T[0]=xt;
	T[1]=yt;
	T[2]=zt;
	translate(T);

	return;
}


void molecule::rotate(double *R)
{
	int ind_atm;
	double xtmp,ytmp,ztmp;

	for(ind_atm=0;ind_atm<num_atm;ind_atm++)
	{
		xtmp=R[0]*mol[ind_atm].x+R[1]*mol[ind_atm].y+R[2]*mol[ind_atm].z;
		ytmp=R[3]*mol[ind_atm].x+R[4]*mol[ind_atm].y+R[5]*mol[ind_atm].z;
		ztmp=R[6]*mol[ind_atm].x+R[7]*mol[ind_atm].y+R[8]*mol[ind_atm].z;
		mol[ind_atm].x=xtmp;
		mol[ind_atm].y=ytmp;
		mol[ind_atm].z=ztmp;
	}
	fix();

	return;
}
void molecule::rotate_x(double teta)
{
	double R[9];

	R[0]=1;			R[1]=0;			R[2]=0;
	R[3]=0; 		R[4]=cos(teta);		R[5]=sin(teta);
	R[6]=0; 		R[7]=-sin(teta);	R[8]=cos(teta);
	rotate(R);
}
void molecule::rotate_y(double teta)
{
	double R[9];

	R[0]=cos(teta);		R[1]=0;			R[2]=-sin(teta);
	R[3]=0; 		R[4]=1;			R[5]=0;
	R[6]=sin(teta); 	R[7]=0;			R[8]=cos(teta);
	rotate(R);
}
void molecule::rotate_z(double teta){
	double R[9];

	R[0]=cos(teta);		R[1]=sin(teta);		R[2]=0;
	R[3]=-sin(teta); 	R[4]=cos(teta);		R[5]=0;
	R[6]=0; 		R[7]=0;			R[8]=1;
	rotate(R);
}
void molecule::rotate_axis(double x, double y, double z, double teta){
	double norm_axis,R[9];
	double xp_x,xp_y,xp_z,yp_x,yp_y,yp_z,norm;

	norm_axis=sqrt(x*x+y*y+z*z);
	if((norm_axis>1+1e-10)||(norm_axis<1-1e-10)){
		cerr<<"WARNING: Rotation axis is not a versor norm = "<<norm_axis<<endl;
	}
	xp_x=0.0;
	xp_y=z;
	xp_z=-y;
	norm=sqrt(xp_x*xp_x+xp_y*xp_y+xp_z*xp_z);
	xp_x=xp_x/norm;
	xp_y=xp_y/norm;
	xp_z=xp_z/norm;
	yp_x=y*xp_z-z*xp_y;
	yp_y=-x*xp_z+z*xp_x;
	yp_z=x*xp_y-y*xp_x;
	norm=sqrt(yp_x*yp_x+yp_y*yp_y+yp_z*yp_z);
	yp_x=yp_x/norm;
	yp_y=yp_y/norm;
	yp_z=yp_z/norm;
	R[0]=xp_x;	R[1]=xp_y;		R[2]=xp_z;
	R[3]=yp_x; 	R[4]=yp_y;		R[5]=yp_z;
	R[6]=x; 	R[7]=y;			R[8]=z;
	rotate(R);
	rotate_z(teta);
	R[0]=xp_x;	R[1]=yp_x;		R[2]=x;
	R[3]=xp_y; 	R[4]=yp_y;		R[5]=y;
	R[6]=xp_z; 	R[7]=yp_z;		R[8]=z;
	rotate(R);
	return;
}

double molecule::rmsd(const molecule &A){
	double rmsd;
	int ind_atm;

	rmsd=0.0;
	if(cmp(A)!=1){
		//cerr<<"WARNING: Molecules are not identical (cmp = "<<cmp(A)<<")"<<endl;
		if ((cmp(A)==1)||(cmp(A)==2)||(cmp(A)==3))exit(1);
	}
	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		rmsd+=( (mol[ind_atm].x-A.mol[ind_atm].x)*(mol[ind_atm].x-A.mol[ind_atm].x) +
			(mol[ind_atm].y-A.mol[ind_atm].y)*(mol[ind_atm].y-A.mol[ind_atm].y) +
			(mol[ind_atm].z-A.mol[ind_atm].z)*(mol[ind_atm].z-A.mol[ind_atm].z) );
	}
	rmsd=sqrt(rmsd/num_atm);
	return rmsd;
}

int molecule::cmp(const molecule &A){
	int ind;
	
	//cout<<num_atm<<" "<<A.num_atm<<endl;
	//cout<<num_res<<" "<<A.num_res<<endl;
	//cout<<num_chn<<" "<<A.num_chn<<endl;
	if((num_atm!=A.num_atm)||	//Same number of atoms, residues, chains
	  (num_res!=A.num_res)||
	  (num_chn!=A.num_chn)){
		return 1;
	}
	for(ind=0;ind<num_res;ind++){	//Same number of atom in each residue
		if(res_end[ind]!=A.res_end[ind])return 2;
	}
	for(ind=0;ind<num_chn;ind++){	//Same number of atom in each chain
		if(chn_end[ind]!=A.chn_end[ind])return 3;
	}
	for(ind=0;ind<num_atm;ind++){	//Same atoms
		if(! ( (mol[ind]).cmp(A.mol[ind]) ) )return 4;
	}

	return 0;
}

//MOLECULE PROPRIETIES
double molecule::mass() {
	double mass=0.0;
	int ind_atm;
	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		mass+=mol[ind_atm].mass();
	}
	return mass;
}
double molecule::centermass(double *bar) {
	double mass=0.0;
	int ind_atm;
	bar[0]=bar[1]=bar[2]=0.0;
	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		bar[0]+=(mol[ind_atm].mass()*mol[ind_atm].x);
		bar[1]+=(mol[ind_atm].mass()*mol[ind_atm].y);
		bar[2]+=(mol[ind_atm].mass()*mol[ind_atm].z);
		mass+=mol[ind_atm].mass();
	}
	bar[0]=bar[0]/mass;
	bar[1]=bar[1]/mass;
	bar[2]=bar[2]/mass;
	return mass;
}
int molecule::inertia(double *I, double x0, double y0, double z0){
	int ind_atm;
	I[0]=I[1]=I[2]=I[3]=I[4]=I[5]=I[6]=I[7]=I[8]=0.0;

	for(ind_atm=0;ind_atm<num_atm;ind_atm++){
		I[0]+=(mol[ind_atm].mass())* (
			(mol[ind_atm].y-y0)*(mol[ind_atm].y-y0) +
			(mol[ind_atm].z-z0)*(mol[ind_atm].z-z0) );
		I[4]+=(mol[ind_atm].mass())* (
			(mol[ind_atm].x-x0)*(mol[ind_atm].x-x0) +
			(mol[ind_atm].z-z0)*(mol[ind_atm].z-z0) );
		I[8]+=(mol[ind_atm].mass())* (
			(mol[ind_atm].x-x0)*(mol[ind_atm].x-x0) +
			(mol[ind_atm].y-y0)*(mol[ind_atm].y-y0) );
		I[1]-=(mol[ind_atm].mass())* (
			(mol[ind_atm].x-x0)*(mol[ind_atm].y-y0) );
		I[2]-=(mol[ind_atm].mass())* (
			(mol[ind_atm].x-x0)*(mol[ind_atm].z-z0) );
		I[5]-=(mol[ind_atm].mass())* (
			(mol[ind_atm].y-y0)*(mol[ind_atm].z-z0) );
	}
	I[3]=I[1];
	I[6]=I[2];
	I[7]=I[5];
	return 0;
}

//FRIEND OPERATORS 
atom& molecule::operator[](int ind_atm)
{
	return mol[ind_atm];
}
ostream &operator<<(ostream &stream, molecule &obj)
{
	int ind_atm,ind_chn,str_atm,end_atm,ind_res;

	obj.fix();
	stream<<"REMARK "<<obj.name<<"\n";
	stream<<"REMARK #CHAINS = "<<obj.num_chn<<" #RESIDUES = "<<obj.num_res
		<<" #ATOMS = "<<obj.num_atm<<"\n";
	stream<<"REMARK MOLECULE CENTER x = "<<obj.xcenter<<" y = "<<obj.ycenter<<" z = "<<obj.zcenter<<"\n";
	stream<<"REMARK MOLECULE LIMITS xmin = "<<obj.xmin<<" ymin = "<<obj.ymin<<" zmin = "<<obj.zmin<<"\n";
	stream<<"REMARK MOLECULE LIMITS xmax = "<<obj.xmax<<" ymax = "<<obj.ymax<<" zmmax= "<<obj.zmax<<"\n";
	stream<<"REMARK MOLECULE LIMITS rmax = "<<obj.rmax<<"\n";
	stream<<"REMARK MOLECULE CHARGE = "<<obj.chg<<"\n";
	for(ind_chn=0;ind_chn<obj.num_chn;ind_chn++)
	{
		if(ind_chn==0)str_atm=0;
		else str_atm=(obj.chn_end[ind_chn-1]);
		end_atm=obj.chn_end[ind_chn]-1;
		stream<<"REMARK CHAIN "<<obj.mol[str_atm].chn<<" ("<<ind_chn<<")"
			<<" STARTS AT ATOM "<<obj.mol[str_atm].num<<" ("<<str_atm<<")"
			<<" ENDS AT ATOM "<<obj.mol[end_atm].num<<" ("<<end_atm<<")\n";
	}
	for(ind_chn=0;ind_chn<obj.num_chn;ind_chn++)
	{
		if(ind_chn==0){
			str_atm=0;
			ind_res=0;
		}
		else{
			str_atm=(obj.chn_end[ind_chn-1]);
		}
		end_atm=obj.chn_end[ind_chn];
		for(ind_atm=str_atm;ind_atm<end_atm;ind_atm++)
		{
			if(ind_atm==obj.res_end[ind_res]){
				stream<<"REMARK RESIDUE "<<obj.mol[ind_atm-1].res<<obj.mol[ind_atm-1].res_num
					<<" CHARGE = "<<obj.res_chg[ind_res];
				switch (obj.res_ter[ind_res]){
					case -1:
						stream<<" N-Terminal"<<endl;
						break;
					case 0:
						stream<<" Residue in between"<<endl;
						break;
					case 1:
						stream<<" C-Terminal"<<endl;
						break;
					default:
						cerr<<"ERROR WRITING THE MOLECULE"<<endl;
						exit(1);
						break;
				}
				ind_res++;
			}
			stream<<obj.mol[ind_atm]<<"\n";
		}
		stream<<"REMARK RESIDUE "<<obj.mol[ind_atm-1].res<<obj.mol[ind_atm-1].res_num
			<<" CHARGE = "<<obj.res_chg[ind_res]<<endl;
		stream<<"TER\n";
	}
	stream<<"END\n";

	return stream;
}

//Destructor
molecule::~molecule()
{
#ifdef DEBUG	
	cerr<<"DISTRUTTORE MOLECOLA\n";
#endif
	delete[] mol;
	delete[] chn_end;
	delete[] res_end;

	return;
}

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

atom::atom(){
	num=0;			//Serial Number
	res_num=0;		//Numero residuo
	strcpy(name,"NULL");	//Name	
	strcpy(res,"NULL");	//Residue name
	chn=' ';		//Identificativo di catena
	x=0; y=0; z=0;		//Coordinate cartesiane
	r=0; t=0;		//Coordinate cilindriche
	rad=0;			//Raggio 
	chg=0;			//Carica
	return;
}
atom::atom(int NUM,int RES_NUM,
		char *NAME, char *RES,
		char CHN,
		double X,double Y,double Z,
		double R,double T,
		double RAD,double CHG){
	num=NUM;		//Serial Number
	res_num=RES_NUM;	//Numero residuo
	strncpy(name,NAME,4);	//Name	
	name[5]='\n';
	strncpy(res,RES,4);	//Residue name
	res[5]='\n';
	chn=CHN;		//Identificativo di catena
	x=X; y=Y; z=Z;		//Coordinate cartesiane
	r=R; t=T;		//Coordinate cilindriche
	rad=RAD;		//Raggio 
	chg=CHG;		//Carica
	return;
}

int atom::clear()
{
	num=0;			//Serial number           
	strcpy(name,"NULL");	//Atom type               
	strcpy(res,"NULL");	//Residue name            
	res_num=0;		//Residue serial number   
	chn=' ';		//Chain identifier        
	x=0; y=0; z=0;		//Cartesian coordinates   
	r=0; t=0;		//Cilindrical coordinates 
	chg=0;			//Charge                  
	rad=0;			//Radius                  

	return 0;
}

atom atom::readpdb(char *line)
{
	char stmp[201];


	if((strncmp(line,"ATOM",4)!=0)||(strlen(line)<54))
	{
		cerr<<"ERROR THIS IS NOT AN ATOM LINE\nLINE: "<<line<<endl;
		exit(1);
	}
	else
	{
		strncpy(stmp,line+5,6);//Atom number
		stmp[6]='\0';
		num=atoi(stmp);
		strncpy(name,line+12,4);//Atom name
		name[4]='\0';
		strncpy(res,line+17,4);//Res name
		res[4]='\0';
		chn=line[21];	//Chain
		strncpy(stmp,line+22,4);//Res number
		stmp[4]='\0';
		res_num=atoi(stmp);
		strncpy(stmp,line+30,8);//x
		stmp[8]='\0';
		x=atof(stmp);
		strncpy(stmp,line+38,8);//y
		stmp[8]='\0';
		y=atof(stmp);
		strncpy(stmp,line+46,8);//z
		stmp[8]='\0';
		z=atof(stmp);
		if(strlen(line)>59)
		{
			strncpy(stmp,line+54,6);
			stmp[6]='\0';
			rad=atof(stmp);
		}
		if(strlen(line)>65)
		{
			strncpy(stmp,line+60,6);
			stmp[6]='\0';
			chg=atof(stmp);
		}
	}
	crt2cln(x,y,r,t);	//Compute cilindrical coordinates

	return *this;
}

int cmpname(char *name,char *tmpl)
{
	char stmp1[10],stmp2[10];
	int i,j;

	j=0;
	for(i=0;i<strlen(tmpl);i++)
	{
		if(tmpl[i]!=' ')stmp1[j++]=tmpl[i];
	}
	stmp1[j]='\0';
	j=0;
	for(i=0;i<strlen(name);i++)
	{
		if(name[i]!=' ')stmp2[j++]=name[i];
	}
	stmp2[j]='\0';
	if(strcmp(stmp1,stmp2)==0)return 1;
	else return 0;
}

int atom::cmpatomname(char *atmname)
{
	return cmpname(name,atmname);
}
int atom::cmpresname(char *resname)
{
	return cmpname(res,resname);
}

int atom::isbackbone()
{
	if(
	   (cmpatomname("CA"))||
	   (cmpatomname("C"))||
	   (cmpatomname("N"))||
	   (cmpatomname("O"))||
	   (cmpatomname("H"))||
	   (cmpatomname("HA"))||
	   ((cmpatomname("3HD")) && (cmpresname("PRO")) ) ||
	   ((cmpatomname("2HD")) && (cmpresname("PRO")) ) ||
	   ( cmpresname("GLY") )
	   )
		return 1;
	else
		return 0;
}

int atom::cmp(atom atm)
{
	//cout<<(*this)<<endl<<atm<<endl;
	if (!(cmpatomname(atm.name)))return 0;
	if (strcmp(res,atm.res)!=0)return 0;
	if ((x<(atm.x-0.2))||(x>(atm.x+0.2)))return 0;
	if ((y<(atm.y-0.2))||(y>(atm.y+0.2)))return 0;
	if ((z<(atm.z-0.2))||(z>(atm.z+0.2)))return 0;
	if ((chg<(atm.chg-0.001))||(chg>(atm.chg+0.001)))return 0;

	return 1;
}

double atom::distance(atom atm){
	return sqrt((x-atm.x)*(x-atm.x)+(y-atm.y)*(y-atm.y)+(z-atm.z)*(z-atm.z));
}

double atom::mass(){
	int i=0;
	while(!isalpha(name[i]))i++;
	switch (name[i]){
		case 'O':return 16;
		case 'N':return 14;
		case 'C':return 12;
		case 'S':return 32;
		case 'H':return 1;
		default:
			 cerr<<"WARNING: mass unknown for atom "<<name<<endl;
			 return 0;
	}
}

//FRIEND FUNCTIONS
double distance(atom atm1, atom atm2)
{
	return sqrt(
			(atm1.x-atm2.x)*(atm1.x-atm2.x)+
			(atm1.y-atm2.y)*(atm1.y-atm2.y)+
			(atm1.z-atm2.z)*(atm1.z-atm2.z));
}

ostream &operator<<(ostream &stream, atom atm)
{
	int num,res_num;

	num=(atm.num%100000);
	res_num=(atm.res_num%10000);
	stream<<setiosflags(ios::fixed);
	stream<<setprecision(3);
	if((atm.chn>=65)&&(atm.chn<=90))
	{
		stream<<"ATOM  "<<setw(5)<<num<<" "<<setw(4)<<atm.name<<" "<<
			setw(4)<<atm.res<<setw(1)<<atm.chn<<setw(4)<<res_num<<
			"    "<<setw(8)<<atm.x<<setw(8)<<atm.y<<setw(8)<<atm.z
			<<setw(6)<<atm.rad<<setw(6)<<atm.chg;
	}
	else
	{
		stream<<"ATOM  "<<setw(5)<<num<<" "<<setw(4)<<atm.name<<" "<<
			setw(4)<<atm.res<<setw(1)<<' '<<setw(4)<<res_num<<
			"    "<<setw(8)<<atm.x<<setw(8)<<atm.y<<setw(8)<<atm.z
			<<setw(6)<<atm.rad<<setw(6)<<atm.chg;
	}
	stream<<resetiosflags(ios::fixed);

	return stream;
}

atom::~atom()
{
	return;
}

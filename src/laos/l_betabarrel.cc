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

#define MAX_LEN_SHEET 100	

int printhelp(char *prm_name);

int main(int argc,char *argv[])
{
	int ind_arg;
	char *prm_name;
	char stmp[201];

	time_t time_start,time_end;

	char fileout[201];
	char filein[201];
	ofstream file_out;
	ifstream file_in;

	int N,S;	//N = numero sheet; S = Pendenza sheet
	int len_sheet;	//Numero aminoacidi nel singolo sheet
	double a,b;	//a = Distanza CA lungo lo sheet;b = Distanza sheet 
	double CaCb;	//Distanza CA - CB
	double R,alpha;	//R = Raggio canale;alpha = angolo asse canale/sheet
	double da;	//Scostamento lungo asse sheet da introdurre affinche' la
			//congiungente tra Ca di beta-sheet differenti sia ortogonale
	double trasl;	//Variabili utilizzate per traslare correttamente i foglietti
	int nres;	//beta lungo l'asse del barrel
	bool seq;	//Sequenza aminoacidica letta da file ?
	int itrs;	//Variabile utilizzata per indirizzare correttamente il primo
			//aminoacido del secondo foglietto
	int jtrs;	//Variabile utilizzata per il corretto posizionamento delle
			//catene laterali
	int n_seqA,n_seqB;	//Numero sequenze lette da file
	typedef struct {	//Sequenze lette da file
		char res_name[4];
		int res_num;
	} seq_sheet_type;
	seq_sheet_type seq_sheetA[MAX_LEN_SHEET];
	seq_sheet_type seq_sheetB[MAX_LEN_SHEET];

	int iz,it,n;
	char line[201];
	double t;
	int i,n_read;		
	atom atm_tmp;
	double CaC,CaN,CN,CO;	//Lunghezze legami
	double d1,d2;		//Distanza CaC e CaN lungo l'asse CaCa
	double Caxs,Naxs;	//Distanze di N e C dall'asse CaCa
	double phi,psi,x,y,z;

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
	strcpy(fileout,prm_name);
	strcat(fileout,".out.pdb");
	strcpy(filein,"NULL");

//Parameter Initialization
	N=0;
	S=0;
	len_sheet=0;
	seq=false;
	//cout<<"Distance between CA along the sheet [3.3 A] = ";
	//cin>>a;
	a=3.3;
	//cout<<"Distance between sheets [4.4 A] = ";
	//cin>>b;
	b=4.4;
	//cout<<"Distance between CA and CB [1.5 A] = ";
	//cin>>CaCb;
	CaCb=1.5;
	//cout<<"Distance Ca-C along the CaCa axes = ";
	//cin>>d1;
	//cout<<"Distance Ca-N along the CaCa axes = ";
	//cin>>d2;
	d1=d2=1.35;
	//Inizializzazione parametri
	//Lunghezze legami da parm99.dat; AMBER FORCE FIELD LEAPRC.FF03
	CaC=1.522;	//Lungezza legame CA-C
	CaN=1.449;	//Lunghezza legame CA-N
	CO=1.229;	//Lunghezza legame C=0
	Caxs=sqrt((CaC*CaC) - (d1*d1));
	Naxs=sqrt((CaN*CaN) - (d2*d2));
	CN=sqrt((Caxs+Naxs)*(Caxs+Naxs)+(a-d1-d2)*(a-d1-d2));
	file_out<<"REMARK Caxs = "<<Caxs<<" Naxs = "<<Naxs<<"\n";
	file_out<<"REMARK C-N bond length = "<<CN<<"\n";

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
	      	  else if(strcmp(argv[ind_arg]+1,"i")==0)//Input File - Aminoacid Sequences
	      	  {
			  seq=true;
	      		  strcpy(filein,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"o")==0)//Output File
	      	  {
	      		  strcpy(fileout,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"N")==0)//Number of Beta-Sheet
	      	  {
	      		  strcpy(stmp,argv[ind_arg+1]);
			  N=atoi(stmp);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"S")==0)//Sheer Number
	      	  {
	      		  strcpy(stmp,argv[ind_arg+1]);
			  S=atoi(stmp);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"len_sheet")==0)//Sheer Number
	      	  {
	      		  strcpy(stmp,argv[ind_arg+1]);
			  len_sheet=atoi(stmp);
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
	if(seq){
  		file_in.open(filein);
		if(!file_in)
		{
			cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filein<<endl;
			return 1;
		}
	}
  	file_out.open(fileout);
	if(!file_out)
	{
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}

//Parameter control
	if ((N<=0)||((N%2)!=0)||(((len_sheet)<=0)&&(!seq)))
	{
		cerr<<"N = "<<N<<" a positive even number is required\n";
		cerr<<"len_sheet = "<<len_sheet<<" a positive number is required\n";
		cerr<<"USAGE: laos [-i NULL] -N N -S S [-len_sheet 0] [-o l_betabarrel.out.pdb]\n";
		return 1;
	}
	if ((S<N)||(S>(2*N)))
	{
		cerr<<"S = "<<S<<" S:  "<<N<<" <= S <= "<<2*N<<"\n";
		cerr<<"USAGE: laos [-i NULL] -N N -S S [-len_sheet 0] [-o l_betabarrel.out.pdb]\n";
		return 1;
	}

//Derivate parameters
	if(seq){
		n_seqA=n_seqB=i=0;
		while(!file_in.eof()){
			file_in.getline(line,100);
			if((line[0]!='#')&&(strlen(line)!=0)){
				if(strncmp(line,"TER",3)==0){
					n_seqA=i;
					i=0;
				}else{
					if(n_seqA==0){
						n_read=sscanf(line,"%s %i",
							(seq_sheetA[i]).res_name,&((seq_sheetA[i]).res_num));
					}else{
						n_read=sscanf(line,"%s %i",
							(seq_sheetB[i]).res_name,&((seq_sheetB[i]).res_num));
					}
					if(n_read!=2){
						cerr<<"ERROR READING FILE "<<filein<<endl;
						return 1;
					}
					i++;
				}
			}
		}
		n_seqB=i;

		//START DEBUG
		//cout<<n_seqA<<" "<<n_seqB<<endl;
		//for(i=0;i<n_seqA;i++){
		//		cout<<seq_sheetA[i].res_name<<" "<<seq_sheetA[i].res_num<<endl;
		//}
		//for(i=0;i<n_seqB;i++){
		//		cout<<seq_sheetB[i].res_name<<" "<<seq_sheetB[i].res_num<<endl;
		//}
		//END DEBUG

		if((n_seqA!=n_seqB)||(n_seqA<=0)){
			cerr<<"ERROR READING FILE "<<filein<<endl;
			return 1;
		}
		if((len_sheet!=0)&&(len_sheet!=n_seqA))cerr<<"WARNING -len_sheet IGNORED"<<endl;
		len_sheet=n_seqA;
	}

//Derivate parameters
	alpha=atan((S*a)/(N*b));
	R=(b)/(2*sin(PI/N)*cos(alpha));
	file_out<<"REMARK Radius of the beta barrel = "<<R<<"\n";
	//da=b*tan(alpha)-a;
	da=b*tan(alpha);
	file_out<<"REMARK Beta-sheet translation = "<<da<<"\n";
	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
	if((len_sheet%2)==0)itrs=1;
	else itrs=0;
//START-UP	END
//---------------------------------------------------------
 

//---------------------------------------------------------
//MAKING BETA-BARREL 	BEGIN
	atm_tmp.clear();
	atm_tmp.num=1;
	atm_tmp.res_num=1;
	for(it=(((int)(N/2))-1);it>=0;it--)	{
	//I due cicli servono per produrre coppie di foglietti antiparalleli con
	//la numerazione corretta
		phi=0.0;
		psi=0.0;
		n=0;
		trasl=da*it*2;	//Traslazione lungo la direzione Ca-Ca tra i foglietti beta
		nres=0;		//Numero residui traslati
		while(trasl>(a-0.01)){
			trasl-=a;
			nres++;
		}
		cout<<"Traslazione reale = "<<da*(it)*2;
		cout<<" Indice = "<<it*2<<" Traslazione = "<<trasl<<" Numero residui = "<<nres<<endl;
		//if((nres%2)!=0)n++;
		//Primo foglietto
		for(iz=(len_sheet-1);iz>=0;iz--)
		{
			//Carbonio Alpha
		      	t=(it*2)*((2*PI)/N) + (iz*(a*sin(alpha)/R)) + (-trasl*sin(alpha)/R);
		      	atm_tmp.x=R*cos(t);
		      	atm_tmp.y=R*sin(t);
		       	atm_tmp.z=((iz*a)-trasl)*cos(alpha);
			strcpy(atm_tmp.name,"CA");
			if(seq){
				//cout<<n<<"["<<n_seqA<<"]"<<endl;
				strcpy(atm_tmp.res,seq_sheetA[n].res_name);
				atm_tmp.res_num=seq_sheetA[n].res_num;
				atm_tmp.chn=('A'+it);

			}else strcpy(atm_tmp.res,"ALA");
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;

			//Gruppo aminico
			strcpy(atm_tmp.name,"N");
			phi+=((+180.0/360.0)*2.0*PI);
			x=pow(-1.0,(nres+1)%2)*Naxs*sin(phi);
			y=pow(-1.0,(nres+1)%2)*Naxs*cos(phi);
			z=+d2;
			atm_tmp.x+=(cos(t)*x-sin(t)*cos(alpha)*y-sin(t)*sin(alpha)*z);
			atm_tmp.y+=(sin(t)*x+cos(t)*cos(alpha)*y+cos(t)*sin(alpha)*z);
			atm_tmp.z+=(-sin(alpha)*y+cos(alpha)*z);
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;

			//Catena Laterale
			if(strncmp(atm_tmp.res,"GLY",3)!=0){
				strcpy(atm_tmp.name,"CB");
		      		atm_tmp.x=(R-CaCb*pow(-1.0,(n%2)))*cos(t);
		      		atm_tmp.y=(R-CaCb*pow(-1.0,(n%2)))*sin(t);
		       		atm_tmp.z=((iz*a)-trasl)*cos(alpha);
		      		file_out<<atm_tmp<<"\n";
		      		atm_tmp.num++;
			}
			n++;
			
			//Gruppo corbosillico
			strcpy(atm_tmp.name,"C");
		      	atm_tmp.x=R*cos(t);
		      	atm_tmp.y=R*sin(t);
		       	atm_tmp.z=((iz*a)-trasl)*cos(alpha);
			psi+=((-180.0/360.0)*2.0*PI);
			x=pow(-1.0,(nres+1)%2)*Caxs*sin(psi);
			y=pow(-1.0,(nres+1)%2)*Caxs*cos(psi);
			z=-d1;
			atm_tmp.x+=(cos(t)*x-sin(t)*cos(alpha)*y-sin(t)*sin(alpha)*z);
			atm_tmp.y+=(sin(t)*x+cos(t)*cos(alpha)*y+cos(t)*sin(alpha)*z);
			atm_tmp.z+=(-sin(alpha)*y+cos(alpha)*z);
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;
			strcpy(atm_tmp.name,"O");
		      	atm_tmp.x=R*cos(t);
		      	atm_tmp.y=R*sin(t);
		        atm_tmp.z=((iz*a)-trasl)*cos(alpha);
			x=pow(-1.0,(nres+1)%2)*(Caxs+CO)*sin(psi);
			y=pow(-1.0,(nres+1)%2)*(Caxs+CO)*cos(psi);
			z=-d1;
			atm_tmp.x+=(cos(t)*x-sin(t)*cos(alpha)*y-sin(t)*sin(alpha)*z);
			atm_tmp.y+=(sin(t)*x+cos(t)*cos(alpha)*y+cos(t)*sin(alpha)*z);
			atm_tmp.z+=(-sin(alpha)*y+cos(alpha)*z);
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;

			if(!seq)atm_tmp.res_num++;
		}
		file_out<<"TER\n";

		phi=0.0;
		psi=0.0;
		n=0;
		trasl=da*((it*2)+1);	//Traslazione lungo la direzione Ca-Ca tra i foglietti beta
		nres=0;			//Numero residui traslati
		while(trasl>(a-0.01)){
			trasl-=a;
			nres++;
		}
		cout<<"Traslazione reale = "<<da*((it*2)+1);
		cout<<" Indice = "<<(it*2)+1<<" Traslazione = "<<trasl<<" Numero residui = "<<nres<<endl;
		//if((nres%2)!=0)n++;
		if((nres%2)!=0)jtrs=1;
		else jtrs=0;
		//Secondo foglietto
		for(iz=0;iz<len_sheet;iz++)
		{
			//Carbonio Alpha
		      	t=(1+(it*2))*((2*PI)/N) + (iz*(a*sin(alpha)/R)) + (-trasl*sin(alpha)/R);
		      	atm_tmp.x=R*cos(t);
		      	atm_tmp.y=R*sin(t);
		       	atm_tmp.z=((iz*a)-trasl)*cos(alpha);
			strcpy(atm_tmp.name,"CA");
			if(seq){
				//cout<<n<<"["<<n_seqA<<"]"<<endl;
				strcpy(atm_tmp.res,seq_sheetB[n].res_name);
				atm_tmp.res_num=seq_sheetB[n].res_num;

			}else strcpy(atm_tmp.res,"ALA");
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;
		      
			//Gruppo aminico
			strcpy(atm_tmp.name,"N");
			phi+=((+180.0/360.0)*2.0*PI);
			x=pow(-1.0,(nres+itrs)%2)*Naxs*sin(phi);
			y=pow(-1.0,(nres+itrs)%2)*Naxs*cos(phi);
			z=-d2;
			atm_tmp.x+=(cos(t)*x-sin(t)*cos(alpha)*y-sin(t)*sin(alpha)*z);
			atm_tmp.y+=(sin(t)*x+cos(t)*cos(alpha)*y+cos(t)*sin(alpha)*z);
			atm_tmp.z+=(-sin(alpha)*y+cos(alpha)*z);
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;

			//Catena Laterale
			if(strncmp(atm_tmp.res,"GLY",3)!=0){
				strcpy(atm_tmp.name,"CB");
		      		atm_tmp.x=(R-CaCb*pow(-1.0,(n+jtrs%2)))*cos(t);
		      		atm_tmp.y=(R-CaCb*pow(-1.0,(n+jtrs%2)))*sin(t);
		       		atm_tmp.z=((iz*a)-trasl)*cos(alpha);
		      		file_out<<atm_tmp<<"\n";
		      		atm_tmp.num++;
			}
			n++;
			
			//Gruppo corbosillico
			strcpy(atm_tmp.name,"C");
		      	atm_tmp.x=R*cos(t);
		      	atm_tmp.y=R*sin(t);
		       	atm_tmp.z=((iz*a)-trasl)*cos(alpha);
			psi+=((-180.0/360.0)*2.0*PI);
			x=pow(-1.0,(nres+itrs)%2)*Caxs*sin(psi);
			y=pow(-1.0,(nres+itrs)%2)*Caxs*cos(psi);
			z=+d1;
			atm_tmp.x+=(cos(t)*x-sin(t)*cos(alpha)*y-sin(t)*sin(alpha)*z);
			atm_tmp.y+=(sin(t)*x+cos(t)*cos(alpha)*y+cos(t)*sin(alpha)*z);
			atm_tmp.z+=(-sin(alpha)*y+cos(alpha)*z);
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;
			strcpy(atm_tmp.name,"O");
		      	atm_tmp.x=R*cos(t);
		      	atm_tmp.y=R*sin(t);
		       	atm_tmp.z=((iz*a)-trasl)*cos(alpha);
			x=pow(-1.0,(nres+itrs)%2)*(Caxs+CO)*sin(psi);
			y=pow(-1.0,(nres)+itrs%2)*(Caxs+CO)*cos(psi);
			z=+d1;
			atm_tmp.x+=(cos(t)*x-sin(t)*cos(alpha)*y-sin(t)*sin(alpha)*z);
			atm_tmp.y+=(sin(t)*x+cos(t)*cos(alpha)*y+cos(t)*sin(alpha)*z);
			atm_tmp.z+=(-sin(alpha)*y+cos(alpha)*z);
		      	file_out<<atm_tmp<<"\n";
		      	atm_tmp.num++;

		      	if(!seq)atm_tmp.res_num++;
		}
		file_out<<"TER\n";
	}
	file_out<<"END\n";
//MAKING BETA-BARREL 	END
//---------------------------------------------------------

	time(&time_end);
	cout<<endl<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;

	file_out.close();

  	return 0;
}

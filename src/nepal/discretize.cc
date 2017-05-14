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
#include <lib/molecule.h>
#include <src/nepal/functions.h>


#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

////////////////////////////////////////////////////////////////////////////////
//GRID
//Discretizza la molecola sulla griglia
int grid(int nr, int nt,int nz,double zmin, double zmax, double rmax,//Dimensioni griglia
		double rad_probe,//Raggio di probe utilizzato nella discretizzazione
		double zmem_min, double zmem_max,//Limiti regione di membrana
		double epsH2O, double epsMOL,//Costanti dielettriche relative
		double phimem, double phiterm, double *phibnd,
		double *Dk, double *Dcl,
		double rad_k, double rad_cl,
		double zminpaine, double zmaxpaine,
		molecule &chn,//Molecola da discretizzare
		int disc_chg,//Metodo discretizzazione
		int PNP,// = 0 Soluzione Poisson = 1 Soluzione PNP
		double kext, double kint, double clext, double clint,//Condizioni al contorno
		double *chg_mol,//Carica discretizzata sulla gliglia
		int *inside_end, int *outside_begin,//Limiti molecola
		int *is_sol,//Elementi in soluzione
		double *sum_phi,
		double *Acs_k,double  *bcs_k,double *Acs_cl,double  *bcs_cl,
		double *A_phi,double *bcs_phi,
		int &min_cell4atom,//Minimo Numero di celle usate per discretizzare un atomo
		int &max_cell4atom,//Massimo Numero di celle usate per discretizzare un atomo
		double &distancemed,//Distanza media tra centro reale atomo e centro carica discretizzata
		double &distancemax//Distanza media tra centro reale atomo e centro carica discretizzata
		)

{
	unsigned long int nrntnz; //Dimensioni griglia
	int nrnt;
	double dr,dt,dz;
	int izmem_min,izmem_max;//Limiti membrana

	double *diff_k,*diff_cl;//Coefficienti di diffusione nella griglia3D
	int flag[chn.num_atm];//Variabile booleana
	double Rpr_atm[chn.num_atm];//Proiezioni raggi atomici sulle fette
	double t,z;//Coordinate
	double R2,deltaz2,deltat,deltatmax,Rpr2,rsen2,rint_med;//Valori temporanei
	double inside,outside,distancetmp;
	int in,out;

	int ir,it,iz,ind_atm;//Indici
	unsigned long int iv;

	double PI2;

        //Inizializzazione
	nrntnz=nr*nt*nz;
	nrnt=nr*nt;
	dr=rmax/nr;
	dt=(2*PI)/nt;
	dz=(zmax-zmin)/nz;
  	izmem_min=(int)round((zmem_min-zmin)/dz);
  	izmem_min++;	//Numero primo elemento di membrana
  	izmem_max=(int)round((zmem_max-zmin)/dz); //Numero ultimo elemento di membrana
	PI2=PI/2;

	diff_k = new double [nrntnz];
	diff_cl = new double [nrntnz];

//-------------------------------------------------
//RICERCA CONFINI
//-------------------------------------------------
	//Aggiunta strato di esclusione
	for(ind_atm=0;ind_atm<chn.num_atm;ind_atm++)chn[ind_atm].rad+=rad_probe;
	
	//La discretizzazione avviene fetta per fetta.....
	for(iz=0;iz<nz;iz++)
	{
		rint_med=0.0;
		z=zmin+(iz*dz)+(dz/2);	//Ci si posiziona a meta' della fetta
		//Si selezionano solo gli atomi che intersecano la fetta
		for(ind_atm=0;ind_atm<chn.num_atm;ind_atm++)
		{
			if( 
		 		((chn[ind_atm].z)>(z-chn[ind_atm].rad)) &&
		 		((chn[ind_atm].z)<(z+chn[ind_atm].rad)) 
			  )
			{
			//Atomi che intersecano la fetta
				flag[ind_atm]=1;
				R2=(chn[ind_atm].rad)*(chn[ind_atm].rad);
				deltaz2=(chn[ind_atm].z-z)*(chn[ind_atm].z-z);
				//Proiezione del raggio sulla meta' della fetta
				if((R2-deltaz2)<=0)
				{
			cerr<<"WARNING: NEGATIVE SQRT IN THE DISCRETIZATION ALGORITHM \n";
			cerr<<"\trad atom = "<<chn[ind_atm].rad<<" rad^2 = "<<R2<<" (R2)\n";
			cerr<<"\tz slice = "<<z<<"\n";
			cerr<<"\tz atom = "<<chn[ind_atm].z<<"\n";
			cerr<<"\t(z atom - z slice)^2 = "<<deltaz2<<" (deltaz2)\n";
			cerr<<"\tR2 - deltaz2 = "<<R2-deltaz2<<"\n";
			cerr<<"\tSetting R2 = deltaz2\n";
			Rpr_atm[ind_atm]=0.0;		
				}
				else
				{
			Rpr_atm[ind_atm]=sqrt(R2-deltaz2);		
				}
			}
			else
			{
			//Atomi che non intersecano la fetta
				flag[ind_atm]=0;
			}
		}
		//.... e spicchio per spicchio
		for(it=0;it<nt;it++)
		{
			t=(dt/2)+it*dt;		//Ci si posiziona a meta' dello spicchio
			inside=rmax;
			outside=0;
			for(ind_atm=0;ind_atm<chn.num_atm;ind_atm++)
			{
				//Si considerano solo gli atomi che intersecano la fetta
				if(flag[ind_atm]==1)
				{
					deltat=dif_rad(chn[ind_atm].t,t);
					//Si considerano solo gli atomi nel semipiano corretto
					if((deltat<(PI2))&&(deltat>(-PI2)))
					{
					//deltatmax e' lo scostamento a cui la proiezione dell'atomo sulla
					//fetta e' tangente allo spicchio
						if(chn[ind_atm].r<1e-5){
							deltatmax=PI2;
						}
						else{
							double pappa=Rpr_atm[ind_atm]/chn[ind_atm].r;
							if(pappa>1.00){
								pappa=1.00;
							}
							deltatmax=asin(pappa);
							//deltatmax=asin(Rpr_atm[ind_atm]/chn[ind_atm].r);
						}
						
					//Si considerano solo gli atomi che possono intersecare la direzione 
					//dello spicchio
						if (abs(deltat)<deltatmax)
						{
					//A questo punto sono stati selezionati atomi le cui proiezioni sulla
					//meta' fetta intersecano la meta' dello spicchio
							Rpr2=(Rpr_atm[ind_atm])*(Rpr_atm[ind_atm]);
							rsen2=(chn[ind_atm].r*sin(deltat))*
								(chn[ind_atm].r*sin(deltat));
							if (rsen2>Rpr2)
							{
								cerr<<"ERROR IN THE GRID SUBROUTINE\n";
								exit(1);
							}
							distancetmp=(chn[ind_atm].r*cos(deltat))
								-sqrt(Rpr2-rsen2);
						//E' la distanza a cui la direzione dello spicchio interseca
						//l'atomo _ verso il centro del cilindro
							if(distancetmp<inside)inside=distancetmp;
							distancetmp=(chn[ind_atm].r*cos(deltat))
								+sqrt(Rpr2-rsen2);
						//E' la distanza a cui la direzione dello spicchio interseca
						//l'atomo _ verso l'esterno del cilindro
							if(distancetmp>outside)outside=distancetmp;
						}
					}
				}
			}
			
			in=(int)round(inside/dr);
			if(in>nr){
				cerr<<"ERROR in the grid subroutine: internal bound outside boundaries\n";
				exit(1);
			}
			out=(int)round(outside/dr);
			if(out>nr){
				cerr<<"WARNING in the grid subroutine: external bound outside boundaries\n";
				exit(1);
			}
			if(in<2){
				cerr<<"WARNING in the grid subroutine: The channel hole is too small\n";
				cerr<<"dr = "<<dr<<" dt = "<<dt<<" dz = "<<dz<<endl;
				cerr<<"iz = "<<iz<<" it = "<<it<<endl;
				cerr<<"z = "<<zmin+(iz*dz)+(dz/2.0)<<endl;
				cerr<<"inside = "<<inside<<endl;
				if(PNP!=0)exit(1); //Continue only if solving the Poisson equation
			}
			if((nr-out)<2){
				cerr<<"ERROR in the grid subroutine: external bound too close to boundaries\n";
				exit(1);
			}
			//Se out==0 non si e' intersecato alcuno atomo
			//---> inside == rmax ---> in == nr
			if((out==0)&&(in!=nr)){
				cerr<<"ERROR in the grid subroutine\n";
				exit(1);
			}
			//Non si e' intersecato alcun atomo ---> soluzione estesa fino al confine
			if(out==0) out=nr;
			//Introduzione membrana
			if((iz>=(izmem_min-1))&&(iz<=(izmem_max-1))){
				out=nr;
			}

//A questo punto sono stati definiti:
//	in   = Ultimo elemento interno di soluzione (numerazione parte da 1)
//		= nr se non ha intersecato niente o l'intersezione avviene proprio
//		     nell'ultimo elemento
//	out  = Ultimo elemento di molecola (numerazione parte da 1) 
//		= nr se siamo in membrana
//		= nr se non ha intersecato niente
//
			inside_end[iz*nt+it]=in;
			outside_begin[iz*nt+it]=out;
			rint_med+=(dr*(inside_end[iz*nt+it]));
		}//Fine primo ciclo in it

		rint_med=(rint_med/(double)nt)+rad_probe;//Raggio medio della fetta iz
		z=iz*dz+(dz/2.0)+zmin;
		//DEBUG
		//if (((zminpaine==0.0)&&(zmaxpaine==0.0))||((z>zminpaine)&&(z<zmaxpaine)))
		//	cerr<<z<<" "<<zminpaine<<" "<<zmaxpaine<<" "<<Dk[iz]<<"\n";
		if(Dk[iz]<0)//Modalita' PaineSherr
		{
			if((rint_med!=rmax)&&
			   (((zminpaine==0.0)&&(zmaxpaine==0.0))||((z>zminpaine)&&(z<zmaxpaine)))){
				Dk[iz]=Dk[iz] / (
						0.64309+
						0.00044*exp(rad_k/(rint_med*0.06894))+
						0.35647*exp(rad_k/(rint_med*0.19409))
							  );
			}
			Dk[iz]=-Dk[iz];
		}
		if(Dcl[iz]<0)//Modalita' PaineSherr
		{
			if((rint_med!=rmax)&&
			   (((zminpaine==0.0)&&(zmaxpaine==0.0))||((z>zminpaine)&&(z<zmaxpaine)))){
				Dcl[iz]=Dcl[iz] / (
						0.64309+
						0.00044*exp(rad_cl/(rint_med*0.06894))+
						0.35647*exp(rad_cl/(rint_med*0.19409))
							    );
			}
			Dcl[iz]=-Dcl[iz];
		}

		//Secondo ciclo in it
		//Memorizzazione costanti dielettriche e coefficienti di diffusione
		//	diff_k/diff_cl	Coefficienti di diffusione
		//	sum_phi		Costante dielettrica
		for(it=0;it<nt;it++) 	
		{
			iv=nrnt*iz+nr*it;//Posizione iv a inizio spicchio
			for(ir=0;ir<inside_end[iz*nt+it];ir++)
			{
				sum_phi[iv]=epsH2O;
				diff_k[iv]=Dk[iz];
				diff_cl[iv]=Dcl[iz];
				is_sol[iv]=1;
				iv++;
			}
			for(ir=inside_end[iz*nt+it];ir<outside_begin[iz*nt+it];ir++)
			{
				sum_phi[iv]=epsMOL;
				diff_k[iv]=0;
				diff_cl[iv]=0;
				is_sol[iv]=0;
				iv++;
			}
			for(ir=outside_begin[iz*nt+it];ir<nr;ir++)
			{
				sum_phi[iv]=epsH2O;
				diff_k[iv]=Dk[iz];
				diff_cl[iv]=Dcl[iz];
				is_sol[iv]=1;
				iv++;
			}
		}//Fine secondo ciclo in it
	}//Fine ciclo in iz

	//Ripristino raggi atomici originali
	for(ind_atm=0;ind_atm<chn.num_atm;ind_atm++)chn[ind_atm].rad-=rad_probe;

//-------------------------------------------------
//DISCRETIZZAZIONE CARICA
	setchg(nr,nt,nz,zmin,zmax,rmax,//Dimensioni griglia
		disc_chg,//Modalita' discretizzazione della carica
		chn,//Molecola
		chg_mol,//Vettore in cui e' restituita la carica discretizzata sulla griglia
		min_cell4atom,max_cell4atom,distancemed,distancemax);//Qualita' discretizzazione
//-------------------------------------------------

//-------------------------------------------------
//COSTRUZIONE STRUTTURE DATI
//In questa fase alcuni vettori sono utilizzati con un significato differente dal solito:
//
//	sum_phi contiene le costanti dielettriche delle varie celle
//	Conterra' la somma al denominatore dei coeffificienti dell'equazioni 
//	di Poisson al termine di mkDSP
//
//	kodd e clodd	Contengono i coefficienti di diffusione delle celle
//
//	keven	Contiene la carica della molecola
//	
//	cleven	E' utilizzato per memorizzazioni temporanee durante 
//	la discretizzazione della carica
//	
	mkDSP(nr,nt,nz,dr,dt,dz,phiterm,phimem,phibnd,sum_phi,chg_mol,A_phi,bcs_phi);//Strutture dati Poisson;
	if(PNP!=0){//Strutture dati Nernst
		mkDSN(nr,nt,nz,dr,dt,dz,inside_end,outside_begin,diff_k,kext,kint,Acs_k,bcs_k);
		mkDSN(nr,nt,nz,dr,dt,dz,inside_end,outside_begin,diff_cl,clext,clint,Acs_cl,bcs_cl);
	}
//-------------------------------------------------

	delete[] diff_k;
	delete[] diff_cl;

	return 0;
}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//SETCHG
//Discretizzazione carica con modalita' distribuzione sferica
int setchg(int nr, int nt,int nz,double zmin, double zmax, double rmax,//Dimensioni griglia
		int disc_chg,//Metodo discretizzazione
		molecule &chn,//Molecola da discretizzare
		double *chg_mol,//Carica discretizzata sulla gliglia
		int &min_cell4atom,//Minimo Numero di celle usate per discretizzare un atomo
		int &max_cell4atom,//Massimo Numero di celle usate per discretizzare un atomo
		double &distancemed,//Distanza media tra centro reale atomo e centro carica discretizzata
		double &distancemax//Distanza media tra centro reale atomo e centro carica discretizzata
		)
{
	unsigned long int nrntnz; //Dimensione griglia
	int nrnt;	//Dimensioni griglia
	double dr,dt,dz;//Step griglia
	double r,t,z;	//Posizione in griglia
	double *chgtmp;	//Vettore utilizzato in fase di calcolo per discretizzare la carica
	int ir,it,iz,ind_atm;	//Indici
	unsigned long int iv;
	int itmp,flag_bg;

	double x_atm,y_atm,z_atm,r_atm,t_atm,rad_atm,chg_atm;
	int ir_min,ir_max,iz_min,iz_max,ir_atm,it_atm,iz_atm;
	double x_cell,y_cell,z_cell;
	double x_center,y_center,z_center;
	double distancetmp,sum;
	int ind_atmdistmax;


        //Inizializzazione
	nrntnz=nr*nt*nz;
	nrnt=nr*nt;
	dr=rmax/nr;
	dt=(2*PI)/nt;
	dz=(zmax-zmin)/nz;
	min_cell4atom=nr*nt*nz;
	max_cell4atom=0;
	distancemax=distancemed=0.0;
	ind_atmdistmax=-1;
	
	//Allocazione
	chgtmp = new double [nrntnz];

	for(ind_atm=0;ind_atm<chn.num_atm;ind_atm++){
		itmp=0;	//Contatore numero elementi di griglia su cui si divide la carica
		flag_bg=0;//Flag per riconoscere quando il centro elemento e' lontano
			  //dall'atomo piu' del raggio
		x_atm=chn.mol[ind_atm].x;	//Valori Atomo
		y_atm=chn.mol[ind_atm].y;
		z_atm=chn.mol[ind_atm].z;
		r_atm=chn.mol[ind_atm].r;
		t_atm=chn.mol[ind_atm].t;
		rad_atm=chn.mol[ind_atm].rad;
		chg_atm=chn.mol[ind_atm].chg;
		//Definizione coordinate di griglia atomo
		//Il cast a intero effettua una operazione di troncamento
		//-->Ottengo gli indici dell'elemento di griglia che contiene l'atomo
		ir_atm=(int)floor((r_atm)/dr);
		it_atm=(int)((t_atm+PI)/dt);
		if(it_atm==nt)it_atm=(nt-1);//For atom with t==PI
		if(it_atm==-1)it_atm=0;//For atom with t==-PI
		iz_atm=(int)floor((z_atm-zmin)/dz);
		//Definizione anello intersecato dall'atomo
		//Il cast a intero effettua una operazione di troncamento
		//-->Ottengo gli indici dell'elemento di griglia che il centro dell'atomo
		//   a cui e' aggiunto o tolto il suo raggio
		ir_min=(int)((r_atm-rad_atm)/dr);
		ir_max=(int)((r_atm+rad_atm)/dr);
		iz_min=(int)((z_atm-rad_atm-zmin)/dz);
		iz_max=(int)((z_atm+rad_atm-zmin)/dz);
		if((ir_max>=nr)||(iz_max>=nz)){
			cerr<<"ERROR IN CHARGE DISCRETIZATION: ATOM OUTSIDE BOUNDARIES\n";
			
			cout << "ind_atm: " << ind_atm << endl;
			cout << "x_atm: " << x_atm << endl;
			cout << "y_atm: " << y_atm << endl;
			cout << "z_atm: " << z_atm << endl;
			
			cout << "ir_max: " << ir_max << endl;
			cout << "nr: " << nr << endl;
			cout << "iz_max: " << iz_max << endl;
			cout << "nz: " << nz << endl;
			
			exit(1);
		}
		//Ciclo per individuare celle intersecate dall'atomo
		for(iz=iz_min;iz<(iz_max+1);iz++)
		{
			z_cell=(dz/2.0)+iz*dz+zmin;
			for(it=0;it<nt;it++)
			{
				for(ir=ir_min;ir<(ir_max+1);ir++)
				{
					iv=(nrnt*iz)+(nr*it)+ir;
					x_cell=((dr/2)+ir*dr)*cos((dt/2)+it*dt);
					y_cell=((dr/2)+ir*dr)*sin((dt/2)+it*dt);
					distancetmp=sqrt((x_atm-x_cell)*(x_atm-x_cell)
							+(y_atm-y_cell)*(y_atm-y_cell)
							+(z_atm-z_cell)*(z_atm-z_cell)
							);
					//Il centro dell'elemento di cella dista dall'atomo
					//meno del raggio di quest'ultimo
					//		OR 
					//L'elemento di cella contiene l'atomo
					//---->L'elemento acquistera' parte della carica atomica
					if(distancetmp <= rad_atm)
					{
				switch(disc_chg)
				{
					case 1:
						chgtmp[iv]=distancetmp;
						break;
					case 0:
						chgtmp[iv]=volume_cell(ir,dr,dt,dz);
						break;
				}
				itmp++;
					}
					else
					{
						chgtmp[iv]=0.0;
					}
					if( ((ir==ir_atm)&&(it==it_atm)&&(iz==iz_atm)) &&
						(distancetmp > rad_atm) )
					{
						flag_bg=1;
					}
				}
			}
		}
		//Caso di atomo totalmente contenuto in un elemento di griglia
		//il cui centro dista dalla posizione dell'atomo piu' del raggio
		sum=0.0;
		if(disc_chg==2)itmp=0;
		if(((flag_bg==1)&&(itmp==0))||(disc_chg==2))
		{
			for(iz=iz_min;iz<(iz_max+1);iz++)
			{
				for(it=0;it<nt;it++)
				{
					for(ir=ir_min;ir<(ir_max+1);ir++)
					{
				iv=(nrnt*iz)+(nr*it)+ir;
				if((ir==ir_atm)&&(it==it_atm)&&(iz==iz_atm))
				{
					chgtmp[iv]=1.0;
					sum=1.0;
					itmp++;
				}
				else
				{
					chgtmp[iv]=0.0;
				}
					}
				}
			}
		}
		else
		{
			//Definizione pesi
			for(iz=iz_min;iz<(iz_max+1);iz++)
			{
				for(it=0;it<nt;it++)
				{
					for(ir=ir_min;ir<(ir_max+1);ir++)
					{
			iv=(nrnt*iz)+(nr*it)+ir;
			switch(disc_chg)
			{
				case 1:
			    if(chgtmp[iv]!=0.0)
			    {
				    chgtmp[iv]=(rad_atm-chgtmp[iv]);
			    }
				case 0:
					sum+=chgtmp[iv];
			}
					}
				}
			}
		}

		//Normalizzazione pesi
		for(iz=iz_min;iz<(iz_max+1);iz++)
			for(it=0;it<nt;it++)
				for(ir=ir_min;ir<(ir_max+1);ir++)
				{
					iv=(nrnt*iz)+(nr*it)+ir;
					if(chgtmp[iv]!=0.0)chgtmp[iv]=(chgtmp[iv]/sum);
				}
		//Definizione centro carica discretizzata
		x_center=y_center=z_center=0.0;
		for(iz=iz_min;iz<(iz_max+1);iz++)
		{
			z=iz*dz+(dz/2.0)+zmin;
			for(it=0;it<nt;it++)
			{
				t=it*dt+(dt/2.0);
				for(ir=ir_min;ir<(ir_max+1);ir++)
				{
					iv=(nrnt*iz)+(nr*it)+ir;
					r=ir*dr+(dr/2.0);
					x_center+=(chgtmp[iv]*r*cos(t));
					y_center+=(chgtmp[iv]*r*sin(t));
					z_center+=(chgtmp[iv]*z);
				}
			}
		}
		distancetmp=sqrt(
				((x_center-x_atm)*(x_center-x_atm))+
				((y_center-y_atm)*(y_center-y_atm))+
				((z_center-z_atm)*(z_center-z_atm))
				);

		if(itmp < min_cell4atom)min_cell4atom=itmp;
		if(itmp > max_cell4atom)max_cell4atom=itmp;
		distancemed+=distancetmp;
		if( (distancetmp>distancemax) && ((itmp>1)||(disc_chg==2)) )
		{
			distancemax=distancetmp;
			ind_atmdistmax=ind_atm;
		}

		//Ciclo per assegnare le cariche
		for(iz=iz_min;iz<(iz_max+1);iz++)
		{
			for(it=0;it<nt;it++)
			{
				for(ir=ir_min;ir<(ir_max+1);ir++)
				{
			iv=(nrnt*iz)+(nr*it)+ir;
			//Assegnazione della corretta frazione di carica
			//Dimensionalmente:
			//e2mVA		Fattore di conversione per portare includere eps0
			//		e portare la dimensione a mV*A
			chg_mol[iv]+=( e2mVA * chgtmp[iv] * chg_atm) ;
				}
			}
		}
	}

	distancemed=distancemed/chn.num_atm;

	//cout<<"CHARGE DISCRETIZATION ROUTINE:\n";
	//cout<<"Distance max discretized charge center Vs actual center = "
	//	<<distancemax<<" for atom :\n"<<chn.mol[ind_atmdistmax]<<"\n";
	//cout<<"Distance med discretized charge center Vs actual center = "
	//	<<distancemed<<"\n\n";

	delete[] chgtmp;

	return 0;
}
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//		FUNZIONI AUSILIARIE
//
//GETSLICE
//Restituisce la slice corrispondente ad una determinata z
//La slice restituita e' la piu' alta NON occupata per piu' della meta' da un oggetto
//che si estende fino a z
//Il ciclo for(iz=izinf;iz<izsup;iz++) Copre un oggetto esteso tra zinf e zsup
int getslice(double z, double zmin, double zmax, double dz)
{
	double ztmp;

	if((z>zmax)||(z<zmin))
	{
		cerr<<"ERROR Z OUTSIDE BOUNDARIES\n";
		cerr<<"z = "<<z<<" zmin = "<<zmin<<" zmax = "<<zmax<<"\n";
		exit(1);
	}
	ztmp=z-zmin-(dz/2.0);
	if(ztmp<0)return 0;
	else return ((int)(ztmp/dz)+1);
}

int lin2exp(double phiterm, unsigned long int nrntnz,double *phi, double *k, double *cl)
{
	unsigned long int iv;
	double beta=1.0/phiterm;

	for(iv=0;iv<nrntnz;iv++)
	{
		k[iv]=k[iv]*exp(beta*phi[iv]);
		cl[iv]=cl[iv]*exp(-beta*phi[iv]);
	}

	return 0;
}
int exp2lin(double phiterm, unsigned long int nrntnz,double *phi, double *k, double *cl)
{
	unsigned long int iv;
	double beta=1.0/phiterm;

	for(iv=0;iv<nrntnz;iv++)
	{
		k[iv]=k[iv]*exp(-beta*phi[iv]);
		cl[iv]=cl[iv]*exp(beta*phi[iv]);
	}

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
  
///////////////////////////////////////////////////////////////////////////////
//MKDSP
//Ipotesi:
//	La costante dielettrica di un elemento oltre il confine e' la stessa
//		dell'elemento che vi confina
//	La costante dielettrica sulle superfici di contorno e' la media tra i valori
//		precedentemente assegnati alle varie celle
//
//Costruisce:
//	Vettore sum_phi contiene la somma che compare al denominatore dei coefficienti
//	Matrice A_phi contiene le sestette dei coefficienti
//		Sono diversi da zero anche quelli con elementi esterni ma vengono poi
//		esclusi in fase di prodotto matriciale
//	Vettore bcs_phi
//		Ci finisce la parte dovuta alle condizioni al contorno
int mkDSP(int nr, int nt,int nz, double dr, double dt, double dz,
		double phiterm,double phimem, double *phibnd,
		double *sum_phi,
		double *chg_mol,
		double *A_phi,
		double *bcs_phi
		)
{
	unsigned long int nrntnz;

	int ir,it,iz;				//Indici posizione
	unsigned long int iv,im;		//Indici matrice e vettore
	double here,bk,fw,dx,sx,dw,up,tmpsum;	//variabili usate per memorizzare temporaneamente
	int nrnt,nrnt_nr;		//Variabili di comodo per velocizzare il calcolo
	double dtdz,dzDdt,dr2dtDdz;
	double half;

	//Inizializzazione
	nrntnz=nr*nt*nz;
	nrnt=nr*nt;
	nrnt_nr=nrnt-nr;
	dtdz=dt*dz;
	dzDdt=dz/dt;
	dr2dtDdz=(dr*dr*dt)/dz;
	half=1.0/2.0;
	
	//Inizializzazione 
	iv=0;
	im=0;

	//Fetta iz=0
	iz=0;
		//Primo spicchio it=0
		it=0;
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2.0;	//FW
			dx=(here+sum_phi[iv+nrnt_nr])/2.0;	//DX
			sx=(here+sum_phi[iv+nr])/2.0;	//SX
			dw=here;		//DW
			up=(here+sum_phi[iv+nrnt])/2.0;	//UP
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phimem;	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
				sx=(here+sum_phi[iv+nr])/2;	//SX
				dw=here;		//DW
				up=(here+sum_phi[iv+nrnt])/2;	//UP
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);		//DX
				A_phi[im++]=sx*dzDdt/(ir+half);		//SX
				A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
				//Definizione contributo al termine noto
				bcs_phi[iv]=A_phi[im++]*phimem;	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;	//BK
			fw=here;			//FW
			dx=(here+sum_phi[iv+nrnt_nr])/2;//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=here;			//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//DW+FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
			//Definizione contributo al termine noto
			bcs_phi[iv]+=A_phi[im++]*phimem;	//DW+FW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			iv++;
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;				//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=here;			//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			im++;				//BK
			A_phi[im++]=fw*dtdz*(ir+1);	//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im]=dw*dr2dtDdz*(ir+half);//DW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phimem;//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;	//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv-nr])/2;	//DX
				sx=(here+sum_phi[iv+nr])/2;	//SX
				dw=here;			//DW
				up=(here+sum_phi[iv+nrnt])/2;	//UP
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;		//BK
				A_phi[im++]=fw*dtdz*(ir+1);	//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//SX
				A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
				//Definizione contributo al termine noto
				bcs_phi[iv]=A_phi[im++]*phimem;	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;	//BK
			fw=here;			//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=here;			//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//DW+FW
			A_phi[im++]=dx*dzDdt/(ir+half);		//DX
			A_phi[im++]=sx*dzDdt/(ir+half);		//SX
			A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
			//Definizione contributo al termine noto
			bcs_phi[iv]+=A_phi[im++]*phimem;	//DW+FW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			iv++;
		}
		//Ultimo spicchio it = nt-1
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;				//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv-nrnt_nr])/2;//SX
			dw=here;			//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phimem;	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;	//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv-nr])/2;	//DX
				sx=(here+sum_phi[iv-nrnt_nr])/2;//SX
				dw=here;			//DW
				up=(here+sum_phi[iv+nrnt])/2;	//UP
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);		//DX
				A_phi[im++]=sx*dzDdt/(ir+half);		//SX
				A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
				//Definizione contributo al termine noto
				bcs_phi[iv]=A_phi[im++]*phimem;	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
			dw=here;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//DW+FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im]=dw*dr2dtDdz*(ir+half);	//DW
			//Definizione contributo al termine noto
			bcs_phi[iv]+=A_phi[im++]*phimem;	//DW+FW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			iv++;

	//Fette interne iz=1:nz-2
	for(iz=1;iz<nz-1;iz++)
	{
		//Primo spicchio it=0
		it=0;
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
				sx=(here+sum_phi[iv+nr])/2;	//SX
				dw=(here+sum_phi[iv-nrnt])/2;		//DW
				up=(here+sum_phi[iv+nrnt])/2;	//UP
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//sX
				A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			iv++;
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv-nr])/2;	//DX
				sx=(here+sum_phi[iv+nr])/2;	//SX
				dw=(here+sum_phi[iv-nrnt])/2;		//DW
				up=(here+sum_phi[iv+nrnt])/2;	//UP
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//sX
				A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			iv++;
		}
		//Ultimo spicchio it = nt-1
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv-nr])/2;	//DX
				sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
				dw=(here+sum_phi[iv-nrnt])/2;		//DW
				up=(here+sum_phi[iv+nrnt])/2;	//UP
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//sX
				A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;		//DW
			up=(here+sum_phi[iv+nrnt])/2;	//UP
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			iv++;
	}


	//Fetta iz=nz-1
		//Primo spicchio it=0
		it=0;
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;	//UP
			up=here;		//DW
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			//Definizione contributo al termine noto
			//bcs_phi[iv]=A_phi[im++]*phimem;	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
				sx=(here+sum_phi[iv+nr])/2;	//SX
				dw=(here+sum_phi[iv-nrnt])/2;	//UP
				up=here;		//DW
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//sX
				A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				//Definizione contributo al termine noto
				//bcs_phi[iv]=A_phi[im++]*phimem;	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv+nrnt_nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;	//UP
			up=here;		//DW
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//UP+FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			//Definizione contributo al termine noto
			//bcs_phi[iv]+=A_phi[im++]*phimem;	//UP+FW
			iv++;
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;	//UP
			up=here;		//DW
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			//Definizione contributo al termine noto
			//bcs_phi[iv]=A_phi[im++]*phimem;	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv-nr])/2;	//DX
				sx=(here+sum_phi[iv+nr])/2;	//SX
				dw=(here+sum_phi[iv-nrnt])/2;	//UP
				up=here;		//DW
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//sX
				A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				//Definizione contributo al termine noto
				//bcs_phi[iv]=A_phi[im++]*phimem;	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv+nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;	//UP
			up=here;		//DW
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//UP+FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			//Definizione contributo al termine noto
			//bcs_phi[iv]+=A_phi[im++]*phimem;//UP+FW
			iv++;
		}
		//Ultimo spicchio it = nt-1
			//Primo elemento ir=0
			ir=0;
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=0;			//BK
			fw=(here+sum_phi[iv+1])/2;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;	//UP
			up=here;		//DW
			//Definizione elementi matrice
			im++;	//BK
			A_phi[im++]=fw*dtdz*(ir+1);		//FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
			//Definizione contributo al termine noto
			//bcs_phi[iv]=A_phi[im++]*phimem;	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				//Definizione valore eps nell'elemento in analisi
				here=sum_phi[iv];
				//Definizione valori di eps sulle superfici di contorno
				bk=(here+sum_phi[iv-1])/2;			//BK
				fw=(here+sum_phi[iv+1])/2;	//FW
				dx=(here+sum_phi[iv-nr])/2;	//DX
				sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
				dw=(here+sum_phi[iv-nrnt])/2;	//UP
				up=here;		//DW
				//Definizione elementi matrice
				A_phi[im++]=bk*dtdz*ir;			//BK
				A_phi[im++]=fw*dtdz*(ir+1);		//FW
				A_phi[im++]=dx*dzDdt/(ir+half);	//DX
				A_phi[im++]=sx*dzDdt/(ir+half);	//sX
				A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
				A_phi[im++]=up*dr2dtDdz*(ir+half);	//UP
				//Definizione contributo al termine noto
				//bcs_phi[iv]=A_phi[im++]*phimem;	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			//Definizione valore eps nell'elemento in analisi
			here=sum_phi[iv];
			//Definizione valori di eps sulle superfici di contorno
			bk=(here+sum_phi[iv-1])/2;			//BK
			fw=here;	//FW
			dx=(here+sum_phi[iv-nr])/2;	//DX
			sx=(here+sum_phi[iv-nrnt_nr])/2;	//SX
			dw=(here+sum_phi[iv-nrnt])/2;	//UP
			up=here;		//DW
			//Definizione elementi matrice
			A_phi[im++]=bk*dtdz*ir;			//BK
			A_phi[im]=fw*dtdz*(ir+1);		//FW
			//Definizione contributo al termine noto
			bcs_phi[iv]=A_phi[im++]*phibnd[iz];	//UP+FW
			A_phi[im++]=dx*dzDdt/(ir+half);	//DX
			A_phi[im++]=sx*dzDdt/(ir+half);	//SX
			A_phi[im++]=dw*dr2dtDdz*(ir+half);	//DW
			A_phi[im++]=up*dr2dtDdz*(ir+(half));	//UP
			//Definizione contributo al termine noto
			//bcs_phi[iv]+=A_phi[im++]*phimem;	//UP+FW
			iv++;
			

	if(iv!=nrntnz)
	{
		cerr<<"ERROR in Poisson data structure definition\n";
		exit(1);
	}

	//Normalizzazione degli indici
	im=0;
	iv=0;
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
		tmpsum=A_phi[im]+A_phi[im+1]+A_phi[im+2]+A_phi[im+3]+A_phi[im+4]+A_phi[im+5];
		A_phi[im]=-A_phi[im]/tmpsum;
		im++;
		A_phi[im]=-A_phi[im]/tmpsum;
		im++;
		A_phi[im]=-A_phi[im]/tmpsum;
		im++;
		A_phi[im]=-A_phi[im]/tmpsum;
		im++;
		A_phi[im]=-A_phi[im]/tmpsum;
		im++;
		A_phi[im]=-A_phi[im]/tmpsum;
		im++;
		//inclusione contributo cariche nelle condizioni al contorno
		bcs_phi[iv]=((chg_mol[iv]/phiterm)+bcs_phi[iv])/tmpsum;
		sum_phi[iv]=volume_cell(ir,dr,dt,dz)*mmol2mVA/(tmpsum*phiterm);
		iv++;
			}
		}
	}
	if(iv!=nrntnz)
	{
		cerr<<"ERROR in Poisson data structure definition\n";
		exit(1);
	}

	return 0;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//MKDSN
//Input:
//	D		Vettore contenente il valore dei coeffcienti di diffusione delle varie celle
//				queste informazioni sono in:
//					kodd per k
//					clodd per c
//	Cext		Concentrazione extracellulare
//	Cint		Concentrazione intracellulare
//
//Ipotesi:
//	I coefficienti di diffusione alle superfici di contorno e' la media tra i valori
//		precedentemente assegnati alle varie celle
//	La diffusione verso le basi del cilindro e' sempre permessa e con il coefficiente di
//		diffusione dell'elemento confinante
//	La diffusione verso le pareti laterali e' sempre proibita
//
//Costruisce:
//	Acs matrice che contiene a sestetti il fattore moltiplicativo
//		(indipendente dal potenziale) dei coefficienti
//	bcs vettore contenente il fattore moltiplicativo fisso del termine noto
int mkDSN(int nr, int nt,int nz, double dr, double dt, double dz,
		int *inside_end,int *outside_begin,
		double *D, double Cext, double Cint, double *Acs, double *bcs)
{
	unsigned long int nrntnz;
	int ir,it,iz;				//Indici posizione
	unsigned long int iv,im;		//Indici matrice e vettore
	int end,begin;				//indici inizio fine molecola
	double here,bk,fw,dx,sx,dw,up;	//variabili usate per memorizzare temporaneamente
	int nrnt,nrnt_nr;		//Variabili di comodo per velocizzare il calcolo
	double dtdz,dzDdt,dr2dtDdz;
	double half;

	//Inizializzazione
	nrntnz=nr*nt*nz;
	dtdz=dt*dz;
	dzDdt=dz/dt;
	dr2dtDdz=(dr*dr*dt)/dz;
	nrnt=nr*nt;
	nrnt_nr=nrnt-nr;

	//Inizializzazioni
	iv=0;
	im=0;
	half=1.0/2.0;

	//Fetta iz=0
	iz=0;
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv+nrnt_nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv];		//DW
			Acs[im]=dw*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cint;
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;	
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv+nrnt_nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv];		//DW
			Acs[im]=dw*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cint;
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;

			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv+nrnt_nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv+nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv];		//DW
					Acs[im]=dw*dr2dtDdz*(ir+half);
					//Definizione contributo al termine noto
					bcs[iv]=Acs[im]*Cint;
					im++;
					up=D[iv+nrnt];
					if(up!=0)		//UP
					{
						up=(up+here)/2;
						Acs[im]=up*dr2dtDdz*(ir+half);
					}
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv];		//DW
			Acs[im]=dw*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cint;
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv];		//DW
			Acs[im]=dw*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cint;
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv-nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv+nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv];		//DW
					Acs[im]=dw*dr2dtDdz*(ir+half);
					//Definizione contributo al termine noto
					bcs[iv]=Acs[im]*Cint;
					im++;
					up=D[iv+nrnt];
					if(up!=0)		//UP
					{
						up=(up+here)/2;
						Acs[im]=up*dr2dtDdz*(ir+half);
					}
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
		}
		//Ultimo spicchio it=nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv-nrnt_nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv];		//DW
			Acs[im]=dw*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cint;
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv-nrnt_nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv];		//DW
			Acs[im]=dw*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cint;
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv-nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv-nrnt_nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv];		//DW
					Acs[im]=dw*dr2dtDdz*(ir+half);
					//Definizione contributo al termine noto
					bcs[iv]=Acs[im]*Cint;
					im++;
					up=D[iv+nrnt];
					if(up!=0)		//UP
					{
						up=(up+here)/2;
						Acs[im]=up*dr2dtDdz*(ir+half);
					}
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv];		//DW
				Acs[im]=dw*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cint;
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
//cout<<"iz = "<<iz<<" iv = "<<iv<<" im = "<<im<<"\n";
	//Fette interne iz=1:nz-2
	for(iz=1;iz<nz-1;iz++)
	{
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv+nrnt_nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv+nrnt_nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv+nrnt_nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv+nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv-nrnt];		//DW
					if(dw!=0)	
					{
						dw=(dw+here)/2;
						Acs[im]=dw*dr2dtDdz*(ir+half);
					}
					im++;
					up=D[iv+nrnt];
					if(up!=0)		//UP
					{
						up=(up+here)/2;
						Acs[im]=up*dr2dtDdz*(ir+half);
					}
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv-nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv+nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv-nrnt];		//DW
					if(dw!=0)	
					{
						dw=(dw+here)/2;
						Acs[im]=dw*dr2dtDdz*(ir+half);
					}
					im++;
					up=D[iv+nrnt];
					if(up!=0)		//UP
					{
						up=(up+here)/2;
						Acs[im]=up*dr2dtDdz*(ir+half);
					}
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
		}
		//Ultimo spicchio it=nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv-nrnt_nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv-nrnt_nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv+nrnt];
			if(up!=0)		//UP
			{
				up=(up+here)/2;
				Acs[im]=up*dr2dtDdz*(ir+half);
			}
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv-nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv-nrnt_nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv-nrnt];		//DW
					if(dw!=0)	
					{
						dw=(dw+here)/2;
						Acs[im]=dw*dr2dtDdz*(ir+half);
					}
					im++;
					up=D[iv+nrnt];
					if(up!=0)		//UP
					{
						up=(up+here)/2;
						Acs[im]=up*dr2dtDdz*(ir+half);
					}
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv+nrnt];
				if(up!=0)		//UP
				{
					up=(up+here)/2;
					Acs[im]=up*dr2dtDdz*(ir+half);
				}
				im++;
				iv++;
			}
//cout<<"iz = "<<iz<<" iv = "<<iv<<" im = "<<im<<"\n";
	}
	//Fetta iz=nz-1
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv+nrnt_nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv];	//UP
			Acs[im]=up*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cext;
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv+nrnt_nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv];	//UP
			Acs[im]=up*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cext;
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv+nrnt_nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv+nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv-nrnt];		//DW
					if(dw!=0)	
					{
						dw=(dw+here)/2;
						Acs[im]=dw*dr2dtDdz*(ir+half);
					}
					im++;
					up=D[iv];	//UP
					Acs[im]=up*dr2dtDdz*(ir+half);
					//Definizione contributo al termine noto
					bcs[iv]=Acs[im]*Cext;
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv+nrnt_nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
			}
//cout<<"Primo spicchio iz = "<<iz<<" iv = "<<iv<<" im = "<<im<<"\n";
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv];	//UP
			Acs[im]=up*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cext;
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv+nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv];	//UP
			Acs[im]=up*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cext;
			im++;
			iv++;
			
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}

			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv-nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv+nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv-nrnt];		//DW
					if(dw!=0)	
					{
						dw=(dw+here)/2;
						Acs[im]=dw*dr2dtDdz*(ir+half);
					}
					im++;
					up=D[iv];	//UP
					Acs[im]=up*dr2dtDdz*(ir+half);
					//Definizione contributo al termine noto
					bcs[iv]=Acs[im]*Cext;
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv+nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
			}
		}
//cout<<"Spicchi interni iz = "<<iz<<" iv = "<<iv<<" im = "<<im<<"\n";
		//Ultimo spicchio it=nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=D[iv];
			bk=0;			//BK
			im++;
			fw=D[iv+1];		//FW 
			if(fw!=0)
			{
				fw=(fw+here)/2;	
				Acs[im]=fw*dtdz*(ir+1);
			}
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv-nrnt_nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv];	//UP
			Acs[im]=up*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cext;
			im++;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=D[iv];
			bk=D[iv-1];			//BK
			if(bk!=0)
			{
				bk=(bk+here)/2;	
				Acs[im]=bk*dtdz*ir;
			}
			im++;
			fw=0;		//FW 
			im++;
			dx=D[iv-nr]; 	//DX
			if(dx!=0) 
			{
				dx=(dx+here)/2;
				Acs[im]=dx*dzDdt/(ir+half);
			}
			im++;
			sx=D[iv-nrnt_nr]; 		//SX
			if(sx!=0)
			{
				sx=(sx+here)/2;
				Acs[im]=sx*dzDdt/(ir+half);
			}
			im++;
			dw=D[iv-nrnt];		//DW
			if(dw!=0)	
			{
				dw=(dw+here)/2;
				Acs[im]=dw*dr2dtDdz*(ir+half);
			}
			im++;
			up=D[iv];	//UP
			Acs[im]=up*dr2dtDdz*(ir+half);
			//Definizione contributo al termine noto
			bcs[iv]=Acs[im]*Cext;
			im++;
			iv++;
			//Per saltare gli elementi non parte del sistema
//cout<<"iz = "<<iz<<" iv = "<<iv<<" im = "<<im<<"\n";
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
//cout<<"iz = "<<iz<<" iv = "<<iv<<" im = "<<im<<"\n";
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				here=D[iv];
				bk=0;			//BK
				im++;
				fw=D[iv+1];		//FW 
				if(fw!=0)
				{
					fw=(fw+here)/2;	
					Acs[im]=fw*dtdz*(ir+1);
				}
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=D[iv];
					bk=D[iv-1];			//BK
					if(bk!=0)
					{
						bk=(bk+here)/2;	
						Acs[im]=bk*dtdz*ir;
					}
					im++;
					fw=D[iv+1];		//FW 
					if(fw!=0)
					{
						fw=(fw+here)/2;	
						Acs[im]=fw*dtdz*(ir+1);
					}
					im++;
					dx=D[iv-nr]; 	//DX
					if(dx!=0) 
					{
						dx=(dx+here)/2;
						Acs[im]=dx*dzDdt/(ir+half);
					}
					im++;
					sx=D[iv-nrnt_nr]; 		//SX
					if(sx!=0)
					{
						sx=(sx+here)/2;
						Acs[im]=sx*dzDdt/(ir+half);
					}
					im++;
					dw=D[iv-nrnt];		//DW
					if(dw!=0)	
					{
						dw=(dw+here)/2;
						Acs[im]=dw*dr2dtDdz*(ir+half);
					}
					im++;
					up=D[iv];
					up=D[iv];	//UP
					Acs[im]=up*dr2dtDdz*(ir+half);
					//Definizione contributo al termine noto
					bcs[iv]=Acs[im]*Cext;
					im++;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=D[iv];
				bk=D[iv-1];			//BK
				if(bk!=0)
				{
					bk=(bk+here)/2;	
					Acs[im]=bk*dtdz*ir;
				}
				im++;
				fw=0;		//FW 
				im++;
				dx=D[iv-nr]; 	//DX
				if(dx!=0) 
				{
					dx=(dx+here)/2;
					Acs[im]=dx*dzDdt/(ir+half);
				}
				im++;
				sx=D[iv-nrnt_nr]; 		//SX
				if(sx!=0)
				{
					sx=(sx+here)/2;
					Acs[im]=sx*dzDdt/(ir+half);
				}
				im++;
				dw=D[iv-nrnt];		//DW
				if(dw!=0)	
				{
					dw=(dw+here)/2;
					Acs[im]=dw*dr2dtDdz*(ir+half);
				}
				im++;
				up=D[iv];	//UP
				Acs[im]=up*dr2dtDdz*(ir+half);
				//Definizione contributo al termine noto
				bcs[iv]=Acs[im]*Cext;
				im++;
				iv++;
			}

	if(iv!=nrntnz)
	{
		cerr<<"ERROR in Poisson data structure definition\n";
		exit(1);
	}

	return 0;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//UPDATEDSN
//
//Input:
//Espliciti
//	phi		Vettore del potenziale
//Impliciti
//	Acs_k		Matrici contenenti il fattore moltiplicativo costante dei coefficienti
//	Acs_cl
//	bcs_k		Vettori contenenti il fattore moltiplicativo costante dei termini noti
//	bcs_cl
//
//
//Costruisce:
//	A_k A_cl	Matrici dei coefficienti
//	b_k b_cl 	Vettori dei termini noti
int updateDSN(int nr,int nt, int nz, double dr, double dt, double dz,
		int *inside_end,int *outside_begin,
		double *phi,
		double phimem,double phiterm,
		double *Acs_k, double *Acs_cl, double *bcs_k, double *bcs_cl,
		double *A_k, double *A_cl, double *b_k, double *b_cl
	     )
{
	unsigned long int nrntnz;
	int nrnt,nrnt_nr;

	double beta2;

	int ir,it,iz;				//Indici posizione
	unsigned long int iv,im;		//Indici matrice e vettore
	int end,begin;				//indici inizio fine molecola
	double here,onePbk,onePfw,onePdx,onePsx,onePdw,onePup,sumk;	//variabili usate per memorizzare temporaneamente
	double oneMbk,oneMfw,oneMdx,oneMsx,oneMdw,oneMup,sumcl;	//variabili usate per memorizzare temporaneamente
	double dtdz,dzDdt,dr2dtDdz;

	//Inizializzazione
	nrntnz=nr*nt*nz;
	nrnt=nr*nt;
	nrnt_nr=nr*nt-nr;
	beta2=1.0/(2.0*phiterm);
	dtdz=dt*dz;
	dzDdt=dz/dt;
	dr2dtDdz=(dr*dr*dt)/dz;

	//Inizializzazioni
	iv=0;
	im=0;
	//Fetta iz=0
		//Primo spicchio it=0
		end=inside_end[0];
		begin=outside_begin[0];
			//Primo elemento ir=0
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phimem - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
			b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phimem - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
			b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phimem - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
					b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
			}
		//Spicchi interni it=1:nt-1
		for(it=1;it<(nt-1);it++)
		{
		end=inside_end[it];
		begin=outside_begin[it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phimem - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
			b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phimem - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
			b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phimem - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
					b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
			}
		}
		//Ultimo spicchio it =nt-1
		end=inside_end[it];
		begin=outside_begin[it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phimem - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
			b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phimem - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
			b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phimem - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
					b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phimem - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePdw/sumk);		//BOUNDARY DW
				b_cl[iv]=(bcs_cl[iv])*(oneMdw/sumcl);	
				iv++;
			}
	//Fette interne iz=0:nz-1
	for(iz=1;iz<(nz-1);iz++)
	{
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
			}
		//Spicchi interni it=1:nt-1
		for(it=1;it<(nt-1);it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
			}
		}
		//Ultimo spicchio it =nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(phi[iv+nrnt] - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				iv++;
			}
	}
	//Fetta iz=nz-1
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(0.0 - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
			b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(0.0 - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
			b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(0.0 - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
					b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv+nrnt_nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
			}
		//Spicchi interni it=1:nt-1
		for(it=1;it<(nt-1);it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(0.0 - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
			b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(0.0 - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
			b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(0.0 - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
					b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv+nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
			}
		}
		//Ultimo spicchio it =nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phi[iv];
			onePbk=0;					//BK
			oneMbk=0;
			onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
			oneMfw=2.0 - onePfw;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(0.0 - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			sumk=0.0;
			sumcl=0.0;
			im++;
			A_k[im]=onePfw*Acs_k[im];	//FW
			A_cl[im]=oneMfw*Acs_cl[im];
			sumk+=oneMfw*Acs_k[im];
			sumcl+=onePfw*Acs_cl[im];
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
			b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phi[iv];
			onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
			oneMbk=2.0 - onePbk;
			onePfw=0.0;					//FW
			oneMfw=0.0;
			onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
			oneMdx=2.0 - onePdx;
			onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
			oneMsx=2.0 - onePsx;
			onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
			oneMdw=2.0 - onePdw;
			onePup=1.0 + beta2*(0.0 - here);	//UP
			oneMup=2.0 - onePup;
			A_k[im]=onePbk*Acs_k[im];	//BK
			A_cl[im]=oneMbk*Acs_cl[im];
			sumk=oneMbk*Acs_k[im];
			sumcl=onePbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			sumk+=0.0;
			sumcl+=0.0;
			im++;
			A_k[im]=onePdx*Acs_k[im];	//DX
			A_cl[im]=oneMdx*Acs_cl[im];
			sumk+=oneMdx*Acs_k[im];
			sumcl+=onePdx*Acs_cl[im];
			im++;
			A_k[im]=onePsx*Acs_k[im];	//SX
			A_cl[im]=oneMsx*Acs_cl[im];
			sumk+=oneMsx*Acs_k[im];
			sumcl+=onePsx*Acs_cl[im];
			im++;
			A_k[im]=onePdw*Acs_k[im];	//DW
			A_cl[im]=oneMdw*Acs_cl[im];
			sumk+=oneMdw*Acs_k[im];
			sumcl+=onePdw*Acs_cl[im];
			im++;
			A_k[im]=onePup*Acs_k[im];	//UP
			A_cl[im]=oneMup*Acs_cl[im];
			sumk+=oneMup*Acs_k[im];
			sumcl+=onePup*Acs_cl[im];
			im++;
			A_k[im-6]=-(A_k[im-6])/sumk;
			A_k[im-5]=-(A_k[im-5])/sumk;
			A_k[im-4]=-(A_k[im-4])/sumk;
			A_k[im-3]=-(A_k[im-3])/sumk;
			A_k[im-2]=-(A_k[im-2])/sumk;
			A_k[im-1]=-(A_k[im-1])/sumk;
			A_cl[im-6]=-(A_cl[im-6])/sumcl;
			A_cl[im-5]=-(A_cl[im-5])/sumcl;
			A_cl[im-4]=-(A_cl[im-4])/sumcl;
			A_cl[im-3]=-(A_cl[im-3])/sumcl;
			A_cl[im-2]=-(A_cl[im-2])/sumcl;
			A_cl[im-1]=-(A_cl[im-1])/sumcl;
			b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
			b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phi[iv];
				onePbk=0;					//BK
				oneMbk=0;
				onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
				oneMfw=2.0 - onePfw;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				sumk=0.0;
				sumcl=0.0;
				im++;
				A_k[im]=onePfw*Acs_k[im];	//FW
				A_cl[im]=oneMfw*Acs_cl[im];
				sumk+=oneMfw*Acs_k[im];
				sumcl+=onePfw*Acs_cl[im];
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phi[iv];
					onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
					oneMbk=2.0 - onePbk;
					onePfw=1.0 + beta2*(phi[iv+1] - here);		//FW
					oneMfw=2.0 - onePfw;
					onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
					oneMdx=2.0 - onePdx;
					onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
					oneMsx=2.0 - onePsx;
					onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
					oneMdw=2.0 - onePdw;
					onePup=1.0 + beta2*(0.0 - here);	//UP
					oneMup=2.0 - onePup;
					A_k[im]=onePbk*Acs_k[im];	//BK
					A_cl[im]=oneMbk*Acs_cl[im];
					sumk=oneMbk*Acs_k[im];
					sumcl=onePbk*Acs_cl[im];
					im++;
					A_k[im]=onePfw*Acs_k[im];	//FW
					A_cl[im]=oneMfw*Acs_cl[im];
					sumk+=oneMfw*Acs_k[im];
					sumcl+=onePfw*Acs_cl[im];
					im++;
					A_k[im]=onePdx*Acs_k[im];	//DX
					A_cl[im]=oneMdx*Acs_cl[im];
					sumk+=oneMdx*Acs_k[im];
					sumcl+=onePdx*Acs_cl[im];
					im++;
					A_k[im]=onePsx*Acs_k[im];	//SX
					A_cl[im]=oneMsx*Acs_cl[im];
					sumk+=oneMsx*Acs_k[im];
					sumcl+=onePsx*Acs_cl[im];
					im++;
					A_k[im]=onePdw*Acs_k[im];	//DW
					A_cl[im]=oneMdw*Acs_cl[im];
					sumk+=oneMdw*Acs_k[im];
					sumcl+=onePdw*Acs_cl[im];
					im++;
					A_k[im]=onePup*Acs_k[im];	//UP
					A_cl[im]=oneMup*Acs_cl[im];
					sumk+=oneMup*Acs_k[im];
					sumcl+=onePup*Acs_cl[im];
					im++;
					A_k[im-6]=-(A_k[im-6])/sumk;
					A_k[im-5]=-(A_k[im-5])/sumk;
					A_k[im-4]=-(A_k[im-4])/sumk;
					A_k[im-3]=-(A_k[im-3])/sumk;
					A_k[im-2]=-(A_k[im-2])/sumk;
					A_k[im-1]=-(A_k[im-1])/sumk;
					A_cl[im-6]=-(A_cl[im-6])/sumcl;
					A_cl[im-5]=-(A_cl[im-5])/sumcl;
					A_cl[im-4]=-(A_cl[im-4])/sumcl;
					A_cl[im-3]=-(A_cl[im-3])/sumcl;
					A_cl[im-2]=-(A_cl[im-2])/sumcl;
					A_cl[im-1]=-(A_cl[im-1])/sumcl;
					b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
					b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phi[iv];
				onePbk=1.0 + beta2*(phi[iv-1] - here);		//BK
				oneMbk=2.0 - onePbk;
				onePfw=0.0;					//FW
				oneMfw=0.0;
				onePdx=1.0 + beta2*(phi[iv-nr] - here);	//DX
				oneMdx=2.0 - onePdx;
				onePsx=1.0 + beta2*(phi[iv-nrnt_nr] - here);		//SX
				oneMsx=2.0 - onePsx;
				onePdw=1.0 + beta2*(phi[iv-nrnt] - here);		//DW
				oneMdw=2.0 - onePdw;
				onePup=1.0 + beta2*(0.0 - here);	//UP
				oneMup=2.0 - onePup;
				A_k[im]=onePbk*Acs_k[im];	//BK
				A_cl[im]=oneMbk*Acs_cl[im];
				sumk=oneMbk*Acs_k[im];
				sumcl=onePbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				sumk+=0.0;
				sumcl+=0.0;
				im++;
				A_k[im]=onePdx*Acs_k[im];	//DX
				A_cl[im]=oneMdx*Acs_cl[im];
				sumk+=oneMdx*Acs_k[im];
				sumcl+=onePdx*Acs_cl[im];
				im++;
				A_k[im]=onePsx*Acs_k[im];	//SX
				A_cl[im]=oneMsx*Acs_cl[im];
				sumk+=oneMsx*Acs_k[im];
				sumcl+=onePsx*Acs_cl[im];
				im++;
				A_k[im]=onePdw*Acs_k[im];	//DW
				A_cl[im]=oneMdw*Acs_cl[im];
				sumk+=oneMdw*Acs_k[im];
				sumcl+=onePdw*Acs_cl[im];
				im++;
				A_k[im]=onePup*Acs_k[im];	//UP
				A_cl[im]=oneMup*Acs_cl[im];
				sumk+=oneMup*Acs_k[im];
				sumcl+=onePup*Acs_cl[im];
				im++;
				A_k[im-6]=-(A_k[im-6])/sumk;
				A_k[im-5]=-(A_k[im-5])/sumk;
				A_k[im-4]=-(A_k[im-4])/sumk;
				A_k[im-3]=-(A_k[im-3])/sumk;
				A_k[im-2]=-(A_k[im-2])/sumk;
				A_k[im-1]=-(A_k[im-1])/sumk;
				A_cl[im-6]=-(A_cl[im-6])/sumcl;
				A_cl[im-5]=-(A_cl[im-5])/sumcl;
				A_cl[im-4]=-(A_cl[im-4])/sumcl;
				A_cl[im-3]=-(A_cl[im-3])/sumcl;
				A_cl[im-2]=-(A_cl[im-2])/sumcl;
				A_cl[im-1]=-(A_cl[im-1])/sumcl;
				b_k[iv]=(bcs_k[iv])*(onePup/sumk);		//BOUNDARY UP
				b_cl[iv]=(bcs_cl[iv])*(oneMup/sumcl);	
				iv++;
			}

	if(iv!=nrntnz)
	{
		cerr<<"ERROR in Poisson data structure definition\n";
		exit(1);
	}

	return 0;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//UPDATEDSN_EXP
//
//Input:
//Espliciti
//	phi		Vettore del potenziale
//Impliciti
//	Acs_k		Matrici contenenti il fattore moltiplicativo costante dei coefficienti
//	Acs_cl
//	bcs_k		Vettori contenenti il fattore moltiplicativo costante dei termini noti
//	bcs_cl
//
//
//Costruisce:
//	A_k A_cl	Matrici dei coefficienti
//	b_k b_cl 	Vettori dei termini noti
int updateDSN_exp(int nr, int nt,int nz,double dr, double dt, double dz,
		int *inside_end,int *outside_begin,
		double phiterm, double phimem,
		double *phiexp,
		double *Acs_k, double *Acs_cl, double *bcs_k, double *bcs_cl,
		double *A_k, double *A_cl, double *b_k, double *b_cl
		)
{
	unsigned long int nrntnz;
	int nrnt,nrnt_nr;

	double beta2;
	int ir,it,iz;				//Indici posiziexp
	unsigned long int iv,im;		//Indici matrice e vettore
	int end,begin;				//indici inizio fine molecola
	//expP?? = exp((phi??-here)/phiterm)
	//expM?? = exp(-(phi??-here)/phiterm)
	double here,expPbk,expPfw,expPdx,expPsx,expPdw,expPup,sumk;
		//variabili usate per memorizzare temporaneamente
	double expMbk,expMfw,expMdx,expMsx,expMdw,expMup,sumcl;	
		//variabili usate per memorizzare temporaneamente
	
	nrntnz=nr*nt*nz;
	nrnt=nr*nt;
	nrnt_nr=nr*nt-nr;
	beta2=1.0/(2.0*phiterm);

	iv=0;
	im=0;
	//Fetta iz=0
		//Primo spicchio it=0
		end=inside_end[0];
		begin=outside_begin[0];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (exp(beta2*phimem) * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
			b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (exp(beta2*phimem) * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
			b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv+nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (exp(beta2*phimem) * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (phiexp[iv+nrnt] * here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
					b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
					b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
					sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
					sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
					A_k[im-6]=-A_k[im-6]/sumk;
					A_cl[im-6]=-A_cl[im-6]/sumcl;
					A_k[im-5]=-A_k[im-5]/sumk;
					A_cl[im-5]=-A_cl[im-5]/sumcl;
					A_k[im-4]=-A_k[im-4]/sumk;
					A_cl[im-4]=-A_cl[im-4]/sumcl;
					A_k[im-3]=-A_k[im-3]/sumk;
					A_cl[im-3]=-A_cl[im-3]/sumcl;
					A_k[im-2]=-A_k[im-2]/sumk;
					A_cl[im-2]=-A_cl[im-2]/sumcl;
					A_k[im-1]=-A_k[im-1]/sumk;
					A_cl[im-1]=-A_cl[im-1]/sumcl;
					b_k[iv]=b_k[iv]/sumk;
					b_cl[iv]=b_cl[iv]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
		//Spicchi interni it=1:nt-1
		for(it=1;it<(nt-1);it++)
		{
		end=inside_end[it];
		begin=outside_begin[it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (exp(beta2*phimem) * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
			b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (exp(beta2*phimem) * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
			b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv-nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv+nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (exp(beta2*phimem) * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (phiexp[iv+nrnt] * here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
					b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
					b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
					sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
					sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
					A_k[im-6]=-A_k[im-6]/sumk;
					A_cl[im-6]=-A_cl[im-6]/sumcl;
					A_k[im-5]=-A_k[im-5]/sumk;
					A_cl[im-5]=-A_cl[im-5]/sumcl;
					A_k[im-4]=-A_k[im-4]/sumk;
					A_cl[im-4]=-A_cl[im-4]/sumcl;
					A_k[im-3]=-A_k[im-3]/sumk;
					A_cl[im-3]=-A_cl[im-3]/sumcl;
					A_k[im-2]=-A_k[im-2]/sumk;
					A_cl[im-2]=-A_cl[im-2]/sumcl;
					A_k[im-1]=-A_k[im-1]/sumk;
					A_cl[im-1]=-A_cl[im-1]/sumcl;
					b_k[iv]=b_k[iv]/sumk;
					b_cl[iv]=b_cl[iv]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
		}
		//Ultimo spicchio it =nt-1
		end=inside_end[it];
		begin=outside_begin[it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (exp(beta2*phimem) * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (exp(beta2*phimem) * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv-nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (exp(beta2*phimem) * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (phiexp[iv+nrnt] * here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (exp(beta2*phimem) * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMdw*exp(phimem);	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPdw*exp(-phimem);	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
	//Fette interne iz=0:nz-1
	for(iz=1;iz<(nz-1);iz++)
	{
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv+nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (phiexp[iv-nrnt] * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (phiexp[iv+nrnt] * here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
					sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
					sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
					A_k[im-6]=-A_k[im-6]/sumk;
					A_cl[im-6]=-A_cl[im-6]/sumcl;
					A_k[im-5]=-A_k[im-5]/sumk;
					A_cl[im-5]=-A_cl[im-5]/sumcl;
					A_k[im-4]=-A_k[im-4]/sumk;
					A_cl[im-4]=-A_cl[im-4]/sumcl;
					A_k[im-3]=-A_k[im-3]/sumk;
					A_cl[im-3]=-A_cl[im-3]/sumcl;
					A_k[im-2]=-A_k[im-2]/sumk;
					A_cl[im-2]=-A_cl[im-2]/sumcl;
					A_k[im-1]=-A_k[im-1]/sumk;
					A_cl[im-1]=-A_cl[im-1]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
			}
		//Spicchi interni it=1:nt-1
		for(it=1;it<(nt-1);it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv-nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv+nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (phiexp[iv-nrnt] * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (phiexp[iv+nrnt] * here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
					sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
					sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
					A_k[im-6]=-A_k[im-6]/sumk;
					A_cl[im-6]=-A_cl[im-6]/sumcl;
					A_k[im-5]=-A_k[im-5]/sumk;
					A_cl[im-5]=-A_cl[im-5]/sumcl;
					A_k[im-4]=-A_k[im-4]/sumk;
					A_cl[im-4]=-A_cl[im-4]/sumcl;
					A_k[im-3]=-A_k[im-3]/sumk;
					A_cl[im-3]=-A_cl[im-3]/sumcl;
					A_k[im-2]=-A_k[im-2]/sumk;
					A_cl[im-2]=-A_cl[im-2]/sumcl;
					A_k[im-1]=-A_k[im-1]/sumk;
					A_cl[im-1]=-A_cl[im-1]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
			}
		}
		//Ultimo spicchio it =nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (phiexp[iv+nrnt] * here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv-nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (phiexp[iv-nrnt] * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (phiexp[iv+nrnt] * here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (phiexp[iv+nrnt] * here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
				iv++;
			}
	}
	//Fetta iz=nz-1
		//Primo spicchio it=0
		it=0;
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv+nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (phiexp[iv-nrnt] * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv+nrnt_nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
		//Spicchi interni it=1:nt-1
		for(it=1;it<(nt-1);it++)
		{
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv+nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv-nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv+nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (phiexp[iv-nrnt] * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv+nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
		}
		//Ultimo spicchio it =nt-1
		end=inside_end[iz*nt+it];
		begin=outside_begin[iz*nt+it];
			//Primo elemento ir=0
			ir=0;
			here=phiexp[iv];
			expPbk=0;					//BK
			expMbk=0;
			expPfw= (phiexp[iv+1] * here);		//FW
			expMfw=1.0/expPfw;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (here);	//UP
			expMup=1.0/expPup;
			A_k[im]=0.0;			//BK
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMfw*Acs_k[im];	//FW
			A_cl[im]=expPfw*Acs_cl[im];
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}
			//Ultimo elemento ir=end-1
			here=phiexp[iv];
			expPbk= (phiexp[iv-1] * here);		//BK
			expMbk=1.0/expPbk;
			expPfw=0.0;					//FW
			expMfw=0.0;
			expPdx= (phiexp[iv-nr] * here);	//DX
			expMdx=1.0/expPdx;
			expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
			expMsx=1.0/expPsx;
			expPdw= (phiexp[iv-nrnt] * here);		//DW
			expMdw=1.0/expPdw;
			expPup= (here);	//UP
			expMup=1.0/expPup;
			A_k[im]=expMbk*Acs_k[im];	//BK
			A_cl[im]=expPbk*Acs_cl[im];
			im++;
			A_k[im]=0.0;	//FW
			A_cl[im]=0.0;
			im++;
			A_k[im]=expMdx*Acs_k[im];	//DX
			A_cl[im]=expPdx*Acs_cl[im];
			im++;
			A_k[im]=expMsx*Acs_k[im];	//SX
			A_cl[im]=expPsx*Acs_cl[im];
			im++;
			A_k[im]=expMdw*Acs_k[im];	//DW
			A_cl[im]=expPdw*Acs_cl[im];
			im++;
			A_k[im]=expMup*Acs_k[im];	//UP
			A_cl[im]=expPup*Acs_cl[im];
			im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
			sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
			sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
			A_k[im-6]=-A_k[im-6]/sumk;
			A_cl[im-6]=-A_cl[im-6]/sumcl;
			A_k[im-5]=-A_k[im-5]/sumk;
			A_cl[im-5]=-A_cl[im-5]/sumcl;
			A_k[im-4]=-A_k[im-4]/sumk;
			A_cl[im-4]=-A_cl[im-4]/sumcl;
			A_k[im-3]=-A_k[im-3]/sumk;
			A_cl[im-3]=-A_cl[im-3]/sumcl;
			A_k[im-2]=-A_k[im-2]/sumk;
			A_cl[im-2]=-A_cl[im-2]/sumcl;
			A_k[im-1]=-A_k[im-1]/sumk;
			A_cl[im-1]=-A_cl[im-1]/sumcl;
			b_k[iv]=b_k[iv]/sumk;
			b_cl[iv]=b_cl[iv]/sumcl;
			iv++;
			//Per saltare gli elementi non parte del sistema
			for(ir=end;ir<begin;ir++)
			{
				iv++;
				im+=6;
			}
			//Elementi di sistema esterni all'interruziexp
			if (begin<nr)
			{
			//Primo elemento ir=begin
				ir=begin;
				here=phiexp[iv];
				expPbk=0;					//BK
				expMbk=0;
				expPfw= (phiexp[iv+1] * here);		//FW
				expMfw=1.0/expPfw;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=0.0;			//BK
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMfw*Acs_k[im];	//FW
				A_cl[im]=expPfw*Acs_cl[im];
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					here=phiexp[iv];
					expPbk= (phiexp[iv-1] * here);		//BK
					expMbk=1.0/expPbk;
					expPfw= (phiexp[iv+1] * here);		//FW
					expMfw=1.0/expPfw;
					expPdx= (phiexp[iv-nr] * here);	//DX
					expMdx=1.0/expPdx;
					expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
					expMsx=1.0/expPsx;
					expPdw= (phiexp[iv-nrnt] * here);		//DW
					expMdw=1.0/expPdw;
					expPup= (here);	//UP
					expMup=1.0/expPup;
					A_k[im]=expMbk*Acs_k[im];	//BK
					A_cl[im]=expPbk*Acs_cl[im];
					im++;
					A_k[im]=expMfw*Acs_k[im];	//FW
					A_cl[im]=expPfw*Acs_cl[im];
					im++;
					A_k[im]=expMdx*Acs_k[im];	//DX
					A_cl[im]=expPdx*Acs_cl[im];
					im++;
					A_k[im]=expMsx*Acs_k[im];	//SX
					A_cl[im]=expPsx*Acs_cl[im];
					im++;
					A_k[im]=expMdw*Acs_k[im];	//DW
					A_cl[im]=expPdw*Acs_cl[im];
					im++;
					A_k[im]=expMup*Acs_k[im];	//UP
					A_cl[im]=expPup*Acs_cl[im];
					im++;
					b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
					b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
					sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
					sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
					A_k[im-6]=-A_k[im-6]/sumk;
					A_cl[im-6]=-A_cl[im-6]/sumcl;
					A_k[im-5]=-A_k[im-5]/sumk;
					A_cl[im-5]=-A_cl[im-5]/sumcl;
					A_k[im-4]=-A_k[im-4]/sumk;
					A_cl[im-4]=-A_cl[im-4]/sumcl;
					A_k[im-3]=-A_k[im-3]/sumk;
					A_cl[im-3]=-A_cl[im-3]/sumcl;
					A_k[im-2]=-A_k[im-2]/sumk;
					A_cl[im-2]=-A_cl[im-2]/sumcl;
					A_k[im-1]=-A_k[im-1]/sumk;
					A_cl[im-1]=-A_cl[im-1]/sumcl;
					b_k[iv]=b_k[iv]/sumk;
					b_cl[iv]=b_cl[iv]/sumcl;
					iv++;
				}
				//Ultimo elemento ir=nr-1
				here=phiexp[iv];
				expPbk= (phiexp[iv-1] * here);		//BK
				expMbk=1.0/expPbk;
				expPfw=0.0;					//FW
				expMfw=0.0;
				expPdx= (phiexp[iv-nr] * here);	//DX
				expMdx=1.0/expPdx;
				expPsx= (phiexp[iv-nrnt_nr] * here);		//SX
				expMsx=1.0/expPsx;
				expPdw= (phiexp[iv-nrnt] * here);		//DW
				expMdw=1.0/expPdw;
				expPup= (here);	//UP
				expMup=1.0/expPup;
				A_k[im]=expMbk*Acs_k[im];	//BK
				A_cl[im]=expPbk*Acs_cl[im];
				im++;
				A_k[im]=0.0;	//FW
				A_cl[im]=0.0;
				im++;
				A_k[im]=expMdx*Acs_k[im];	//DX
				A_cl[im]=expPdx*Acs_cl[im];
				im++;
				A_k[im]=expMsx*Acs_k[im];	//SX
				A_cl[im]=expPsx*Acs_cl[im];
				im++;
				A_k[im]=expMdw*Acs_k[im];	//DW
				A_cl[im]=expPdw*Acs_cl[im];
				im++;
				A_k[im]=expMup*Acs_k[im];	//UP
				A_cl[im]=expPup*Acs_cl[im];
				im++;
				b_k[iv]=bcs_k[iv]*expMup;	//Termine Noto
				b_cl[iv]=bcs_cl[iv]*expPup;	//Termine Noto
				sumk=A_k[im-6]+A_k[im-5]+A_k[im-4]+A_k[im-3]+A_k[im-2]+A_k[im-1];
				sumcl=A_cl[im-6]+A_cl[im-5]+A_cl[im-4]+A_cl[im-3]+A_cl[im-2]+A_cl[im-1];
				A_k[im-6]=-A_k[im-6]/sumk;
				A_cl[im-6]=-A_cl[im-6]/sumcl;
				A_k[im-5]=-A_k[im-5]/sumk;
				A_cl[im-5]=-A_cl[im-5]/sumcl;
				A_k[im-4]=-A_k[im-4]/sumk;
				A_cl[im-4]=-A_cl[im-4]/sumcl;
				A_k[im-3]=-A_k[im-3]/sumk;
				A_cl[im-3]=-A_cl[im-3]/sumcl;
				A_k[im-2]=-A_k[im-2]/sumk;
				A_cl[im-2]=-A_cl[im-2]/sumcl;
				A_k[im-1]=-A_k[im-1]/sumk;
				A_cl[im-1]=-A_cl[im-1]/sumcl;
				b_k[iv]=b_k[iv]/sumk;
				b_cl[iv]=b_cl[iv]/sumcl;
				iv++;
			}

	if(iv!=nrntnz)
	{
		cerr<<"ERROR in Nernst data structure definition\n";
		exit(1);
	}

	return 0;
}
///////////////////////////////////////////////////////////////////////////////

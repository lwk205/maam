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
#include <src/nepal/global.h>
#include <lib/etc.h>
#include <lib/molecule.h>

///////////////////////////////////////////////////////////////////////////////
//GET_FLUX_LIN
//Calcolo dei flussi ionici, per risoluzioni di Nernst nel campo lineare
//CONVENZIONE DEI SEGNI
//IL FLUSSO E' SEMPRE POSITIVO QUANDO USCENTE
//IL FLUSSO ASSIALE QUANDO VERSO L'ESTERNO
//LA CORRENTE QUANDO VERSO L'ESTERNO
//
//Input:
//	double* D 	Coefficiente di diffusione 
//	double* c	Concentrazione ionica
//	double* phi	Potenziale elettrostatico
//	double cint	Concentrazione dal lato intracellulare 
//	double cext	Concentrazione dal lato intracellulare 
//	int val		Segno della carica
///////////////////////////////////////////////////////////////////////////////
int get_flux(double *D, double *c, double cint, double cext, double *phi, int val,
		double *Jbk,double *Jfw,double *Jsx,double *Jdx,double *Jdw,double *Jup)
{
	int ir,it,iz,idx,isx,end,begin;
	double betav,betav2,Dslice,Ddw,Dup,drdz;
	double Sbk,Sfw,Sdwup,dstep;
	unsigned long int iv;

	//Inizializzazioni
	drdz=dr*dz;
	betav=((double)val)/phiterm;
	betav2=betav/2.0;
	if(lin_exp==1)
		for(iv=0;iv<nrntnz;iv++)c[iv]=c[iv]*exp(betav*phi[iv]);

	iv=0;
	for(iz=0;iz<nz;iz++)
	{
		Dslice=D[iz];
		if(iz==0) Ddw=D[iz];
		else Ddw=(D[iz]+D[iz-1])/2.0;
		if(iz==(nz-1)) Dup=D[iz];
		else Dup=(D[iz]+D[iz+1])/2.0;

		for(it=0;it<nt;it++)
		{
			end=inside_end[iz*nt+it];
			begin=outside_begin[iz*nt+it];

			for(ir=0;ir<nr;ir++)
			{
	if((ir<begin)&&(ir>=end))
	{
		iv++;
	}
	else
	{
		Sbk=(2*PI*(ir)*dr*dz)/nt;
		Sfw=(2*PI*(ir+1)*dr*dz)/nt;
		Sdwup=(PI*dr*dr*(2.0*ir+1))/nt;
		dstep=2*PI*(ir+0.5)*dr/nt;

		//BK
		if((ir==0)||(ir==begin)) Jbk[iv]=0;
		else 	if(lin_exp==0) 	Jbk[iv]= is_sol[iv-1] * Sbk * (-Dslice) * ( ((c[iv-1]-c[iv])/dr) +
			      		betav * ((c[iv-1]+c[iv])/2.0) * ((phi[iv-1]-phi[iv])/dr) );
			else	Jbk[iv]= is_sol[iv-1] * Sbk * (-Dslice) * exp (-betav2*(phi[iv-1]+phi[iv]))
				* ((c[iv-1]-c[iv])/dr); 
		
		//FW
		if((ir==(nr-1))||(ir==(end-1))) Jfw[iv]=0;
		else 	if(lin_exp==0)	Jfw[iv]= is_sol[iv+1] * Sfw * (-Dslice) * ( ((c[iv+1]-c[iv])/dr) +
			      		betav * ((c[iv+1]+c[iv])/2.0) * ((phi[iv+1]-phi[iv])/dr) );
			else	Jfw[iv]= is_sol[iv+1] * Sfw * (-Dslice) * exp (-betav2*(phi[iv+1]+phi[iv])) 
				* ((c[iv+1]-c[iv])/dr); 

		//DX
		if(it==0) idx=iv+nrnt_nr;
		else idx=iv-nr;
		if(lin_exp==0) 	Jdx[iv]= is_sol[idx] * drdz * (-Dslice) * ( ((c[idx]-c[iv])/dstep) +
			 	betav * ((c[idx]+c[iv])/2.0) * ((phi[idx]-phi[iv])/dstep) );
		else	Jdx[iv]= is_sol[idx] * drdz * (-Dslice) 
			* exp (-betav2*(phi[idx]+phi[iv])) * ((c[idx]-c[iv])/dstep); 
		
		//SX
		if(it==(nt-1)) isx=iv-nrnt_nr;
		else isx=iv+nr;
		if(lin_exp==0)	Jsx[iv]= is_sol[isx] * drdz * (-Dslice) * ( ((c[isx]-c[iv])/dstep) +
			 	betav * ((c[isx]+c[iv])/2.0) * ((phi[isx]-phi[iv])/dstep) );
		else	Jsx[iv]= is_sol[isx] * drdz * (-Dslice) * exp (-betav2*(phi[isx]+phi[iv])) 
			* ((c[isx]-c[iv])/dstep); 

		//DW
		if(iz==0) 	if(lin_exp==0)  Jdw[iv]= Sdwup * (-Ddw) * ( ((cint-c[iv])/dz) +
				   	 	betav * ((cint+c[iv])/2.0) * ((phimem-phi[iv])/dz) );
				else	Jdw[iv]= Sdwup * (-Ddw) * exp (-betav2*(phimem+phi[iv])) 
					* ((cint*exp(betav*phimem)-c[iv])/dz); 
		else 	if(lin_exp==0)	Jdw[iv]= is_sol[iv-nrnt] * Sdwup * (-Ddw) * ( ((c[iv-nrnt]-c[iv])/dz)
					+ betav * ((c[iv-nrnt]+c[iv])/2.0) * ((phi[iv-nrnt]-phi[iv])/dz) );
			else	Jdw[iv]= is_sol[iv-nrnt] * Sdwup * (-Ddw) *
				exp (-betav2*(phi[iv-nrnt]+phi[iv])) * ((c[iv-nrnt]-c[iv])/dz); 

		//UP
		if(iz==(nz-1)) 	if(lin_exp==0)	Jup[iv]= Sdwup * (-Dup) * ( ((cext-c[iv])/dz) +
				        	betav * ((cext+c[iv])/2.0) * ((0.0-phi[iv])/dz) );
				else 	Jup[iv]= Sdwup * (-Dup) * exp (-betav2*(0.0+phi[iv])) 
					* ((cext-c[iv])/dz); 
		else 	if(lin_exp==0)Jup[iv]= is_sol[iv+nrnt] * Sdwup * (-Dup) * ( ((c[iv+nrnt]-c[iv])/dz) 
				      + betav * ((c[iv+nrnt]+c[iv])/2.0) * ((phi[iv+nrnt]-phi[iv])/dz) );
		 	else  	Jup[iv]= is_sol[iv+nrnt] * Sdwup * (-Dup) *
				exp (-betav2*(phi[iv+nrnt]+phi[iv])) * ((c[iv+nrnt]-c[iv])/dz); 

		iv++;
	}
			}
		}
	}

	if(lin_exp==1)
		for(iv=0;iv<nrntnz;iv++)c[iv]=c[iv]*exp(-betav*phi[iv]);

	return 0;
}

int fluxes(double* k, double *cl, double *phi)
{
	unsigned long int iv,num_sol;
	double *Jkbk,*Jkfw,*Jkdx,*Jksx,*Jkdw,*Jkup;
	double *Jclbk,*Jclfw,*Jcldx,*Jclsx,*Jcldw,*Jclup;
	double *Jslicek, *Jslicecl;
	double Jerrk_max,Jerrk_med,Jerrcl_max,Jerrcl_med,Jerrtmp;
	int iz_errk,it_errk,ir_errk;
	int iz_errcl,it_errcl,ir_errcl;
	double Jaxisk_med,Jaxisk_max,Jaxisk_min,rmsdk;	
	int iz_axisk_max,iz_axisk_min;
	double Jaxiscl_med,Jaxiscl_max,Jaxiscl_min,rmsdcl;	
	int iz_axiscl_max,iz_axiscl_min;
	double Ik,Icl,Itot;
	int iz,it,ir;
	double z,t,r,x,y,J_z,J_t,J_r,J_x,J_y;
	double Jk_zmax,Jk_tmax,Jk_rmax;
	double Jcl_zmax,Jcl_tmax,Jcl_rmax;

	//Allocazione & Inizializzaione memoria
	Jslicek = new double [nz];
	Jslicecl = new double [nz];
	Jkbk = new double [nrntnz];
	Jkfw = new double [nrntnz];
	Jkdx = new double [nrntnz];
	Jksx = new double [nrntnz];
	Jkdw = new double [nrntnz];
	Jkup = new double [nrntnz];
	Jclbk = new double [nrntnz];
	Jclfw = new double [nrntnz];
	Jcldx = new double [nrntnz];
	Jclsx = new double [nrntnz];
	Jcldw = new double [nrntnz];
	Jclup = new double [nrntnz];
	for(iv=0;iv<nrntnz;iv++)
	{
		Jkbk[iv]=Jkfw[iv]=Jkdx[iv]=Jksx[iv]=Jkdw[iv]=Jkup[iv]=0.0;
		Jclbk[iv]=Jclfw[iv]=Jcldx[iv]=Jclsx[iv]=Jcldw[iv]=Jclup[iv]=0.0;
	}
	for(iz=0;iz<nz;iz++)
	{
		Jslicek[iz]=Jslicecl[iz]=0.0;
	}
	iz_errk=it_errk=ir_errk=-1;
	iz_errcl=it_errcl=ir_errcl=-1;

	//Calcolo campi vettoriali di flusso
	get_flux(Dk,k,kint,kext,phi,+1,Jkbk,Jkfw,Jkdx,Jksx,Jkdw,Jkup);
	get_flux(Dcl,cl,clint,clext,phi,-1,Jclbk,Jclfw,Jcldx,Jclsx,Jcldw,Jclup);

	//Estrazione parametri dai dati di flusso
	num_sol=0;
	Jerrk_max=Jerrcl_max=Jerrk_med=Jerrcl_med=0.0;
	Jk_zmax=Jk_tmax=Jk_rmax=Jcl_zmax=Jcl_tmax=Jcl_rmax=0.0;
	for(iv=0;iv<nrntnz;iv++)
	{
		if(is_sol[iv])
		{
			iz=iv/nrnt;	//Calcolo coordinate
			it=(iv%nrnt)/nr;
			ir=iv-(iz*nrnt)-(it*nr);
			num_sol++;	//Incremento numero elementi in soluzione
			//Errore massimo e medio
			Jerrtmp=abs(Jkbk[iv]+Jkfw[iv]+Jkdx[iv]+Jksx[iv]+Jkdw[iv]+Jkup[iv]);
			Jerrk_med+=Jerrtmp;
			if(Jerrtmp>Jerrk_max)
			{
				Jerrk_max=Jerrtmp;
				iz_errk=iz;
				it_errk=it;
				ir_errk=ir;
			}
			Jerrtmp=abs(Jclbk[iv]+Jclfw[iv]+Jcldx[iv]+Jclsx[iv]+Jcldw[iv]+Jclup[iv]);
			Jerrcl_med+=Jerrtmp;
			if(Jerrtmp>Jerrcl_max)
			{
				Jerrcl_max=Jerrtmp;
				iz_errcl=iz;
				it_errcl=it;
				ir_errcl=ir;
			}
			//Flusso assiale
			Jslicek[iz]+=Jkup[iv];
			Jslicecl[iz]+=Jclup[iv];
			//Flussi massimi
			if(abs(Jkfw[iv]-Jkbk[iv])>Jk_rmax) Jk_rmax=abs(Jkfw[iv]-Jkbk[iv]);
			if(abs(Jclfw[iv]-Jclbk[iv])>Jcl_rmax) Jcl_rmax=abs(Jclfw[iv]-Jclbk[iv]);
			if(abs(Jksx[iv]-Jkdx[iv])>Jk_tmax) Jk_tmax=abs(Jksx[iv]-Jkdx[iv]);
			if(abs(Jclsx[iv]-Jcldx[iv])>Jcl_tmax) Jcl_tmax=abs(Jclsx[iv]-Jcldx[iv]);
			if(abs(Jkup[iv]-Jkdw[iv])>Jk_zmax) Jk_zmax=abs(Jkup[iv]-Jkdw[iv]);
			if(abs(Jclup[iv]-Jcldw[iv])>Jcl_zmax) Jcl_zmax=abs(Jclup[iv]-Jcldw[iv]);
		}
	}
	if(Jk_rmax==0.0)Jk_rmax=1.0;
	if(Jk_tmax==0.0)Jk_tmax=1.0;
	if(Jk_zmax==0.0)Jk_zmax=1.0;
	if(Jcl_rmax==0.0)Jcl_rmax=1.0;
	if(Jcl_tmax==0.0)Jcl_tmax=1.0;
	if(Jcl_zmax==0.0)Jcl_zmax=1.0;
	//Errore medio
	Jerrk_med=Jerrk_med/num_sol;
	Jerrcl_med=Jerrcl_med/num_sol;
	Jaxisk_max=Jaxiscl_max=-1e-40;
	Jaxisk_min=Jaxiscl_min=1e-40;
	iz_axisk_max=iz_axisk_min=iz_axiscl_max=iz_axiscl_min=-1;
	Jaxisk_med=Jaxiscl_med=0.0;
	//FLusso assiale medio, massimo, minimo
	for(iz=0;iz<nz;iz++)
	{
		Jaxisk_med+=Jslicek[iz];
		if(Jslicek[iz]>Jaxisk_max)
		{
			Jaxisk_max=Jslicek[iz];
			iz_axisk_max=iz;
		}
		if(Jslicek[iz]<Jaxisk_min)
		{
			Jaxisk_min=Jslicek[iz];
			iz_axisk_min=iz;
		}
		Jaxiscl_med+=Jslicecl[iz];
		if(Jslicecl[iz]>Jaxiscl_max)
		{
			Jaxiscl_max=Jslicecl[iz];
			iz_axiscl_max=iz;
		}
		if(Jslicecl[iz]<Jaxiscl_min)
		{
			Jaxiscl_min=Jslicecl[iz];
			iz_axiscl_min=iz;
		}
	}
	//FLusso medio
	Jaxisk_med=Jaxisk_med/nz;
	Jaxiscl_med=Jaxiscl_med/nz;
	rmsdk=rmsdcl=0.0;
	for(iz=0;iz<nz;iz++)
	{
		rmsdk+=((Jaxisk_med-Jslicek[iz])*(Jaxisk_med-Jslicek[iz]));
		rmsdcl+=((Jaxiscl_med-Jslicecl[iz])*(Jaxiscl_med-Jslicecl[iz]));
	}
	rmsdk=sqrt(rmsdk/nz);
	rmsdcl=sqrt(rmsdcl/nz);
	//Correnti
	Ik=mmolns2pA*Jaxisk_med;
	Icl=-mmolns2pA*Jaxiscl_med;
	Itot=Ik+Icl;

	//Output: Valori medii ed errori
	cout<<setw(25)<<"\t"<<setw(15)<<"CATIONS "<<" "<<setw(15)<<""<<"\t"<<setw(15)<<"ANIONS"<<"\n";
	cout<<setw(25)<<"MAX ERROR [A^3mmol/ns]\t"<<setw(15)<<Jerrk_max<<" "<<setw(5)
		<<iz_errk<<setw(5)<<it_errk<<setw(5)<<ir_errk<<"\t"
		<<setw(15)<<Jerrcl_max<<setw(5)<<iz_errcl<<setw(5)<<it_errcl<<setw(5)<<ir_errcl<<"\n";
	cout<<setw(25)<<"MEAN ERROR [A^3mmol/ns]\t"<<setw(15)<<Jerrk_med<<setw(15)<<"\t"
		<<setw(15)<<Jerrcl_med<<"\n";
	cout<<setw(25)<<"RMSD FLUX [A^3mmol/ns]\t"<<setw(15)<<rmsdk<<setw(15)<<"\t"<<setw(15)<<rmsdcl<<"\n";
	cout<<setw(25)<<"MAX FLUX [A^3mmol/ns]\t"<<setw(15)<<Jaxisk_max<<" "
		<<setw(5)<<" "<<setw(5)<<" "<<setw(5)<<iz_axisk_max<<"\t"
		<<setw(15)<<Jaxiscl_max<<setw(5)<<" "<<setw(5)<<" "<<setw(5)<<iz_axiscl_max<<"\n";
	cout<<setw(25)<<"MIN FLUX [A^3mmol/ns]\t"<<setw(15)<<Jaxisk_min<<" "
		<<setw(5)<<" "<<setw(5)<<" "<<setw(5)<<iz_axisk_min<<"\t"
		<<setw(15)<<Jaxiscl_min<<setw(5)<<" "<<setw(5)<<" "<<setw(5)<<iz_axiscl_min<<"\n";
	cout<<setw(25)<<"AXIAL FLUX [A^3mmol/ns]\t"<<setw(15)<<Jaxisk_med<<setw(15)<<"\t"
		<<setw(15)<<Jaxiscl_med<<"\n\n";
	cout<<setw(25)<<"CURRENT [pA]\t"<<setw(15)<<Ik<<setw(5)<<""<<setw(5)<<" (+) "<<setw(5)<<""<<"\t"
		<<setw(15)<<Icl<<" = "<<setw(10)<<Itot<<" [pA]\n";
	cout<<setw(50)<<""<<setw(5)<<"   (/)"<<setw(5)<<""<<"\t"
		<<setw(15)<<""<<" = "<<setw(10)<<Ik/Icl<<"\n";
	cout<<"#"<<setw(14)<<"Vmembrane [mV]"<<"\t"<<setw(15)<<" I [pA]"<<"\t"<<"\n";
	cout<<setw(15)<<phimem<<"\t"<<setw(15)<<Itot<<"\n";



	//Output: Flussi Su File
	if (strcmp(fileflx,"NEPAL.flx")!=0)
	{
		file_flx.open(fileflx);
		if (!file_flx)
		{
			cerr<<"ERROR IT IS NOT POSSIBLE TO OPEN "<<fileflx<<"\n";
		  	exit(1);
		}
		file_flx<<"#FLUX FIELD\n";
		file_flx<<"#Usage: splot 'NEPAL.flx' index 1:nrntnz  using 1:2:3 (using 4:5:6 anion) w l\n";
		iv=0;
		for(iz=0;iz<nz;iz++)
		{
			z=zmin+(dz/2.0)+iz*dz;
			for(it=0;it<nt;it++)
			{
				t=(dt/2.0)+it*dt;
				for(ir=0;ir<nr;ir++)
				{
					r=(dr/2.0)+ir*dr;
					if(is_sol[iv])
					{
				x=r*cos(t);
				y=r*sin(t);
				//Punto inizio vettori
				file_flx<<setw(10)<<x<<"  "<<setw(10)<<y<<"  "<<setw(10)<<z<<"  ";
				file_flx<<setw(10)<<x<<"  "<<setw(10)<<y<<"  "<<setw(10)<<z<<"\n";
				J_z=Jkup[iv]-Jkdw[iv];
				J_z=(dz/2.0)*(J_z/Jk_zmax);
				J_t=Jksx[iv]-Jkdx[iv];
				J_t=(dr/2.0)*dt*(J_t/Jk_tmax);
				J_r=Jkfw[iv]-Jkbk[iv];
				J_r=(dr/2.0)*(J_r/Jk_rmax);
				J_x=J_r*cos(t)-J_t*sin(t);
				J_y=J_r*sin(t)+J_t*cos(t);
				//Punto fine vettore: campo vettoriale Jk
				file_flx<<setw(10)<<x+J_x<<"  "<<setw(10)<<y+J_y<<"  "<<setw(10)<<z+J_z<<" ";
				J_z=Jclup[iv]-Jcldw[iv];
				J_z=(dz/2.0)*(J_z/Jcl_zmax);
				J_t=Jclsx[iv]-Jcldx[iv];
				J_t=(dr/2.0)*dt*(J_t/Jcl_tmax);
				J_r=Jclfw[iv]-Jclbk[iv];
				J_r=(dr/2.0)*(J_r/Jcl_rmax);
				J_x=J_r*cos(t)-J_t*sin(t);
				J_y=J_r*sin(t)+J_t*cos(t);
				//Punto fine vettore: campo vettoriale Jcl
				file_flx<<setw(10)<<x+J_x<<"  "<<setw(10)<<y+J_y<<"  "<<setw(10)<<z+J_z<<"\n";
				file_flx<<"\n\n";
					}
					iv++;
				}
			}
		}
		file_flx<<"#FLUXES THROUGH THE SLICES\n";
		file_flx<<"#Usage: plot 'NEPAL.flx' index 0 using 1:2 (1:3 anion) w l\n";
		file_flx<<"#"<<setw(9)<<"z"<<"  "<<setw(10)<<"Jslicek"<<"  "<<setw(10)<<"Jslicecl"<<"\n";
		for(iz=0;iz<nz;iz++)
		{
		   z=zmin+(dz/2.0)+iz*dz;
		   file_flx<<setw(10)<<z<<"  "<<setw(10)<<Jslicek[iz]<<"  "<<setw(10)<<Jslicecl[iz]<<"\n";
		}
		file_flx<<"\n\n";
		file_flx.close();
	}

	delete[] Jslicek;
	delete[] Jslicecl;
	delete[] Jkbk;
	delete[] Jkfw;
	delete[] Jkdx;
	delete[] Jksx;
	delete[] Jkdw;
	delete[] Jkup;
	delete[] Jclbk;
	delete[] Jclfw;
	delete[] Jcldx;
	delete[] Jclsx;
	delete[] Jcldw;
	delete[] Jclup;
	return 0;
}

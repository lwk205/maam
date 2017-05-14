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
#include <src/nepal/global.h>
#include <src/nepal/functions.h>

////////////////////////////////////////////////////////////////////////////////
//SHOWPRM
//Esegue l'output dei parametri di esecuzione
int showprm(const molecule &chn,const int min_cell4atom, const int max_cell4atom)
{
  double chgmol;
  unsigned long int iv,num_sol;
  int iz;

  cout<<"INPUT FILES\n";
  cout<<"file.pdb = "<<filepdb<<"\n";
  cout<<"file.str = "<<filestr<<"\n";
  cout<<"\n";
  cout<<"OUTPUT FILES\n";
  cout<<"file.grd = "<<filegrd<<"\n";
  cout<<"file.dat = "<<filedat<<"\n";
  cout<<"file.flx = "<<fileflx<<"\n";
  cout<<"file.end = "<<fileend<<"\n";
  cout<<"\n";
  cout<<"INPUT/OUTPUT CONTROL\n";
  cout<<"Verbose level : "<<verbose<<"\n";
  cout<<"\n";
  cout<<"GRID PARAMETERS\n";
  cout<<"Number of grid elements = "<<nrntnz<<"\n";
  num_sol=0; for(iv=0;iv<nrntnz;iv++)if(is_sol[iv])num_sol++;
  cout<<"Number of solvent grid elements = "<<num_sol<<"\n";
  cout<<"nstep r = "<<nr<<"\tnstep teta = "<<nt<<"\tnstep z = "<<nz<<"\n";
  cout<<"rmax = "<<rmax<<" [A]"<<"\tDelta r = "<<dr<<" [A]"<<"\n";
  cout<<"Delta teta = "<<dt<<" [rad]"<<" Arc max = "<<dr*dt*nr<<" [A]\n";
  cout<<"zmin = "<<zmin<<" [A]"<<"\tzmax = "<<zmax<<" [A]"<<"\tDelta z = "<<dz<<" [A]"<<"\n";
  switch(disc_chg){
	  case 0:
		  cout<<"Charge discretizazion method : uniform on atom volume\n";
	 	  break;
	  case 1:
		  cout<<"Charge discretizazion method : distance dependent\n";
	 	  break;
	  case 2:
		  cout<<"Charge discretizazion method : concentrated in 1 grid element\n";
	 	  break;
  }
  cout<<"Min grid element volume = "<<volume_cell(0,dr,dt,dz)<<" [A^3] = "
	  <<volume_cell(0,dr,dt,dz)/((4.0*PI/3.0))<<" [1A sphere]\n";
  cout<<"Max grid element volume = "<<volume_cell((nr-1),dr,dt,dz)<<" [A^3] = "
	  <<volume_cell((nr-1),dr,dt,dz)/((4.0*PI/3.0))<<" [1A sphere]\n";
  cout<<"Min number of cells used to discretize an atom = "<<min_cell4atom<<"\n";
  cout<<"Max number of cells used to discretize an atom = "<<max_cell4atom<<"\n";
  cout<<"\n";
  cout<<"MEMBRANE PARAMETERS\n";
  cout<<"Membrane lower limit = "<<zmem_min<<" [A]"<<"\tizmem_min = "<<izmem_min<<"\n";
  cout<<"Membrane upper limit = "<<zmem_max<<" [A]"<<"\tizmem_max = "<<izmem_max<<"\n";
  cout<<"\n";
  cout<<"BOUNDARY CONDITIONS\n";
  cout<<"phi[0] = "<<phibnd[0]<<" [Vt]\tphi["<<(nz/2)<<"] = "
	  <<phibnd[(nz/2)]<<" [Vt]\tphi["<<nz-1<<"] = "<<phibnd[nz-1]<<" [Vt]\n";
  cout<<"kext = "<<kext<<" [mmol]"<<"\tkint = "<<kint<<" [mmol]\n";
  cout<<"clext = "<<clext<<" [mmol]"<<"\tclint = "<<clint<<" [mmol]\n";
  cout<<"\n";
  cout<<"PHYSICAL PARAMETERS\n";
  cout<<"Termic potential = "<<phiterm<<" [mV]\n";
  cout<<"epsH2O = "<<epsH2O<<"\tepsMOL = "<<epsMOL<<"\n";
  cout<<"Exclusion layer thickness = "<<rad_probe<<" [A]\n";
  cout<<"Cation Radius = "<<rad_k<<" [A]     Anion Radius = "<<rad_cl<<" [A]\n";
  cout<<"Boundary for Paine-Sheer (if used) = "<<zminpaine<<" "<<zmaxpaine<<"\n";
  cout<<"Dk = "<<setw(10)<<Dk[0]<<" [A^2/ns]"<<"\t\t"<<setw(10)<<zmin<<" [A]";
  for(iz=1;iz<nz;iz++){
	  if (Dk[iz]!=Dk[iz-1]){
		  cout<<"\t"<<setw(10)<<zmin+(iz)*dz<<" [A]"<<"\n";
  		  cout<<"Dk = "<<setw(10)<<Dk[iz]<<" [A^2/ns]"<<"\t\t"<<setw(10)<<zmin+iz*dz<<" [A]";
	  }
  }
  cout<<"\t"<<setw(10)<<zmax<<" [A]"<<"\n";
  cout<<"Dcl = "<<setw(10)<<Dcl[0]<<" [A^2/ns]"<<"\t\t"<<setw(10)<<zmin<<" [A]";
  for(iz=1;iz<nz;iz++)
  {
	  if (Dcl[iz]!=Dcl[iz-1])
	  {
		  cout<<"\t"<<setw(10)<<zmin+(iz)*dz<<" [A]"<<"\n";
  		  cout<<"Dcl = "<<setw(10)<<Dcl[iz]<<" [A^2/ns]"<<"\t\t"<<setw(10)<<zmin+iz*dz<<" [A]";
	  }
  }
  cout<<"\t"<<setw(10)<<zmax<<" [A]"<<"\n";
  cout<<"\n";
  cout<<"MOLECULE PARAMETERS\n";
  cout<<"Number of atoms = "<<chn.num_atm<<"\n";
  chgmol=0;
  for(iv=0;iv<nrntnz;iv++)chgmol+=chg_mol[iv];
  cout<<"Total charge = "<<chn.chg<<" [e] (discretized on the grid = "
	  <<chgmol/e2mVA<<" [e])\n";
  cout<<"Maximum radius = "<<chn.rmax<<"\n";
  cout<<"zmin = "<<chn.zmin<<" [A] zmax = "<<chn.zmax<<" [A]\n";
  cout<<"\n";
  cout<<"ITERATION CONTROL\n";
  cout<<"Poisson-Nernst-Planck equations will be solved\n";
  if(lin_exp==0) cout<<"Nernst equation will be solved in a linear space\n";
  else cout<<"Nernst equation will be solved in an exponential space\n";
  cout<<"Max iterations = "<<max_it<<"\n";
  cout<<"Tollerance = "<<tol<<"\n";
  cout<<"Max iteration Poisson = "<<max_itP<<"\n";
  cout<<"Tollerance Poisson = "<<tol_P<<"\n";
  cout<<"Max iteration Nernst = "<<max_itN<<"\n";
  cout<<"Tollerance Nernst = "<<tol_N<<"\n";
  cout<<"Restart files written every "<<it_restart<<" steps\n";
  switch(smooth){
	  case 0:
		  cout<<"No smoothing (wP = "<<wP<<" and wN = "<<wN<<")\n";
		  break;
	  case 1:
		  cout<<"Smoothing by wP = "<<wP<<" and  wN = "<<wN<<"\n";
		  break;
	  case 2:
		  cout<<"Smoothing by wP and wN defined according to max allowed deviation at the first iteration\n";
		  break;
	  case 3:
		  cout<<"Smoothing by wP = "<<wP<<" and  wN = "<<wN<<" and iteration updating\n";
		  break;
	  default:
		  cerr<<"ERROR: wrong kind of smoothing: "<<smooth<<"\n";
		  exit(1);
  }
  cout<<"\n";
  cout.flush();

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
//WRITE_VOL[CHG/CHN]
//
//Esegue l'output per utilizzo in VMD
////////////////////////////////////////////////////////////////////////////////
void start_outvol(ofstream &file_out, double &xmax, int &nx, double &dx, int new_nz=nz, int iz_start=0){
  xmax=rmax*((sqrt(2.0))/2.0);
  nx=(2*((int)(nr*((sqrt(2.0))/2.0))) + 1);
  dx=(2.0*xmax)/(nx-1);
  file_out<<setiosflags(ios::fixed);
  file_out<<setprecision(3);
  file_out<<"object 1 class gridpositions counts "<<new_nz<<" "<<nx<<" "<<nx<<endl;
  file_out<<"origin "<<-xmax<<" "<<-xmax<<" "<<zmin+(iz_start*dz)+(dz/2.0)<<endl;
  file_out<<"delta 0.000000e+00 0.000000e+00 "<<dz<<endl;
  file_out<<"delta 0.000000e+00 "<<dx<<" 0.000000e+00"<<endl;
  file_out<<"delta "<<dx<<" 0.000000e+00 0.000000e+00"<<endl;
  file_out<<"object 2 class gridconnections counts "<<new_nz<<" "<<nx<<" "<<nx<<endl;
  file_out<<"object 3 class array type double rank 0 items "<<nx*nx*new_nz<<" data follows"<<endl;
  return;
}
void cln2crt(double x, double y, int iz,  unsigned long int &iv){
  double r_grd,t_grd;
  int ir_cln, it_cln;
  crt2cln(x,y,r_grd,t_grd);
  if(t_grd < 0.0)t_grd+=2*PI;
  ir_cln=(int)(r_grd/dr);
  it_cln=(int)(t_grd/dt);
  iv=ir_cln+it_cln*nr+iz*nrnt;
  return;
}
void end_outvol(ofstream &file_out){
  file_out<<endl;
  file_out<<"attribute \"dep\" string \"positions\""<<endl;
  file_out<<"object \"regular positions regular connections\" class field"<<endl;
  file_out<<"component \"positions\" value 1"<<endl;
  file_out<<"component \"connections\" value 2"<<endl;
  file_out<<"component \"data\" value 3"<<endl;
  file_out<<endl;
  return;
}

int write_volchg(ofstream &file_out){
  double xmax,dx,x,y;
  int nx,ix,iy,iz,dummy_cnt;
  unsigned long int iv;
	
  start_outvol(file_out,xmax,nx,dx);
  dummy_cnt=1;
  for(iz=0;iz<nz;iz++){
  	for(iy=0;iy<nx;iy++){
  		y=-xmax+(iy*dx);
  		for(ix=0;ix<nx;ix++){
  			x=-xmax+(ix*dx);
			cln2crt(x,y,iz,iv);
			//Output carica unita' di misura cariche parziali [e/1000]
  			file_out<<setw(8)<<1000*(chg_mol[iv]/e2mVA)<<"\t";
  			if((dummy_cnt%3)==0)file_out<<endl;
  			dummy_cnt++;
  		}
  	}
  }
  end_outvol(file_out);
  file_out.flush();
  return 0;
}

int write_volchn(ofstream &file_out_profile,ofstream &file_out_radius){
  double xmax,dx,x,y,z,t,rint_med,rint,rint_tmp,rext_med,rext,rext_tmp,rad,dtmp;
  int nx,ix,iy,iz,it,dummy_cnt;
	
  bool started_chn,ended_chn;
  int iz_start,iz_end;

  started_chn=ended_chn=false;
  rint_med=rext_med=0.0;
  rint=dr*nr;
  rext=0.0;
  for(iz=0;iz<nz;iz++){
	z=zmin+(dz/2.0)+iz*dz;
	rint_med=0.0;
	rint=dr*nr;
	for(it=0;it<nt;it++)
	{
		t=(dt/2.0)+it*dt;
		rint_tmp=dr*(inside_end[iz*nt+it]);
		rint_med+=rint_tmp;
		if(rint_tmp<rint)rint=rint_tmp;
	}
	rint_med=rint_med/(double)nt;
	if((!started_chn)&&(rint<(rmax-dr))){
		iz_start=iz;
		started_chn=true;
	}
	if((started_chn)&&(rint>(rmax-dr))&&(!ended_chn)){
		iz_end=iz;
		ended_chn=true;
	}
  }

  dummy_cnt=1;
  start_outvol(file_out_profile,xmax,nx,dx,iz_end-iz_start,iz_start);
  start_outvol(file_out_radius,xmax,nx,dx,iz_end-iz_start,iz_start);
  for(iz=iz_start;iz<iz_end;iz++){
  //start_outvol(file_out_profile,xmax,nx,dx);
  //start_outvol(file_out_radius,xmax,nx,dx);
  //for(iz=0;iz<nz;iz++){
	z=zmin+(dz/2.0)+iz*dz;
	rint_med=rext_med=0.0;
	rint=dr*nr;
	rext=0.0;
	for(it=0;it<nt;it++)
	{
		t=(dt/2.0)+it*dt;
		rint_tmp=dr*(inside_end[iz*nt+it]);
		rext_tmp=dr*(outside_begin[iz*nt+it]);
		rint_med+=rint_tmp;
		rext_med+=rext_tmp;
		if(rint_tmp<rint)rint=rint_tmp;
		if(rext_tmp>rext)rext=rext_tmp;
	}
	rint_med=rint_med/(double)nt;
	rext_med=rext_med/(double)nt;
	cout<<"rint = "<<rint<<" rext = "<<rext<<" rmax = "<<rmax<<endl;
  	for(iy=0;iy<nx;iy++){
  		y=-xmax+(iy*dx);
  		for(ix=0;ix<nx;ix++){
  			x=-xmax+(ix*dx);
			rad=sqrt(x*x+y*y);
			//if((rad<rint)||(rad>rext_med)||(rint==rext))dtmp=0;
			if((rad<rint))dtmp=0;
			else dtmp=1;
  			file_out_profile<<setw(8)<<dtmp<<"\t";
  			file_out_radius<<setw(8)<<rint<<"\t";
  			if((dummy_cnt%3)==0){
				file_out_profile<<endl;
				file_out_radius<<endl;
			}
  			dummy_cnt++;
  		}
  	}
  }
  end_outvol(file_out_profile);
  end_outvol(file_out_radius);
  file_out_profile.flush();
  file_out_radius.flush();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//WRITEMOLECULE
//
//Esegue l'output della molecola discretizzata sulla griglia
////////////////////////////////////////////////////////////////////////////////
int writemolecule()
{
	int it,iz;
	double t,z,rint_tmp,rext_tmp;
	double rint,rext,rint_med,rext_med,xint,yint,xext,yext;
	double r_circle;//Raggio circonferenza con area equivalente a quella della sezione
	double Dk_med,Dk_circle;//Diffusione secondo Paine-Scherr

	file_grd<<setiosflags(ios::fixed);
	file_grd<<setprecision(4);
	file_grd<<"#Internal-External radius of molecule: "<<filepdb
		<<" discretized on a grid : "<<nr<<"x"<<nt<<"x"<<nz<<endl;
	file_grd<<"#Usage internal radius: plot '"<<filegrd
		<<"' index 0 using 1:(2 for min value)(4 for mead values) w l"<<endl;
	file_grd<<"#Usage external radius: plot '"<<filegrd
		<<"' index 0 using 1:(3 for min value)(5 for mead values) w l"<<endl;
	file_grd<<"#Usage diffusion coefficient: plot '"<<filegrd
		<<"' index 0 using 1:(6 for Dk[iz])(7 for Dcl[iz]) w l"<<endl;
	file_grd<<"#"<<setw(9)<<"z"<<setw(10)<<"rint"<<setw(10)<<"rext"
	        <<setw(10)<<"rint_med"<<setw(10)<<"rext_med"
	        <<setw(10)<<"Dk"<<setw(10)<<"Dcl"<<endl;
	for(iz=0;iz<nz;iz++)
	{
		z=zmin+(dz/2.0)+iz*dz;
		rint_med=rext_med=r_circle=0.0;
		rint=dr*nr;
		rext=0;
		for(it=0;it<nt;it++)
		{
			t=(dt/2.0)+it*dt;
			rint_tmp=dr*(inside_end[iz*nt+it]);
			rext_tmp=dr*(outside_begin[iz*nt+it]);
			rint_med+=rint_tmp;
			rext_med+=rext_tmp;
			//Aggiunge area spicchio
			r_circle+=(PI*rint_tmp*rint_tmp)/(double)nt;
			if(rint_tmp<rint)rint=rint_tmp;
			if(rext_tmp>rext)rext=rext_tmp;
		}
		rint_med=rint_med/(double)nt;
		rext_med=rext_med/(double)nt;
		r_circle=sqrt(r_circle/PI);
		file_grd<<setw(10)<<z<<setw(10)<<rint<<setw(10)<<rext
		        <<setw(10)<<rint_med<<setw(10)<<rext_med
			<<setw(10)<<Dk[iz]<<setw(10)<<Dcl[iz]<<endl;
			
	}
	file_grd<<"\n";
	file_grd<<"\n";

	file_grd<<"#Molecule "<<filepdb
		<<" discretized on a grid : "<<nr<<"x"<<nt<<"x"<<nz<<"\n";
	file_grd<<"#Usage: splot '"<<filegrd
		<<"' index 1 using 1:2:5 w l,'' using 3:4:5 w l\n";
	file_grd<<"#"<<setw(9)<<"xint"<<"\t"<<setw(10)<<"yint"<<"\t"<<setw(10)
		<<"xext"<<"\t"<<setw(10)<<"yext"<<"\t"<<setw(10)<<"z"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		z=zmin+(dz/2.0)+iz*dz;
		for(it=0;it<nt;it++)
		{
			t=(dt/2.0)+it*dt;
			rint=dr*(inside_end[iz*nt+it]);
			rext=dr*(outside_begin[iz*nt+it]);
			xint=rint*cos(t);
			yint=rint*sin(t);
			xext=rext*cos(t);
			yext=rext*sin(t);
			file_grd<<setw(10)<<xint<<"\t"<<setw(10)<<yint<<"\t"
				<<setw(10)<<xext<<"\t"<<setw(10)<<yext<<"\t"
				<<setw(10)<<z<<"\n";
		}
		t=(dt/2.0);
		rint=dr*(inside_end[iz*nt]);
		rext=dr*(outside_begin[iz*nt]);
		xint=rint*cos(t);
		yint=rint*sin(t);
		xext=rext*cos(t);
		yext=rext*sin(t);
		file_grd<<setw(10)<<xint<<"\t"<<setw(10)<<yint<<"\t"
			<<setw(10)<<xext<<"\t"<<setw(10)<<yext<<"\t"
			<<setw(10)<<z<<"\n";
		file_grd<<"\n";
	}

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
//WRITEDATA
//
//Esegue l'output di potenziale e concentrazioni ioniche
////////////////////////////////////////////////////////////////////////////////
int writedata(double *phi, double *k, double *cl)
{
	int ir,it,iz;
	double r,t,z,x,y,r_med,phi_med,k_med,cl_med;
	unsigned long int iv,iv_start,iv_end;

	file_dat.open(filedat);
 	if (!file_dat){
	  	cerr<<"ERROR IT IS NOT POSSIBLE TO OPEN "<<filedat<<"\n";
	  	exit(1);
	}

//Dati lungo asse centrale
	file_dat<<"#Data along the channel axis ir=0 it=0 iz=0:nz-1\n";
	file_dat<<"#Usage-Potential: plot 'tmp.dat' index 0 using 1:2 w l\n";
	file_dat<<"#Usage-K concentration: plot 'tmp.dat' index 0 using 1:3 w l\n";
	file_dat<<"#Usage-Cl concentration: plot 'tmp.dat' index 0 using 1:4 w l\n";
	file_dat<<"#z\t\tphi [mV]\t\tk [mmol]\t\tcl [mmol]\n";
	ir=0;
	it=0;
	for(iz=0;iz<nz;iz++)
	{
		z=zmin+(dz/2.0)+iz*dz;
		iv=(iz*nrnt)+(it*nr)+ir;
		file_dat<<setw(10)<<z<<'\t'<<setw(10)<<phi[iv]<<'\t'<<setw(10)
			<<k[iv]<<'\t'<<setw(10)<<cl[iv]<<'\n';

	}
	file_dat<<"\n\n";

//Dati alle diverse fette	
	if(verbose>1){
		for(iz=0;iz<nz;iz++){
			file_dat<<"#Data on the different slices\n";
			file_dat<<"#Usage-Potential: splot 'tmp.dat' index iz+1 using 1:2:3 w l\n";
			file_dat<<"#Usage-K concentration: splot 'tmp.dat' index iz+1 using 1:2:4 w l\n";
			file_dat<<"#Usage-Cl concentration: splot 'tmp.dat' index iz+1 using 1:2:5 w l\n";
			file_dat<<"#"<<setw(9)<<"x"<<setw(12)<<"y"
				<<setw(12)<<"phi[iv]"<<setw(12)<<"k[iv]"<<setw(12)<<"cl[iv]"<<"\n";
			iv=iz*nrnt;
			it=0;
			ir=0;
			iv_start=nrnt*iz;
			iv_end=nrnt*(iz+1);
			for(iv=iv_start;iv<iv_end;iv++)
			{
				t=(dt/2.0)+it*dt;
				r=(dr/2.0)+ir*dr;
				x=r*cos(t);
				y=r*sin(t);
				file_dat<<setw(10)<<x<<"  "<<setw(10)<<y<<"  "
					<<setw(10)<<phi[iv]<<"  "<<setw(10)<<k[iv]<<"  "<<setw(10)<<cl[iv]<<"\n";
				ir++;
				if(ir==nr)
				{
					ir=0;
					it++;
				}
			}
			file_dat<<"\n\n";
		}
	}

//Dati mediati di superficie interna	
	file_dat<<"#Mean data\n";
	file_dat<<"#It gives on a surface with radius equal to protein mean radius\n";
	file_dat<<"#the mean potential and concetrations at the boundaries with the protein\n";
	file_dat<<"#Usage-Potential: splot 'file_dat' index ind using 1:2:3:4 w pm3d\n";
	file_dat<<"#Usage-K concentration: splot 'file_dat' index ind using 1:2:3:5 w pm3d\n";
	file_dat<<"#Usage-Cl concentration:  splot 'file_dat' index ind using 1:2:3:6 w pm3d\n";
	file_dat<<"#x\t\ty\t\tz\t\tphi [mV]\t\tk [mmol]\t\tcl [mmol]\n";
	for(iz=0;iz<nz;iz++)
	{
		z=zmin+(dz/2.0)+iz*dz;
		//Calcolo raggio medio e dati medii
		r_med=0;
		phi_med=0;
		k_med=0;
		cl_med=0;
		for(it=0;it<nt;it++)
		{
			ir=inside_end[iz*nt+it];
			r_med+=(dr*ir);
			iv=(iz*nrnt)+(it*nr)+(ir-1);
			phi_med+=phi[iv];
			k_med+=k[iv];
			cl_med+=cl[iv];
		}
		r_med=r_med/nt;
		phi_med=phi_med/nt;
		k_med=k_med/nt;
		cl_med=cl_med/nt;
		for(it=0;it<(nt/2);it++)
		{
			t=(dt/2.0)+it*dt;
			x=r_med*cos(t);
			y=r_med*sin(t);
			file_dat<<setw(10)<<x<<"\t"<<setw(10)<<y<<"\t"<<setw(10)<<z<<"\t"
				<<setw(10)<<phi_med<<"\t"<<setw(10)<<k_med<<'\t'<<setw(10)<<cl_med<<"\n";
		}
		file_dat<<"\n";
	}
	file_dat<<"\n\n";

	file_dat.close();
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
//WRITEBACKUP
//
//Esegue l'output dei dati per una successiva ripresa del calcolo
////////////////////////////////////////////////////////////////////////////////
int writebackup(double *phi, double *k, double *cl){
	unsigned long int iv;
	int iz;

  	for(iv=0;iv<nrntnz;iv++)phi[iv]=(phi[iv]*phiterm);
  	for(iz=0;iz<nz;iz++)phibnd[iz]=(phibnd[iz]*phiterm);
  	phimem=phimem*phiterm;
	fluxes(k,cl,phi);
	writedata(phi,k,cl);
  	file_end.open(fileend,ios::binary);
  	if (!file_end){
  	  cerr<<"ERROR it is not possible to write "<<fileend<<"\n";
  	  return 1;
  	}
	writerestart(phi,k,cl);
	file_end.close();
  	for(iv=0;iv<nrntnz;iv++)phi[iv]=(phi[iv]/phiterm);
  	for(iz=0;iz<nz;iz++)phibnd[iz]=(phibnd[iz]/phiterm);
  	phimem=phimem/phiterm;
	return 1;
}
////////////////////////////////////////////////////////////////////////////////
//WRITERESTART
//
//Esegue l'output dei dati per una successiva ripresa del calcolo
////////////////////////////////////////////////////////////////////////////////
int writerestart(double *phi, double *k, double *cl)
{
	//Inserisco nel file binario la dimensione della griglia per avere poi un
	//controllo in fase di lettura
	file_end.write((char *)&nr,sizeof(int));
	file_end.write((char *)&nt,sizeof(int));
	file_end.write((char *)&nz,sizeof(int));
	//Scrittura nel file binario di potenziale e concentrazioni ioniche
	file_end.write((char *)phi,nrntnz*sizeof(double));
	file_end.write((char *)k,nrntnz*sizeof(double));
	file_end.write((char *)cl,nrntnz*sizeof(double));
	//Parametri griglia
	file_end.write((char *)&zmin,sizeof(double));
	file_end.write((char *)&zmax,sizeof(double));
	file_end.write((char *)&rmax,sizeof(double));
	//Parametri fisici
	file_end.write((char *)&rad_k,sizeof(double));
	file_end.write((char *)&rad_cl,sizeof(double));
	file_end.write((char *)Dk,nz*sizeof(double));
	file_end.write((char *)Dcl,nz*sizeof(double));

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
//OUTPUTDS
//
//Esegue l'output di matrici coefficienti e vettori termini noti
///////////////////////////////////////////////////////////////////////////////
int outputDS()
{
	int iz,it,ir;
	unsigned long int iv,im;

	im=0;
	cout<<"A_PHI\n";
	cout<<setw(5)<<"iz"<<setw(5)<<"it"<<setw(5)<<"ir"
		<<setw(8)<<"bk"<<setw(8)<<"fw"<<setw(8)<<"dx"
		<<setw(8)<<"sx"<<setw(8)<<"dw"<<setw(8)<<"up"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
					<<setw(8)<<A_phi[im]<<setw(8)<<A_phi[im+1]<<setw(8)<<A_phi[im+2]
					<<setw(8)<<A_phi[im+3]<<setw(8)<<A_phi[im+4]<<setw(8)<<A_phi[im+5]<<"\n";
				im=im+6;
			}
		}
	}
	cout<<"\n";
	im=0;
	cout<<"ACS_K and A_K\n";
	cout<<setw(5)<<"iz"<<setw(5)<<"it"<<setw(5)<<"ir"
		<<setw(8)<<"bk"<<setw(8)<<"fw"<<setw(8)<<"dx"
		<<setw(8)<<"sx"<<setw(8)<<"dw"<<setw(8)<<"up"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
					<<setw(8)<<Acs_k[im]<<setw(8)<<Acs_k[im+1]<<setw(8)<<Acs_k[im+2]
					<<setw(8)<<Acs_k[im+3]<<setw(8)<<Acs_k[im+4]<<setw(8)<<Acs_k[im+5]<<"\n";
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
					<<setw(8)<<A_k[im]<<setw(8)<<A_k[im+1]<<setw(8)<<A_k[im+2]
					<<setw(8)<<A_k[im+3]<<setw(8)<<A_k[im+4]<<setw(8)<<A_k[im+5]<<"\n";
				im=im+6;
			}
		}
	}
	cout<<"\n";
	im=0;
	cout<<"ACS_CL and A_CL\n";
	cout<<setw(5)<<"iz"<<setw(5)<<"it"<<setw(5)<<"ir"
		<<setw(8)<<"bk"<<setw(8)<<"fw"<<setw(8)<<"dx"
		<<setw(8)<<"sx"<<setw(8)<<"dw"<<setw(8)<<"up"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
				<<setw(8)<<Acs_cl[im]<<setw(8)<<Acs_cl[im+1]<<setw(8)<<Acs_cl[im+2]
				<<setw(8)<<Acs_cl[im+3]<<setw(8)<<Acs_cl[im+4]<<setw(8)<<Acs_cl[im+5]<<"\n";
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
				<<setw(8)<<A_cl[im]<<setw(8)<<A_cl[im+1]<<setw(8)<<A_cl[im+2]
				<<setw(8)<<A_cl[im+3]<<setw(8)<<A_cl[im+4]<<setw(8)<<A_cl[im+5]<<"\n";
				im=im+6;
			}
		}
	}
	cout<<"\n";
	iv=0;
	cout<<"BCS_PHI - B_PHI\n";
	cout<<setw(5)<<"iz"<<setw(5)<<"it"<<setw(5)<<"ir"
		<<setw(8)<<"bcs"<<setw(8)<<"b"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
				<<setw(8)<<bcs_phi[iv]<<setw(8)<<b_phi[iv]<<"\n";
				iv++;
			}
		}
	}
	cout<<"\n";
	iv=0;
	cout<<"BCS_K - B_K and BCS_CL - B_CL\n";
	cout<<setw(5)<<"iz"<<setw(5)<<"it"<<setw(5)<<"ir"
		<<setw(10)<<"bcs"<<setw(10)<<"b"<<setw(10)<<"bcs"<<setw(10)<<"b"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
				<<setw(10)<<bcs_k[iv]<<setw(10)<<b_k[iv]
				<<setw(10)<<bcs_cl[iv]<<setw(10)<<b_cl[iv]<<"\n";
				iv++;
			}
		}
	}
	cout<<"\n";

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
//OUTPUTRESULTS
//
//Esegue l'output del risultato di un passo di iterazione
///////////////////////////////////////////////////////////////////////////////
int outputresults(double *phi, double *k, double *cl)
{
	int iz,ir,it;
	unsigned long int iv;

	cout<<"RESULT OUTPUT\n";
	for(iz=nz-1;iz>=0;iz--)
	{
	        cout<<"----------------------------------------------------------------------\n";
	        cout<<"iz = "<<iz<<" PHI\n";
	        for(it=0;it<nt;it++)
	        {
	      	  cout<<setw(8)<<"";
	      	  for(ir=0;ir<nr;ir++)
	      	  {
	      		  iv=(nrnt*iz)+(nr*it)+ir;
	      		  cout<<setw(8)<<phi[iv];
	      	  }
	      	  cout<<setw(8)<<phibnd[iz]<<"\n";
	        }
	        cout<<"iz = "<<iz<<" K\n";
	        for(it=0;it<nt;it++)
	        {
	      	  cout<<setw(8)<<"";
	      	  for(ir=0;ir<nr;ir++)
	      	  {
	      		  iv=(nrnt*iz)+(nr*it)+ir;
	      		  cout<<setw(8)<<k[iv];
	      	  }
	      	  cout<<"\n";
	        }
	        cout<<"iz = "<<iz<<" CL\n";
	        for(it=0;it<nt;it++)
	        {
	      	  cout<<setw(8)<<"";
	      	  for(ir=0;ir<nr;ir++)
	      	  {
	      		  iv=(nrnt*iz)+(nr*it)+ir;
	      		  cout<<setw(8)<<cl[iv];
	      	  }
	      	  cout<<"\n";
	        }
	        cout<<"----------------------------------------------------------------------\n";
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
//OUTPUTVECTOR
//
//Esegue l'output di un singolo vettore
///////////////////////////////////////////////////////////////////////////////
int outputvector(double *v)
{
	unsigned long int iv;
	int iz,it,ir;

	iv=0;
	for(iz=0;iz<nz;iz++)
	{
	        cout<<"----------------------------------------------------------------------\n";
	        for(it=0;it<nt;it++)
	        {
	      	  cout<<setw(8)<<"";
	      	  for(ir=0;ir<nr;ir++)
	      	  {
	      		  iv=(nrnt*iz)+(nr*it)+ir;
	      		  cout<<setw(8)<<v[iv];
	      	  }
	      	  cout<<"\n";
	        }
	        cout<<"----------------------------------------------------------------------\n";
	}
	cout<<"\n";

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
//OUTPUTMATRIX
//
//Esegue l'output di una singola matrice
///////////////////////////////////////////////////////////////////////////////
int outputmatrix(double *A)
{
	unsigned long int im;
	int iz,it,ir;

	im=0;
	cout<<"A\n";
	cout<<setw(5)<<"iz"<<setw(5)<<"it"<<setw(5)<<"ir"
		<<setw(8)<<"bk"<<setw(8)<<"fw"<<setw(8)<<"dx"
		<<setw(8)<<"sx"<<setw(8)<<"dw"<<setw(8)<<"up"<<"\n";
	for(iz=0;iz<nz;iz++)
	{
		for(it=0;it<nt;it++)
		{
			for(ir=0;ir<nr;ir++)
			{
				cout<<setw(5)<<iz<<setw(5)<<it<<setw(5)<<ir
					<<setw(8)<<A[im]<<setw(8)<<A[im+1]<<setw(8)<<A[im+2]
					<<setw(8)<<A[im+3]<<setw(8)<<A[im+4]<<setw(8)<<A[im+5]<<"\n";
				im=im+6;
			}
		}
	}
	cout<<"\n";

	return 0;
}



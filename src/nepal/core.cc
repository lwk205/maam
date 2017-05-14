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
//M4V_P(double *vout, double* m, double *vin)
//
//
// Esegue il prodotto matrice*vettore nel caso di griglia cilindrica regolare
// (usato per l'equazione di Poisson) 
// Ipotesi alla base del funzionamento:
// 	Il coefficiente che moltiplica l'elemento omologo e' 1
// 	La matrice contiene a gruppi di 6 elementi successivi 
// 	i coeffcienti delle 6 celle adiacenti
// 	Il vettore contiene in linea tutte le celle
// 	L'ordine degli indici e': ir, it, iz (ir e' il piu' veloce)
// 	La presenza dei confini e' contenuta nell'implementazione 
// 	del prodotto non nei termini della matrice
// 	(cioe' non bisogna moltiplicare ad esempio per le celle sotto quando iz = 0)
int m4v_P(double *vout, double* m, double *vin)
{
	int ir,it,iz;
	unsigned long int iv,im;

	iv=0;
	im=0;
	//Prima fetta iz=0
		//Primo spicchio it=0
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
		}
		//Ultimo spicchi it=nt-1
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
	//Fette interne iz=1:nz-2
	for(iz=1;iz<nz-1;iz++)
	{
		//Primo spicchio it=0
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			//NUCLEO DELLA ROUTINE E' QUESTO IL CICLO ESEGUITO PIU' VOLTE
			  
			//Versione standard
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Versione sperimentale
			//iv_bk=iv;
			//im_bk=im;
			//for(ir=1;ir<nr-1;ir++)
			//{
			//	vout[iv]=vin[iv];
			//	vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			//	vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			//	vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			//	vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			//	im+=4;
			//	iv++;
			//}
			//iv=iv_bk;
			//im=im_bk;
			//for(ir=1;ir<nr-1;ir++)
			//{
			//	im+=4;
			//	vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			//	iv++;
			//}
			//iv=iv_bk;
			//im=im_bk;
			//for(ir=1;ir<nr-1;ir++)
			//{
			//	im+=5;
			//	vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			//	iv++;
			//}

			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
		}
		//Ultimo spicchi it=nt-1
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
	}
	//Ultima fetta iz=nz-1
		//Primo spicchio it=0
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
		//Spicchi interni it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
			//Ultimo elemento ir=nr-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
		}
		//Ultimo spicchi it=nt-1
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Elementi interni ir=1:nr-2
			for(ir=1;ir<nr-1;ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
			vout[iv]=vin[iv];
			//Ultimo elemento ir=nr-1
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;

			
	//if(iv!=nrntnz)
	//{
	//	cerr<<"ERROR in Poisson matrix moltiplication\n";
	//	exit(1);
	//}

	return 0;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//M4V_N(double *vout, double* m, double *vin)
// Esegue il prodotto matrice*vettore nel caso di griglia irregolare
// (equazione di Nernst-Planck)
// Ipotesi alla base del funzionamento:
// 	La matrice contiene a gruppi di 6 elementi successivi 
// 	i coefficienti delle 6 celle adiacenti
// 	Il vettore contiene in linea tutte le celle
// 	L'ordine degli indici e': ir, it, iz (ir e' il piu' veloce)
// 	I vettori inside_end e ouside_begin contengo rispettivamente
// 		numero elementi prima della prima interruzione
// 			nel caso non ci sia interruzione = nr
// 		numero primo elemento nuovamente parte del sistema
// 			nel caso il sistema non riprende = nr
// 			nel caso non ci sia interruzione = nr
// 	ir = 0 e' sempre parte del sistema
// 	ir = 1 e' sempre parte del sistema
// 		La prima ipotesi e' necessaria perche' il canale sia aperto
// 		La seconda ipotesi semplifica il calcolo evitando 
// 		l'eccezione di un elemento
// 		contemporaneamente primo ed ultimo
// 	Se ir = nr-1 e' parte del sistema lo e' anche ir=nr-2
// 		Questo semplifica il calcolo per la stessa ragione 
// 		esposta al punto precedente
// 		si evita cioe' l'eccezione di un elemento sia primo 
// 		che ultimo della soluzione 
int m4v_N(double *vout, double* m, double *vin)
{
	int ir,it,iz,iznt,end,begin;
	unsigned long int iv,im;

	iv=0;
	im=0;
	//Prima fetta iz=0
		//Primo spicchio it=0
		end=inside_end[0];
		begin=outside_begin[0];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
					im++;	//DW
					vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
		//Spicchi it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
		end=inside_end[it];
		begin=outside_begin[it];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
					im++;	//DW
					vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
		}
		//Spicchio it=nt-1
		end=inside_end[it];
		begin=outside_begin[it];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			im++;	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
					im++;	//DW
					vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				im++;	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
	//Fette iz=1:nz-2
	for(iz=1;iz<nz-1;iz++)
	{
		//Primo spicchio it=0
		iznt=iz*nt;
		end=inside_end[iznt];
		begin=outside_begin[iznt];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
					vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
					vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
		//Spicchi it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			end=inside_end[iznt+it];
			begin=outside_begin[iznt+it];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
			//Per saltare gli elementi non parte del sistema
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
					vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
					vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
		}
		//Spicchio it=nt-1
			end=inside_end[iznt+it];
			begin=outside_begin[iznt+it];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
					vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
					vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				vout[iv]+=((m[im++])*(vin[iv+nrnt]));	//UP
				iv++;
			}
	}
	//Fetta iz=nz-1
		//Primo spicchio it=0
		iznt=iz*nt;
		end=inside_end[iznt];
		begin=outside_begin[iznt];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
					vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
					im++;	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv+nrnt_nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
		//Spicchi it=1:nt-2
		for(it=1;it<nt-1;it++)
		{
			end=inside_end[iznt+it];
			begin=outside_begin[iznt+it];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
					vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
					im++;	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv+nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
		}
		//Spicchio it=nt-1
			end=inside_end[iznt+it];
			begin=outside_begin[iznt+it];
			//Primo elemento ir=0
			vout[iv]=vin[iv];
			im++;	//BK
			vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Elementi interni ir=1:end-2
			for(ir=1;ir<(end-1);ir++)
			{
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}
			//Ultimo elemento ir=end-1
			vout[iv]=vin[iv];
			vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
			im++;	//FW
			vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
			vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
			vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
			im++;	//UP
			iv++;
			//Per saltare gli elementi non parte del sistema
			iv+=(begin-end);
			im+=(6*(begin-end));
			//Elementi di sistema esterni all'interruzione
			if (begin<nr)
			{
				//Primo elemento ir=begin
				vout[iv]=vin[iv];
				im++;	//BK
				vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
				//Elementi interni ir=begin+1:nr-2
				for(ir=begin+1;ir<(nr-1);ir++)
				{
					vout[iv]=vin[iv];
					vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
					vout[iv]+=((m[im++])*(vin[iv+1]));	//FW
					vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
					vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
					vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
					im++;	//UP
					iv++;
				}
				//Ultimo elemento ir=nr-1
				vout[iv]=vin[iv];
				vout[iv]+=((m[im++])*(vin[iv-1]));	//BK
				im++;	//FW
				vout[iv]+=((m[im++])*(vin[iv-nr]));	//DX
				vout[iv]+=((m[im++])*(vin[iv-nrnt_nr]));	//SX
				vout[iv]+=((m[im++])*(vin[iv-nrnt]));	//DW
				im++;	//UP
				iv++;
			}


	//if (iv!=nrntnz)
	//{
	//	cerr<<"ERRORE in matrix vector moltiplication iv = "<<iv<<"\n";
	//	exit(1);
	//}
	return 0;
}
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//CONJUGATE GRADIENT
//
//	x		Vettore incognito
//	tol		Tolleranza raggiunta
//	it		Numero iterazioni eseguite
//
//	A		Matrice sistema
//	b 		vettore termini noti
//	x0		Soluzione di primo tentativo
//	max_iter	Numero massimo iterazioni permesse
//	tol_aim		Tolleranza obbiettivo
//
//	Funzionamento algoritmo
//
//	r = b - A*x0
//	resid = norm(r)/norm(b)
//	Controlla che il residuo non sia gia' al di sotto della tolleranza desiderata
//	
//	Prima iterazione
//	rh1=<r,r>
//	p = r
//	q = A*p
//	alpha = rh1 / <p,q>
//	x = x + alpha * p
//	r = r - alpha * q
//	Testa convergenza
//
//	Iterazione fino a convergenza o numero massimo di iterazioni
//		rh2 = rh1
//		rh1 = <r,r>
//		beta = rh1 / rh2
//		p = r + beta * p
//		q = A*p
//		alpha = rh1 / <p,q>
//		x = x + alpha * p
//		r = r - alpha * q
//		Testa convergenza
int CG_P(double *x, double &tol, int &it, double *A, double *b, double *x0, int max_iter, double tol_aim)
{
	 
  unsigned long int iv;
  double normr,normb,rh1,rh2,alpha,beta,tmp;
  double *p, *r, *q; 
  int reached;

  //Allocazione memoria
  p=new double [nrntnz];
  r=new double [nrntnz];
  q=new double [nrntnz];

  //Norma termine noto
        normb=0;
	for(iv=0;iv<nrntnz;iv++) 
	{
		normb+=(b[iv])*(b[iv]);
	}
	normb=sqrt(normb);

  //Definizione vettore residuo iniziale r = b - A*x0
  	m4v_P(r,A,x0);	//r=A*x0
	for(iv=0;iv<nrntnz;iv++) 
	{
		r[iv]=b[iv] - r[iv];
	}

  //In modo da evitare una divisione per zero
  	if (normb == 0.0) normb = 1;
  
  //Controllo di non essere gia' al di sotto della tolleranza
        normr=0;
	for(iv=0;iv<nrntnz;iv++) 
	{
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;
  	if (tol <= tol_aim) 
  	{
    		it = 0;
		for(iv=0;iv<nrntnz;iv++) //x = x0
		{
			x[iv]=x0[iv];
		}
		delete[] p;
		delete[] r;
		delete[] q;
    		return 0;
  	}

  	it=1;
  	rh1=0;		//rh1=<r,r> & p = r
	for(iv=0;iv<nrntnz;iv++) 
	{
		rh1+=(r[iv])*(r[iv]);
		p[iv] = r [iv];
	}
  	m4v_P(q,A,p);		//q = A*p
  	tmp=0;		//alpha = rh1/<p,q>
	for(iv=0;iv<nrntnz;iv++) 
	{
		tmp+=(p[iv])*(q[iv]);
	}
  	alpha=rh1/tmp;
	for(iv=0;iv<nrntnz;iv++) //x = x + alpha * p &  r = r - alpha * q
	{
		x[iv]=x[iv] + alpha*p[iv];
		r[iv]=r[iv] - alpha*q[iv];
	}
        normr=0;	//Test convergenza
	for(iv=0;iv<nrntnz;iv++) 
	{
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;
  	if (tol <= tol_aim) reached=1;
  	else reached=0;

  //Ciclo iterativo
	  while((it<=max_iter)&&(!reached))
	  {
		rh2=rh1;
	  	rh1=0;		//rh1=<r,r>
		for(iv=0;iv<nrntnz;iv++) 
		{
			rh1+=(r[iv])*(r[iv]);
			p[iv] = r [iv];
		}
		beta=rh1/rh2;
		for(iv=0;iv<nrntnz;iv++) //p = r + beta * p
		{
			p[iv]=r[iv] + beta*p[iv];
		}
	  	m4v_P(q,A,p);		//q = A*p
	  	tmp=0;		//alpha = rh1/<p,q>
		for(iv=0;iv<nrntnz;iv++) 
		{
			tmp+=(p[iv])*(q[iv]);
		}
	  	alpha=rh1/tmp;
		for(iv=0;iv<nrntnz;iv++) //x = x + alpha * p &  r = r - alpha * q
		{
			x[iv]=x[iv] + alpha*p[iv];
			r[iv]=r[iv] - alpha*q[iv];
		}
	        normr=0;	//Test convergenza
		for(iv=0;iv<nrntnz;iv++) 
		{
			normr+=(r[iv])*(r[iv]);
		}
		normr=sqrt(normr);
	  	tol = normr/normb;
	  	if (tol <= tol_aim) reached=1;
	  	else reached=0;
		it++;			//Incremento numero iterazioni
	  }
  

	delete[] p;
	delete[] r;
	delete[] q;
  	
	if(reached)
	{
		return 0;
	}
	else 
	{
	  	cerr<<"Convergence not reached in iterative solution\n";
	  	cerr<<"# iteration = "<<it<<" Tollerance = "<<tol<<"\n";
		return 1;
	}
}

int CG_N(double *x, double &tol, int &it, double *A, double *b, double *x0, int max_iter, double tol_aim)
{
	 
  unsigned long int iv;
  double normr,normb,rh1,rh2,alpha,beta,tmp;
  double *p, *r, *q; 
  int reached;

  //Allocazione memoria
  p=new double [nrntnz];
  r=new double [nrntnz];
  q=new double [nrntnz];

  //Norma termine noto
        normb=0;
	for(iv=0;iv<nrntnz;iv++) 
	{
		normb+=(b[iv])*(b[iv]);
		p[iv]=r[iv]=q[iv]=0.0;
	}
	normb=sqrt(normb);

  //Definizione vettore residuo iniziale r = b - A*x0
  	m4v_N(r,A,x0);	//r=A*x0
	for(iv=0;iv<nrntnz;iv++) 
	{
		r[iv]=b[iv] - r[iv];
	}

  //In modo da evitare una divisione per zero
  	if (normb == 0.0) normb = 1;
  
  //Controllo di non essere gia' al di sotto della tolleranza
        normr=0;
	for(iv=0;iv<nrntnz;iv++) 
	{
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;
  	if (tol <= tol_aim) 
  	{
    		it = 0;
		for(iv=0;iv<nrntnz;iv++) //x = x0
		{
			x[iv]=x0[iv];
		}
		delete[] p;
		delete[] r;
		delete[] q;
    		return 0;
  	}

  	it=1;
  	rh1=0;		//rh1=<r,r> & p = r
	for(iv=0;iv<nrntnz;iv++) 
	{
		rh1+=(r[iv])*(r[iv]);
		p[iv] = r [iv];
	}
  	m4v_N(q,A,p);		//q = A*p
  	tmp=0;		//alpha = rh1/<p,q>
	for(iv=0;iv<nrntnz;iv++) 
	{
		tmp+=(p[iv])*(q[iv]);
	}
  	alpha=rh1/tmp;
	for(iv=0;iv<nrntnz;iv++) //x = x + alpha * p &  r = r - alpha * q
	{
		x[iv]=x[iv] + alpha*p[iv];
		r[iv]=r[iv] - alpha*q[iv];
	}
        normr=0;	//Test convergenza
	for(iv=0;iv<nrntnz;iv++) 
	{
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;
  	if (tol <= tol_aim) reached=1;
  	else reached=0;

  //Ciclo iterativo
	  while((it<=max_iter)&&(!reached))
	  {
		rh2=rh1;
	  	rh1=0;		//rh1=<r,r>
		for(iv=0;iv<nrntnz;iv++) 
		{
			rh1+=(r[iv])*(r[iv]);
			p[iv] = r [iv];
		}
		beta=rh1/rh2;
		for(iv=0;iv<nrntnz;iv++) //p = r + beta * p
		{
			p[iv]=r[iv] + beta*p[iv];
		}
	  	m4v_N(q,A,p);		//q = A*p
	  	tmp=0;		//alpha = rh1/<p,q>
		for(iv=0;iv<nrntnz;iv++) 
		{
			tmp+=(p[iv])*(q[iv]);
		}
	  	alpha=rh1/tmp;
		for(iv=0;iv<nrntnz;iv++) //x = x + alpha * p &  r = r - alpha * q
		{
			x[iv]=x[iv] + alpha*p[iv];
			r[iv]=r[iv] - alpha*q[iv];
		}
	        normr=0;	//Test convergenza
		for(iv=0;iv<nrntnz;iv++) 
		{
			normr+=(r[iv])*(r[iv]);
		}
		normr=sqrt(normr);
	  	tol = normr/normb;
		/*
		switch(debug_Nernst)
		{
			case 2:
				writedata(x,r,r);
			case 1:
				cout<<"\t\t"<<setw(5)<<it<<setw(20)<<normr<<setw(20)<<tol<<"\n";
				cin.get();
				break;
		}
		*/
	  	if (tol <= tol_aim) reached=1;
	  	else reached=0;
		it++;			//Incremento numero iterazioni
	  }
  

	delete[] p;
	delete[] r;
	delete[] q;
  	
	if(reached)
	{
		return 0;
	}
	else 
	{
	  	cerr<<"Convergence not reached in the Nernst iterative solution\n";
	  	cerr<<"# iteration = "<<it<<" Tollerance = "<<tol<<"\n";
		return 1;
	}
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//BICONJUGATE GRADIENT STABLE
//
//	x		Vettore incognito
//	tol		Tolleranza raggiunta
//	it		Numero iterazioni eseguite
//
//	A		Matrice sistema
//	b 		vettore termini noti
//	x0		Soluzione di primo tentativo
//	max_iter	Numero massimo iterazioni permesse
//	tol_aim		Tolleranza obbiettivo
//
//	Funzionamento algoritmo
//
//	r = b - A*x0
//	rtilde = r
//
//	tol = norm(r) / norm(b)
//	Controlla che il residuo non sia gia' al di sotto della tolleranza desiderata
//
//	Prima iterazione
//	rh1= <rtilde,r>
//	Se rh1 == 0 l'algoritmo e' fallito
//	p = r
//	v = A*p
//	alpha = rh1 / <rtilde,v>
//	s = r - alpha*v
//	tol = norm(s) / norm(b)
//	Testa convergenza, nel caso x = x + alpha * p
//	t = A*s
//	omega = <t,s> / <t,t>
//	x = x + alpha*p + omega*s
//	r = s - omega * t
//	rh2 = rh1
//	tol = norm(r) / norm(b)
//	Testa convergenza
//	Testa omega != 0, altrimenti algoritmo fallito
//
//      Iterazione fino a convergenza o numero massimo iterazioni
//		rh1= <rtilde,r>
//		Se rh1 == 0 l'algoritmo e' fallito
//		beta = (rh1/rh2) * (alpha/omega)
//		p = r + beta * (p - omega*v)
//		v = A*p
//		alpha = rh1 / <rtilde,v>
//		s = r - alpha*v
//		tol = norm(s) / norm(b)
//		Testa convergenza, nel caso x = x + alpha * p
//		t = A*s
//		omega = <t,s> / <t,t>
//		x = x + alpha*p + omega*s
//		r = s - omega * t
//		rh2 = rh1
//		tol = norm(r) / norm(b)
//		Testa convergenza
//		Testa omega != 0, altrimenti algoritmo fallito
int BiCGSTAB_P(double *x,double &tol,int &it,double *A,double *b,double *x0,int max_iter,double tol_aim)
{
  double normb,normr,norms,rh1,rh2,alpha,beta,omega;
  double *r, *rtilde, *p, *s, *t, *v; 
  unsigned long int iv;
  double tmp,tmp2;

  r = new double [nrntnz];
  rtilde = new double [nrntnz];
  p = new double [nrntnz];
  s = new double [nrntnz];
  t = new double [nrntnz];
  v = new double [nrntnz];

  //Norma termine noto normb=norm(b)
  //+Inizializzazione
        it=0;
  	normb=0;
	for(iv=0;iv<nrntnz;iv++)
	{
		normb+=(b[iv])*(b[iv]);
		x[iv]=x0[iv];
	}
	normb=sqrt(normb);
  //In modo da evitare una divisione per zero
	if (normb == 0.0) normb = 1;

  //Definizione vettore residui iniziali
	m4v_P(r,A,x0);			//r=A*x0

  	normr=0;
	for(iv=0;iv<nrntnz;iv++)	//r=b-A*x0 & rtilde=r & norm(r)
	{
		r[iv]=b[iv] - r[iv];
		rtilde[iv]=r[iv];
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;

  //Controllo di non essere gia' al di sotto della tolleranza
	if (tol <= tol_aim) 
	{
  		delete[] r;
  		delete[] rtilde;
  		delete[] p;
  		delete[] s;
  		delete[] t;
  		delete[] v;
		return 0;
	}

  //Prima iterazione

	rh1=0;
	for(iv=0;iv<nrntnz;iv++)	//rh1= <rtilde,r> & p = r
	{
		rh1+=(rtilde[iv])*(r[iv]);
		p[iv]=r[iv];
	}
	m4v_P(v,A,p);			//v = A*p
	tmp=0;
	for(iv=0;iv<nrntnz;iv++)	//tmp= <rtilde,v>
	{
		tmp+=(rtilde[iv])*(v[iv]);
	}
	alpha=rh1/tmp;			//alpha = rh1 / <rtilde,v>
  	norms=0; 			
	for(iv=0;iv<nrntnz;iv++)	//s = r - alpha*v & norm(s)
	{
		s[iv]=r[iv] - alpha*v[iv];
		norms+=(r[iv])*(r[iv]);
	}
	norms=sqrt(norms);
  	tol = norms/normb;
	it++;
	if (tol <= tol_aim) 			//Test convergenza
	{
		for(iv=0;iv<nrntnz;iv++)	//x = x + alpha*p
		{
			x[iv]=x[iv] + alpha*p[iv];
		}
  		delete[] r;
  		delete[] rtilde;
  		delete[] p;
  		delete[] s;
  		delete[] t;
  		delete[] v;
		return 0;
	}

	m4v_P(t,A,s);			//t = A*s
	tmp=0;
	tmp2=0;
	for(iv=0;iv<nrntnz;iv++)	//tmp= <t,s> & tmp2=<t,t>
	{
		tmp+=(t[iv])*(s[iv]);
		tmp2+=(t[iv])*(t[iv]);
	}
	omega= tmp/tmp2;		//omega = <t,s> / <t,t>
  	normr=0;
	for(iv=0;iv<nrntnz;iv++) 	//x = x + alpha*p + omega*s & norm(r)
	{
		x[iv]=x[iv] + alpha*p[iv] + omega*s[iv];
		r[iv]=s[iv] - omega*t[iv];
		normr+=(r[iv])*(r[iv]);
	}
	rh2=rh1;			//rh2 = rh1
	normr=sqrt(normr);
  	tol = normr/normb;
	it++;
	if (tol <= tol_aim) 		//Test convergenza
	{
  		delete[] r;
  		delete[] rtilde;
  		delete[] p;
  		delete[] s;
  		delete[] t;
  		delete[] v;
		return 0;
	}
	if(omega==0)
	{
		cerr<<"ERROR in the Poisson iterative procedure - BiCGSTAB routine\n";
		cerr<<"omega parameter equal to zero\n";
		exit(1);
	}

  //Iterazione
	for(it=2;it<max_iter;it++)
	{
		rh1=0;
		for(iv=0;iv<nrntnz;iv++)	//rh1= <rtilde,r>
		{
			rh1+=(rtilde[iv])*(r[iv]);
		}
		if(rh1==0)
		{
			cerr<<"ERROR in the Poisson iterative procedure - BiCGSTAB routine\n";
			cerr<<"rh1 parameter equal to zero\n";
			exit(1);
		}
		beta= (rh1/rh2) * (alpha/omega);//beta = (rh1/rh2) * (alpha/omega)
		for(iv=0;iv<nrntnz;iv++) 	//p = r + beta * (p - omega*v)
		{
			p[iv]=r[iv] + beta*(p[iv] - omega*v[iv]);
		}
		m4v_P(v,A,p);			//v = A*p
		tmp=0;
		for(iv=0;iv<nrntnz;iv++)	//tmp= <rtilde,v>
		{
			tmp+=(rtilde[iv])*(v[iv]);
		}
		alpha=rh1/tmp;	//alpha = rh1 / <rtilde,v>
  		norms=0;
		for(iv=0;iv<nrntnz;iv++)	//s = r - alpha*v & norms
		{
			s[iv]=r[iv] - alpha*v[iv];
			norms+=(s[iv])*(s[iv]);
		}
		norms=sqrt(norms);
  		tol = norms/normb;
		if (tol <= tol_aim) 		//Test convergenza
		{
			for(iv=0;iv<nrntnz;iv++)	//x = x + alpha*p
			{
				x[iv]=x[iv] + alpha*p[iv];
			}
  			delete[] r;
  			delete[] rtilde;
  			delete[] p;
  			delete[] s;
  			delete[] t;
  			delete[] v;
			return 0;
		}

		m4v_P(t,A,s);			//t = A*s
		tmp=0;
		tmp2=0;
		for(iv=0;iv<nrntnz;iv++)	//tmp= <t,s> & tmp2=<t,t>
		{
			tmp+=(t[iv])*(s[iv]);
			tmp2+=(t[iv])*(t[iv]);
		}
		omega= tmp/tmp2;		//omega = <t,s> / <t,t>
  		normr=0;
		for(iv=0;iv<nrntnz;iv++) 	//x = x + alpha*p + omega*s & r = s - omega*t & norm(r)
		{
			x[iv]=x[iv] + alpha*p[iv] + omega*s[iv];
			r[iv]=s[iv] - omega*t[iv];
			normr+=(r[iv])*(r[iv]);
		}
		rh2=rh1;			//rh2 = rh1
		normr=sqrt(normr);
  		tol = normr/normb;
		if (tol <= tol_aim) 		//Test convergenza
		{
  			delete[] r;
  			delete[] rtilde;
  			delete[] p;
  			delete[] s;
  			delete[] t;
  			delete[] v;
			return 0;
		}
		if(omega==0)
		{
			cerr<<"ERROR in the Poisson iterative procedure - BiCGSTAB routine\n";
			cerr<<"omega parameter equal to zero\n";
			exit(1);
		}
	}

  cerr<<"Convergence not reached in the Poisson iterative solution - BiCGSTAB routine\n";
  cerr<<"# iteration = "<<it<<" Tollerance = "<<tol<<"\n";
  delete[] r;
  delete[] rtilde;
  delete[] p;
  delete[] s;
  delete[] t;
  delete[] v;

  return 1;
}
int BiCGSTAB_N(double *x,double &tol,int &it,double *A,double *b,double *x0,int max_iter,double tol_aim)
{
  double normb,normr,norms,rh1,rh2,alpha,beta,omega;
  double *r, *rtilde, *p, *s, *t, *v; 
  unsigned long int iv;
  double tmp,tmp2;

  r = new double [nrntnz];
  rtilde = new double [nrntnz];
  p = new double [nrntnz];
  s = new double [nrntnz];
  t = new double [nrntnz];
  v = new double [nrntnz];

  //Norma termine noto normb=norm(b)
  //+Inizializzazione
  	it=0;
  	normb=0;
	for(iv=0;iv<nrntnz;iv++)
	{
		normb+=(b[iv])*(b[iv]);
		r[iv]=rtilde[iv]=p[iv]=s[iv]=t[iv]=v[iv]=0.0;
		x[iv]=x0[iv];
	}
	normb=sqrt(normb);
  //In modo da evitare una divisione per zero
	if (normb == 0.0) normb = 1;

  //Definizione vettore residui iniziali
	m4v_N(r,A,x0);			//r=A*x0
	
  	normr=0;
	for(iv=0;iv<nrntnz;iv++)	//r=b-A*x0 & rtilde=r & Calcolo norm(r)
	{
		r[iv]=b[iv] - r[iv];
		rtilde[iv]=r[iv];
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;

  //Controllo convergenza	
	if (tol <= tol_aim) 
	{
  		delete[] r;
  		delete[] rtilde;
  		delete[] p;
  		delete[] s;
  		delete[] t;
  		delete[] v;
		return 0;
	}

  //Prima iterazione
	rh1=0;
	for(iv=0;iv<nrntnz;iv++)	//rh1= <rtilde,r> & p = r
	{
		rh1+=(rtilde[iv])*(r[iv]);
		p[iv]=r[iv];
	}
	m4v_N(v,A,p);			//v = A*p
	tmp=0;
	for(iv=0;iv<nrntnz;iv++)	//tmp= <rtilde,v>
	{
		tmp+=(rtilde[iv])*(v[iv]);
	}
	alpha=rh1/tmp;			//alpha = rh1 / <rtilde,v>
  	norms=0; 			//tol = norm(s) / norm(b)
	for(iv=0;iv<nrntnz;iv++)	//s = r - alpha*v & norm(s)
	{
		s[iv]=r[iv] - alpha*v[iv];
		norms+=(s[iv])*(s[iv]);
	}
	norms=sqrt(norms);
  	tol = norms/normb;
	it++;
  //Controllo convergenza	
	if (tol <= tol_aim) 			//Test convergenza
	{
		for(iv=0;iv<nrntnz;iv++)	//x = x + alpha*p
		{
			x[iv]=x[iv] + alpha*p[iv];
		}
  		delete[] r;
  		delete[] rtilde;
  		delete[] p;
  		delete[] s;
  		delete[] t;
  		delete[] v;
		return 0;
	}

	m4v_N(t,A,s);			//t = A*s
	tmp=0;
	tmp2=0;
	for(iv=0;iv<nrntnz;iv++)	//tmp= <t,s> & tmp2=<t,t>
	{
		tmp+=(t[iv])*(s[iv]);
		tmp2+=(t[iv])*(t[iv]);
	}
	omega= tmp/tmp2;		//omega = <t,s> / <t,t>
  	normr=0;
	for(iv=0;iv<nrntnz;iv++) 	//x = x + alpha*p + omega*s & r = s - omega * t & norm(r)
	{
		x[iv]=x[iv] + alpha*p[iv] + omega*s[iv];
		r[iv]=s[iv] - omega*t[iv];
		normr+=(r[iv])*(r[iv]);
	}
	rh2=rh1;			//rh2 = rh1
	normr=sqrt(normr);
  	tol = normr/normb;
  //Controllo convergenza	
	if (tol <= tol_aim) 		//Test convergenza
	{
  		delete[] r;
  		delete[] rtilde;
  		delete[] p;
  		delete[] s;
  		delete[] t;
  		delete[] v;
		return 0;
	}
	if(omega==0)
	{
		cerr<<"ERROR in the Nernst iterative procedure - BiCGSTAB routine\n";
		cerr<<"omega parameter equal to zero\n";
		exit(1);
	}

  //Iterazione
	for(it=2;it<max_iter;it++)
	{
		rh1=0;
		for(iv=0;iv<nrntnz;iv++)	//rh1= <rtilde,r>
		{
			rh1+=(rtilde[iv])*(r[iv]);
		}
		if(rh1==0)
		{
			cerr<<"ERROR in the Nernst iterative procedure - BiCGSTAB routine\n";
			cerr<<"rh1 parameter equal to zero\n";
			exit(1);
		}
		beta= (rh1/rh2) * (alpha/omega);//beta = (rh1/rh2) * (alpha/omega)
		for(iv=0;iv<nrntnz;iv++) 	//p = r + beta * (p - omega*v)
		{
			p[iv]=r[iv] + beta*(p[iv] - omega*v[iv]);
		}
		m4v_N(v,A,p);			//v = A*p
		tmp=0;
		for(iv=0;iv<nrntnz;iv++)	//tmp= <rtilde,v>
		{
			tmp+=(rtilde[iv])*(v[iv]);
		}
		alpha=rh1/tmp;	//alpha = rh1 / <rtilde,v>
  		norms=0; 
		for(iv=0;iv<nrntnz;iv++)	//s = r - alpha*v & norm(s)
		{
			s[iv]=r[iv] - alpha*v[iv];
			norms+=(s[iv])*(s[iv]);
		}
		norms=sqrt(norms);
  		tol = norms/normb;
  		//Controllo convergenza	
		if (tol <= tol_aim) 		//Test convergenza
		{
			for(iv=0;iv<nrntnz;iv++)	//x = x + alpha*p
			{
				x[iv]=x[iv] + alpha*p[iv];
			}
  			delete[] r;
  			delete[] rtilde;
  			delete[] p;
  			delete[] s;
  			delete[] t;
  			delete[] v;
			return 0;
		}

		m4v_N(t,A,s);			//t = A*s
		tmp=0;
		tmp2=0;
		for(iv=0;iv<nrntnz;iv++)	//tmp= <t,s> & tmp2=<t,t>
		{
			tmp+=(t[iv])*(s[iv]);
			tmp2+=(t[iv])*(t[iv]);
		}
		omega= tmp/tmp2;		//omega = <t,s> / <t,t>
  		normr=0;
		for(iv=0;iv<nrntnz;iv++) 	//x = x + alpha*p + omega*s & r = s - omega*t & norm(r)
		{
			x[iv]=x[iv] + alpha*p[iv] + omega*s[iv];
			r[iv]=s[iv] - omega*t[iv];
			normr+=(r[iv])*(r[iv]);
		}
		rh2=rh1;			//rh2 = rh1
		normr=sqrt(normr);
  		tol = normr/normb;
  		//Controllo convergenza	
		if (tol <= tol_aim) 		//Test convergenza
		{
  			delete[] r;
  			delete[] rtilde;
  			delete[] p;
  			delete[] s;
  			delete[] t;
  			delete[] v;
			return 0;
		}
  		//Output debug	
		if(omega==0)
		{
			cerr<<"ERROR in the Nernst iterative procedure - BiCGSTAB routine\n";
			cerr<<"omega parameter equal to zero\n";
			exit(1);
		}
	}

  cerr<<"Convergence not reached in the Nernst iterative solution - BiCGSTAB routine\n";
  cerr<<"# iteration = "<<it<<" Tollerance = "<<tol<<"\n";
  delete[] r;
  delete[] rtilde;
  delete[] p;
  delete[] s;
  delete[] t;
  delete[] v;
  return 1;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//CONJUGATE GRADIENT SQUARED METHOD
//
//	x		Vettore incognito
//	tol		Tolleranza raggiunta
//	it		Numero iterazioni eseguite
//
//	A		Matrice sistema
//	b 		vettore termini noti
//	x0		Soluzione di primo tentativo
//	max_iter	Numero massimo iterazioni permesse
//	tol_aim		Tolleranza obbiettivo
//
//	Funzionamento algoritmo
//
//
//	Prima iterazione
//	r = b - A*x0
//	rtilde = r
//	rh1= <rtilde,r>
//	u = r
//	p = r
//	v = A * p
//	alpha = rh1 / <rtilde,v>
//	q = u - alpha * v
//	utilde = u + q
//	x = x + alpha * utilde
//	q = A * utilde
//	r = r - alpha * q
//	Controllo convergenza
//	
//	Iterazione fino a convergenza o numero massimo di iterazioni
//		rh2=rh1
//		rh1= <rtilde,r>
//		Se rh1 == 0 il metodo e' fallito
//		u = r + beta * q
//		p = u + beta * (q + beta * p)
//		v = A * p
//		alpha = rh1 / <rtilde,v>
//		q = u - alpha * v
//		utilde = u + q
//		x = x + alpha * utilde
//		q = A * utilde
//		r = r - alpha * q
//		Controllo convergenza
int CGS_N(double *x,double &tol,int &it,double *A,double *b,double *x0,int max_iter,double tol_aim)
{
  double normb,normr,rh1,rh2,alpha,beta;
  double *r, *rtilde, *p, *q, *u, *utilde, *v; 
  unsigned long int iv;
  double tmp;

  r = new double [nrntnz];
  rtilde = new double [nrntnz];
  p = new double [nrntnz];
  q = new double [nrntnz];
  u = new double [nrntnz];
  utilde = new double [nrntnz];
  v = new double [nrntnz];

  //Output debug	
  /*
	if(debug_Nernst)
	{
		cout<<"\t\tNERNST ITERATION DEBUG\n";
		cout<<"\t\t"<<setw(5)<<"#it"<<setw(20)<<"Res. Norm"<<setw(20)<<"Tollerance"<<"\n";
	}
*/

  //Norma termine noto normb=norm(b)
  	normb=0;
	for(iv=0;iv<nrntnz;iv++)
	{
		normb+=(b[iv])*(b[iv]);
	}
	normb=sqrt(normb);

  //Definizione vettore residui iniziali
	m4v_N(r,A,x0);			//r=A*x0
	

	for(iv=0;iv<nrntnz;iv++)	//r=b-A*x0 & rtilde=r
	{
		r[iv]=b[iv] - r[iv];
		rtilde[iv]=r[iv];
	}

  //In modo da evitare una divisione per zero
	if (normb == 0.0) normb = 1;
  
  //Controllo di non essere gia' al di sotto della tolleranza
  	normr=0;
	for(iv=0;iv<nrntnz;iv++)
	{
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;
	for(iv=0;iv<nrntnz;iv++)	//x = x0
	{
		x[iv]=x0[iv];
	}
  //Controllo convergenza	
	it = 0;
	if (tol <= tol_aim) 
	{
		return 0;
	}
  //Output debug	
  /*
	switch(debug_Nernst)
	{
		case 2:
			writedata(x,r,p);
		case 1:
			cout<<"\t\t"<<setw(5)<<it<<setw(20)<<normr<<setw(20)<<tol<<"\n";
			cin.get();
			break;
	}
	*/

  //Prima iterazione
	rh1=0;
	for(iv=0;iv<nrntnz;iv++)	//rh1= <rtilde,r> & u = r & p = r
	{
		rh1+=(rtilde[iv])*(r[iv]);
		u[iv]=r[iv];
		p[iv]=r[iv];
	}
	m4v_N(v,A,p);			//v = A*p
	tmp=0;
	for(iv=0;iv<nrntnz;iv++)	//tmp= <rtilde,v>
	{
		tmp+=(rtilde[iv])*(v[iv]);
	}
	alpha=rh1/tmp;			//alpha = rh1 / <rtilde,v>
	for(iv=0;iv<nrntnz;iv++)	//q = u - alpha*v & utilde = u + q & x = x + alpha * utilde
	{
		q[iv]=u[iv] - alpha*v[iv];
		utilde[iv] = u[iv] + q[iv];
		x[iv] = x[iv] + alpha*utilde[iv];
	}
	m4v_N(q,A,utilde);			//q = A*utilde
	for(iv=0;iv<nrntnz;iv++)	//r = r - alpha * q
	{
		r[iv]=r[iv] - alpha*q[iv];
	}
  	normr=0;
	for(iv=0;iv<nrntnz;iv++)
	{
		normr+=(r[iv])*(r[iv]);
	}
	normr=sqrt(normr);
  	tol = normr/normb;
  //Controllo convergenza	
  	it=1;
	if (tol <= tol_aim) 			//Test convergenza
	{
		return 0;
	}

  //Iterazione	
	for(it=1;it<max_iter;it++)
	{
		rh2=rh1;
		rh1=0;
		for(iv=0;iv<nrntnz;iv++)	//rh1= <rtilde,r>
		{
			rh1+=(rtilde[iv])*(r[iv]);
		}
		if(rh1==0)
		{
			cerr<<"ERROR in the Nernst iterative solution: rh1 == 0\n";
			exit(1);
		}
		beta=rh1/rh2;
		for(iv=0;iv<nrntnz;iv++)	//u = r + beta * q & p = u + beta*(q+beta*p) 
		{
			u[iv]=r[iv] + beta*q[iv];
			p[iv]=u[iv] + beta*(q[iv]+ beta*p[iv]);
		}
		m4v_N(v,A,p);			//v = A*p
		tmp=0;
		for(iv=0;iv<nrntnz;iv++)	//tmp= <rtilde,v>
		{
			tmp+=(rtilde[iv])*(v[iv]);
		}
		alpha=rh1/tmp;			//alpha = rh1 / <rtilde,v>
		for(iv=0;iv<nrntnz;iv++)	//q = u - alpha*v & utilde = u + q & x = x + alpha * utilde
		{
			q[iv]=u[iv] - alpha*v[iv];
			utilde[iv] = u[iv] + q[iv];
			x[iv] = x[iv] + alpha*utilde[iv];
		}
		m4v_N(q,A,utilde);			//q = A*utilde
		for(iv=0;iv<nrntnz;iv++)	//r = r - alpha * q
		{
			r[iv]=r[iv] - alpha*q[iv];
		}
	  	normr=0;
		for(iv=0;iv<nrntnz;iv++)
		{
			normr+=(r[iv])*(r[iv]);
		}
		normr=sqrt(normr);
	  	tol = normr/normb;
	  //Controllo convergenza	
		if (tol <= tol_aim) 			//Test convergenza
		{
			return 0;
		}
  		//Output debug	
		/*
		switch(debug_Nernst)
		{
			case 2:
			writedata(x,r,p);
			case 1:
				cout<<"\t\t"<<setw(5)<<it<<setw(20)<<normr<<setw(20)<<tol<<"\n";
				cin.get();
			break;
		}
		*/
	}

  cerr<<"Convergence not reached in the Nernst iterative solution - CGS_N routine\n";
  cerr<<"# iteration = "<<it<<" Tollerance = "<<tol<<"\n";
  //exit(1);
  return 1;
}
////////////////////////////////////////////////////////////////////////////////

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

//CRT2CLN
//Passaggio da coordinate cartesiane a coordinate cilindriche
void crt2cln(double x,double y,double &r,double &t)
{
	r=sqrt(x*x+y*y);
	t=atan2(y,x);

	return;
}

//DEG2RAD
//Trasformazione di un angolo in radianti
//Restituisce un valore nell'intervallo ]-PI:PI]
double deg2rad(double teta)
{
	teta=((2*PI)*teta)/360.0;
	if(teta>PI)teta=teta-(2*PI);
	return teta;
}

//DIF_RAD
//return = t1 - t2
//The value is returned in the interval [-pi;pi[
double dif_rad(double t1,double t2)
{
	double deltat;

	deltat=t1-t2;
	if (deltat<-PI) deltat=deltat+(2*PI);
	if (deltat>=PI) deltat=deltat-(2*PI);
	if ((deltat<-PI)||(deltat>PI))
	{
		cerr<<"ERROR IN dif_rad SUBROUTINE\n";
		exit(1);
	}

	return deltat;
}

//SEARCHCHAR
//Position of a character in string
//-1 if the character is not found
int searchchar(char *line,char c)
{
	int i=0;
	while((line[i])&&(line[i]!=c))i++;
	if (line[i])return i;
	else return -1;
}
//SEARCHNCHAR
//Position of the first character not equal to c
//-1 if all the character are equal to c
int searchNchar(char *line,char c)
{
	int i=0;
	while((line[i])&&(line[i]==c))i++;
	if (line[i])return i;
	else return -1;
}

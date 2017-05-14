#ifndef INCLUDE_COMMON_H
#define INCLUDE_COMMON_H

//Config header
#if HAVE_CONFIG_H
#include <config.h>
#endif

//Standard libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>
#include <ctime>
#include <csignal>

//Floating point control
#if HAVE_FENV_H
#include <fenv.h>
#endif

//Physical constants
const double PI=3.14;		//Pi greca
const double INF=1.7e308;	//Infinite
const double zero=1.7e-308;	//Infinitesimal

const double e_chg=1.609e-19;		//Elementary charge [C]
const double bltz_cnst=1.381e-23;	//Boltzmann constant [J/K]

//Measure Unit Conversion
const double mmol2ionA3=6.0224*1e-7;		//mmol --> ion/A^3
const double mmol2mVA=(0.60225*1.6022)/8.8544;	//mmol --> mV*A/A^3
						//It is used to rapidly compute charges
						//in the Poisson equation
const double mmolns2pA=1.6022*6.0225*1e-5;	//mmol/ns --> pA	

const double e2mVA=(1.6022/8.8544)*1e6;		//From elementary charge to mV*A
						//It is used to include the charge
						//in the Poisson equation, division
						//by eps0 is included

const double bltzONe=8.16796e-5;	//bltz_cnst/e_chg

const double J2cal=0.238846;	//From Joule to calories

const double eA2Deybe=16.022/3.335664;	//From eA to Deybe

//Namespace definition
using namespace std;

#endif //INCLUDE_COMMON_H

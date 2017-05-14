# Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AM_CXXFLAGS = -DDATADIR=\"$(datadir)\" 

bin_PROGRAMS = l_trsl l_rttn l_radchg l_fixnumber \
	       l_amber2gromacs l_dipole l_betabarrel \
	       l_flatbarrel l_expandrad l_resetchg \
	       l_mkcrown l_pdbrq2pqr l_findaxis

l_trsl_SOURCES = l_trsl.cc printhelp.cc
l_trsl_LDADD = -lmolecule
l_trsl_LDFLAGS = -L$(top_srcdir)/lib
l_trsl_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a

l_rttn_SOURCES = l_rttn.cc printhelp.cc
l_rttn_LDADD = -lmolecule
l_rttn_LDFLAGS = -L$(top_srcdir)/lib
l_rttn_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libetc.a

l_radchg_SOURCES = l_radchg.cc printhelp.cc
l_radchg_LDADD = -lmolecule
l_radchg_LDFLAGS = -L$(top_srcdir)/lib
l_radchg_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libatom.a \
			$(top_srcdir)/lib/libetc.a

l_fixnumber_SOURCES = l_fixnumber.cc printhelp.cc
l_fixnumber_LDADD = -lmolecule
l_fixnumber_LDFLAGS = -L$(top_srcdir)/lib
l_fixnumber_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a

l_amber2gromacs_SOURCES = l_amber2gromacs.cc printhelp.cc
l_amber2gromacs_LDADD = -lmolecule
l_amber2gromacs_LDFLAGS = -L$(top_srcdir)/lib
l_amber2gromacs_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a

l_dipole_SOURCES = l_dipole.cc printhelp.cc
l_dipole_LDADD = -lmolecule
l_dipole_LDFLAGS = -L$(top_srcdir)/lib
l_dipole_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a

l_betabarrel_SOURCES = l_betabarrel.cc printhelp.cc
l_betabarrel_LDADD = -lmolecule
l_betabarrel_LDFLAGS = -L$(top_srcdir)/lib
l_betabarrel_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a

l_flatbarrel_SOURCES = l_flatbarrel.cc printhelp.cc
l_flatbarrel_LDADD = -lmolecule
l_flatbarrel_LDFLAGS = -L$(top_srcdir)/lib
l_flatbarrel_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a

l_expandrad_SOURCES = l_expandrad.cc printhelp.cc
l_expandrad_LDADD = -lmolecule
l_expandrad_LDFLAGS = -L$(top_srcdir)/lib
l_expandrad_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libetc.a

l_resetchg_SOURCES = l_resetchg.cc printhelp.cc
l_resetchg_LDADD = -lmolecule
l_resetchg_LDFLAGS = -L$(top_srcdir)/lib
l_resetchg_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libatom.a \
			$(top_srcdir)/lib/libetc.a

l_mkcrown_SOURCES = l_mkcrown.cc printhelp.cc
l_mkcrown_LDADD = -lmolecule
l_mkcrown_LDFLAGS = -L$(top_srcdir)/lib
l_mkcrown_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libatom.a \
			$(top_srcdir)/lib/libetc.a

l_pdbrq2pqr_SOURCES = l_pdbrq2pqr.cc printhelp.cc
l_pdbrq2pqr_LDADD = -lmolecule
l_pdbrq2pqr_LDFLAGS = -L$(top_srcdir)/lib
l_pdbrq2pqr_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libatom.a \
			$(top_srcdir)/lib/libetc.a

l_findaxis_SOURCES = l_findaxis.cc printhelp.cc
l_findaxis_LDADD = -lmolecule
l_findaxis_LDFLAGS = -L$(top_srcdir)/lib
l_findaxis_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libatom.a \
			$(top_srcdir)/lib/libetc.a

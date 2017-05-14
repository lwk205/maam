# Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

check_PROGRAMS = libmolecule

libmolecule_SOURCES = libmolecule.cc
libmolecule_LDADD = -lmolecule -letc
libmolecule_LDFLAGS = -L$(top_srcdir)/lib
libmolecule_DEPENDENCIES = $(top_srcdir)/lib/libmolecule.a $(top_srcdir)/lib/libetc.a

TESTS = $(check_PROGRAMS)

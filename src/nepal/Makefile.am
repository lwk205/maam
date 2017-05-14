# Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

noinst_HEADERS = global.h functions.h

if FORCEMOMENTA
bin_PROGRAMS = nepal b_axis b_bind b_maps b_forces b_forcesmomenta b_watervolume b_vol
else
bin_PROGRAMS = nepal b_axis b_bind b_maps b_forces b_watervolume b_vol
endif

nepal_SOURCES = nepal.cc read.cc write.cc discretize.cc core.cc analysis.cc auxiliary.cc
nepal_LDADD = -letc -lmolecule
nepal_LDFLAGS = -L$(top_srcdir)/lib
nepal_DEPENDENCIES = $(top_srcdir)/lib/libetc.a $(top_srcdir)/lib/libmolecule.a \
		     global.h functions.h

b_axis_SOURCES = b_axis.cc  auxiliary.cc discretize.cc write.cc analysis.cc
b_axis_LDADD = -letc -lmolecule 
b_axis_LDFLAGS = -L$(top_srcdir)/lib
b_axis_DEPENDENCIES = $(top_srcdir)/lib/libetc.a $(top_srcdir)/lib/libmolecule.a

b_bind_SOURCES = b_bind.cc auxiliary.cc discretize.cc write.cc analysis.cc
b_bind_LDADD = -letc -lmolecule
b_bind_LDFLAGS = -L$(top_srcdir)/lib
b_bind_DEPENDENCIES = $(top_srcdir)/lib/libetc.a  $(top_srcdir)/lib/libmolecule.a \
		      global.h

b_maps_SOURCES = b_maps.cc auxiliary.cc
b_maps_LDADD = -letc -lmolecule
b_maps_LDFLAGS = -L$(top_srcdir)/lib
b_maps_DEPENDENCIES = $(top_srcdir)/lib/libetc.a  $(top_srcdir)/lib/libmolecule.a \
		      global.h

b_forces_SOURCES = b_forces.cc discretize.cc auxiliary.cc
b_forces_LDADD = -letc -lmolecule
b_forces_LDFLAGS = -L$(top_srcdir)/lib
b_forces_DEPENDENCIES = $(top_srcdir)/lib/libetc.a  $(top_srcdir)/lib/libmolecule.a

b_forcesmomenta_SOURCES = b_forcesmomenta.cc discretize.cc auxiliary.cc
b_forcesmomenta_LDADD = -letc -lmolecule -llapack -lblas -lg2c -lm
b_forcesmomenta_LDFLAGS = -L$(top_srcdir)/lib
b_forcesmomenta_DEPENDENCIES = $(top_srcdir)/lib/libetc.a  $(top_srcdir)/lib/libmolecule.a

b_watervolume_SOURCES = b_watervolume.cc auxiliary.cc
b_watervolume_LDADD = -letc -lmolecule
b_watervolume_LDFLAGS = -L$(top_srcdir)/lib
b_watervolume_DEPENDENCIES = $(top_srcdir)/lib/libetc.a  $(top_srcdir)/lib/libmolecule.a functions.h

b_vol_SOURCES = b_vol.cc auxiliary.cc read.cc discretize.cc
b_vol_LDADD = -letc -lmolecule
b_vol_LDFLAGS = -L$(top_srcdir)/lib
b_vol_DEPENDENCIES = $(top_srcdir)/lib/libetc.a  $(top_srcdir)/lib/libmolecule.a \
		      global.h

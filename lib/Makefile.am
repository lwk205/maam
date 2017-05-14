# Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

noinst_LIBRARIES = libetc.a libatom.a libmolecule.a

libetc_a_SOURCES = etc.cc etc.h

libatom_a_SOURCES = atom.cc atom.h etc.cc etc.h
libatom_a_DEPENDENCIES = libetc.a

libmolecule_a_SOURCES = molecule.cc molecule.h atom.cc atom.h etc.cc etc.h
libmolecule_a_DEPENDENCIES = libetc.a libatom.a

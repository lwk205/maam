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

int main(int argc, char *argv[])
{
	molecule in,out,test;
	double xt,yt,zt;

	ifstream in_pdb,test_pdb;
	ofstream out_pdb;

	in_pdb.open("../dat/test.pdbrq");
	out_pdb.open("../dat/test_tmp.pdbrq");

	in.readpdb(in_pdb);
	in.fix();
	xt=in.xcenter;
	yt=in.ycenter;
	zt=in.zcenter;
	out.copy(in);
	out.translate(-xt,-yt,-zt);
	out.fix();
	out.rotate_x(PI/2); //Rotation 2*PI around each axis
	out.rotate_x(PI/2);
	out.rotate_x(PI/2);
	out.rotate_x(PI/2);
	out.rotate_y(PI/2);
	out.rotate_y(PI/2);
	out.rotate_y(PI/2);
	out.rotate_y(PI/2);
	out.rotate_z(PI/2);
	out.rotate_z(PI/2);
	out.rotate_z(PI/2);
	out.rotate_z(PI/2);
	out.translate(+xt,+yt,+zt);
	out_pdb<<out<<endl;
	out_pdb.close();

	test_pdb.open("../dat/test_tmp.pdbrq");
	test.readpdb(test_pdb);

	test_pdb.close();
	in_pdb.close();
	remove("../dat/test_tmp.pdbrq");

	if((test.cmp(in))) 
	{
		return 0;
	}
	else 
	{
		return 1;
	}
}

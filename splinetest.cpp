#include <cmath>
#include "bsplines.h"
#include <armadillo>

using namespace std;

typedef unsigned int uing;

int pingcount = 0;

void ping() {pingcount++; cout << "ping" << pingcount << endl;}

int main()
{
	string knotsfile = "knots.txt";
	uint order = 4;
	uint nphysknots = 4;
	uint nknots = nphysknots + 2*(order - 1);	
	uint nsplines = nknots - order;
	uint ngrid = 100;

	vector <bspline> splines;
	splines.reserve(nsplines);

	//bspline initspline;
	//initspline.setsource(knotsfile);
	//initspline.setorder(order);
	//initspline.setgridpts(ngrid);
	//initspline.readknots();
	

	bspline initspline(knotsfile,0,order,ngrid);	
	initspline.makegrid();

	for (uint k = 0; k < nsplines; k++)
	{
		splines.push_back(initspline);
		splines[k].setindex(k);
		splines[k].calcspline();
		splines[k].writespline();	
	}
}

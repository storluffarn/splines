
#include <cmath>
#include "bsplines.h"
#include <armadillo>

using namespace std;

typedef unsigned int uing;

//int pingcount = 0;

//void ping() {pingcount++; cout << "ping" << pingcount << endl;}

int main()
{
	string knotfile = "testknots.txt";
	uint order = 4;
	uint gridpts = 1000;
	double tolerance = 10e-10;

	bsplines splines(knotfile,order,gridpts,tolerance);
	
	//splines.knotgrid();
	splines.makegrid();

	splines.printknots();
	splines.printgrid();

	splines.calcsplines();

	splines.writesplines();

}


#include <cmath>
#include <iostream>
#include <vector>
#include "bsplines.h"
#include <armadillo>

double pi = 4*atan(1);
int pingcount = 0;

typedef unsigned int uint;
using namespace std;

void ping() {pingcount++; cout << "ping" << pingcount << endl;}

double sphere(double r)
{
	double V = 4.0 / 3.0 * pi * pow(r,3);

	return V;
}

double spherecharge(double r, double Q0)
{
	double V = sphere(r);
	double Q = Q0 / V;

	return Q;
}

double shellcharge(double r, double rmin, double Q0)
{
	double V = sphere(r) - sphere(rmin);
	double Q = Q0 / V;

	return Q;
}

void calcknots(double rmin, double rmax, double Q0, int steps, int order, string filename)
{
	vector <double> knots;
	knots.reserve(steps + 2*(order - 1));
	
	double dr = (rmax-rmin)/steps;

	for (double r = dr; r < rmax; r += dr)
	{
		double knot = - r * 4 * pi * spherecharge(r,Q0);
		//double knot = - r * 4 * pi * shellcharge(r,rmin,Q0);
		
		cout << knot << endl;
		knots.push_back(knot);
	}

	//for (int k = 0; k < order - 1; k++)
	//{
	//	knots.insert(knots.begin(),knots.front());
	//	knots.insert(knots.end(),knots.back());
	//}

	ofstream fs;
	fs.open(filename);

	for (auto& el : knots)
		fs << el << endl;

	fs.close();
}

void buildmatrix (arma::mat* splinemat, vector<bspline>* splines, int order, int gridpoints)
{
	int eqs = gridpoints - 2 * (order - 1);

	for (int k = 0; k < eqs; k++)
	{
		//(*splines)[k].printvals();
		arma::rowvec splinerow(eqs+2,arma::fill::zeros);

		for (int l = 0; l < order - 1; l++)
		{
			vector <double>* d2s = (*splines)[k+l].getvals();
			//for(auto& el : (*d2s))
			//	cout << el << endl;

			splinerow(k+l) = (*d2s)[order+k-1];
			
			//cout << splinerow(k+l) << endl;
		}
		
		(*splinemat).insert_rows((*splinemat).n_rows,splinerow);
	}
}

//void boundcond (arma::mat* splinemat, int order, int gridpoints)
//{
//	int elements = gridpoints - 2 * (order - 1);
//
//	arma::rowvec toprow(elements,arma::fill::zeros);
//	arma::rowvec botrow(elements,arma::fill::zeros);
//
//	toprow(0) = 1;
//	botrow(botrow.n_cols-1) = 1;
//
//	(*splinemat).insert_rows(0,toprow);
//	(*splinemat).insert_rows((*splinemat).n_rows,botrow);
//}

//void solving(arma::mat* splinemat, arma::mat* sols, bspline* spline, int order, double Q0)
//{
//	vector<double>* knotsptr = spline->getknots();
//
//	int knotpts = spline->getknotpts();
//	int eqs = knotpts - 2 * (order - 1);
//	
//	arma::colvec knots(eqs+2,arma::fill::zeros);
//
//	for (int k = 0; k < eqs; k++)
//		knots(k+1) = (*knotsptr)[k+order-1];
//	
//	knots(0) = 0;
//	knots(knots.n_rows-2) = (Q0);
//	
//	knots.print();
//
//	cout << arma::size(knots) << " " << arma::size((*splinemat)) << endl;
//
//	(*sols) = arma::solve((*splinemat),knots);
//}

int main ()
{	
	int order = 4;
	string knotfile = "knots.txt";
	int gridpoints = 4;	// knotgrid() overrides this
	int physpoints = 4;
	int nsplines = gridpoints - order;

	double rmin = 0; 
	double rmax = 1;
	double totcharge = 1;

	calcknots(rmin,rmax,totcharge,physpoints,order,knotfile);
	
	vector<bspline> splines (nsplines);
	splines.reserve(nsplines);
	arma::mat splinemat;

	bspline initspline(knotfile,0,order,gridpoints);
	initspline.knotgrid();

	for (int k = 0; k < nsplines; k++)
		splines.push_back(initspline);

	for (auto& el : splines)
	{	
		auto i = &el - &splines[0];
		
		//el.printknots();
		//spline.printgrid();
		el.setindex(i);
		el.calcspline();
		el.writespline();
	}
		
	buildmatrix(&splinemat, &splines, order, gridpoints);	
	//boundcond(&splinemat, &spline, order);
	//splinemat.print();
	//splinemat.save("splinematrix.txt",arma::raw_ascii);

	//arma::mat sols;
	//solving(&splinemat, &sols, &spline,order,totcharge);
	
	//sols.print();
}

























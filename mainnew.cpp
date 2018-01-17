
#include <cmath>
#include <iostream>
#include <vector>
#include "bsplines.h"
#include <armadillo>
#include <boost/math/special_functions/laguerre.hpp>

double pi = 4*atan(1);
int pingcount = 0;
//double hbar = 1.0545718e-34;
//double me = 9.10938356-31;
//double e = 1.60217662e-19;
//double e0 = 8.85418781e-12;
double a0 = 5.29e-11;

typedef unsigned int uint;
using namespace std;

void ping() {pingcount++; cout << "ping" << pingcount << endl;}

constexpr uint64_t factorial(uint64_t n)
{ 
	    return n == 0 ? 1 : n*factorial(n-1);
}

double sphere(double r)
{
	double V = 4.0 / 3.0 * pi * pow(r,3);

	return V;
}

double shell(double rmax, double rmin)
{
	double V = sphere(rmax) - sphere(rmin);

	return V;
}

double hpsi(double r, int n, int l)
{
	//double Q = exp(-r/a0)/(sqrt(pi)*pow(a0,1.5));

	double a = 1;

	double t111 = factorial(n-l-1.0);	// 1
	double t112 = 2.0*n*factorial(n+l);	// 2
	double t11 = t111/t112;				// 0.5
	double t1 = sqrt(pow(2.0/(n*a),3)*t11);		// sqrt(8*0.5) = 2
	double t2 = exp(-r/(n*a))*pow(2.0*r/(n*a),l);	// 0
	double t3 = pow(boost::math::laguerre(n-l-1,2.0*r/(n*a)),2*l+1);

	double psi = t1*t2*t3;

	//cout << t111 << " " << t112 << " " << t11 << " " << t1 << " " << t2 << " " << t3 << endl;
	
	if (r == 0)
		psi = 0;

	return psi;
}

void calcpotsphere (vector <double>* realpots, uint gridpts, double rmax, double Q0, vector <double>* grid)
{
	realpots->reserve(gridpts);

	for (auto& el : (*grid))
	{
		double pot;

		if (el <= rmax)
			pot = Q0/(2*rmax)*(3-pow(el/rmax,2));
		else
			pot = Q0/el;
	
		realpots->push_back(pot);
	}
}

void calcpotshell (vector <double>* realpots, uint gridpts, double rmin, double rmax, double Q0, vector <double>* grid)
{
	realpots->reserve(gridpts);

	for (auto& el : (*grid))
	{
		double pot;

		if (el < rmin)
			pot = 2*pi*Q0 * ( pow(rmax,2) - pow(rmin,2) );
		else if (el >= rmin && el <= rmax)
			pot = 2*pi*Q0*pow(rmax,2) - 4*pi*Q0/3*(pow(el,2)/2 + pow(rmin,3)/el);
		else
			pot = Q0/el;
		
		cout << el << " " << pot << endl;	
		realpots->push_back(pot);
	}
}

void calchpot (vector <double>* realpots, uint gridpts, int n, int l, vector <double>* grid)
{
	realpots->reserve(gridpts);

	for (auto& el : (*grid))
	{
		double pot = pow(hpsi(el,n,l),2);

		pot *= 4*pi;
	
		//cout << pot << endl;
		realpots->push_back(pot);
	}
}

void calcknots(double rmin, double rmax, double rend, int n, int l, double Q0, int gridpts, int order, string knotfile, vector <double>* printgrid, vector<double>* dens, bool state)
{
	vector <double> knots;
	knots.reserve(gridpts + 2*(order - 1));
	dens->reserve(gridpts + 2*(order - 1));

	double dr = (rend)/gridpts;

	for (int k = 0; k <= gridpts; k++)
		printgrid->push_back(k*dr);

	if(state)
	{
		for (auto& r : (*printgrid))
		{
			double den = -4*pi*r*pow(hpsi(r,n,l),2);

			//cout << r << " " << den << endl;
			knots.push_back(r);
			dens->push_back(den);
		}
	}
	else
	{
		for (auto& r : (*printgrid))
		{
			double den = 0;
	
			if (rmin == 0)
			{
				double rho = Q0/sphere(rmax);
				
				if (r <= rmax)
					den = - r * 4 * pi * rho;
				else
					den = 0;
			}
			else
			{
				double rho = Q0/shell(rmax,rmin);
				
				if (r < rmin)
					den = 0;
				else if (r >= rmin && r <= rmax )
					den = - r * 4 * pi * rho;
				else
					den = 0;
			}			
		
			//cout << r << " " << den << endl;

			knots.push_back(r);
			dens->push_back(den);
		}
	}

	ofstream knotstream;
	knotstream.open(knotfile);

	for (auto& el : knots)
		knotstream << el << endl;
	
	knotstream.close();
}

void buildmatrix (arma::mat* splinemat, arma::colvec* dens, vector<double>* densvec, bsplines* splines, uint order, uint gridpts)
{
	uint eqs = gridpts;
	
	for (uint k = 0; k < eqs; k++)
	{
		//arma::rowvec splinerow(eqs+2,arma::fill::zeros);

		//for (uint l = 0; l < order - 1; l++)
		for (uint l = 1; l < eqs+2; l++)
		{
			//splinerow(k+l) = (*splines->getd2(order-1,k+l))[k];
			(*splinemat)(k,l) = (*splines->getd2(order-1,l))[k];
		}

		//splinemat->insert_rows((*splinemat).n_rows,splinerow);
	}
	
	dens->zeros(eqs+2);

	for (uint k = 0; k < eqs; k++)
		(*dens)(k+1) = (*densvec)[k];
}

void boundcond (arma::mat* splinemat, arma::colvec* dens, uint gridpts, double Q0)
{
	uint el = gridpts + 2;

	arma::mat tmat (el,el,arma::fill::zeros);

	tmat(0,0) = 1;
	tmat(tmat.n_cols-1,tmat.n_rows-1) = 1;

	tmat.submat(1,0,tmat.n_rows-2,tmat.n_cols-1) = (*splinemat);

	(*splinemat) = tmat;
	
	(*dens)(0) = 0;
	(*dens)((*dens).n_rows-1) = (Q0);
}

void solving (arma::mat* splinemat, arma::colvec* dens, arma::colvec* sols)
{
	// dealing with r = 0
	//splinemat->shed_col(0);
	//splinemat->shed_row(0);
	//dens->shed_row(0);
	//(*splinemat).shed_col((*splinemat).n_cols-1);
	//(*splinemat).shed_row((*splinemat).n_rows-1);
	//physpts.shed_row(physpts.n_rows-1);

	

	// LU-decomposition for increased performance	
	
	//cout << "matrices before solving " << endl;
	//splinemat->print();
	//cout << endl;
	//cout << "density: " << endl;
	//dens->print();
	//cout << endl;
	//cout << "sizes: " << splinemat->n_rows << "x" << splinemat->n_cols << " " << dens->size() << endl;

	arma::mat U,L;
	//lu(U,L,splinemat->submat(1,1,splinemat->n_rows-1.,splinemat->n_cols-1));
	lu(U,L,(*splinemat));
	arma::mat LUsplinemat = U*L;
	
	//physpts.print();
	//(*splinemat).print();
	//LUsplinemat.print();

	(*sols) = arma::solve(LUsplinemat,(*dens));

	//sols->insert_rows(0,0);
}


void estpot (vector <double>* estpots, arma::colvec* sols, bsplines* splines, uint gridpts, uint order)
{
	estpots->reserve(gridpts);
	vector<double>* grid = splines->getgrid();

	for (auto& el : (*grid))
	{
		uint k = &el - &(*grid)[0];
		double accpot = 0; 

		for (uint	l = 1; l < gridpts + 2; l++)
		{
			double x = (*sols)(l)*(*splines->getvals(order-1,l))[k];
				
			accpot += x;
		}

		estpots->push_back(accpot/el);
	}

	estpots->insert(estpots->begin(),0);
}


int main ()
{	
	string knotfile = "knots.dat";
	string outfile = "output.csv";
	uint order = 4;
	uint gridpts = 500;
	double tolerance = 5e-9;
	vector <double> printgrid;
	vector <double> dens;

	double rmin = 0.0; 
	double rmax = 1.0;
	double rend = 4*rmax;
	double Q0 = 1;
	bool state = 0;
	int n = 1;
	int l = 0;

	calcknots(rmin,rmax,rend,n,l,Q0,gridpts,order,knotfile,&printgrid,&dens,state);

	for (auto& el : printgrid)
		cout << el << endl;

	bsplines splines(knotfile,order,gridpts,tolerance);	
	splines.knotgrid();		

	// check norm
	//vector <double> diffgrid;
	//for (uint k = 0; k < printgrid.size() - 1; k++)
	//{

	//}

	//splines.printknots();

	//arma::mat splinemat;
	arma::mat splinemat(gridpts,gridpts + order - 2,arma::fill::zeros);
	arma::colvec densvec;
	vector <double> estpots;
	vector <double> realpots;

	splines.calcsplines();
	splines.writesplines();

	buildmatrix(&splinemat, &densvec, &dens, &splines, order, gridpts);	
	
	//cout << "splinemat after build function: " << endl;
	//splinemat.print();
	//cout << "densvec after build function: " << endl;
	//densvec.print();
	
	boundcond(&splinemat, &densvec, gridpts, Q0);

	//cout << "splinemat after boundary conds: " << endl;
	//splinemat.print();
	//cout << "densvec after boundary conds: " << endl; 
	//densvec.print();
	
	//splinemat.save("splinematrix.txt",arma::raw_ascii);

	//cout << "into solve: " << endl; 

	//splinemat.print();
	//cout << rcond(splinemat) << endl;
	//densvec.print();

	arma::colvec sols;
	solving(&splinemat, &densvec, &sols);

	//cout << "solutions are: " << endl;	
	//sols.print();

	estpot(&estpots,&sols,&splines,gridpts,order);

	ofstream outstream;
	outstream.open(outfile);
	
	if (!state)
	{
		if (rmin == 0)
			calcpotsphere(&realpots,gridpts,rmax,Q0,&printgrid);
		else
			calcpotshell(&realpots,gridpts,rmin,rmax,Q0,&printgrid);
		
		for (uint k = 0; k < estpots.size(); k++)
			outstream << estpots[k] << "," << realpots[k] << endl;
	}		
	else
	{
		calchpot(&realpots,gridpts,n,l,&printgrid);
		
		for (uint k = 0; k < estpots.size(); k++)
			outstream << estpots[k] << "," << realpots[k] << endl;
	}


	outstream.close();
}

























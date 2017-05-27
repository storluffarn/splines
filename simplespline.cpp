
#include <iostream>
#include <vector>
#include <algorithm>	// for minmax
#include <fstream>

using namespace std;

class bspline
{
	public:
	vector<double> val;
	vector<double> d1;
	vector<double> d2;

	double calcspline(int,int,double,vector<double>*);
	double calcsplined1(int,int,double,vector<double>*);
	double calcsplined2(int,int,double,vector<double>*);

	void calcstuff(int,int,vector<double>*,vector<double>*);
};

double bspline::calcspline(int i, int k, double x, vector <double>* tlist)
{
	double ti		=	(*tlist)[i];
	double tip1		=	(*tlist)[i+1];
	double tik		=	(*tlist)[i+k];
	double tikn1	=	(*tlist)[i+k-1];
	
	if (k == 1)
	{
		if (x >= ti && x < tip1)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		double b1 = calcspline(i,k-1,x,tlist);
		double term1 = (x - ti) / (tikn1 - ti) * b1;
		
		double b2 = calcspline(i+1,k-1,x,tlist);
		double term2 = (tik - x) / (tik - tip1) * b2;
		
		if(isnan(term1) || isnan(term2))
		{
			cout << "error: recursion unexplectedly stoped at k = " << k << endl;
			cout << "first term: " << term1 << endl; 
			cout << "t_i+k-1: " << tikn1 << " t_i: " << ti << endl;
			cout << "second term: " << term2 << endl; 
			cout << "t_i+k: " << tik << " t_i+1: " << tip1 << endl;
			cout << "quotient " << (x - ti) / (tikn1 - ti) << endl;
			cout << "from term 1 iteration: " << b1 << endl;
			cout << "from term 2 iteration: " << b2 << endl;
			return 0;
		}

				
		double output = term1 + term2;

		//cout << "outputting: " << output << endl;
				
		return output;
	}
}

double bspline::calcsplined1(int i, int k, double x, vector <double>* tlist)
{
	double ti		=	(*tlist)[i];
	double tip1		=	(*tlist)[i+1];
	double tik		=	(*tlist)[i+k];
	double tikn1	=	(*tlist)[i+k-1];
	
	double term1 = 1/(tikn1 - ti)*calcspline(i,k-1,x,tlist);
	double term2 = 1/(tik - tip1)*calcspline(i+1,k-1,x,tlist);
	double output = (k-1)*(term1 + term2);

	return output;
}
	

double bspline::calcsplined2(int i, int k, double x, vector <double>* tlist)
{
	double ti		=	(*tlist)[i];
	double tip1		=	(*tlist)[i+1];
	double tip2		=	(*tlist)[i+2];
	double tik		=	(*tlist)[i+k];
	double tikn1	=	(*tlist)[i+k-1];
	double tikn2	=	(*tlist)[i+k-2];

	double term1 = 1/( (tikn1 - ti) * (tikn2 - ti) )  *  calcspline(i,k-2,x,tlist);
	double term2 = 1/( (tikn1 - ti) * (tikn1-tip1) )  *  calcspline(i,k-2,x,tlist);
	double term3 = 1/( (tik - tip1) * (tikn1-tip1) )  *  calcspline(i,k-2,x,tlist);
	double term4 = 1/( (tik - tip1) * (tik-tip2) )  *  calcspline(i,k-2,x,tlist);

	double output = (k-1)*(k-2)*(term1 - term2 - term3 + term4);

	return output;
}

void bspline::calcstuff(int i, int k, vector<double>* grid, vector<double>* tlist)
{
	for (auto& el : (*grid))
	{
		val.push_back(calcspline(i,k,el,tlist));
		d1.push_back(calcsplined1(i,k,el,tlist));
		d2.push_back(calcsplined2(i,k,el,tlist));
	}
}

int main ()
{
	vector <double> tlist = {1,1,1,1,1.5,4,5,5,5,5};
	
	auto minmax = minmax_element(begin(tlist), end(tlist));
	double min = *(minmax.first);
	double max = *(minmax.second);

	
	double gridpoints = 100;
	double diff = max - min;
	vector <double> grid;
	for (double m = 0; m < gridpoints; m++)
		grid.push_back(1+diff*m/gridpoints);
	
	for (auto& el : grid)
		cout << el << " ";
	cout << endl;
	
	int k = 4;
	int i = 4;
	
	bspline spline;
	
	for (auto& el : grid)
		 spline.val.push_back(spline.calcspline(i,k,el,&tlist));

	ofstream fs;
	fs.open("simplespline.dat");

	for (auto& el : spline.val)
		fs << el << endl;
	//cout << endl;
	
	fs.close();
}
























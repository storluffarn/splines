
//
// Header file for calculating splines
//
// use		x spline class to create a spline
//			x getval() to get pointer to vals
//			x getd1() to get pointer to first order derivatives
//			x getd2() to get pointer to second order derivatives
//			x calcstuff(i,k,&grid,&knots) to calculate splines
//			x makegrid() to st up grid
//			x readknots() reads knot points from a file
//
// For some reason beyond my comperehension this header
// uses 1-indexation rather than 0-indexation for the splines.
// 			


#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>	// for minmax
#include <fstream>
#include <functional>

using namespace std;

class bspline
{
	string source;
	int index;
	int order;
	int knotpts;
	int gridpts;
	
	vector<double> knots;
	vector<double> grid;
	
	vector<double> vals;
	vector<double> d1;
	vector<double> d2;

	double calcsplinevals(int,int,double);
	double calcsplined1(int,int,double);
	double calcsplined2(int,int,double);
	
	public:

	void readknots();
	void readknotsnorm();
	void makegrid();
	void makegridnorm();
	void knotgrid();
	void calcspline();
	void writespline();
	void printknots();
	void printgrid();
	void printvals();
	void printd1();
	void printd2();
	void clearall(){vals.clear();d1.clear();d2.clear();}
	int getknotpts(){return knotpts;}
	int getgridpts(){return gridpts;}
	double getval(double x){return vals[x];}
	double getd1(double x){return d1[x];}
	double getd2(double x){return d2[x];}
	int getindex(){return index;}
	vector <double>* getknots(){return &knots;}
	vector <double>* getgrid(){return &grid;}
	vector <double>* getvals(){return &vals;}
	vector <double>* getd1(){return &d1;}
	vector <double>* getd2(){return &d2;}
	void setindex(int i){index = i;}
	void setorder(int i){order = i;}
	void setsource(string filename){source = filename;}

	//bspline (string isource, int iindex, int iorder, int igridpts)
	//	: source(isource), index(iindex), order(iorder), gridpts(igridpts)
	//{
	//}
};

double bspline::calcsplinevals(int i, int k, double x)
{
	double ti		=	knots[i];
	double tip1		=	knots[i+1];
	double tik		=	knots[i+k];
	double tikn1	=	knots[i+k-1];
	
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
		double b1 = calcsplinevals(i,k-1,x);
		double term1 = (x - ti) / (tikn1 - ti) * b1;
		
		double b2 = calcsplinevals(i+1,k-1,x);
		double term2 = (tik - x) / (tik - tip1) * b2;

		if(isnan(term1))
			term1 = 0;
		if(isnan(term2))
			term2 = 0;
				
		double output = term1 + term2;

		return output;
	}
}

double bspline::calcsplined1(int i, int k, double x) 
{
	double ti		=	knots[i];
	double tip1		=	knots[i+1];
	double tik		=	knots[i+k];
	double tikn1	=	knots[i+k-1];
	
	double b1 = calcsplinevals(i,k-1,x);
	double b2 = calcsplinevals(i+1,k-1,x);
	double term1 = 1.0/(tikn1 - ti)*b1;
	double term2 = 1.0/(tik - tip1)*b2;	

	if(isnan(term1))
		term1 = 0;
	if(isnan(term2))
		term2 = 0;

	double output = (k-1)*(term1 - term2);

	return output;
}

double bspline::calcsplined2(int i, int k, double x)
{
	double ti		=	knots[i];
	double tip1		=	knots[i+1];
	double tip2		=	knots[i+2];
	double tik		=	knots[i+k];
	double tikn1	=	knots[i+k-1];
	double tikn2	=	knots[i+k-2];
	
	double b1 = calcsplinevals(i,k-2,x);
	double b2 = calcsplinevals(i+1,k-2,x);
	double b3 = b2;
	double b4 = calcsplinevals(i+2,k-2,x);
	
	double term1 = 1.0/( (tikn1 - ti) * (tikn2 - ti) )    *  b1;
	double term2 = 1.0/( (tikn1 - ti) * (tikn1 - tip1) )  *  b2;
	double term3 = 1.0/( (tik - tip1) * (tikn1 - tip1) )  *  b3;
	double term4 = 1.0/( (tik - tip1) * (tik-tip2) )      *  b4;
	
	if(isnan(term1))
		term1 = 0;
	if(isnan(term2))
		term2 = 0;
	if(isnan(term3))
		term3 = 0;
	if(isnan(term4))
		term4 = 0;

	double output = (k-1)*(k-2)*(term1 - term2 - term3 + term4);

	return output;
}

void bspline::calcspline()
{
	for (auto& el : grid)
	{
		vals.push_back(calcsplinevals(index,order,el));
		d1.push_back(calcsplined1(index,order,el));
		d2.push_back(calcsplined2(index,order,el));
	}
}

void bspline::readknots()
{
	double x;
	knots.reserve(knotpts + 2*(order - 1));

	ifstream readknots(source);
	while (readknots >> x)
		knots.push_back(x);

	for (int k = 0; k < order - 1; k++)
	{   
		knots.insert(knots.begin(),knots.front());
		knots.insert(knots.end(),knots.back());
	}

	knotpts = knots.size();

	//printknots();
}

void bspline::readknotsnorm()
{
	double x;
	knots.reserve(knotpts + 2*(order - 1));

	ifstream readknots(source);
	while (readknots >> x)
		knots.push_back(x);
	
	auto minmax = minmax_element(begin(knots), end(knots));
	double min = *(minmax.first);
	double max = *(minmax.second);

	for (auto& el : knots)
		el = (el - min) / (max-min);
	
	for (int k = 0; k < order - 1; k++)
	{   
		knots.insert(knots.begin(),knots.front());
		knots.insert(knots.end(),knots.back());
	}

	knotpts = knots.size(); 
}

void bspline::makegrid ()
{
	grid.reserve(gridpts);
	
	auto minmax = minmax_element(begin(knots), end(knots));
	double min = *(minmax.first);
	double max = *(minmax.second);

	double diff = max - min;

	for (double m = 0; m < gridpts; m++)
		grid.push_back(1+diff*m/gridpts);

	gridpts = grid.size();
}

void bspline::makegridnorm ()
{
	grid.reserve(gridpts);

	for (double m = 0; m < gridpts; m++)
		grid.push_back(m/gridpts);

	gridpts = grid.size();
}

void bspline::knotgrid()
{
	grid = knots;
	gridpts = knotpts;
}

void bspline::printknots()
{
	cout << "content in knot vector: " << endl;
	for (auto& el : knots)
		cout << el << " ";
	cout << endl;
}

void bspline::printgrid()
{
	cout << "content in grid vector: " << endl;
	for (auto& el : grid)
		cout << el << " ";
	cout << endl;
}

void bspline::printvals()
{
	cout << "content in vals vector: " << endl;
	for (auto& el : vals)
		cout << el << " ";
	cout << endl;
}

void bspline::writespline()
{
	ostringstream namestream;
	
	namestream << "B(" << fixed << setprecision(1) << index << "," << fixed << setprecision(1) << order << ").csv";
	string filename = namestream.str();
	
	ofstream fs;
	fs.open(filename);

	for (int k = 0; k < gridpts; k++)
		 fs << vals[k] << ","
			<< d1[k] << ","
			<< d2[k] << endl;
	fs.close();
	
	cout << "wrote " << gridpts << " numbers to the file: " << filename << endl;
}




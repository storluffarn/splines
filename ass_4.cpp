
#include <iostream>
#include <vector>
#include <algorithm>	// for minmax
#include <armadillo>


using namespace std;

double pi = 4*atan(1);

class bspline
{
	int num;
	double val;
	double d1;
	double d2;
	
	double calcsplinesub(int,int,double,vector<double>*,int);
	
	public:	
	
	void calcspline(int,int,double,vector<double>*);

	double getval(){return val;}
	double getd1(){return d1;}
	double getd2(){return d2;}

	bspline(int inum, int ii, int ik, double ir, vector <double>* tlist)
		: num(inum)
	{
		val = calcsplinesub(ii,ik,ir,tlist,1);
		d1 = calcsplinesub(ii,ik,ir,tlist,2);
		d2 = calcsplinesub(ii,ik,ir,tlist,3);
	}
};


double bspline::calcsplinesub(int i, int k, double x, vector <double>* tlist, int mode)
{
	double ti		=	(*tlist)[i];
	double tip1		=	(*tlist)[i+1];
	double tip2		=	(*tlist)[i+2];
	double tik		=	(*tlist)[i+k];
	double tikn1	=	(*tlist)[i+k-1];
	double tikn2	=	(*tlist)[i+k-2];
	double tikp1	=	(*tlist)[i+k+1];
	
	//cout << "k: " <<  k << endl;
	//cout << "ti: " << ti << endl << "ti+1: " << tip1 << endl << "tik: " << tik << endl << "tik+1: " << tikp1 << endl;

	if (ti == tip1)
		return 0;

	if (k == 0)
	{
		//cout << "ti " << ti << " x " << x << " tip1 " << tip1 << endl;
		if (ti == tip1 || ti > tip1)
		{
			cout << "error: bad comparison ti < ti+1 is false" << " ti is: " << ti << " ti+1 is: " << tip1 << endl; 
			return 0;
		}
		else if (x >= ti && x < tip1)
		{
			//cout << "accepted piece with ti: " << ti <<  " x: " << x  << " ti+1: " << tip1 << endl;
			return 1;
		}
		else
		{
		//cout << "rejected piece with: ti: " << ti << " x: " << x << " ti+1: " << tip1 << endl;
			return 0;
		}
	}
	else
	{
		switch (mode)
		{
			case 1: 
			{
				double out = (x - ti) / (tik - ti) * calcsplinesub(i,k-1,x,tlist,1)  +  
							 (tikp1 - x) / (tikp1 - tip1) * calcsplinesub(i+1,k-1,x,tlist,1);

				//cout << out << endl; 
				return out;
			}

			case 2: 
			{
				double out = (k - 1) * (calcsplinesub(i,k-1,x,tlist,1) / (tikn1 - ti)  -  
										calcsplinesub(i+1,k-1,x,tlist,1) / (tik - tip1));

				return out;
			}

			case 3: 
			{	
				double out =  (k - 1) * (k - 2) * (
				calcsplinesub(i,k-2,x,tlist,1)	 /	((tikn1 - ti)   * (tikn2 - ti))		-  
				calcsplinesub(i+1,k-2,x,tlist,1) /	((tikn1 - ti)   * (tikn1 - tip1))	-
				calcsplinesub(i+1,k-2,x,tlist,1) /	((tik   - tip1) * (tikn1 - tip1))	+
				calcsplinesub(i,k-2,x,tlist,1)	 /	((tik   - tip1) * (tik   - tip2))	);
				
				return out;
			}
		}
		
		cout << "error: recursion unexplectedly stoped at number " << k << endl;
		return 0;
	}
}

double spherecharge(double r)	// assume sphere or radius 1
{
	if (r < 0)
		cout << "error: negative radii not allowed" << endl;

	if (r > 1)
		return -4*pi*1.0/(4.0/3.0*pi);
	else 
		return -r*4*pi*1.0/(4.0/3.0*pi*pow(r,3));
}

double shellcharge(double r, double r0)		// assume shell of outer radius 1
{
	if (r < 0)
		cout << "error: negative radii not allowed" << endl;
	
	if (r > 1)
		return -4*pi*1.0/(4.0/3.0*pi);
	else if (r < r0)
		return 0;
	else 
		return -r*4*pi*1.0/(4.0/3.0*pi*r);
}

double testfunc(double r){return r*r;}

double scaling(double x, double min, double max)
{
	return (x - min) / (max - min);
}

int main ()
{
	double order = 4;
	unsigned int points = 10;

	vector <double> tlist;
	for (unsigned int k = 0; k < points; k++)
		tlist.push_back(testfunc(k));
		//tlist.push_back(spherecharge((double)(k)/points));
	
	for (auto& el : tlist)
		cout << el << " ";
	cout << endl;

	auto minmax = minmax_element(begin(tlist), end(tlist));
	double min = *(minmax.first);
	double max = *(minmax.second);

	for (auto& el : tlist)
		el = scaling(el,min,max);

	//cout << min << " " << max << endl;
	//for (auto& el : tlist)
	//	cout << el << " ";
	//cout << endl;
	
	for (int k = 0; k < order - 1; k++)			// ghost points
		{
			tlist.insert(tlist.begin(),tlist.front());
			tlist.push_back(tlist.back()); 
		}
	
	for (auto& el : tlist)
		cout << el << " ";
	cout << endl;
	
	// Sanity check NOTE need k = k to points + 1 in first loop.
	//bspline a = bspline(0,8,4,tlist[9],&tlist);
	//cout << a.getval() << " " << a.getd1() << " " << a.getd2() << endl;

	vector <bspline> splines;
	int m = 0;
	for (unsigned int t = order - 1; t < tlist.size() - (order - 1); t++)	// only physical points hence k=1;
		for (unsigned int i = t; i < t + 3; i++)
		{
			int n = i - order + 1;
			splines.push_back(bspline(m,n,order,tlist[t],&tlist));
			m++;
			//cout << "B(" << n << "," << order << "," << t << ")" << endl;
		}
	
	for (auto& el : splines)
		cout << el.getd1() << " ";
	cout << endl;

//	arma::mat sys(points,points,arma::fill::zeros);
//	
//	m = 0;
//	for (unsigned int k = 0; k < points; k++)
//		for (unsigned int l = k; l < k + order - 1; l++)
//		{	
//			sys(k,l) = splines[m].getval();
//			m++;
//		}
//
//	sys.print();
}
























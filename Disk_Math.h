#define PI 3.14159265358979323846

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex>
#include<iostream>
#include<cstdlib>
#include<time.h>
#include<algorithm>
#include<iostream>

using namespace std;  


double Fact(double n)					//Factorial, invalid for n>30, for it leaks out.
{
	double s=1;
	if((double)(n)<0.01&&(double)(n)>-0.01)
		return 1.;
	if(n<0.0)
		return 0.;
	for(double i=n;i>0;i=i-1.0)
	{
		s=s*i;
	}
	
	return s;
}
double Fact(double n, double l)		//overloading Factorial n(n-1)...(n-l+1)
{
	double s=1;
	if((double)(n)<0.01&&(double)(n)>-0.01||l==0.)
		return 1.;
	if(n<0.0||n<l)
		return 0.;
	for(double i=0;i<l;i=i+1.0)
	{
		s=s*(n-i);
	}
	
	return s;
}
double C(double up, double low)			//Number of combinations
{
	if(low<0)
		return 0;
	if(up<low)
		return 0;
	if(up==low&&up==0)
		return 1;
	double multi=1;
	for(double i=up;i>up-low;i--)
	{
		multi=multi*i;
	}
	multi=multi/Fact(low);
	return multi;
	
}

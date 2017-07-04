#include<Polynomial.h>
double Nm(int n, int m) //Normalization factor
{
        return pow(-1,n)*sqrt(Fact(n)/(2.*PI*pow(2,m)*Fact(n+m)));
}
/**
class Jas(): public Polynomial	//Jasfactor for variant z
{
	void Jas(Cdouble z*, int n);
};
void Jas::
**/

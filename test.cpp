//#include <Polynomial.h>
//#include "Polynomial.h"
//#include "LA1.h"
#include "Monca.h"

using namespace std;
void test()  
{  
    Polynomial p1;  
	Cdouble c1(3.,0.);
	Cdouble c2(4.,0);
	Cdouble d1(3.,4);
	p1.NewTerm(c1,0);
	p1.NewTerm(c2,1);
    Polynomial p2;	 
    p2.NewTerm(d1,0);  //cout<<p2;

    cout<<"("<<p1<<") * ("<<p2<<") = "<<p1*p2<<endl;  
}  

int main()
{
	Monca experiment(6,3,1);
	qmnumber qm[6];
	qm[0].n=0; qm[0].m=0;
	qm[1].n=0; qm[1].m=1;
	qm[2].n=0; qm[2].m=2;
	qm[3].n=0; qm[3].m=3;
	qm[4].n=1; qm[4].m=-1;
	qm[5].n=1; qm[5].m=0;
	cout<<"6,3,1   "<<experiment.Metrop(2000000,qm)<<endl;
}

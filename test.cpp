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
	//test();
	Monca experiment(6,3,1);
	cout<<experiment.Metrop(500000)<<endl;
	system("Pause");
}

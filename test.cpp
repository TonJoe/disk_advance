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
//    p1.NewTerm(6.5,2);  
//	p1.NewTerm(2,3);
	p1.NewTerm(c1,0);
	p1.NewTerm(c2,1);
    //p1.NewTerm(2.1,3);  
	//cout<<p1.Eval(d1 )<<"       "<<endl;
    Polynomial p2;	
  //  p2.NewTerm(1,2);  
  //  p2.NewTerm(1,3);  
    p2.NewTerm(d1,0);  //cout<<p2;

    cout<<"("<<p1<<") * ("<<p2<<") = "<<p1*p2<<endl;  
    //cout<<"F(x=2) = "<<(p1+p2).Eval(2)<<endl;  
    //cout<<"("<<p1<<")+("<<p2<<") = "<<p1 + p2<<endl;
	//Polynomial p2(p1);
	//cout<<"p1="<<p1<<endl;
	//cout<<p1.Deriv()<<endl;
	//cout<<p1.Deriv(1)<<endl;
	//cout<<p1.Deriv(2)<<endl;
	

}  

int main()
{
	//test();
	Monca experiment(6,2,1);
	//experiment.Metrop(1000);
	cout<<experiment.Metrop(500000)<<endl;
	system("Pause");
}

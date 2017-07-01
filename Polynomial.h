#define PI 3.14159265358979323846
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex>
#include<iostream>
#include<cstdlib>
#include<time.h>
#include<algorithm>
using namespace std;
typedef complex<double> Cdouble;

class Polynomial;
class Term	//Term in a polynomial.
{
	friend Polynomial;
public:
	double coef;
	int exp;
};
class Polynomial
{
	friend ostream & operator<<(ostream &o,const Polynomial & poly);  
public:
	Polynomial();
	Polynomial(const Polynomial & poly);
	~Polynomial();
	Polynomial operator+(const Polynomial & poly); //Polynomial addition.
	Polynomial operator*(const Polynomial & poly);	
	Cdouble Eval(Cdouble z);	//Evaluate the polynomial;
	void NewTerm(double coef, int exp);
private:
	void InsertTerm(const Term & term);
private:
	Term *termArray;	//Term Array
	int capacity;		//Capacity of the polynomial
	int terms;			//Number of Non-zero terms
};

Polynomial::Polynomial()	//Initializaing a polynomial with 10 terms.
{
	this->terms=0;
	this->capacity=10;
	termArray= new Term[this->capacity];
}

Polynomial::Polynomial(const Polynomial & b)  
{  
    this->terms=0;  
    this->capacity=b.capacity;  
    termArray = new Term[this->capacity];  
    for(int i=0;i<b.terms;i++)
	{  
        NewTerm(b.termArray[i].coef,b.termArray[i].exp);  
    }  
}  
  
Polynomial::~Polynomial()  
{  
    delete [] termArray;  
}  

Polynomial Polynomial::operator+(const Polynomial & b)  
{  
    Polynomial c;  
    int aPos=0;  
    int bPos=0;  
    while(aPos<terms && bPos<b.terms){  
        if(termArray[aPos].exp == b.termArray[bPos].exp){  
            float coef=termArray[aPos].coef+b.termArray[bPos].coef;  
            if(coef)c.NewTerm(coef,termArray[aPos].exp);  
            aPos++;bPos++;  
        }else if(termArray[bPos].exp < b.termArray[bPos].exp){  
            c.NewTerm(b.termArray[bPos].coef,b.termArray[bPos].exp);  
            bPos++;  
        }else{  
            c.NewTerm(termArray[aPos].coef,termArray[aPos].exp);  
            aPos++;  
        }  
    }  
    while (aPos < terms){  
        c.NewTerm(termArray[aPos].coef,termArray[aPos].exp);  
        aPos++;  
    }  
    while (bPos < b.terms){  
        c.NewTerm(b.termArray[bPos].coef,b.termArray[bPos].exp);  
        bPos++;  
    }  
    return c;  
}  
Polynomial Polynomial::operator*(const Polynomial & b)  
{  
    Polynomial c;  
    for(int i=0; i<terms; i++){  
        for(int j=0; j<b.terms; j++){  
            float coef = termArray[i].coef*b.termArray[j].coef;  
            int exp = termArray[i].exp + b.termArray[j].exp;  
            c.NewTerm(coef,exp);  
        }  
    }  
    return c;  
}
void Polynomial::NewTerm(double coef, int exp)  
{  
    if(terms == capacity){  
        capacity *= 2;  
        Term *tmp = new Term[capacity];  
        copy(termArray,termArray+terms,tmp);  
        delete [] termArray;  
        termArray = tmp;  
    }  
    Term ATerm;  
    ATerm.coef=coef;ATerm.exp=exp;  
    InsertTerm(ATerm);  
} 
void Polynomial::InsertTerm(const Term & term)  
{  
    int i;  
    for(i=0; i<terms && term.exp<termArray[i].exp; i++){  
    }  
    if(term.exp == termArray[i].exp){  
        termArray[i].coef += term.coef;  
        if(!termArray[i].coef){  
            for(int j=i; j<terms-1; j++)  
                termArray[j]= termArray[j+1];  
            terms--;  
        }  
    }else{  
        for(int j=terms-1; j>=i;j--)  
            termArray[j+1]=termArray[j];  
        termArray[i] = term;  
        terms++;  
    }  
}  
Cdouble Polynomial::Eval(Cdouble z)  
{  
    Cdouble res=polar(0.,0.);  
    for(int i=0;i<terms; i++){  
        res =res+ termArray[i].coef * pow(z,termArray[i].exp);  
    }  
    return res;  
}

#include<Disk_Math.h>

using namespace std;
typedef complex<double> Cdouble;
class Monca;
class Polynomial;
class Term	//Term in a polynomial.
{
	friend Polynomial;
public:
	Cdouble coef;
	int exp;
};
class Polynomial
{
	friend ostream & operator<<(ostream &o,const Polynomial & poly);
	friend Monca;
	//friend Polynomial Der(const Polynomial & poly, int n);
public:
	Polynomial();
	Polynomial(const Polynomial & poly);
	~Polynomial();
	Polynomial& operator=(const Polynomial & poly);
	Polynomial operator+(const Polynomial & poly); //Polynomial addition.
	Polynomial operator*(const Polynomial & poly);
	
	void Clear1();
	void Clear0();
	Polynomial Deriv();
	Polynomial Deriv(int n);
	void Jas(Cdouble *z, int length, int i);	//This is the methodd to construct a Jas factor from a series of coordinates z*, length 'length', the center coordinate z_i. J_i=\sum(z_i-z_k),k!=j
	
	//void operator=(const Polynomial & poly);
	Cdouble Eval(Cdouble z);	//Evaluate the polynomial;
	void NewTerm(Cdouble coef, int exp);
protected:
	void InsertTerm(const Term & term);
protected:
	Term *termArray;	//Term Array
	int capacity;		//Capacity of the polynomial
	int terms;			//Number of Non-zero terms
};

Polynomial::Polynomial()	//Initializaing a polynomial with 5 terms.
{
	this->terms=0;
	this->capacity=5;
	termArray= new Term[this->capacity];
	for(int i=0;i<capacity;i++)
	{
		termArray[i].coef=polar(0.,0.);
		termArray[i].exp=0;

	}
}

Polynomial::Polynomial(const Polynomial & b)  //Initializaing a polynomial copying from b
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
void Polynomial::Clear1()
{
	//if(termArray!=NULL)
		
	delete [] termArray;
	
	this->terms=0;
	this->capacity=5;
	this->termArray = new Term[capacity];
	for(int i=0;i<capacity;i++)
	{
		termArray[i].coef=polar(0.,0.);
		termArray[i].exp=0;

	}
	Cdouble c(1.,0.);
	this->NewTerm(c,0);	
}
void Polynomial::Clear0()
{
	//if(termArray!=NULL)
		
	delete [] termArray;
	
	this->terms=0;
	this->capacity=5;
	this->termArray = new Term[capacity];
	for(int i=0;i<capacity;i++)
	{
		termArray[i].coef=polar(0.,0.);
		termArray[i].exp=0;

	}
}




Polynomial Polynomial::Deriv()	//Derivative
{
	Polynomial c;
	//int aPos=0;
	int bPos=0;
	while(bPos<this->terms)
	{
		Cdouble bexp(this->termArray[bPos].exp,0);
		Cdouble coef(this->termArray[bPos].coef*bexp); cout<<coef<<norm(coef)<<"  ";
		if(norm(coef)>0.0000001)
		{
			c.NewTerm(coef,this->termArray[bPos].exp-1);
			cout<<"A new term created";
			cout<<c<<endl;
		}
		bPos++;
	}
	return c;
}

Polynomial Polynomial::Deriv(int n)		//n'th order derivative
{
	Polynomial c;
	//int aPos=0;
	int bPos=0;
	while(bPos<this->terms)
	{
		Cdouble bexp(this->termArray[bPos].exp,0.);
		
		if(this->termArray[bPos].exp>0)		//positive exponent
		{
			if(this->termArray[bPos].exp>=n)
			{
				double fact=Fact(this->termArray[bPos].exp,double(n));
				Cdouble coef(fact,0.);
				coef=coef*this->termArray[bPos].coef;
				if(norm(coef)>0.0000001)
					c.NewTerm(coef,this->termArray[bPos].exp-n);
				bPos++;
			}
			else	//Derivative order> exp, term 0. Next term
				bPos++;			
		}
		else							//negative exponent
		{
			double fact=Fact(abs(this->termArray[bPos].exp)+double(n-1),double(n))*pow(-1.,n);
			Cdouble coef(fact,0.);
			coef=coef*this->termArray[bPos].coef;
			c.NewTerm(coef,this->termArray[bPos].exp-n);
			bPos++;
		}
	}
	return c;
}

void Polynomial::Jas(Cdouble* z, int length, int i)
{
	Clear1();
	cout<<"this"<<*this<<endl;
	Cdouble coef(1.,0);
	for(int n=0;n<length;n++)
	{
		if(n!=i)
		{
			Polynomial c;
			c.Clear1();
			c.NewTerm(coef,1);
			c.NewTerm(-z[n],0);
			cout<<"term"<<c<<endl;
			(*this)=(*this)*c;
			cout<<"  now"<<(*this)<<endl;
		}
	}
}

/******Basic Functions*******/
Polynomial& Polynomial::operator=(const Polynomial & poly)
{
	if(termArray!=NULL)
	{
		delete [] termArray;
	}
	this->capacity=poly.capacity;
	this->terms=poly.terms;
	this->termArray = new Term[this->capacity];
	for(int i=0;i<terms;i++)
	{
		this->termArray[i].coef=poly.termArray[i].coef;
		this->termArray[i].exp=poly.termArray[i].exp;
	}
	return *this;
}

Polynomial Polynomial::operator+(const Polynomial & b)  
{  
    Polynomial c;  
    int aPos=0;  
    int bPos=0;  
    while(aPos<terms && bPos<b.terms){  
        if(termArray[aPos].exp == b.termArray[bPos].exp){  
            Cdouble coef=termArray[aPos].coef+b.termArray[bPos].coef;  
            if(norm(coef)>0.00001)c.NewTerm(coef,termArray[aPos].exp);  
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
            Cdouble coef = termArray[i].coef*b.termArray[j].coef;  
            int exp = termArray[i].exp + b.termArray[j].exp;  
            c.NewTerm(coef,exp);  
        }  
    }  
    return c;  
}

void Polynomial::NewTerm(Cdouble coef, int exp)  
{  
    if(terms >= capacity){  
        capacity *= 2;  
        Term *tmp = new Term[capacity];
		for(int i=0;i<capacity;i++)
		{
			tmp[i].coef=polar(0.,0.);
			tmp[i].exp=0;
		}
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
        if(termArray[i].coef==polar(0.,0.)){  
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

ostream & operator<<(ostream & o,const Polynomial & poly)  
{  
    for(int i=0;i<poly.capacity-1;i++){  
        o<<poly.termArray[i].coef<<"x^"<<poly.termArray[i].exp<<" + ";  
    }  
    o<<poly.termArray[poly.capacity-1].coef<<"x^"<<poly.termArray[poly.capacity-1].exp;  
    return o;  
}

/********************************/
/***Other stuff
/********************************/

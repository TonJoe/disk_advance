#include<Polynomial.h>
/**
class Monca
{
public:
	Polynomial *LY;			//Landaul level
	Polynomial *JAS;		//Jastrow Factor
	Polynomial *ynm;		//Projected landau level function
public:
	//Monca();
	//~Monca();
	
	void Build();
	Cdouble CF_Wave(Cdouble z*);
	void Metrop(int N_p, int N, int P, int steps);
	
	void G_Metrop(int n_p, int n, int p, int steps);	//Ground state Metropolis. particle number, filling factor n, flux number p,(v=n/(2pn+1)), Metropolis steps
	void Ex_Metrop(int n_p, int n, int p, int steps);	//Exciton Metropolis
	Cdouble *z;
	int n_p,n,p;
	double RN;
};
void Monca::Build()
{
	JAS=new Polynomial[n_p];
	ynm=new Polynomial[n_p];
	for(int i=0;i<n_p;i++)
	{
		JAS[i].Clear1();
		JAS[i].Jas(z,n_p,i);
		ynm[i].Clear0();
	}
	int i=0;
	for(int lvl=0;lvl<n;lvl++)
	{
		for(int m=-lvl;m<-lvl+n_p/n&&i<n_p;m++)
		{
			for(int k=0;k<lvl;k++)
			{
				Polynomial tmp1,tmp2;
				Cdouble c1(double(k+m),0.), c2(pow(-1.,k)*C(lvl+m,lvl-k)/Fact(k),0.);
				tmp1.NewTerm(c1,k+m);
				tmp2.NewTerm(c2,0);
				ynm[i]=ynm[i]+tmp2*(tmp1*JAS[i]).Deriv(k);
			}
			i++;
		}
	}
}

Cdouble CF_Wave(Cdouble z*)	//Calculate determinant wave-funtion
{
	Cdouble product(1.0,0.0);
	Cdouble sum(0.,0.);
	for(int j=0;j<n_p;j++)
		{
			for(int i=0;i<n_p;i++)
			{
				product=product*ynm[(i+j)%n_p].Eval(z[i]);
			}
			sum=sum+product;
		}
	///	for
}

void Monca::Metrop(int N_p, int N, int P, int steps)
{
	n_p=N_p;
	n=N;
	p=P;
	Cdouble r[n_p];
	
	for(int i=0;i<n;i++)
	{
		z[i]=polar((double)(double(rand())/double(RAND_MAX))*RN,(double)(double(rand())/double(RAND_MAX))*2.*PI);
		//cout<<z[i]<<endl;
	}
	for(int st=0;st<steps;st++)
	{
		
		
		for(int j=0;j<n;j++)	//Copy to a buffer
		{
			r[j]=z[j];
		}		
		r[st%n_p]=r[st%n_p]+polar((double(rand())/double(RAND_MAX))*0.2*RN,(double(rand())/double(RAND_MAX))*2.*PI);
		if(norm(CF_Wave[r]/CF_Wave[z])>1.)
	}
}**/


#include<Polynomial.h>
#include<LA1.h>
#include<det.h>
double EllipticE(double x);
class Monca
{
public:
	///Polynomial *LY;			//Landaul level
	Polynomial *JAS;		//Jastrow Factor
	Polynomial *ynm;		//Projected landau level function
public:
	//Monca();
	Monca(int N_p, int N, int P);
	//~Monca();
	
	void DoJas();
	void Build();
	Cdouble CF_Wave(const Cdouble *z);
	double Metrop(int steps);
	double Vee(int n, complex<double> *z);
	double Vbb(int n, double niu);
	double Vbe(int n, double niu, double RN, complex<double> *z);

	
	//void G_Metrop(int n_p, int n, int p, int steps);	//Ground state Metropolis. particle number, filling factor n, flux number p,(v=n/(2pn+1)), Metropolis steps
	//void Ex_Metrop(int n_p, int n, int p, int steps);	//Exciton Metropolis
	
	Cdouble *z;
	int n_p,n,p;
	double niu,RN;
};
Monca::Monca(int N_p, int N, int P)
{
	n_p=N_p;
	n=N;
	p=P;
	niu=double(n)/(2.*n*n_p+1.);
	RN=sqrt(2.0*n_p/niu);
	
	z=new Cdouble[n_p];
	///LY=new Polynomial [n_p*n_p];
	ynm=new Polynomial [n_p*n_p];
	JAS=new Polynomial [n_p];
}
void Monca::DoJas()
{
	if(JAS==NULL)
	{
		cout<<"Void Jastrow."<<endl;
		return;
	}
	delete [] ynm;
	ynm=new Polynomial[n_p*n_p];
	for(int i=0;i<n_p*n_p;i++)
		ynm[i].Clear0();
	int l=0,m=0,t=0;
	int i,j;
	for(i=0;i<n_p;i++)	//i'th row
	{
		for(j=0;j<n_p;j++,t++)	//j'th column
		{
			l=i/(n_p/n);
			m=i%(n_p/n)-l;
			/**Now start calculate wave function**/
			
			for(int k=0;k<=n;k++)
			{
				double c;
				c=pow(-1.,k)*C(l+m,l-k)/Fact(k);
				Cdouble coef(c,0.);
				Polynomial tmp;
				tmp.Clear0();
				tmp.NewTerm(coef,k+m);
				ynm[i*n_p+j]=ynm[i*n_p+j]+(tmp*JAS[j]).Deriv(k);
				
			}
			//cout<<endl<<"t="<<t<<"  ynm("<<l<<","<<m<<")J"<<j<<"="<<ynm[i*n_p+j]<<endl;
		}
	}
}

void Monca::Build()	//This gives single CF_wave_function.
{
	//JAS=new Polynomial[n_p];
	//ynm=new Polynomial[n_p];
	for(int i=0;i<n_p;i++)
	{
		//cout<<"z[i]="<<z[i]<<endl;
		//cout<<n_p<<endl;
		JAS[i].Clear1();
		JAS[i].Jas(z,n_p,i); 
		
		//cout<<JAS[i]<<endl;
		
		
	}
	DoJas();
}

Cdouble Monca::CF_Wave(const Cdouble *z)	//Calculate determinant of wave-funtion: many-body function
{
	Cdouble matrix[n_p*n_p];
	double sum=0.;
	for(int i=0;i<n_p;i++)	// i'th function
	{
		for(int j=0;j<n_p;j++)	//j'th coordinate
		{
			//cout<<ynm[i]<<endl;
			matrix[i*n_p+j]=ynm[i*n_p+j].Eval(z[j]);
		}
	}
	for(int i=0;i<n_p;i++)
	{
		sum=sum+norm(z[i]);
	}
	sum=exp(-0.25*sum);
	Cdouble ep(sum,0.);
	//return polar(1.,0.);
	return cpxdbl_det0(matrix, n_p)*ep;
}

double Monca::Metrop(int steps)
{
	cout<<"The current program only applies for the case p=1, since I didn't make the exponent of Jastrow."<<endl;
	Cdouble r[n_p];
	double Energy=0;
	
	srand (time(NULL));
	for(int i=0;i<n_p;i++)
	{
		z[i]=polar((double)(double(rand())/double(RAND_MAX))*RN,(double)(double(rand())/double(RAND_MAX))*2.*PI);
		cout<<z[i]<<endl;  
		
	}
	/**
	for(int i=0;i<n_p;i++)
	{
		JAS[i].Jas(z,n_p,i); 
		cout<<JAS[i]<<endl;
	}
	DoJas();
	//cout<<CF_Wave(z)<<endl;
	**/
	for(int st=0;st<steps;st++)
	{
		
		Build();
		
		for(int j=0;j<n_p;j++)	//Copy to a buffer
		{
			r[j]=z[j]; //cout<<r[j]<<"  "<<z[j]<<endl;
		}	
		//cout<<norm(CF_Wave(z))<<"  "<<norm(CF_Wave(r))<<endl;
		r[st%n_p]=r[st%n_p]+polar((double(rand())/double(RAND_MAX))*0.15*RN,(double(rand())/double(RAND_MAX))*2.*PI);
		//cout<<norm(CF_Wave(z))<<"  "<<norm(CF_Wave(r))<<endl<<endl;
		if(norm(CF_Wave(r)/CF_Wave(z))>(double(rand())/double(RAND_MAX))&&norm(r[st%n_p])<RN*RN)
		{	
			//cout<<"accept"<<endl;
			z[st%n_p]=r[st%n_p];
		}
		else
		{
			//cout<<"refuse"<<endl;
		}
		Energy=Energy+Vee(n_p, z);
		
	}
	return Energy/steps;
}


/***********************************/
/**Three potentials: e-e, e-b, b-b**/
/***********************************/
double Monca::Vee(int n, complex<double> *z)
{
	double V=0;
	for(int i=0;i<n;i++)
		for(int j=i+1;j<n;j++)
		{
			V=V+1./abs(z[i]-z[j]);
		}
	return V;
}
double Monca::Vbb(int n, double niu)
{
	return n*8./3./PI*sqrt(niu*n/2);
}
double Monca::Vbe(int n, double niu, double RN, complex<double> *z)
{
	double sum=0.0;
	for(int i=0;i<n;i++)
	{
		sum=sum+EllipticE((abs(z[i])/RN)*(abs(z[i])/RN));
	}
	return -sqrt(2.*niu*n)*sum*2/PI;
}

double EllipticE(double x)	//Elliptic function, for the use of the evaluation for some wave functions.
{
	return PI/2.-PI/8.*x-3.*PI/128.*x*x-5.*PI/512.*pow(x,3)-175.*PI/32768.*pow(x,4)-441.*PI/131072.*pow(x,5)-4851*PI/2097152.*pow(x,6)-14157.*PI/8388608.*pow(x,7)-2760615.*PI/2147483648.*pow(x,8)-8690825.*PI/8589934592*pow(x,9)-112285459.*PI/137438953472.*pow(x,10);

}
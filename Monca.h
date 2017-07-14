
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
	Polynomial *tmpJAS;		//temp polynomial, for the buffer in Monte-Carlo
	Polynomial *tmpynm;
	
public:
	//Monca();
	Monca(int N_p, int N, int P);
	//~Monca();
	
	void DoJas(const Polynomial *JAS, Polynomial *ynm);
	void Build(Polynomial *JAS, Polynomial *ynm, const Cdouble *z);
	Cdouble CF_Wave( Polynomial *tynm, const Cdouble *tz);
	double Metrop(int steps);
	double Vee(int n, complex<double> *z);
	double Vbb(int n, double niu);
	double Vbe(int n, double niu, double RN, complex<double> *z);

	
	//void G_Metrop(int n_p, int n, int p, int steps);	//Ground state Metropolis. particle number, filling factor n, flux number p,(v=n/(2pn+1)), Metropolis steps
	//void Ex_Metrop(int n_p, int n, int p, int steps);	//Exciton Metropolis
	
	Cdouble *z,*r;
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
	r=new Cdouble[n_p];
	///LY=new Polynomial [n_p*n_p];
	ynm=new Polynomial [n_p*n_p];	
	tmpynm=new Polynomial [n_p*n_p];
	
	JAS=new Polynomial [n_p];
	tmpJAS=new Polynomial [n_p];
}
void Monca::DoJas(const Polynomial *tJAS, Polynomial *tynm)
{
	if(tJAS==NULL)
	{
		cout<<"Void Jastrow."<<endl;
		return;
	}
	
	if(tynm!=NULL)
	{
	for(int i=0;i<n_p*n_p;i++)
	{
		tynm[i].Clear0();
		//cout<<tynm[i];
	}
	int l=0,m=0;
	int i,j;
	for(i=0;i<n_p;i++)	//i'th row
	{
		for(j=0;j<n_p;j++)	//j'th column
		{
			l=i/(n_p/n);
			m=i%(n_p/n)-l;
			/**Now start calculating wave function**/
			
			for(int k=0;k<=n;k++)
			{
				double c;
				c=pow(-1.,k)*C(l+m,l-k)/Fact(k);
				Cdouble coef(c,0.);
				Polynomial tmp;
				tmp.Clear0();
				tmp.NewTerm(coef,k+m);
				tynm[i*n_p+j]=tynm[i*n_p+j]+(tmp*tJAS[j]).Deriv(k);
				
			}
			//cout<<endl<<"t="<<t<<"  ynm("<<l<<","<<m<<")J"<<j<<"="<<ynm[i*n_p+j]<<endl;
		}
	}
	}
}

void Monca::Build(Polynomial *tJAS, Polynomial *tynm, const Cdouble *tz)	//This gives single CF_wave_function.
{
	//JAS=new Polynomial[n_p];
	//ynm=new Polynomial[n_p];
	for(int i=0;i<n_p;i++)
	{
		//tJAS[i].Clear1();
		tJAS[i].Jas(tz,n_p,i); 
	}
	DoJas(tJAS, tynm);
}

Cdouble Monca::CF_Wave( Polynomial *tynm, const Cdouble *tz)	//Calculate determinant of wave-funtion: many-body function
{
	Cdouble *matrix;
	matrix=new Cdouble[n_p*n_p];
	double sum=0.;
	for(int i=0;i<n_p;i++)	// i'th function
	{
		for(int j=0;j<n_p;j++)	//j'th coordinate
		{
			//cout<<ynm[i]<<endl;
			matrix[i*n_p+j]=tynm[i*n_p+j].Eval(tz[j]);
			//cout<<tynm[i*n_p+j]<<endl;
		}
	}
	for(int i=0;i<n_p;i++)
	{
		sum=sum+norm(tz[i]);
	}
	sum=exp(-0.25*sum);
	Cdouble ep(sum,0.);
	//return polar(1.,0.);
	return cpxdbl_det0(matrix, n_p)*ep;
}

double Monca::Metrop(int steps)
{
	cout<<"The current program only applies for the case p=1, since I didn't make the exponent of Jastrow."<<endl;
	double Energy=0;
	
	srand (time(NULL));
	for(int i=0;i<n_p;i++)
	{
		z[i]=polar((double)(double(rand())/double(RAND_MAX))*RN,(double)(double(rand())/double(RAND_MAX))*2.*PI);
		r[i]=z[i];
		//cout<<z[i]<<endl;  
		
	}
	
	
	for(int st=0;st<steps;st++)
	{	
		//for(int i=0;i<n_p;i++)cout<<"&&&&&&&&&&&&&&&&&"<<z[i]-r[i]<<endl;
		
		for(int j=0;j<n_p;j++)	//Copy to a buffer
		{
			r[j]=z[j]; 
		}	
		
		r[st%n_p]=r[st%n_p]+polar((double(rand())/double(RAND_MAX))*0.15*RN,(double(rand())/double(RAND_MAX))*2.*PI);
		
		Build(JAS, ynm, z);
		Build(tmpJAS, tmpynm, r);
		cout<<RAND_MAX<<" Step:"<<st<<":: "<<norm(CF_Wave(ynm, z))<<"  "<<norm(CF_Wave(tmpynm,r))<<endl;
		if(norm(CF_Wave(tmpynm, r)/CF_Wave(ynm, z))>(double(rand())/double(RAND_MAX))&&norm(r[st%n_p])<RN*RN)
		{	
			cout<<"accept:"<<st<<endl;
			
			z[st%n_p]=r[st%n_p];
			
			//Build(JAS, ynm, z);
		}
		else
		{
			cout<<"refuse:"<<st<<endl;
		}
		
		Energy=Energy+Vee(n_p, z);
	}
	//cout<<"Now error sofar"<<endl;
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
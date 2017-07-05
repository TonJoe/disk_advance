double Nm(int n, int m) //Normalization factor
{
        return pow(-1,n)*sqrt(Fract(n)/(2.*PI*pow(2,m)*Fract(n+m)));
}
class Monca
{
public:
	Polynomial *LY;			//Landaul level
	Polynomial *JAS;		//Jastrow Factor
		Polynomial *ynm;		//Projected landau level function
public:
	//Monca();
	//~Monca();
	
	void Build(int N_P, int N, int P);
	void Symmtrize();
	void G_Metrop(int n_p; int n; int p; int steps);	//Ground state Metropolis. particle number, filling factor n, flux number p,(v=n/(2pn+1)), Metropolis steps
	void Ex_Metrop(int n_p; int n; int p; int steps);	//Exciton Metropolis
	Cdouble *z;
	int n_p,n,p;
	double RN;
};
void Monca::Build(int N_P, int N, int P)
{
	n_p=N_P;
	n=N;
	p=P;
	RN=sqrt
}

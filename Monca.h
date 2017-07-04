class Monca
{
public:
	Polynomial *LY;			//Landaul level
	Polynomial *JAS;		//Jastrow Factor
	Polynomial *ynm;		//Projected landau level function
public:
	Monca();
	~Monca();
	
	void Symmtrize();
	void G_Metrop(int n_p; int n; int p; int steps);	//Ground state Metropolis. particle number, filling factor n, flux number p,(v=n/(2pn+1)), Metropolis steps
	void Ex_Metrop(int n_p; int n; int p; int steps);	//Exciton Metropolis
	
};
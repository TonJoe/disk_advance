//Feb5

#define NR_END 1
#define LOG2 0.693147180559945
#include <complex>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#define TINY    1.5e-300
#define TINYL    1.5e-1000L

#define FREE_ARG char*

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run_time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}


void cpxdbl_ludcmp0(complex<double> **a,int n,double *d)
{
	double veevee[n+1];
	int i,imax,j,k;
	double big,dum,temp;
	complex<double> sum,dum2;

  *d=1.0;
  for(i=0;i<n;i++)
  {
    big=0.0;
    for(j=0;j<n;j++)
      if( (temp=cabs(a[i][j]))>big)       big=temp;
    if(big==0.0)      	{
				nrerror("Singular matrix in routine cpx_ludcmp\n");
			}
    veevee[i]=1.0/big;      /* Save the scaling */
  }
  for(j=0;j<n;j++){
    for(i=0;i<j;i++){
      sum= a[i][j];
      for(k=0;k<i;k++)     sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<n;i++){
      sum=a[i][j];
      for(k=0;k<j;k++)
        sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;

      /* Is the figure of matrix for the pivot better than the best so far? */
      if( (dum=veevee[i]*cabs(sum)) >= big){
        big=dum;
        imax=i;
      }
    }

    /*Do we need to interchange rows? Yes, do so.. */
    if(j!=imax){
      for(k=0;k<n;k++){
        dum2=a[imax][k];
        a[imax][k]=a[j][k];
  a[imax][k]=a[j][k];
        a[j][k]=dum2;
      }
      *d=-(*d);
      veevee[imax]=veevee[j];
    }
    if(cabs(a[j][j])==0.0)     a[j][j]=TINY+_Complex_I*TINY;
    if(j!=n){
      dum2=1.0/a[j][j];
      for(i=j+1;i<n;i++)      a[i][j]=a[i][j]*dum2;
    }
  }
}


void cpxdbl_ludcmp1(complex<double> **a,int n,double *d)
{
	double veevee[n+1];
	int i,imax,j,k;
	double big,dum,temp;
	complex<double> sum,dum2;

  *d=1.0;
  for(i=1;i<=n;i++)
  {
    big=0.0;
    for(j=1;j<=n;j++)
      if( (temp=cabs(a[i][j]))>big)       big=temp;
    if(big==0.0)      	{
				nrerror("Singular matrix in routine cpx_ludcmp\n");
			}
    veevee[i]=1.0/big;      /* Save the scaling */
  }
  for(j=1;j<=n;j++){
    for(i=1;i<j;i++){
      sum= a[i][j];
      for(k=1;k<i;k++)     sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<=n;i++){
      sum=a[i][j];
      for(k=1;k<j;k++)
        sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;

      /* Is the figure of matrix for the pivot better than the best so far? */
      if( (dum=veevee[i]*cabs(sum)) >= big){
        big=dum;
        imax=i;
      }
    }

    /*Do we need to interchange rows? Yes, do so.. */
    if(j!=imax){
      for(k=1;k<=n;k++){
        dum2=a[imax][k];
        a[imax][k]=a[j][k];
		 a[imax][k]=a[j][k];
        a[j][k]=dum2;
      }
      *d=-(*d);
      veevee[imax]=veevee[j];
    }

	//indx[j]=imax;

    if(cabs(a[j][j])==0.0)     a[j][j]=TINY+_Complex_I*TINY;
    if(j!=n){
      dum2=1.0/a[j][j];
      for(i=j+1;i<=n;i++)      a[i][j]=a[i][j]*dum2;
    }
  }
}

complex<double> cpxdbl_det0(complex<double> **a,int n)
{
  int i;
  double d;
  complex<double> det;

  cpxdbl_ludcmp0(a,n,&d);
  det=d;
  for(i=0;i<n;i++){
    det=det*a[i][i];
  }
return det;
}

complex<double> cpxdbl_det1(complex<double> **a,int n)
{
  int i;
  double d;
  complex<double> det;

  cpxdbl_ludcmp1(a,n,&d);
  det=d;
  for(i=1;i<=n;i++){
    det=det*a[i][i];
  }
return det;
}


double  *dvector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if(!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;
}

 double  *ldvector(int nl,int nh)
/* allocate a  double vector with subscript range v[nl..nh] */
{
   double *v;
  v=( double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof( double)));
  if(!v) nrerror("allocation failure in ldvector()");
  return v-nl+NR_END;
}

int  *ivector(int nl,int nh)
/* allocate a int vector with subscript range v[nl..nh] */
{
  int *v;
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if(!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

  *lvector(int nl,int nh)
/* allocate a  vector with subscript range v[nl..nh] */
{
   *v;
  v=( *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof()));
  if(!v) nrerror("allocation failure in lvector()");
  return v-nl+NR_END;
}

double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /*allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl])    nrerror("allocation failure 2 in dmatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}


double ***dtensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
    int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    double ***t;

    t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double **)));
    if(!t)    nrerror("allocation failure 1 in cpx_tensor()");
    t+=NR_END;
    t-=nrl;

    t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double *)));
    if(!t[nrl])   nrerror("allocation failure 2 in cpx_tensor()");
    t[nrl]+=NR_END;
    t[nrl]-=ncl;

    t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
    if(!t[nrl][ncl])    nrerror("allocation failure 3 in the cpx_tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1;j<=nch;j++)
        t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++)
    {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++)
            t[i][j]=t[i][j-1]+ndep;
    }
    
    return t;
}

void free_dtensor(double ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}


 double **ldmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

  /*allocate pointers to rows */
  m=( double **) malloc((size_t)((nrow+NR_END)*sizeof( double*)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=( double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof( double)));
  if (!m[nrl])    nrerror("allocation failure 2 in dmatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}

short unsigned int **suimatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  short unsigned int **m;

  /*allocate pointers to rows */
  m=(  short unsigned int **) malloc((size_t)((nrow+NR_END)*sizeof(short unsigned int*)));
  if (!m) nrerror("allocation failure 1 in suimatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(short unsigned int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(short unsigned int)));
  if (!m[nrl])    nrerror("allocation failure 2 in suimatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;
  return m;
}

uint8_t **ui8matrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  uint8_t **m;

  /*allocate pointers to rows */
  m=(  uint8_t **) malloc((size_t)((nrow+NR_END)*sizeof(uint8_t*)));
  if (!m) nrerror("allocation failure 1 in suimatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(uint8_t *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(uint8_t)));
  if (!m[nrl])    nrerror("allocation failure 2 in suimatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;
  return m;
}

int8_t **i8matrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int8_t **m;

  /*allocate pointers to rows */
  m=(  int8_t **) malloc((size_t)((nrow+NR_END)*sizeof(int8_t*)));
  if (!m) nrerror("allocation failure 1 in suimatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(int8_t *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int8_t)));
  if (!m[nrl])    nrerror("allocation failure 2 in suimatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;
  return m;
}


complex<double> **cmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  complex<double> **m;

  /*allocate pointers to rows */
  m=(complex<double> **) malloc((size_t)((nrow+NR_END)*sizeof(complex<double> *)));
  if (!m) nrerror("allocation failure 1 in cmatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(complex<double> *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(complex<double>)));
  if (!m[nrl])    nrerror("allocation failure 2 in dmatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}

complex<double> ***ctensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
    int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    complex<double> ***t;

    t=(complex<double> ***) malloc((size_t)((nrow+NR_END)*sizeof(complex<double>**)));
    if(!t)    nrerror("allocation failure 1 in cpx_tensor()");
    t+=NR_END;
    t-=nrl;

    t[nrl]=(complex<double> **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(complex<double>*)));
    if(!t[nrl])   nrerror("allocation failure 2 in cpx_tensor()");
    t[nrl]+=NR_END;
    t[nrl]-=ncl;

    t[nrl][ncl]=(complex<double> *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(complex<double>)));
    if(!t[nrl][ncl])    nrerror("allocation failure 3 in the cpx_tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1;j<=nch;j++)
        t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++)
    {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++)
            t[i][j]=t[i][j-1]+ndep;
    }
    
    return t;
}

void free_ctensor(complex<double> ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}


 complex<double> **lcmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
   complex<double> **m;

  /*allocate pointers to rows */
  m=( complex<double> **) malloc((size_t)((nrow+NR_END)*sizeof( complex<double> *)));
  if (!m) nrerror("allocation failure 1 in cmatrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=( complex<double> *)malloc((size_t)((nrow*ncol+NR_END)*sizeof( complex<double>)));
  if (!m[nrl])    nrerror("allocation failure 2 in dmatrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}

int **imatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /*allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl])    nrerror("allocation failure 2 in matrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}

short int **simatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  short int **m;

  /*allocate pointers to rows */
  m=(short int **) malloc((size_t)((nrow+NR_END)*sizeof(short int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(short int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(short int)));
  if (!m[nrl])    nrerror("allocation failure 2 in matrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}

complex<double> *cvector(int nl,int nh)
{
  complex<double> *v;
  v=(complex<double>*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(complex<double>)));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}

 complex<double> *lcvector(int nl,int nh)
{
   complex<double> *v;
  v=( complex<double>*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof( complex<double>)));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}


short int *sivector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	short int *v;
	v=(short int*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(short int)));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}


bool *boolvector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	bool *v;
	v=(bool*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(bool)));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}


uint8_t *ui8vector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	uint8_t *v;
	v=(uint8_t*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(uint8_t)));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}

int8_t *i8vector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	int8_t *v;
	v=(int8_t*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int8_t)));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}


float **fmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /*allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m +=NR_END;
  m -= nrl;
  
  /*allocate rows and set pointers to them */
  m[nrl]=(float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl])    nrerror("allocation failure 2 in matrix()");
  m[nrl] +=NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)     m[i]=m[i-1]+ncol;

  /*return pointer to array of pointers to row */
  return m;
}



float *fvector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  float *v;
  v=(float*)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float )));
  if(!v) nrerror("allocation failure in cvector()");
  return v-nl+NR_END;
}



void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_ldmatrix( double **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_cmatrix(complex<double> **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_lcmatrix( complex<double> **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


void free_simatrix(short int **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}


void free_suimatrix(short unsigned int **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_fmatrix(float **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_dvector(double *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ldvector( double *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(complex<double> *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_lcvector( complex<double> *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_ivector(int *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_sivector(short int *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void free_suivector(short unsigned int *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_fvector(float *v,int nl,int nh)
/* free a  double vector allocated with ldvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void cpxl_ludcmp( complex<double> **a,int n,int *indx, double *d)
{
	 double veevee[n+1];
	int i,imax,j,k;
	 double big,dum,temp;
	 complex<double> sum,dum2;

  *d=1.0L;
  for(i=1;i<=n;i++)
  {
    big=0.0L;
    for(j=1;j<=n;j++)
		if( (temp=cabsl(a[i][j]))>big)       big=temp;
    if(big==0.0L)      	{
				nrerror("Singular matrix in routine cpx_ludcmp\n");
			}
    veevee[i]=1.0L/big;      /* Save the scaling */
  }
  for(j=1;j<=n;j++){
    for(i=1;i<j;i++){
      sum= a[i][j];
      for(k=1;k<i;k++)     sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0L;
    for(i=j;i<=n;i++){
      sum=a[i][j];
      for(k=1;k<j;k++)
        sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;

      /* Is the figure of matrix for the pivot better than the best so far? */
      if( (dum=veevee[i]*cabsl(sum)) >= big){
        big=dum;
        imax=i;
      }
    }

    /*Do we need to interchange rows? Yes, do so.. */
    if(j!=imax){
      for(k=1;k<=n;k++){
        dum2=a[imax][k];
        a[imax][k]=a[j][k];
		 a[imax][k]=a[j][k];
        a[j][k]=dum2;
      }
      *d=-(*d);
      veevee[imax]=veevee[j];
    }

	indx[j]=imax;

    if(cabsl(a[j][j])==0.0L)     a[j][j]=TINYL+_Complex_I*TINYL;
    if(j!=n){
      dum2=1.0/a[j][j];
      for(i=j+1;i<=n;i++)      a[i][j]=a[i][j]*dum2;
    }
  }
}
void cpxl_lubksb( complex<double> **a,int n,int *indx, complex<double> b[])
{
     complex<double> sum;
    int i,ii=0,ip,j;
    
    for(i=1;i<=n;i++)
    {
        ip=indx[i];sum=b[ip];b[ip]=b[i];
        if(ii)
        {
           for(j=ii;j<=i-1;j++) sum-=a[i][j]*b[j];
        }
        else if(cabsl(sum)) ii=i;
        b[i]=sum;
    }

    for(i=n;i>=1;i--)
    {
        sum=b[i];
        for(j=i+1;j<=n;j++) sum-=a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}
void inv_cpxl_matrix( complex<double> **a, complex<double> **inv_a,int n, complex<double> *determinant)
{
     complex<double>  *col,det;
     double d;
    int i,j,*indx;

    indx=ivector(1,n);col=lcvector(1,n);

    cpxl_ludcmp(a,n,indx,&d);
	  
	det=d; for(i=1;i<=n;i++) det=det*a[i][i]; *determinant=det;

    for(j=1;j<=n;j++)
    {
        for(i=1;i<=n;i++) col[i]=0.0L;
        
        col[j]=1.0L;
        cpxl_lubksb(a,n,indx,col);
        
        for(i=1;i<=n;i++) inv_a[i][j]=col[i];
    }

    free_ivector(indx,1,n);free_lcvector(col,1,n);
}



void cpx_ludcmp(complex<double> **a,int n,int *indx,double *d)
{
	double veevee[n+1];
	int i,imax,j,k;
	double big,dum,temp;
	complex<double> sum,dum2;

  *d=1.0;
  for(i=1;i<=n;i++)
  {
    big=0.0;
    for(j=1;j<=n;j++)
      if( (temp=cabs(a[i][j]))>big)       big=temp;
    if(big==0.0)      	{
				nrerror("Singular matrix in routine cpx_ludcmp\n");
			}
    veevee[i]=1.0/big;      /* Save the scaling */
  }
  for(j=1;j<=n;j++){
    for(i=1;i<j;i++){
      sum= a[i][j];
      for(k=1;k<i;k++)     sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<=n;i++){
      sum=a[i][j];
      for(k=1;k<j;k++)
        sum=sum-a[i][k]*a[k][j];
      a[i][j]=sum;

      /* Is the figure of matrix for the pivot better than the best so far? */
      if( (dum=veevee[i]*cabs(sum)) >= big){
        big=dum;
        imax=i;
      }
    }

    /*Do we need to interchange rows? Yes, do so.. */
    if(j!=imax){
      for(k=1;k<=n;k++){
        dum2=a[imax][k];
        a[imax][k]=a[j][k];
		 a[imax][k]=a[j][k];
        a[j][k]=dum2;
      }
      *d=-(*d);
      veevee[imax]=veevee[j];
    }

	indx[j]=imax;

    if(cabs(a[j][j])==0.0)     a[j][j]=TINY+_Complex_I*TINY;
    if(j!=n){
      dum2=1.0/a[j][j];
      for(i=j+1;i<=n;i++)      a[i][j]=a[i][j]*dum2;
    }
  }
}
void cpx_lubksb(complex<double> **a,int n,int *indx,complex<double> b[])
{
    complex<double> sum;
    int i,ii=0,ip,j;
    
    for(i=1;i<=n;i++)
    {
        ip=indx[i];sum=b[ip];b[ip]=b[i];
        if(ii)
        {
           for(j=ii;j<=i-1;j++) sum-=a[i][j]*b[j];
        }
        else if(cabs(sum)) ii=i;
        b[i]=sum;
    }

    for(i=n;i>=1;i--)
    {
        sum=b[i];
        for(j=i+1;j<=n;j++) sum-=a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}
void inv_cpx_matrix(complex<double> **a,complex<double> **inv_a,int n,complex<double> *determinant)
{
    complex<double>  *col,det;
    double d;
    int i,j,*indx;

    indx=ivector(1,n);col=cvector(1,n);

    cpx_ludcmp(a,n,indx,&d);
	  
	det=d; for(i=1;i<=n;i++) det=det*a[i][i]; *determinant=det;

    for(j=1;j<=n;j++)
    {
        for(i=1;i<=n;i++) col[i]=0.0;
        
        col[j]=1.0;
        cpx_lubksb(a,n,indx,col);
        
        for(i=1;i<=n;i++) inv_a[i][j]=col[i];
    }

    free_ivector(indx,1,n);free_cvector(col,1,n);
}






double DDOT(int n,double *dx, double *dy)
{
	double dotprod=0;
	int i, m;
	m = n % 5;
	if ( m != 0 ) {
	for ( i = 1 ; i <= m ; i++ )
	dotprod = dotprod + dx[i] * dy[i];
	if ( n < 5 )
		return dotprod;
	}
	for ( i = m + 1 ; i <= n ; i = i + 5 )
	dotprod = dotprod + dx[i] * dy[i] + dx[i+1] * dy[i+1] + dx[i+2] * dy[i+2] + dx[i+3] * dy[i+3] + dx[i+4] * dy[i+4];
	return dotprod;
}


void DAX(int n, double da, double *dx,double *dy)
{
	int ix, iy, i, m;
	m = n % 4;
	if ( m != 0 ) 
	{
		for ( i = 1 ; i <= m ; i++ )
		{
			dy[i] = da * dx[i];
		}
		if ( n < 4 )
		{
			return;
		}
	}
	for ( i = m + 1 ; i <= n ; i = i + 4 ) 
	{
		dy[i] = da * dx[i];
		dy[i+1] = da * dx[i+1];
		dy[i+2] = da * dx[i+2];
		dy[i+3] = da * dx[i+3];
	}
	return;
}


void DAXPY(int n,double da,double *dx,double *dy)
{
	//dy=dy+sa*dx
	int ix, iy, i, m;
	m = n % 4;
	if ( m != 0 ) 
	{
		for ( i = 1 ; i <= m ; i++ )
		dy[i] = dy[i] + da * dx[i];
		if ( n < 4 )
		{
			return;
		}
	}
	for ( i = m + 1 ; i <= n ; i = i + 4 ) 
	{
		dy[i] = dy[i] + da * dx[i];
		dy[i+1] = dy[i+1] + da * dx[i+1];
		dy[i+2] = dy[i+2] + da * dx[i+2];
		dy[i+3] = dy[i+3] + da * dx[i+3];
	}
	return;
}


void DGEVM(int dima,int dimb, double **A, double *b, double *out)
{
	int i1;
	for(i1=1;i1<=dima;i1++)
	{
		out[i1]=0;
	}
	for(i1=1;i1<=dimb;i1++)
	{
		DAXPY(dima,b[i1],A[i1],out);
	}
}

void DGEMV(int dima,int dimb, double **A, double *b, double *out)
{
	int i1;
	for(i1=1;i1<=dima;i1++)
	{
		out[i1]=DDOT(dimb,A[i1],b);
	}
}

void DspMV( int nnz, int dimout, int **indx, double *A, double *b, double *out)
{
	 int i1;
	for(i1=1;i1<=dimout;i1++)
	{
		out[i1]=0;
	}
	for(i1=1;i1<=nnz;i1++)
	{
		out[indx[i1][1]]+=A[i1]*b[indx[i1][2]];
	}
}


void DGEVM_SELF(int dima,int dimb, double **A, double *b)
{
	int i1;
	DAX(dima,b[1],A[1],A[1]);
	for(i1=2;i1<=dimb;i1++)
	{
		DAXPY(dima,b[i1],A[i1],A[1]);
	}
}



complex<double> CDOT(int n,complex<double> *cx, complex<double> *cy)
{
	complex<double> dotprod=0;
	int i, m;
	m = n % 5;
	if ( m != 0 ) {
	for ( i = 1 ; i <= m ; i++ )
	dotprod = dotprod + cx[i] * cy[i];
	if ( n < 5 )
		return dotprod;
	}
	for ( i = m + 1 ; i <= n ; i = i + 5 )
	dotprod = dotprod + ~cx[i] * cy[i] + ~cx[i+1] * cy[i+1] + ~cx[i+2] * cy[i+2] + ~cx[i+3] * cy[i+3] + ~cx[i+4] * cy[i+4];
	return dotprod;
}

complex<double> CDOTD(int n,complex<double> *cx, double *dy)
{
	complex<double> dotprod=0;
	int i, m;
	m = n % 5;
	if ( m != 0 ) {
	for ( i = 1 ; i <= m ; i++ )
	dotprod = dotprod + cx[i] * dy[i];
	if ( n < 5 )
		return dotprod;
	}
	for ( i = m + 1 ; i <= n ; i = i + 5 )
	dotprod = dotprod + cx[i] * dy[i] + cx[i+1] * dy[i+1] + cx[i+2] * dy[i+2] + cx[i+3] * dy[i+3] + cx[i+4] * dy[i+4];
	return dotprod;
}


void CAX(int n, complex<double> ca, complex<double> *cx,complex<double> *cy)
{
	int ix, iy, i, m;
	m = n % 4;
	if ( m != 0 ) 
	{
		for ( i = 1 ; i <= m ; i++ )
		{
			cy[i] = ca * cx[i];
		}
		if ( n < 4 )
		{
			return;
		}
	}
	for ( i = m + 1 ; i <= n ; i = i + 4 ) 
	{
		cy[i] = ca * cx[i];
		cy[i+1] = ca * cx[i+1];
		cy[i+2] = ca * cx[i+2];
		cy[i+3] = ca * cx[i+3];
	}
	return;
}


void CAXPY(int n,complex<double> ca,complex<double> *cx,complex<double> *cy)
{
	//dy=dy+sa*dx
	int ix, iy, i, m;
	m = n % 4;
	if ( m != 0 ) 
	{
		for ( i = 1 ; i <= m ; i++ )
		cy[i] = cy[i] + ca * cx[i];
		if ( n < 4 )
		{
			return;
		}
	}
	for ( i = m + 1 ; i <= n ; i = i + 4 ) 
	{
		cy[i] = cy[i] + ca * cx[i];
		cy[i+1] = cy[i+1] + ca * cx[i+1];
		cy[i+2] = cy[i+2] + ca * cx[i+2];
		cy[i+3] = cy[i+3] + ca * cx[i+3];
	}
	return;
}

void CADXPCY(int n,complex<double> ca,double *dx,complex<double> *cy)
{
	//dy=dy+sa*dx
	int ix, iy, i, m;
	m = n % 4;
	if ( m != 0 ) 
	{
		for ( i = 1 ; i <= m ; i++ )
		cy[i] = cy[i] + ca * dx[i];
		if ( n < 4 )
		{
			return;
		}
	}
	for ( i = m + 1 ; i <= n ; i = i + 4 ) 
	{
		cy[i] = cy[i] + ca * dx[i];
		cy[i+1] = cy[i+1] + ca * dx[i+1];
		cy[i+2] = cy[i+2] + ca * dx[i+2];
		cy[i+3] = cy[i+3] + ca * dx[i+3];
	}
	return;
}

void CGEVM(int cima,int cimb, complex<double> **cA, complex<double> *cx, complex<double> *cy)
{
	int i1;
	for(i1=1;i1<=cima;i1++)
	{
		cy[i1]=0;
	}
	for(i1=1;i1<=cimb;i1++)
	{
		CAXPY(cima,cx[i1],cA[i1],cy);
	}
}

void GEDMCV(int dima,int cimb, double **dA, complex<double> *cx, complex<double> *cy)
{
	int i1;
	for(i1=1;i1<=dima;i1++)
	{
		cy[i1]=CDOTD(cimb,cx,dA[i1]);
	}
}

void CGEMV_SELF(int cima,int cimb, complex<double> **cA, complex<double> *cx)
{
	int i1;
	CAX(cima,cx[1],cA[1],cA[1]);
	for(i1=2;i1<=cimb;i1++)
	{
		CAXPY(cima,cx[i1],cA[i1],cA[1]);
	}
}


//Bisection Search
//assuming increasing order in each row and Lexicographical order vertically
int BISEARCH(int sizA, int siza, short int **A, short int *a)
{
	int i1,i2,Max,current,lo=1,hi=sizA;

	if(sizA==0) return 0;

	for(i2=1;i2<=siza;i2++)
	{
		if(A[1][i2]!=a[i2])
		{
			break;
		}
		
		if(i2==siza)
		{
			return 1;
		}
	}
	
	if(sizA==1) return 0;
	
	for(i2=1;i2<=siza;i2++)
	{
		if(A[sizA][i2]!=a[i2])
		{
			break;
		}
		if(i2==siza)return sizA;
	}
	
	current=sizA/2;
	Max=(int)(1+log(sizA)/LOG2);
	for(i1=1;i1<=Max;i1++)
	{
		for(i2=1;i2<=siza;i2++)
		{
			if(A[current][i2]>a[i2])
			{
				hi=current;
				current=(current+lo)/2;
//				goto NEXTITER;
				break;
			}
			if(A[current][i2]<a[i2])
			{
				lo=current;
				current=(current+hi)/2;
//				goto NEXTITER;
				break;
			}
			if(i2==siza)return current;
		}
//		NEXTITER:
	}
	
	return 0;
}


//Bisection Search
//assuming increasing order in each row and Reverse Lexicographical order vertically
int BISEARCH_REV(int sizA, int siza, short int **A, short int *a)
{
	int i1,i2,Max,current,lo=1,hi=sizA;

	if(sizA==0) return 0;

	for(i2=1;i2<=siza;i2++)
	{
		if(A[1][i2]!=a[i2])
		{
			break;
		}
		
		if(i2==siza)
		{
			return 1;
		}
	}
	
	if(sizA==1) return 0;
	
	for(i2=1;i2<=siza;i2++)
	{
		if(A[sizA][i2]!=a[i2])
		{
			break;
		}
		if(i2==siza)return sizA;
	}
	
	current=sizA/2;
	Max=(int)(1+log(sizA)/LOG2);
	for(i1=1;i1<=Max;i1++)
	{
		for(i2=1;i2<=siza;i2++)
		{
			if(A[current][i2]<a[i2])
			{
				hi=current;
				current=(current+lo)/2;
//				goto NEXTITER;
				break;
			}
			if(A[current][i2]>a[i2])
			{
				lo=current;
				current=(current+hi)/2;
//				goto NEXTITER;
				break;
			}
			if(i2==siza)return current;
		}
//		NEXTITER:
	}
	
	return 0;
}


 int LBISEARCH_UINT8( int sizA, int siza, uint8_t **A, uint8_t *a)
{
	 int i1,i2,Max,current,lo=1,hi=sizA;
	if(sizA==0) return 0;

	for(i2=1;i2<=siza;i2++)
	{
		if(A[1][i2]!=a[i2])
		{
			break;
		}
		
		if(i2==siza)
		{
			return 1;
		}
	}
	if(sizA==1) return 0;
	for(i2=1;i2<=siza;i2++)
	{
		if(A[sizA][i2]!=a[i2])
		{
			break;
		}
		if(i2==siza)return sizA;
	}
	
	current=sizA/2;
	Max=(int)(1+log(sizA)/LOG2);
	for(i1=1;i1<=Max;i1++)
	{
		for(i2=1;i2<=siza;i2++)
		{
			if(A[current][i2]>a[i2])
			{
				hi=current;
				current=(current+lo)/2;
//				goto NEXTITER;
				break;
			}
			if(A[current][i2]<a[i2])
			{
				lo=current;
				current=(current+hi)/2;
//				goto NEXTITER;
				break;
			}
			if(i2==siza)return current;
		}
//		NEXTITER:
	}
	
	return 0;
}

 int LBISEARCH_UINT8_REV( int sizA, int siza, uint8_t **A, uint8_t *a)
{
	 int i1,i2,Max,current,lo=1,hi=sizA;
	if(sizA==0) return 0;

	for(i2=1;i2<=siza;i2++)
	{
		if(A[1][i2]!=a[i2])
		{
			break;
		}
		
		if(i2==siza)
		{
			return 1;
		}
	}
	if(sizA==1) return 0;
	for(i2=1;i2<=siza;i2++)
	{
		if(A[sizA][i2]!=a[i2])
		{
			break;
		}
		if(i2==siza)return sizA;
	}
	
	current=sizA/2;
	Max=(int)(1+log(sizA)/LOG2);
	for(i1=1;i1<=Max;i1++)
	{
		for(i2=1;i2<=siza;i2++)
		{
			if(A[current][i2]<a[i2])
			{
				hi=current;
				current=(current+lo)/2;
//				goto NEXTITER;
				break;
			}
			if(A[current][i2]>a[i2])
			{
				lo=current;
				current=(current+hi)/2;
//				goto NEXTITER;
				break;
			}
			if(i2==siza)return current;
		}
//		NEXTITER:
	}
	
	return 0;
}


//BUBBLE SORT a short integer vector resulting in increasing order output
//does not return the sign of the permutation
void BUBBLESORT(int n,short int *v)
{
	int i1,nsw,sw;
	
	while (n>0)
	{
		nsw=0;
		for(i1=1;i1<n;i1++)
		{
			if(v[i1+1]<v[i1])
			{
				sw=v[i1+1];
				v[i1+1]=v[i1];
				v[i1]=sw;
				nsw=i1+1;
			}
		}
		n=nsw;
//		printf("%d\n",nsw);
	}
}



//returns decreasing order result;
void BUBBLESORT_UINT8_DEC(int n,uint8_t *v)
{
	short int i1,nsw,sw;
	while (n>0)
	{
		nsw=0;
		for(i1=1;i1<n;i1++)
		{
			if(v[i1+1]>v[i1])
			{
				sw=v[i1+1];
				v[i1+1]=v[i1];
				v[i1]=sw;
				nsw=i1+1;
			}
		}
		n=nsw;
	}
}

//BUBBLE SORT a short integer vector resulting in increasing order output
//outputs sign of the permutation
short int BUBBLESORT_PERM(int n,short int *v)
{
	int i1,nsw,sw,parity=1;
	while (n>0)
	{
		nsw=0;
		for(i1=1;i1<n;i1++)
		{
			if(v[i1+1]<v[i1])
			{
				sw=v[i1+1];
				v[i1+1]=v[i1];
				v[i1]=sw;
				nsw=i1+1;
				parity=parity*-1;
			}
		}
		n=nsw;
//		printf("%d\n",nsw);
	}
	return parity;
}



//fliplr a short int vector
void FLIPLR(int n,short int *a)
{
	int i1,m,sw;
	if(n%2==0)
	{
		m=n/2;
	}
	if(n%2==1)
	{
		m=(n-1)/2;
	}
	
	for(i1=1;i1<=m;i1++)
	{
		sw=a[i1];
		a[i1]=a[n-i1+1];
		a[n-i1+1]=sw;
	}
}


void FLIPLR_UINT8(int n,uint8_t *a)
{
	int i1,m,sw;
	if(n%2==0)
	{
		m=n/2;
	}
	if(n%2==1)
	{
		m=(n-1)/2;
	}
	for(i1=1;i1<=m;i1++)
	{
		sw=a[i1];
		a[i1]=a[n-i1+1];
		a[n-i1+1]=sw;
	}
}

void FLIPUD_UINT8( int dim,int n,uint8_t **a)
{
	unsigned  int i1,m;
	short int i2,sw;
	if(dim%2==0)
	{
		m=dim/2;
	}
	else
	{
		m=(dim-1)/2;
	}
	for(i1=1;i1<=m;i1++)
	{
		for(i2=1;i2<=n;i2++)
		{
			sw=a[i1][i2];
			a[i1][i2]=a[dim-i1+1][i2];
			a[dim-i1+1][i2]=sw;
		}
	}
}



//fliplr a short int vector and flip all the signs
void FLIPLRWITHSIGN(int n,short int *a)
{
	int i1,m,sw;
	if(n%2==0)
	{
		m=n/2;
	}
	if(n%2==1)
	{
		m=(n-1)/2;
		a[m+1]*=-1;
	}
	
	for(i1=1;i1<=m;i1++)
	{
		sw=a[i1];
		a[i1]=a[n-i1+1]*-1;
		a[n-i1+1]=sw*-1;
	}
}

double factorial(short unsigned int a)
{
	int i1;
	double d=1.0;
	for(i1=1;i1<=a;i1++)
	{
		d=d*((double)i1);
	}
	return d;
}

void removeMultipleLines(char *S)
{
	char name[500];
	sprintf(name,"sort %s | uniq > bbbbbbbbbb",S);
	system(name);
	sprintf(name," mv bbbbbbbbbb %s",S);
	system(name);	
}

double CLBGRD(float j1,float j2, float m1, float m2, float j, float m)
{
	if ((m!=(m1+m2)) || (j<fabs(j1-j2)) || (j>fabs(j1+j2)) || (j1<0) || (j2<0) || (j2<fabs(m2)) || (j1<fabs(m1)) || (fabs(m)>j))
	{
		return 0;
	}
	int i1,k,Mloop,mloop;
	double d1,d2,signchoice;
	if( ( ( ( j1+j2-(int)(j1+j2+0.1) ) - ( j-(int)(j+0.1) ) ) != 0  ) || ( (j1+m1-(int)(j1+m1+0.1))!=0 ) || ( (j2+m2-(int)(j2+m2+0.1))!=0 ))
	{
//	printf("%f %f\n",( j1+j2-(int)(j1+j2+0.1) ),( j-(int)(j+0.1) ));
		return 0;
	}
	double *fact;
	fact=dvector(0,(int)(2+j+j1+j2));
	fact[0]=1.0;
	for(i1=1;i1<=(2+j+j1+j2);i1++)
	{
		fact[i1]=(double)i1*fact[i1-1];
	}
	
	signchoice=1.0;
	if(m<0) 
	{
		signchoice=pow(-1.0,j-j1-j2);
		m1=-m1;
		m2=-m2;
	}
	if(j1<j2)
	{
		d2=j1; j1=j2; j2=d2;
		d2=m1; m1=m2; m2=d2;
		signchoice=signchoice*pow(-1.0,j1+j2-j);
	}
	
	
	d1=1.0;
	//evaluate (j-j2+j1)!/(1+j1+j2+j)!
	for(i1=(int)(j-j2+j1+1.0);i1<=(int)(1.0+j1+j2+j);i1++)
	{
		d1=d1*(double)i1;
//		printf("%f\n",d1);
	}
	d1=(2.0*j+1.0)*fact[(int)(j-j1+j2)]*fact[(int)(j1+j2-j)]/d1;
	
	d1=sqrt(d1);
	
	d1=d1*sqrt(fact[(int)(j+m)]*fact[(int)(j-m)]*fact[(int)(j1-m1)]*fact[(int)(j2-m2)]*fact[(int)(j1+m1)]*fact[(int)(j2+m2)]);
	
	Mloop=j1+j2-j;
	if(Mloop>(j1-m1)) Mloop=(j1-m1);
	if(Mloop>(j2+m2)) Mloop=(j2+m2);
	
	mloop=0;
	if(mloop<(-j+j1+m2)) mloop=-((j-j1-m2));
	if(mloop<(-j+j2-m1)) mloop=-((j-j2+m2));
	
	d2=0;
	for(k=mloop;k<=Mloop;k++)
	{
		d2=d2+pow(-1.0,k)/fact[k]/fact[(int)(j1+j2-j-k)]/fact[(int)(j1-m1-k)]/fact[(int)(j2+m2-k)]/fact[(int)(j-j2+m1+k)]/fact[(int)(j-j1-m2+k)];
//		printf("%f,%f\n",j-j1-m2+k,pow(-1.0,k)/fact[k]/fact[(int)(j1+j2-j-k)]/fact[(int)(j1-m1-k)]/fact[(int)(j2+m2-k)]/fact[(int)(j-j2+m1+k)]/fact[(int)(j-j1-m2+k)]);
	}	
	d1=d1*d2;
	free_dvector(fact,0,(int)(2+j1+j2+j));
	return d1*signchoice;
}

double CLBGRD_FACT(double *fact,float j1,float j2, float m1, float m2, float j, float m)
{
	if ((m!=(m1+m2)) || (j<fabs(j1-j2)) || (j>fabs(j1+j2)) || (j1<0) || (j2<0) || (j2<fabs(m2)) || (j1<fabs(m1)) || (fabs(m)>j))
	{
		return 0;
	}
	int i1,k,Mloop,mloop,i2;
	 double d1,d2,signchoice,a[7],b[7],d3;
	if( ( ( ( j1+j2-(int)(j1+j2+0.1) ) - ( j-(int)(j+0.1) ) ) != 0  ) || ( (j1+m1-(int)(j1+m1+0.1))!=0 ) || ( (j2+m2-(int)(j2+m2+0.1))!=0 ))
	{
//	printf("%f %f\n",( j1+j2-(int)(j1+j2+0.1) ),( j-(int)(j+0.1) ));
		return 0;
	}

	signchoice=1.0;
	if(m<0) 
	{
		signchoice=pow(-1.0,j-j1-j2);
		m1=-m1;
		m2=-m2;
	}
	if(j1<j2)
	{
		d2=j1; j1=j2; j2=d2;
		d2=m1; m1=m2; m2=d2;
		signchoice=signchoice*pow(-1.0,j1+j2-j);
	}
	
	Mloop=j1+j2-j;
	if(Mloop>(j1-m1)) Mloop=(j1-m1);
	if(Mloop>(j2+m2)) Mloop=(j2+m2);
	
	mloop=0;
	if(mloop<(-j+j1+m2)) mloop=-((j-j1-m2));
	if(mloop<(-j+j2-m1)) mloop=-((j-j2+m2));
	
	
	a[1]=sqrt(fact[(int)(j+m)]);	
	a[2]=sqrt(fact[(int)(j-m)]);
	a[3]=sqrt(fact[(int)(j1-m1)]);
	a[4]=sqrt(fact[(int)(j2-m2)]);
	a[5]=sqrt(fact[(int)(j1+m1)]);
	a[6]=sqrt(fact[(int)(j2+m2)]);
	
	for(i1=1;i1<=6;i1++)
	{
		for(i2=(i1);i2<=5;i2++)
		{
			if(a[i2]<a[i2+1])
			{
				d3=a[i2]; a[i2]=a[i2+1]; a[i2+1]=d3;
			}
		}
	}
	d2=0;
	for(k=mloop;k<=Mloop;k++)
	{
	 
		b[1]=fact[k]; b[2]=fact[(int)(j1+j2-j-k)]; b[3]=fact[(int)(j1-m1-k)]; b[4]=fact[(int)(j2+m2-k)]; b[5]=fact[(int)(j-j2+m1+k)]; b[6]=fact[(int)(j-j1-m2+k)];
		for(i1=1;i1<=6;i1++)
		{
			for(i2=(i1);i2<=5;i2++)
			{
			if(b[i2]<b[i2+1])
				{
					d3=b[i2]; b[i2]=b[i2+1]; b[i2+1]=d3;
				}
			}
		}
		
		d2=d2+pow(-1.0,k)*(a[1]/b[1])*(a[2]/b[2])*(a[3]/b[3])*(a[4]/b[4])*(a[5]/b[5])*(a[6]/b[6]);
	}
	d1=d2;	
	for(i1=(int)(j-j2+j1+1.0);i1<=(int)(1.0+j1+j2+j);i1++)
	{
		d1=d1/sqrt((double)i1);
	}
	d1=d1*sqrt((2.0*j+1.0)*fact[(int)(j-j1+j2)]*fact[(int)(j1+j2-j)]);	
	return d1*signchoice;
}

double CLBGRD_LFACT(double *lfact,float j1,float j2, float m1, float m2, float j, float m)
{
	if ((m!=(m1+m2)) || (j<fabs(j1-j2)) || (j>fabs(j1+j2)) || (j1<0) || (j2<0) || (j2<fabs(m2)) || (j1<fabs(m1)) || (fabs(m)>j))
	{
		return 0;
	}
	int i1,k,Mloop,mloop,i2;
	 double d1,d2,signchoice,a[7],b[7],d3;
	if( ( ( ( j1+j2-(int)(j1+j2+0.1) ) - ( j-(int)(j+0.1) ) ) != 0  ) || ( (j1+m1-(int)(j1+m1+0.1))!=0 ) || ( (j2+m2-(int)(j2+m2+0.1))!=0 ))
	{
//	printf("%f %f\n",( j1+j2-(int)(j1+j2+0.1) ),( j-(int)(j+0.1) ));
		return 0;
	}

	signchoice=1.0;
	if(m<0) 
	{
		signchoice=pow(-1.0,j-j1-j2);
		m1=-m1;
		m2=-m2;
	}
	if(j1<j2)
	{
		d2=j1; j1=j2; j2=d2;
		d2=m1; m1=m2; m2=d2;
		signchoice=signchoice*pow(-1.0,j1+j2-j);
	}
	
	Mloop=j1+j2-j;
	if(Mloop>(j1-m1)) Mloop=(j1-m1);
	if(Mloop>(j2+m2)) Mloop=(j2+m2);
	
	mloop=0;
	if(mloop<(-j+j1+m2)) mloop=-((j-j1-m2));
	if(mloop<(-j+j2-m1)) mloop=-((j-j2+m2));
	
	
	a[1]=0.5*lfact[(int)(j+m)];	
	a[2]=0.5*lfact[(int)(j-m)];
	a[3]=0.5*lfact[(int)(j1-m1)];
	a[4]=0.5*lfact[(int)(j2-m2)];
	a[5]=0.5*lfact[(int)(j1+m1)];
	a[6]=0.5*lfact[(int)(j2+m2)];

	for(i1=1;i1<=6;i1++)
	{
		for(i2=(i1);i2<=5;i2++)
		{
			if(a[i2]<a[i2+1])
			{
				d3=a[i2]; a[i2]=a[i2+1]; a[i2+1]=d3;
			}
		}
	}
	d2=0;
	for(k=mloop;k<=Mloop;k++)
	{
	 
		b[1]=lfact[k]; b[2]=lfact[(int)(j1+j2-j-k)]; b[3]=lfact[(int)(j1-m1-k)]; b[4]=lfact[(int)(j2+m2-k)]; b[5]=lfact[(int)(j-j2+m1+k)]; b[6]=lfact[(int)(j-j1-m2+k)];
		for(i1=1;i1<=6;i1++)
		{
			for(i2=(i1);i2<=5;i2++)
			{
			if(b[i2]<b[i2+1])
				{
					d3=b[i2]; b[i2]=b[i2+1]; b[i2+1]=d3;
				}
			}
		}
		d2=d2+pow(-1.0,k)*exp((a[1]-b[1])+(a[2]-b[2])+(a[3]-b[3])+(a[4]-b[4])+(a[5]-b[5])+(a[6]-b[6]));
	}
	d1=d2;	
	for(i1=(int)(j-j2+j1+1.0);i1<=(int)(1.0+j1+j2+j);i1++)
	{
		d1=d1/sqrt((double)i1);
	}
	d1=d1*sqrt((2.0*j+1.0))*exp((lfact[(int)(j-j1+j2)]+lfact[(int)(j1+j2-j)])/2.0);
	return d1*signchoice;
}


void prepareFactorial(double *fact,int nnn)
{
	if(nnn>100)
	{
		printf("large factorials are unreliable\n");
	}
	
	int i1;
	fact[0]=1;
	for(i1=1;i1<=nnn;i1++)
	{
		fact[i1]=fact[i1-1]*((double)i1);
	}
}

void preparelogFactorial(double *lfact,int nnn)
{
	if(nnn>400)
	{
		printf("large factorials are unreliable\n");
	}
	
	int i1;
	lfact[0]=0;
	for(i1=1;i1<=nnn;i1++)
	{
//		lfact[i1]=lfact[i1-1]+log(((double)i1));
		lfact[i1]=lgamma(i1+1);
	}
}


//cpseudop should be a vector with memory available from 0 to MAXm
void CoulombPseudoPotential(int LL, int Q2, int MAXm, double *cpseudop)
{
	int m;
	float l,j,Q,m2,m1;
	Q=(float)Q2/2.0;
	l=Q+LL;
	double v,d;
	double fact[100];
	
	prepareFactorial(fact,99);
	
	for(m=0;m<=MAXm;m++)
	{
		v=0;
		for(m1=-l;m1<=l;m1++)
		{
			for(m2=-l;m2<=l;m2++)
			{
				for(j=fabs(m1-m2);j<=2.0*l;j++)
				{
					d=CLBGRD_FACT(fact,l,j,m1,m2-m1,l,m2)*CLBGRD_FACT(fact,l,j,Q,0,l,Q);
					d=pow(-1.0,(j+m2-m1))*CLBGRD_FACT(fact,l,l,m1,-m1,2.0*l-m,0)*CLBGRD_FACT(fact,l,l,m2,-m2,2.0*l-m,0)*d*d;
					v=v+d;
				}
			}
		}	
		cpseudop[m]=v/sqrt(Q);
	}	
}


void SORTLEXICOGRAPHIC_UINT8( int dim, int n,uint8_t **A)
{
	 int i1,i2,nsw,sw;
	bool flag;
	while (dim>0)
	{
		nsw=0;
		for(i1=1;i1<dim;i1++)
		{
			flag=0;
			for(i2=1;i2<=n;i2++)
			{
				if(A[i1+1][i2]<A[i1][i2])
				{
					flag=1;
					break;
				}
				if(A[i1+1][i2]>A[i1][i2])
				{
					flag=0;
					break;
				}
			}
			if(flag==1)
			{
				for(i2=1;i2<=n;i2++)
				{
					sw=A[i1+1][i2];
					A[i1+1][i2]=A[i1][i2];
					A[i1][i2]=sw;
				}
				nsw=i1+1;
			}
		}
		dim=nsw;
//		printf("%d\n",nsw);
	}
}

void SORTLEXICOGRAPHIC_UINT8_SINGLE( int dim, int n,uint8_t **A)
{
	 int i1,i2,nsw,sw;
	bool flag;

	for(i1=dim-1;i1>=1;i1--)
	{
		flag=0;
		for(i2=1;i2<=n;i2++)
		{
			if(A[i1+1][i2]<A[i1][i2])
			{
				flag=1;
				break;
			}
			if(A[i1+1][i2]>A[i1][i2])
			{
				flag=0;
				break;
			}
		}
		if(flag==1)
		{
			for(i2=1;i2<=n;i2++)
			{
				sw=A[i1+1][i2];
				A[i1+1][i2]=A[i1][i2];
				A[i1][i2]=sw;
			}
			nsw=i1+1;
		}
	}
}

void SORTLEXICOGRAPHIC_UINT8_SINGLE_REV( int dim, int n,uint8_t **A)
{
	 int i1,i2,nsw,sw;
	bool flag;

	for(i1=dim-1;i1>=1;i1--)
	{
		flag=0;
		for(i2=1;i2<=n;i2++)
		{
			if(A[i1+1][i2]>A[i1][i2])
			{
				flag=1;
				break;
			}
			if(A[i1+1][i2]<A[i1][i2])
			{
				flag=0;
				break;
			}
		}
		if(flag==1)
		{
			for(i2=1;i2<=n;i2++)
			{
				sw=A[i1+1][i2];
				A[i1+1][i2]=A[i1][i2];
				A[i1][i2]=sw;
			}
			nsw=i1+1;
		}
	}
}



//written in reverse order: bubble sort from lower end;
void SORTLEXICOGRAPHIC_UINT8_REV( int dim, int n,uint8_t **A)
{
	 int i1,i2,nsw,sw;
	bool flag;
	while (dim>0)
	{
		nsw=0;
		for(i1=1;i1<dim;i1++)
		{
			flag=0;
			for(i2=1;i2<=n;i2++)
			{
				if(A[i1+1][i2]<A[i1][i2])
				{
					flag=1;
					break;
				}
				if(A[i1+1][i2]>A[i1][i2])
				{
					flag=0;
					break;
				}
			}
			if(flag==1)
			{
				for(i2=1;i2<=n;i2++)
				{
					sw=A[i1+1][i2];
					A[i1+1][i2]=A[i1][i2];
					A[i1][i2]=sw;
				}
				nsw=i1+1;
			}
		}
		dim=nsw;
//		printf("%d\n",nsw);
	}
}



void SORTLEXICOGRAPHIC(int dim, int n,short int **A)
{
	int i1,i2,nsw,sw;
	bool flag;
	while (dim>0)
	{
		nsw=0;
		for(i1=1;i1<dim;i1++)
		{
			flag=0;
			for(i2=1;i2<=n;i2++)
			{
				if(A[i1+1][i2]<A[i1][i2])
				{
					flag=1;
					break;
				}
				if(A[i1+1][i2]>A[i1][i2])
				{
					flag=0;
					break;
				}
			}
			if(flag==1)
			{
				for(i2=1;i2<=n;i2++)
				{
					sw=A[i1+1][i2];
					A[i1+1][i2]=A[i1][i2];
					A[i1][i2]=sw;
				}
				nsw=i1+1;
			}
		}
		dim=nsw;
//		printf("%d\n",nsw);
	}
}


//perm returns the permutation that was applied. So applying Anew[perm[i1]] will not give Aold[perm[i1]]
//Aold[perm[i1]]=Anew[i1]
void SORTLEXICOGRAPHIC_PERM(int dim, int n,short int **A,int *perm)
{
	int i1,i2,nsw,sw;
	bool flag;
	
	for(i1=1;i1<=dim;i1++)
	{
		perm[i1]=i1;
	}
	while (dim>0)
	{
		nsw=0;
		for(i1=1;i1<dim;i1++)
		{
			flag=0;
			for(i2=1;i2<=n;i2++)
			{
				if(A[i1+1][i2]<A[i1][i2])
				{
					flag=1;
					break;
				}
				if(A[i1+1][i2]>A[i1][i2])
				{
					flag=0;
					break;
				}
			}
			if(flag==1)
			{
				for(i2=1;i2<=n;i2++)
				{
					sw=A[i1+1][i2];
					A[i1+1][i2]=A[i1][i2];
					A[i1][i2]=sw;
				}
				nsw=i1+1;
				sw=perm[i1+1];
				perm[i1+1]=perm[i1];
				perm[i1]=sw;
			}
		}
		dim=nsw;
//		printf("%d\n",nsw);
	}
}


bool ISLEXICOGRAPHIC(int dim, int n,short int **A)
{
	int i1,i2;
	for(i1=1;i1<dim;i1++)
	{
		for(i2=1;i2<=n;i2++)
		{
			if(A[i1+1][i2]<A[i1][i2])
			{
				return 0;
			}
			if(A[i1+1][i2]>A[i1][i2])
			{
				break;
			}
		}
	}
	return 1;
}

int iRoundd(double d)
{
	if(d<0)
	{
		return ((int)(d-0.5));
	}

	else if(d>0)
	{
		return ((int)(d+0.5));
	}
	
	else
	{
		return (0);
	}

}

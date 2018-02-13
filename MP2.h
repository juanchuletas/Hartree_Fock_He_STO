//----------------------------------------
//Hartree-Fock program
//Using Slater Type Orbitals (STOs)
//Juan José García
//October 2017
//Universidad Autónoma Metropolitana
//-----------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//---------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
static double  get_sum(int a,int b,int c,int d,double vec[],double J[])
{

	int ijkl=0;
	static	double sum;
        for(int i=0;i<nbasis;i++)
        {
                for(int j=0;j<nbasis;j++)
                {
                        //ij = i*nbasis+j;
                        for(int k=0;k<nbasis;k++)
                        {
                                for(int l=0;l<nbasis;l++)
                                {
                                        
					sum = sum + J[ijkl]*vec[i+nbasis*a]*vec[j+nbasis*b]*vec[k+nbasis*c]*vec[l+nbasis*d];
					ijkl++;
					//printf("%d  %d %d %d\n",i+nbasis*a,j+nbasis*b,k+nbasis*c,l+nbasis*d);

                                }
                        }
                }
	}
	return sum;
}

double get_corr_MP2(double J[], double c[], double e[])
{
  double fin_sum = 0.0;
  double fin_sum2 = 0.0,denom;
  for(int a=0;a<N/2;a++) 
    for(int b=0;b<N/2;b++) 
      for(int r=N/2;r<nbasis;r++) 
        for(int s=N/2;s<nbasis;s++) {
	  denom = (e[a]+e[b]-e[r]-e[s]);
	  fin_sum2 += get_sum(a,b,r,s,c,J)*(2.f*get_sum(a,b,r,s,c,J) - get_sum(r,b,s,a,c,J))/denom; 
        }
  	printf("sum %lf\n",fin_sum2);
  return fin_sum2;
}

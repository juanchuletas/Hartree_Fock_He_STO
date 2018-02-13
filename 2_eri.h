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
double rad1J(int qnn[],double alpha[],double j3,int i,int j,int k,int l)
{
	int m,n;
	double a,b,get;

	a = alpha[i]+alpha[j];
	b = alpha[k]+alpha[l];
	m = (qnn[i]+qnn[j]-j3-1.0);
	n = (qnn[k]+qnn[l]+j3);

	get = radial_1(m,n,a,b);

	return get;

}
double rad2J(int qnn[],double alpha[],double j3,int i,int j,int k,int l)
{
	int m,n;
	double a,b,get;

	a = alpha[i]+alpha[j];
	b = alpha[k]+alpha[l];
	n = qnn[k]+qnn[l]-j3-1.0;
	m = qnn[i]+qnn[j]+j3;

	get = radial_2(m,n,a,b);
	return get;

}
double rad1K(int qnn[],double alpha[],double j3,int i,int l,int k,int j)
{
	int m,n;
	double a,b,get;
	
	a = alpha[i]+alpha[l];
	b = alpha[k]+alpha[j];
	m = qnn[i]+qnn[l]-j3-1.0;
	n = qnn[k]+qnn[j]+j3;

	get = radial_1(m,n,a,b);
	return get;

}
double rad2K(int qnn[],double alpha[],double j3,int i,int l,int k,int j)
{
	int m,n;
	double a,b,get;

	a = alpha[i]+alpha[l];
	b = alpha[k]+alpha[j];
	n = qnn[k]+qnn[j]-j3-1.0;
	m = qnn[i]+qnn[l]+j3;

	get = radial_2(m,n,a,b);
	return get;
}
double a_K(int qnn[],double alpha[],double norm[],int i,int l,int k,int j)
{
        double result1,result2,j3=0.0;
                                        
	result1 = (rad1K(qnn,alpha,j3,i,l,k,j) + rad2K(qnn,alpha,j3,i,l,k,j));

                                        
	result2 = norm[i]*norm[j]*norm[k]*norm[l]*16.0*PI*PI*result1;

	return result2;

}

//----------------------------------------------------------------------------
double a_J(int qnn[],double alpha[],double norm[],int i,int j,int k,int l)
{
	double result1,result2,j3=0.0;
					
	result1 = ((rad1J(qnn,alpha,j3,i,j,k,l)+rad2J(qnn,alpha,j3,i,j,k,l)));
	result2= norm[i]*norm[j]*norm[k]*norm[l]*16.0*PI*PI*result1;

	return result2;

}
void get_G(double G[],double P[],double alpha[],int qnn[],double norm[])
{
	int ia,ib,iab,iabcd,icd,ic,id;
	int i,j,k,l,ij,kl,ijkl;
	double suma,suma2;


	iabcd=0;
	suma2 = 0.0;
	for(i=0;i<nbasis;i++)
	{
		for(j=0;j<nbasis;j++)
		{
			ij = j*nbasis+i;
			suma = 0.0;
			//G[ij]=0.0;
			for(k=0;k<nbasis;k++)
			{
				for(l=0;l<nbasis;l++)
				{
					kl = l*nbasis+k;

					//printf("J %lf    K %lf\n",a_J(qnn,alpha,norm,i,j,k,l),a_K(qnn,alpha,norm,i,j,k,l));

					suma = suma + P[kl]*(a_J(qnn,alpha,norm,i,j,k,l) - 0.5*a_K(qnn,alpha,norm,i,l,k,j));

					G[ij] = suma;



				}
			}

		}
	}
}

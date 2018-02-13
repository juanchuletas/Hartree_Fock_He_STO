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

#define nbasis 2
//---------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
double get_corr(double J[],double K[],double hcore[],double c[])
{
	      
        int kl,ij,ijkl,lk,ji;
	int iab,iba,ab;
        double J11,J22,h11,h22,K12,J12;
        J11 = 0.0;
        J22 = 0.0;
	J12 = 0.0;
	K12 = 0.0;
        ijkl = 0;
	int pl;
        for(int i=0;i<nbasis;i++)
        {       
                for(int j=0;j<nbasis;j++)
                {
			pl = i*nbasis+j;
                        for(int k=0;k<nbasis;k++)
                        {       
                                for(int l=0;l<nbasis;l++)
                                {       
					K12 += c[i]*c[i]*c[l]*c[k]*K[ijkl];
	
                                        J11 += c[i]*c[j]*c[k]*c[l]*J[ijkl];

                                        J12 += c[ji]*c[lk]*J[ijkl];

                                        J22 += c[ij]*c[kl]*J[ijkl];


                                        ijkl++;
                                
                                }
                        }
                }
        }

	ab = 0;
	h11 = 0.0;
	h22 = 0.0;
	for(int i=0;i<nbasis;i++)
        {
                for(int j=0;j<nbasis;j++)
                {
                        iab=nbasis+j;
                        iba=nbasis+i;
                        h11 += hcore[ab]*c[i]*c[j];
                        h22 += hcore[ab]*c[iba]*c[iab];
                        ab++;

                }

        }

	
	printf("h11 = %lf\n",h11);
	printf("h22 = %lf\n",h22);
	printf("J11 = %lf\n",J11);
	printf("J12 = %lf\n",J12);
	printf("J22 = %lf\n",J22);
	printf("K12 = %lf\n",K12);


	double e_cero,E_corr,delta,e_2,e_1;

	delta = 0.5*(2.0*h22 - 2.0*h11 + 3.0*J11 + J22);

	E_corr = delta - sqrt(pow(delta,(double)2) + pow(K12,(double)2));

	e_2 = h22 + 2.0*K12 - K12;

	printf("e2 = %lf\n",e_2);
	printf("e1 = %lf\n",h11+J11);


	return E_corr;





}

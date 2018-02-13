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
                        for(int k=0;k<nbasis;k++)
                        {       
                                for(int l=0;l<nbasis;l++)
                                {       
	
                                        J11 += c[i]*c[j]*c[k]*c[l]*J[ijkl];

                                        J12 += c[i]*c[nbasis+j]*c[nbasis+k]*c[l]*J[ijkl];

                                        J22 += c[nbasis+i]*c[nbasis+j]*c[nbasis+k]*c[nbasis+l]*J[ijkl];

                                        K12 += c[i]*c[nbasis+j]*c[k]*c[nbasis+l]*K[ijkl];

					//printf("J= %lf  K= %lf\n",J[ijkl],K[ijkl]);


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

	
	printf("-----------------------------------------------------------\n");
	printf("-----------------DATA FROM CI------------------------------\n");
	printf("h11 = %lf\n",h11);
	printf("h22 = %lf\n",h22);
	printf("J11 = %lf\n",J11);
	printf("J12 = %lf\n",J12);
	printf("J22 = %lf\n",J22);
	printf("K12 = %lf\n",K12);


	double e_cero,E_corr,delta,e2,e1;

	e2 = h22 + 2.0*J12 - K12;
	e1 = h11 + J11;

	delta = 0.5*(2.0*(e2-e1) + J11 + J22 -4.0*J12 + 2.0*K12);

	E_corr = delta - sqrt(pow(delta,(double)2) + pow(K12,(double)2));

	printf("Orb(1) energy by CI  = %lf\n",e1);
	printf("Orb(2) energy by CI  = %lf\n",e2);
        printf("-----------------------------------------------------------\n");


	return E_corr;
	

}

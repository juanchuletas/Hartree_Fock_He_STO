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

#define nbasis 3
#define N 2
#include "factorial.h" /* <------------ Header for compute Factorials--------------*/
#include "Clebsch_Gordan.h"/* <-------Header for Clebsch-Gordan Coeffs------------*/
#include "gamma_function.h"
#include "2_eri.h"
#include "Electron_Repulsion_Integrals.h"
#include "CI_4_index_transformation.h" /* <---- Only for CI and two basis functions --------  */

#include "MP2.h"
#define Z 2.0
#define PI 3.14159
#define ABS(a)     (((a) < 0) ? -(a) : (a))
//------------------------------------------------------------------------------------
//Diagonalization and tranformation module
//-------------------------------------------------------------------------------------
extern void dsyev(char*, char*, int*, double*, int*, double*, double*, int*, int*);

void print_matrix( char* desc, int m, int n, double* a, int lda )
{
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( "  %lf  ", a[i*n+j] );
                printf( "\n" );
        }
}
// INPUT:
// h ---> Matrix to be diagonalized i.e. Fock Matrix
// s ---> Overlap Matrix; It is neccesary to obtain de Transformation Matrix;
// n ---> dimension of Matrix
//
// OUPUT
//
// e ---> eigenvalues
// v ---> eigenvectors

void diag (int n, double *h, double *s, double *e, double *v)
{
  int lwork, i, j, k, nn, ij, info;
  static double small = 1.0e-10;
  static char *V = "V";
  static char *U = "U";
  double *work, *aux, *uno;

  lwork=3*n;
  work = (double *) malloc( lwork * sizeof (double));

  aux = (double *) malloc( n * n * sizeof (double));
  for ( ij = 0; ij < n*n; ij++ )
  {
    aux[ij] = s[ij];
  }

  dsyev( V, U, &n, aux, &n, e, work, &lwork, &info ); 

  if ( info !=0 ) 
  {
    fprintf (stderr, "Overlap matrix diagonalization failed\n");
    exit (1);
  }
// print_matrix( "S matrix ", nbasis, nbasis, aux, nbasis );
 //printf("\n");

 nn = 0;
  for ( i = 0; i < n; i++ ) 
  {
    if ( e[i] > small) 
    {
      for ( j = 0; j < n; j++ ) 
      {
        aux[j+nn*n] = aux[j+i*n] / sqrt(e[i]);
      }
      nn++;
    }
  }
  if ( nn < n ) 
  {
     fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
  }

 //print_matrix( "X matrix ", nbasis, nbasis, aux, nbasis );
 //printf("\n");

  for ( i = 0; i < n; i++) 
  {
    for ( j = 0; j < nn; j++) 
    {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < n; k++) 
      {
        v[i+j*n] = v[i+j*n] + h[i+k*n] * aux[k+j*n];
      }
    }
  }
  uno = (double *) malloc( nn * nn * sizeof (double));
  for ( i = 0; i < nn; i++ ) 
  {
    for ( j = 0; j < nn; j++ ) 
    {
      uno[i+j*nn] = 0.0;
      for  ( k = 0; k < n; k++) 
      {
        uno[i+j*nn] = uno[i+j*nn] + aux[k+i*n] * v[k+j*n];
      }
    }
  }
  //print_matrix( "F\' matrix ", nbasis, nbasis, uno, nbasis );
  //printf("\n");
   /* ---------------- Diagonalization of transformed F or H matrix -------------*/
  dsyev("V", "U", &nn, uno, &nn, e, work, &lwork, &info );

  if ( info !=0 ) 
  {
    fprintf (stderr, "Fock matrix diagonalization failed\n");
    exit (1);
  }


  for ( i = 0; i < n; i++ ) 
  {
    for ( j = 0; j < nn; j++ ) 
    {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < n; k++ ) 
      {
        v[i+j*n] = v[i+j*n] + aux[i+k*n] * uno[k+j*nn];
      }
    }
  }
	//print_matrix( "eigenvalues ", 1, nbasis, e, 1 );
  	//printf("\n");
  free(uno); free(aux); free(work);
  return;
}
//----------------------------------------------------------------------------------
// End of diag & transf module
//----------------------------------------------------------------------------------
void reduced_density_matrix(int n, double *rho, double *d)
{
  int lwork, i, j, k, nn, ij, info;
  static double small = 1.0e-10;
  static char *V = "V";
  static char *U = "U";
  double *work, *aux, *uno;

  lwork=3*n;
  work = (double *) malloc( lwork * sizeof (double));

  aux = (double *) malloc( n * n * sizeof (double));
  for ( ij = 0; ij < n*n; ij++ )
  {
    aux[ij] = rho[ij];
  }


  dsyev( V, U, &n, aux, &n, d, work, &lwork, &info ); 

  if ( info !=0 ) 
  {
    fprintf (stderr, "Density matrix diagonalization failed\n");
    exit (1);
  }
}
void norm_array(double norm[],int qnn[],double alpha[])
{
	int i,ex;
	double al,potAlph;
	for(i=0;i<nbasis;i++)
	{
		al  = 2.0*alpha[i];
		ex  = qnn[i]+qnn[i]+1;
		potAlph = pow(al,(double)ex);
		norm[i] = sqrt(potAlph/(4.0*PI*fact((qnn[i]+qnn[i]))));
	}

}
void overlap(double norm[],double s[],int qnn[],double alpha[])
{
	int i,j,ij,qntot;
	double alphij;
	ij = 0;
	for(i=0;i<nbasis;i++)
	{
		for(j=0;j<nbasis;j++)
		{
			alphij = alpha[i]+alpha[j];
			qntot = qnn[i]+qnn[j];

			s[ij]=4.0*PI*norm[i]*norm[j]*fact(qntot)/((pow(alphij,qntot+1)));

			ij++;

		}
	}

}
void h_Core(double norm[],double pot[],int qnn[],double alpha[],int qnl[],double kin[],double hcore[])
{
	int i,j,ij,qntot;
	double kij1,kij2,kij3;
	double alphij,kinNij;
	ij=0;
	for(i=0;i<nbasis;i++)
	{
		for(j=0;j<nbasis;j++)
		{

			alphij = alpha[i]+alpha[j];
			qntot = qnn[i]+qnn[j];

			kinNij = 2.0*PI*norm[i]*norm[j]; 
			kij1 = (2.0*qnn[j]*alpha[j]*fact(qntot-1.0))/(pow(alphij,qntot));
			kij2  =pow(alpha[j],(double)2)*fact(qntot)/(pow(alphij,qntot+1.0));
			kij3 = (qnl[j]*(qnl[j]+1.0)-qnn[j]*(qnn[j]-1.0))*fact(qntot-2.0)/(pow(alphij,qntot-1.0));
			kin[ij]= kinNij*(kij1-kij2+kij3);

			pot[ij]=(-4.0*Z*PI*norm[i]*norm[j])*fact(qntot-1.0)/((pow(alphij,qntot)));
			hcore[ij] = kin[ij]+pot[ij];

			ij++;

		}
	}
}

void run_SCF(double hcore[],double s[],double norm[],double alpha[],int qnn[])
{

	int i,j,k;
	int ia,ib,ij,ab,a,b,iba;
	int iter,iab,flag;
	double e_new,e_old,energy_cero,E_corr;
	double *c,*F,*G,*e,*v,*P,*rho,*cf,*d,*J,*K;
	double ener_old,delta_E,ener_new;

	F = (double *)malloc(pow(nbasis,2)*sizeof(double));
	J = (double *)malloc(pow(nbasis,4)*sizeof(double));
	K = (double *)malloc(pow(nbasis,4)*sizeof(double));
	P = (double *)malloc(pow(nbasis,2)*sizeof(double));
	G = (double *)malloc(pow(nbasis,2)*sizeof(double));
	c = (double *) malloc ( nbasis* sizeof (double));
	v =  (double *) malloc ( pow(nbasis,2) * sizeof (double) );
	d =  (double *) malloc ( pow(nbasis,2) * sizeof (double) );
	rho =  (double *) malloc ( pow(nbasis,2) * sizeof (double) );
	cf =  (double *) malloc ( pow(nbasis,2) * sizeof (double) );
	e = (double *) malloc (nbasis*sizeof(double));


	for ( ia=0; ia<nbasis; ia++)
        {
		c[ia] = 0.0;
        }
		for(ia=0;ia<nbasis;ia++)
                {
                        for(ib=0;ib<nbasis;ib++)
                        {
                                iab = ib*nbasis+ia;
                                P[iab]= 2.0*c[ia]*c[ib]; /*  First Density Matrix P=0.0 */
                        }
                }


	/*----------------------------------- The SCF starts here -----------------------------------------*/
	printf("-----------------------------------------------------------\n");
	printf("-------------------SCF MODULE------------------------------\n");
	iter=0;  
	delta_E=100000.0;
	ener_old  = 0.0;
	flag=0;
	do
	{
		
		
		if(nbasis<=2)
		{
		get_G_1(G,P,J,K,alpha,qnn,norm); // <---  With this call you are going to get the G(ij), J(ijkl) and K(ijkl) arrays
		}
		else
		{
			get_G(G,P,alpha,qnn,norm);
		}

		//print_matrix( "Density matrix ",nbasis, nbasis, P, nbasis );
  		//printf("\n");
		//print_matrix( "G matrix ",nbasis, nbasis, G, nbasis );
  		//printf("\n");
		for(i=0;i<nbasis;i++)
		{
			for(j=0;j<nbasis;j++)
			{
				ij=j*nbasis+i;

				F[ij] = hcore[ij] + G[ij];

			}
		}
		//print_matrix( "F matrix ", nbasis, nbasis, F, nbasis);
  		//printf("\n");
		diag(nbasis,F,s,e,v); // <---- Here is where you diagonalize the Fock matrix and get the new transformed expansion coeffs
		//print_matrix( "Eigenvalues of F matrix ", 1, nbasis, e, 1 );
  		//printf("\n");
 		for (ia = 0; ia<nbasis; ia++ )
        	{
                	        c[ia] = v[ia]; /* Store Coeffs v[ia] in c[ia]*/
                        	//printf("C=%lf\n",c[ia]);
        	}

		for(ia=0;ia<nbasis;ia++)
                {
                        for(ib=0;ib<nbasis;ib++)
                        {
                                iab = ib*nbasis+ia;

                                P[iab]= 2.0*c[ia]*c[ib]; /*  Density Matrix P */
                                //printf("%lf\n",P[iab]);
                        }
                }
		ener_new = 0.0;
                for(a=0;a<nbasis;a++)
                {
                        for(b=0;b<nbasis;b++)
                        {
				ab=b*nbasis+a;


                                ener_new = ener_new + 0.5*hcore[ab]*P[ab];
                                ener_new = ener_new + 0.5*F[ab]*P[ab];
                        }
                }

 		printf("Energy=%lf  orb=%lf   iteration %d\n",ener_new,e[0],iter);
                delta_E = ABS(ener_old-ener_new);
		ener_old = ener_new;

		iter++;

	}while(delta_E>0.000000001);
	/*----------------------------------- The SCF ends here -----------------------------------------*/

		printf("-----------------------------------------------------------\n");
		printf("--------------Hartree-Fock Final Results-------------------\n");
		print_matrix( "Final Eigenvectors ",nbasis, nbasis, v, nbasis );
		printf("-----------------------------------------------------------\n");
		printf("Hatree-Fock Final Energy = %lf\n",ener_new);
		printf("-----------------------------------------------------------\n");
		printf("Occupied Orbital Energy = %lf\n",e[0]);
		printf("Unoccupied Orbital Energy = %lf\n",e[nbasis-1]);
		printf("-----------------------------------------------------------\n");
		printf("-----------------------------------------------------------\n");
		for (a = 0; a<nbasis; a++ )
        	{
			for(b=0;b<nbasis;b++)
			{
				ab = a*nbasis+b;
                	        cf[ab] = v[ab];
				printf("%lf\n",cf[ab]); //cf means final expansion coeffs
			}
        	}

		printf("-----------------------------------------------------------\n");
		iab = 0;
		for(int i=0;i<nbasis;i++)
                {
                        for(int j=0; j<nbasis;j++)
                        {
                                rho[iab]= cf[i]*cf[j];
				rho[iab] = rho[iab]*cf[nbasis+i]*cf[nbasis+j];
				//printf("%d %d",cf[i],cf[i]);
				iab++;
			     
			
                        }
                }
		print_matrix( "Final Density matrix ",nbasis, nbasis, rho, nbasis );
  		printf("\n");
		printf("-----------------------------------------------------------\n");

		reduced_density_matrix(nbasis,rho,d);

		print_matrix( "First Order Reduced Density Matrix ", 1, nbasis, d, 1 );
		printf("-----------------------------------------------------------\n");

		//-------- At this point you will use de CI_4_index_transformatio header file to compute 
		//-------- the Correlation Energy
		if(nbasis<=2)
		{
		E_corr = get_corr(J,K,hcore,cf); /* <--Be sure that you are sending the correct expansion coeffss */
		printf("E_corr = %lf\n",E_corr);
		printf("E_0  = %lf\n",E_corr+ener_new);
		printf("-----------------------------------------------------------\n");
		}
		else 
		{
			printf("CI only for two basis functions\n");
		}


		printf("Ecorr = %lf\n",get_corr_MP2(J,cf,e));

		free(F),free(J),free(K),free(G),free(c),free(rho),free(d),free(P),free(e);

}





int main ()
{
	double alpha[10] = {2.91093,1.45363,1.45,0.9,0.625,0.4050,4.0};
	int qnn[3]={1,1,2};
	int qnl[3]={0,0,0};
	int qnm[3]={0,0,0};
	double *s,*pot,*kin,*hcore,*norm,*G,*P;

	kin = (double *)malloc (nbasis*nbasis* sizeof(double));
	pot = (double *)malloc (nbasis*nbasis* sizeof(double));
	s = (double *)malloc (nbasis*nbasis* sizeof(double));
	hcore = (double *)malloc (nbasis*nbasis* sizeof(double));

	extern void diag (int , double *, double *, double *, double *);

	norm = (double *)malloc(nbasis* sizeof(double));


	norm_array(norm,qnn,alpha);

	overlap(norm,s,qnn,alpha);

	h_Core(norm,pot,qnn,alpha,qnl,kin,hcore);

	run_SCF(hcore,s,norm,alpha,qnn);



	//free(hcore),free(s),free(pot),free(kin);


	return 0;


}

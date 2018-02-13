//------------------------------
//Clebsch_Gordan coefficients 
//Juan José García Miranda
//Octubre 2017
//Universidad Autónoma Metropolita
//-------------------------------
#ifndef _Clebsch_Gordan_H_
#define _Clebsch_Gordan_H_ 1
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159

static double clebsch_gordan(double j1, double j2,double j3,double m1, double m2, double m3)
{
	int k;
	static double coef;
	double sum,term;

	if(j3<fabs(j1-j2) && (j3>j1+j2) && fabs(m1)>j1 && fabs(m2)>j2 && fabs(m3)>j3)
	{
		coef = 0.0;
	}
	else
	{

		coef = sqrt((j3+j3+1)/fact(j1+j2+j3+1));
		coef = coef*sqrt(fact(j1+j2-j3)*fact(j2+j3-j1)*fact(j3+j1-j2));
		coef = coef*sqrt(fact(j1+m1)*fact(j1-m1)*fact(j2+m2)*fact(j2-m2)*fact(j3+m3)*fact(j3-m3));
		for(k=0;k<99;k++)
		{
			if(j1+j2-j3-k<0) 
			{
				continue;
			}
				if(j3-j1-m2+k<0)
				{
					continue;
				}
					if(j3-j2+m1+k<0)
					{
						continue;
					}
						if(j1-m1-k<0)
						{
							continue;
						}
							if(j2+m2-k<0)
							{
								continue;
							}
			term=fact(j1+j2-j3-k)*fact(j3-j1-m2+k)*fact(j3-j2+m1+k)*fact(j1-m1-k)*fact(j2+m2-k)*fact(k);

			if (k%2==1)
			{
				term = -term;
			}
			sum = sum + 1.0/term;
		}

		coef  = coef*sum;
	}

	return coef;
}
#endif

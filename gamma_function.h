//------------------------------
//Gamma Function
//Juan José García Miranda
//Octubre 2017
//Universidad Autónoma Metropolita
//-------------------------------


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159


double radial_1(int m, int  n, double a, double b)
{
	int i;
	double term,result;
	result = 0.0;
	term  = (fact(n)/pow(b,n+1.0))*(fact(m)/pow(a,m+1.0));

	for(i=0;i<=n;i++)
	{
		result = result + (fact(n)*fact(n+m-i))/((fact(n-i)*pow(b,i+1.0))*pow((a+b),(m+n)-i+1.0));
	}

	return (term-result);
}
double radial_2(int m,int n, double a, double b)
{
	int i;
	double result;
	result = 0.0;
	for(i=0;i<=n;i++)
	{
		result = result + (fact(n)*fact(m+n-i))/(fact(n-i)*pow(b,i+1.0)*pow(a+b,(m+n-i)+1.0));
	}
	return result;
}


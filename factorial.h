#ifndef _factorial_H_
#define _factorial_H_ 1
#include <stdio.h>

double fact(double n)
{
        int i;
        double result = 1.0;

        for (i=1;i<=n;i++)
        {
                result = result*i;
        }
        return result;
}
#endif


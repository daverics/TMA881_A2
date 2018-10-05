#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<pthread.h>
#include<unistd.h>
#include<complex.h>





int main(){

double _Complex a = 2.0 + 2.0*_Complex_I;
double _Complex b = 2.0 + 2.0*_Complex_I;

double _Complex * sum = a*b;

printf("%d", creal(*sum));




  return 0;
}

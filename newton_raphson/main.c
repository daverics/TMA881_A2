#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<pthread.h>
#include<unistd.h>
#include<complex.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
//#include "concatString.h"
#include "concatDiffString.h"


#define PI 3.14159265358979323846
#define NANO 0.000000001

struct cpl_root {
  double complex number;
  int root;
  int num_of_ite;
};

static long int numb_degree;
static long int numb_pixel;
static long int numb_threads;
double complex * roots;
int* block;
static struct cpl_root ** start_space;
char ** colours;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

// BEGINtest
void * ppm_writer(void * arg){
pthread_mutex_lock(&mutex);
  FILE * fp;
  fp = fopen("newton_attractors_xd.ppm","w");

  char n_pix[5];
  sprintf(n_pix, "%d ", numb_pixel);

  char **writeString = malloc(sizeof(char)*4);
  writeString[0] = (char*) malloc(sizeof(char)*4);
  writeString[1] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeString[2] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeString[3] = (char*) malloc(sizeof(char)*6);

  printf("s");

  writeString[0] = "P3\n";
  writeString[1] = n_pix;
  writeString[2] = n_pix;
  writeString[3] = "\n10\n";


  char * stringent = concatDiffString(4,writeString);
  //"P3\n1000 1000\n10\n";
  fwrite(stringent,sizeof(char),strlen(stringent),fp);
  for (size_t i = 0; i < numb_pixel; i++) {
    for (size_t j = 0; j < numb_pixel; j++) {
       fwrite(colours[start_space[i][j].root],sizeof(char),strlen(colours[start_space[i][j].root]),fp);
    }

  }

  fclose(fp);
   pthread_mutex_unlock(&mutex);
}


// SLUTtest




double complex pow_cpx(double complex a, int d){
  double complex a_tmp=a;
  for (size_t i = 1; i < d; i++) {
    a=a*a_tmp;
  }
  return a;
}
double complex newton_raphson_iteration(double complex startpoint,int d){
  double complex numerator1, numerator2,denominator,help,final_numerator, nextpoint;
  numerator1 = pow_cpx(startpoint,d);
  complex double unity = 1 + 0*I;
  numerator1=numerator1-unity;
  help = pow_cpx(startpoint,d-1);
  numerator2=conj(help);
  denominator=numerator2*help;
  final_numerator = numerator1*numerator2;
  __real__ final_numerator = creal(final_numerator)/(creal(denominator)*d);
  __imag__ final_numerator = cimag(final_numerator)/(creal(denominator)*d);
  nextpoint = startpoint - final_numerator;

  return nextpoint;



}

struct cpl_root newton_raphson_method(struct cpl_root start_point,int degree,double complex * roots , int num_ite){
  double complex start = start_point.number;
  double complex unity = 1+0*I;
  double complex test;
  for (size_t i = 0; i < num_ite; i++) {
    start = newton_raphson_iteration(start,degree);
    test = pow_cpx(start,degree)-unity;
      test = test*conj(test);
      double test2 =(double) creal(test);
      double test3 =(double) cimag(test);
      test2 = abs(test2);
      test3 = abs(test3);
      double test4 = cabs(test4);

    if ((test2 > 10000000000) || (test3 >10000000000)){
      printf("FAIL");
      start_point.root = 9;
      start_point.num_of_ite=i;
      break;
    } else if (test4<0.001) {

      for (size_t k = 0; k < degree; k++) {
        if (cabs(start-roots[k]) < 0.001){
          start_point.root = k+1;
          start_point.num_of_ite=i;
          start_point.number = start;
          return start_point;
        }
      }

    }

    }


    start_point.num_of_ite=1000000;
    start_point.root = 9;
    return start_point;


}

void * newton_f(void * arg){
  int * block_n = (int*) arg;
  for (size_t ix = (*block_n); ix < *(block_n+1); ix++) {
    for (size_t jx = 0; jx < numb_pixel; jx++) {
    start_space[ix][jx]= newton_raphson_method(start_space[ix][jx],numb_degree,roots,50000);
    }
  }
  pthread_exit(0);
}




struct cpl_root ** generateStartspace(){
  struct cpl_root * asentries = (struct cpl_root*) malloc(sizeof(struct cpl_root) * numb_pixel*numb_pixel);
  struct cpl_root ** start_space = (struct cpl_root**) malloc(sizeof(struct cpl_root*) * numb_pixel);
  for ( size_t ix = 0, jx = 0; ix < numb_pixel; ++ix, jx+=numb_pixel ){
    start_space[ix] = asentries + jx;
  }
for ( int ix = 0; ix < numb_pixel; ++ix ){
  for ( int jx = 0; jx < numb_pixel; ++jx ){
    double test4 =(double) (4*jx)/(numb_pixel-1);
    double test5 =(double) (4*ix)/(numb_pixel-1);
    __real__ (start_space[ix][jx].number) =(double)-2+test4;
    __imag__ (start_space[ix][jx].number) =(double)-2+test5;
  }
}
return start_space;
}



double complex * generateRoots(int degree){
  double complex * roots = (double complex*) malloc(sizeof(double complex)*degree);
  for (int ix = 0; ix < degree; ix++) {
    __real__ roots[ix] = cos(ix*2*PI/(degree));
    __imag__ roots[ix] = sin(ix*2*PI/(degree));
  }
  return roots;
}





int main(int argc,char * argv[]){
  typedef struct timespec TS;
  TS tstart, tend;
  timespec_get(&tstart, TIME_UTC);









  if (argc < 4) {
    printf("Wrong number of arguments\n");
    return 0;
  }

  if (argv[1][0] == '-' && argv[1][1] == 't') {
    numb_threads = strtol((argv[1]+2),NULL,10);
  }

  if (argv[2][0] == '-' && argv[2][1] == 'l'){
    numb_pixel = strtol((argv[2]+2),NULL,10);
  }

  numb_degree = strtol(argv[3],NULL,10);

  start_space = generateStartspace();
  roots = generateRoots(numb_degree);

  block = (int*) malloc(sizeof(int)*(numb_threads+1));
  for (size_t ix = 0; ix < numb_threads+1; ix++) {
    block[ix] = ix*numb_pixel/numb_threads;
  }

  pthread_t threads[numb_threads];
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  for (size_t ix = 0; ix < numb_threads; ix++) {
    pthread_create(&threads[ix], &attr, newton_f, &block[ix]);
  }

  for (size_t ix = 0; ix < numb_threads; ix++) {
    pthread_join(threads[ix],NULL);
  }





    colours = (char**) malloc(sizeof(char*)*10);
    for (size_t i = 0; i < 10; i++) {
      colours[i] = (char *) malloc(sizeof(char)*6);
    }
    colours[0] = "5 1 1 ";
    colours[1] = "1 5 1 ";
    colours[2] = "1 1 5 ";
    colours[3] = "1 10 5 ";
    colours[4] = "10 5 1 ";
    colours[5] = "1 5 5 ";
    colours[6] = "5 1 5 ";
    colours[7] = "5 5 1 ";
    colours[8] = "5 5 5 ";
    colours[9] = "9 5 1 ";


  pthread_t painter;


  pthread_create(&painter, &attr, ppm_writer, NULL);
  pthread_join(painter,NULL);



  timespec_get(&tend, TIME_UTC);
  double recTime = tend.tv_sec - tstart.tv_sec + NANO*(tend.tv_nsec - tstart.tv_nsec);

    printf("Time: %f \n", recTime);



  return 0;
}

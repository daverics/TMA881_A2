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

// Declare global variables that need to be accessed by all threads

static long int numb_degree;
static long int numb_pixel;
static long int numb_threads;
double complex * roots; // Contains the roots we look for
int* block;

static struct cpl_root ** start_space; // Corresponds to int ** results;
char * item_done; // To check if items are done
char ** colours; // Will contain colour rgb strings
char ** colourIter;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


double complex somePow(double complex a, int d) {
     double complex out = a;
     if (d == 2) {
          out = a*a;
     } else if (d == 3) {
          out = a*a*a;
     } else if (d == 4) {
          double complex b = a*a;
          out = b*b;
     } else if (d == 5) {
          double complex b = a*a;
          out = b*b*a;
     } else if (d == 6) {
          double complex b = a*a*a;
          out = b*b;
     } else if (d == 7) {
          double complex b = a*a;
          double complex c = b*b;
          out = a*b*c;
     } else if (d == 8) {
          double complex b = a*a;
          double complex c = b*b;
          out = c*c;
     } else if (d == 0) {
          out = 1 + 0.0*I;
     }
     return out;
}

void * ppm_writer2(void * arg){
  char n_pix[5];
  sprintf(n_pix, "%d ", numb_pixel);

  char **writeString = malloc(sizeof(char*)*4);
  writeString[0] = (char*) malloc(sizeof(char)*4);
  writeString[1] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeString[2] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeString[3] = (char*) malloc(sizeof(char)*6);

  writeString[0] = "P3\n";
  writeString[1] = n_pix;
  writeString[2] = n_pix;
  writeString[3] = "\n10\n";

  char * stringent = concatDiffString(4,writeString);

  char **writeStringConv = malloc(sizeof(char*)*4);
  writeStringConv[0] = (char*) malloc(sizeof(char)*4);
  writeStringConv[1] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeStringConv[2] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeStringConv[3] = (char*) malloc(sizeof(char)*7);

  writeStringConv[0] = "P3\n";
  writeStringConv[1] = n_pix;
  writeStringConv[2] = n_pix;
  writeStringConv[3] = "\n100\n";

  char * stringentConv = concatDiffString(4,writeStringConv);

  FILE * fp1;
  FILE * fp2;
  fp1 = fopen("newton_attractors_x7.ppm","w");
  fp2 = fopen("newton_convergence_xd.ppm","w");
  fwrite(stringent,sizeof(char),strlen(stringent),fp1);
  fwrite(stringentConv,sizeof(char),strlen(stringentConv),fp2);

  size_t numWrite = 50;

  char ** str1 = (char **) malloc(sizeof(char*)*numWrite);
  char ** str2 = (char **) malloc(sizeof(char*)*numWrite);

  bool a = true;
  for (size_t i = 0; i < numb_pixel; i++) {
    for (size_t j = 0; j < numb_pixel; j +=numWrite) {
      a = true;

      while (a){
        if (start_space[i][j+numWrite - 1].root != 0){
          a = false;
        }
      }
      //printf("aa\n");


      //printf("aa\n");
      for (size_t jx = j; jx < j + numWrite; jx++) {
        str1[jx] = colours[start_space[i][jx].root];
        str2[jx] = colourIter[start_space[i][jx].num_of_ite];
      }

      //printf("aa\n");

      char * out1 = concatDiffString(numWrite,str1);
      //printf("aa\n");

      char * out2 = concatDiffString(numWrite,str2);
      //printf("aa\n");

      fwrite(out1,sizeof(char),strlen(out1),fp1);
      fwrite(out2,sizeof(char),strlen(out2),fp2);
      //printf("aa\n");


      //fwrite(colours[start_space[i][j].root],sizeof(char),strlen(colours[start_space[i][j].root]),fp1);
      //fwrite(colourIter[start_space[i][j].num_of_ite],sizeof(char),strlen(colourIter[start_space[i][j].num_of_ite]),fp2);


    }
  }
  fclose(fp1);
  fclose(fp2);
  pthread_exit(0);
}


  //
  // for (size_t i = 0; i < numb_pixel; i++) {
  //   for (size_t j = 0; j < numb_pixel; j++) {
  //      fwrite(colours[start_space[i][j].root],sizeof(char),strlen(colours[start_space[i][j].root]),fp);
  //   }
  //
  // }






// BEGINtest


// SLUTtest

// The following block contains Davids updated code. Not entirely done, need to double
// check everything that was changed from start_space to asentries

// Calculates one iteration of Newton


double complex newtonUpdate(double complex a, int d) {
    double complex c = somePow(a, d-1);
    //double complex e = c*c;
    // cabs(e) = creal(c)*creal(c) + cimag(c)*cimag(c)
    double norm = 1/(d*(creal(c)*creal(c) + cimag(c)*cimag(c)));
    double complex b = a - (c*a - 1)*conj(c)*norm;
    return b;
}

// Iterates newton until convergence/not

struct cpl_root newton_raphson_method(struct cpl_root start_point, int degree,
     double complex * roots, int num_ite) {
          double complex x_k = start_point.number; // Set initial point
          double abs_fx_k = 1.0;
          size_t j = 0;
          // function value is not close enough to zero, and real/imag part < 10^10
          while (abs_fx_k > 0.001 && cabs(x_k) > 0.001 &&
          abs(creal(x_k)) < 10000000000 && abs(cimag(x_k)) < 10000000000 ) {
               x_k = newtonUpdate(x_k, degree);
               abs_fx_k = cabs(somePow(x_k, degree) - 1);
               j++;
          }

          if (abs_fx_k < 0.001) {
               for (size_t k = 0; k < degree; k++) {
                    if (cabs(x_k - roots[k]) < 0.001) {
                         start_point.root = k+1;
                         if (j < 12) {
                              start_point.num_of_ite = 1;
                         } /*else if (j < 16) {
                              start_point.num_of_ite = j-7;
                         } else if (j < 18) {
                              start_point.num_of_ite = 10;
                         } else if (j < 20) {
                              start_point.num_of_ite = 11;
                         } else if (j < 22) {
                              start_point.num_of_ite = 12;
                         } else if (j < 24) {
                              start_point.num_of_ite = 13;
                         } else if (j < 26) {
                              start_point.num_of_ite = 14;
                         } else if (j < 28) {
                              start_point.num_of_ite = 15;
                         } else if (j < 30) {
                              start_point.num_of_ite = 16;
                         } else if (j < 32) {
                              start_point.num_of_ite = 17;
                         }*/ else {
                              start_point.num_of_ite = 19;
                         }
                         return start_point;
                    }
               }
          } else {
               start_point.root = degree + 1;
               start_point.num_of_ite = 8;
               return start_point;
          }
}

void * newton_f(void * arg){
  int * block_n = (int*) arg;
  for (size_t ix = *(block_n); ix < numb_pixel; ix+=numb_threads) {
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
       start_space[ix][jx].number = -2 + (double) (4*jx)/(numb_pixel-1) + ( - ((double) (4*ix)/(numb_pixel-1)) + 2)*I;
       start_space[ix][jx].root = 0;
  }
}

return start_space;
}



double complex * generateRoots(int degree){
  double complex * roots = (double complex*) malloc(sizeof(double complex)*degree);
  for (int ix = 0; ix < degree; ix++) {
       roots[ix] = cos(ix*2*PI/(degree)) +  sin(ix*2*PI/(degree))*I;
  }

  return roots;
}



int main(int argc,char * argv[]){
  typedef struct timespec TS;
  TS tstart, tend;
  timespec_get(&tstart, TIME_UTC);

  // Define the colours to use in the image //
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

      colourIter = (char**) malloc(sizeof(char*)*20);

      for(size_t i = 0; i < 20; i++) {
           colourIter[i] = (char*) malloc(sizeof(char)*7);
      }



colourIter[0] = "1 1 1  ";
colourIter[1]= "5 5 5  ";
colourIter[2] = "10 10 10  ";
colourIter[3] = "15 15 15  ";
colourIter[4] = "20 20 20  ";
colourIter[5] = "25 25 25  ";
colourIter[6] = "30 30 30  ";
colourIter[7] = "35 35 35  ";
colourIter[8] = "40 40 40  ";
colourIter[9] = "45 45 45  ";
colourIter[10] = "50 50 50  ";
colourIter[11] = "55 55 55  ";
colourIter[12] = "60 60 60  ";
colourIter[13] = "63 63 63  ";
colourIter[14] = "66 66 66  ";
colourIter[15] = "69 69 69  ";
colourIter[16] = "72 72 72  ";
colourIter[17] = "75 75 75  ";
colourIter[18] = "81 81 81  ";
colourIter[19] = "84 84 84  ";

  //------------------------------------------//









  if (argc < 4) {
    printf("Wrong number of arguments\n");
    return 0;
  }

// Extract arguments and store in variables  //
  if (argv[1][0] == '-' && argv[1][1] == 't') {
    numb_threads = strtol((argv[1]+2),NULL,10);
  }

  if (argv[2][0] == '-' && argv[2][1] == 'l'){
    numb_pixel = strtol((argv[2]+2),NULL,10);
  }

  numb_degree = strtol(argv[3],NULL,10);

  // ---------------------------------------//

  start_space = generateStartspace();
  roots = generateRoots(numb_degree);
  block = (int*) malloc(sizeof(int*)*(numb_threads));
  for (size_t ix = 0; ix < numb_threads; ix++) {
    block[ix] = ix;//*numb_pixel/numb_threads;
  }

  pthread_t threads[numb_threads];
  pthread_t painter;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  for (size_t ix = 0; ix < numb_threads; ix++) {
    pthread_create(&threads[ix], &attr, newton_f, &block[ix]);

  }
  pthread_create(&painter, &attr, ppm_writer2, NULL);




  // Here we join the threads


  for (size_t ix = 0; ix < numb_threads; ix++) {
    pthread_join(threads[ix],NULL);
  }

  pthread_join(painter,NULL);

  timespec_get(&tend, TIME_UTC);
  double recTime = tend.tv_sec - tstart.tv_sec + NANO*(tend.tv_nsec - tstart.tv_nsec);
  printf("Time: %f \n", recTime);



  return 0;
}

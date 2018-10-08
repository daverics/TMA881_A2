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
     }
     return out;
}



// BEGINtest
void * ppm_writer(void * arg){
pthread_mutex_lock(&mutex);
  FILE * fp;
  fp = fopen("newton_attractors_x7.ppm","w");


  // Rewrite this, since numb_pixel is string from the beginning
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
  //"P3\n1000 1000\n10\n";
  fwrite(stringent,sizeof(char),strlen(stringent),fp);
  for (size_t i = 0; i < numb_pixel; i++) {
    for (size_t j = 0; j < numb_pixel; j++) {
         //printf("%d\n", start_space[i][j].root);
       fwrite(colours[start_space[i][j].root],sizeof(char),strlen(colours[start_space[i][j].root]),fp);
    }

  }

  fclose(fp);
  fp = fopen("newton_convergence_xd.ppm","w");

  //char n_pix[5];
  //sprintf(n_pix, "%d ", numb_pixel);
/*
  char **writeString = malloc(sizeof(char*)*4);
  writeString[0] = (char*) malloc(sizeof(char)*4);
  writeString[1] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeString[2] = (char*) malloc(sizeof(char)*strlen(n_pix));
  writeString[3] = (char*) malloc(sizeof(char)*6);

  writeString[0] = "P3\n";
  writeString[1] = n_pix;
  writeString[2] = n_pix;
  writeString[3] = "\n100\n";


  char * stringent = concatDiffString(4,writeString);
  //"P3\n1000 1000\n10\n";*/
  writeString[3] = "\n100\n";
  fwrite(stringent,sizeof(char),strlen(stringent),fp);

  for (size_t i = 0; i < numb_pixel; i++) {
    for (size_t j = 0; j < numb_pixel; j++) {
         //printf("%d\n", start_space[i][j].root);
       fwrite(colourIter[start_space[i][j].num_of_ite],sizeof(char),strlen(colourIter[start_space[i][j].num_of_ite]),fp);
    }

  }


  fclose(fp);
   pthread_mutex_unlock(&mutex);
}


// SLUTtest

// The following block contains Davids updated code. Not entirely done, need to double
// check everything that was changed from start_space to asentries

// Calculates one iteration of Newton


double complex newtonUpdate(double complex a, int d) {
     double norm = 1/(d*cabs(somePow(a, d - 1))*cabs(somePow(a, d - 1)));
     double complex b = a - (somePow(a,d) - 1)*somePow(conj(a), d - 1)*norm;
     return b;
}

// Iterates newton until convergence/not

struct cpl_root newton_raphson_method(struct cpl_root start_point, int degree,
     double complex * roots, int num_ite) {
          //printf("insideNewton\n");
          double complex x_k = start_point.number; // Set initial point
          double abs_fx_k = 1.0;
          size_t j = 0;
          // function value is not close enough to zero, and real/imag part < 10^10
          /*printf("%f + %fi\n",creal(x_k),cimag(x_k));
          printf("%f\n", cabs(somePow(x_k, degree) - 1));
          printf("%f\n", cabs(x_k));*/
          while (abs_fx_k > 0.001 && cabs(x_k) > 0.001 &&
          abs(creal(x_k)) < 10000000000 && abs(cimag(x_k)) < 10000000000 ) {
               x_k = newtonUpdate(x_k, degree);
               abs_fx_k = cabs(somePow(x_k, degree) - 1);
               j++;
               //printf("%f + %f i\n",creal(x_k), cimag(x_k) );
          }

          if (abs_fx_k < 0.001) {
               for (size_t k = 0; k < degree; k++) {
                    if (cabs(x_k - roots[k]) < 0.001) {
                         start_point.root = k;
                         if (j < 19) {
                              start_point.num_of_ite = j;
                         } else {
                              start_point.num_of_ite = 19;
                         }
                         //printf("%d\n", start_point.root);
                         return start_point;
                    }
               }
          } else {
               start_point.root = degree + 1;
               start_point.num_of_ite = 19;
               return start_point;
          }
}

void * newton_f(void * arg){
  int * block_n = (int*) arg;
  for (size_t ix = (*block_n); ix < *(block_n+1); ix++) {
    for (size_t jx = 0; jx < numb_pixel; jx++) {
         //printf("Started newton for some block\n");
         // Update this to asentries instead, which is a vector
    start_space[ix][jx]= newton_raphson_method(start_space[ix][jx],numb_degree,roots,50000);
    // asentries[ix + jx] = newton_raphson_method(asentries[ix + jx], numb_degree, roots,50000);
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
       start_space[ix][jx].number = -2 + (double) (4*jx)/(numb_pixel-1) + ( (double) (4*ix)/(numb_pixel-1) - 2)*I;
      /* asentries[ix + jx]
    double test4 =(double) (4*jx)/(numb_pixel-1);
    double test5 =(double) (4*ix)/(numb_pixel-1);
    __real__ (start_space[ix][jx].number) =(double)-2+test4;
    __imag__ (start_space[ix][jx].number) =(double)-2+test5;
    */
  }
}
// This should return asentries instead, change output to struct cpl_root *
//return asentries;
//printf("Generated startspace\n");
return start_space;
}



double complex * generateRoots(int degree){
  double complex * roots = (double complex*) malloc(sizeof(double complex)*degree);
  for (int ix = 0; ix < degree; ix++) {
       roots[ix] = cos(ix*2*PI/(degree)) +  sin(ix*2*PI/(degree))*I;
    //__real__ roots[ix] = cos(ix*2*PI/(degree));
    //__imag__ roots[ix] = sin(ix*2*PI/(degree));
  }
  //printf("Generated roots\n");
  return roots;
}



/*

double complex pow_cpx(double complex a, int d){
  double complex a_tmp=a;
  for (size_t i = 1; i < d; i++) {
    a=a*a_tmp;
  }
  return a;
}

/*

// This function computes one iteration of the newton method

/*
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
*/



/*
  double complex start = start_point.number;
  double complex unity = 1+0*I;
  double complex test;
  // for loop for each iteration, checks convergence etc
  for (size_t i = 0; i < num_ite; i++) {
    start = newton_raphson_iteration(start,degree); // Calculate x_(k+1)
    test = pow_cpx(start,degree)-unity;
    // Not sure what happens here
      test = test*conj(test);
      double test2 =(double) creal(test);
      double test3 =(double) cimag(test);
      test2 = abs(test2);
      test3 = abs(test3);
      double test4 = cabs(test4);

      // Checking if real or imaginary part is large enough to abort
    if ((test2 > 10000000000) || (test3 >10000000000)){
      printf("FAIL");
      start_point.root = 9;
      start_point.num_of_ite=i;
      break;
    } else if (test4<0.001) {

         // If the function value is close enough to zero, check which root is closest
      for (size_t k = 0; k < degree; k++) {
        if (cabs(start-roots[k]) < 0.001){
          start_point.root = k+1; // Set which root it is
          start_point.num_of_ite=i; // Set the number of iterations
          start_point.number = start; // Has no use? Is the current point
          return start_point;
        }
      }
    }
    }
    start_point.num_of_ite=1000000; // If the method doesn't converge in num_of_ite
    start_point.root = 9; // The 'extra' roots
    return start_point;
}

*/







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
  //printf("hej1\n");
// Removed +1. No problem.
  block = (int*) malloc(sizeof(int*)*(numb_threads));
  for (size_t ix = 0; ix < numb_threads + 1; ix++) {
    block[ix] = ix*numb_pixel/numb_threads;
  }
//printf("hej2\n");
  pthread_t threads[numb_threads];
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  for (size_t ix = 0; ix < numb_threads; ix++) {
    pthread_create(&threads[ix], &attr, newton_f, &block[ix]);
    //printf("Created thread %d\n",ix + 1 );
  }

  // Here threads will do stuff //

  // CODE FROM HERE IS NOT TESTED YET
/*
  // Threads do computations
  for ( size_t ix = offset; ix < nmb_work_items; ix += nmb_threads ) {
    int * result = (int*)malloc(sizeof(int)*item_size);
    // compute work item = Do the newton for an item. Then we store it in results
    // In our case, results is startspace
    results[ix] = result;

    pthread_mutex_lock(&item_done_mutex);
    item_done[li] = 1;
    pthread_mutex_unlock(&item_done_mutex);
  }

  // Make the local copy for the write threads

  void *
write_main(
    void * args
    )
{
  char * item_done_loc = (char*)calloc(nmb_items, sizeof(char));


  // Make a local copy of the result and use it to write

  for ( size_t ix = 0; ix < nmb_items; ) {
  pthread_mutex_lock(&item_done_mutex);
  if ( item_done[ix] != 0 )
    memcpy(item_done_loc, item_done, nmb_items*sizeof(char));
  pthread_mutex_unlock(&item_done_mutex);

  if ( item_done_loc[ix] == 0 ) {
    nanosleep(&sleep_timespec, NULL);
    continue;
  }

  for ( ; ix < nmb_items && item_done_loc[ix] != 0; ++ix ) {
    result = results[ix];

    // Here we should first write the roots it converged to, for all newly finished points
    // Then we should repeat it but with the number of iterations

    // write result

    free(result);
  }
}
/// END OF NON-TESTED CODE
*/



  // Here we join the threads

 // printf("hej3\n");
  for (size_t ix = 0; ix < numb_threads; ix++) {
    pthread_join(threads[ix],NULL);
  }







  pthread_t painter;

  //printf("hej4\n");
  pthread_create(&painter, &attr, ppm_writer, NULL);
  pthread_join(painter,NULL);
//printf("hej5\n");


  timespec_get(&tend, TIME_UTC);
  double recTime = tend.tv_sec - tstart.tv_sec + NANO*(tend.tv_nsec - tstart.tv_nsec);

    printf("Time: %f \n", recTime);



  return 0;
}

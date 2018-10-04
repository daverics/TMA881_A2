#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#define TOL 0.001

typedef struct cplx{
  double real;
  double imag;
  int root;
  int iteration;
}cplx;


static long int numb_degree;
static long int numb_pixel;
static long int numb_threads;
int * block;
static cplx ** pixel;
cplx * roots;
char ** colours;


void * ppm_writer(void * arg){
  FILE * fp;
  fp = fopen("newton_attractors_xd.ppm","w");
  fwrite("P3\n10 10\n10",sizeof(char),numb_pixel,fp);
  fclose(fp);
}




cplx cplx_mul(cplx b, cplx c) {
      cplx a;
      a.real = -(b.imag * c.imag) + b.real * c.real;
      a.imag = b.real * c.imag + b.imag * c.real;
      return a;
    }

cplx cplx_pow(cplx a,long int d){
  if (d<1)
    return a;
  d = d-1;
  cplx c = a;
  for (size_t ix = 0; ix < d; ix++) {
    c = cplx_mul(c,a);
  }
  return c;
}

double norm(cplx a){
  double norm = 1.0;
  norm = norm*(a.real*a.real+a.imag*a.imag);
  return norm;
}

double dist(cplx a, cplx b){
  double dist;
  dist = (a.real-b.real)*(a.real-b.real)+(a.imag-b.imag)*(a.imag-b.imag);
  dist = sqrt(dist);
  return dist;
}

cplx cplx_conj(cplx a) {
  a.imag = -a.imag;
  return a;
}


cplx newton_iteration(cplx cord,long int numb_degree){
  cplx a = cord;
  cplx b = cord;
  cplx c;
  cplx d;
  cplx e;



  c = cplx_pow(a,numb_degree);
  b = c;
  c.real = c.real - 1.0;

  d = cplx_pow(a,numb_degree-1);
  double norm1 = numb_degree*norm(d);
  d = cplx_conj(d);

  e = cplx_mul(c,d);
  e.real = (cord.real) - (e.real/norm1);
  e.imag = (cord.imag) - (e.imag/norm1);

  return e;

}

cplx newton(cplx cord, cplx * roots,long int numb_degree){
  int jx = 0;
  while(jx<1000000){
    jx++;
    if (sqrt(norm(cord)) < TOL){
      cord.root = 0;
      cord.iteration = jx;
      //printf("Zero\n");
      return cord;
    }
    if (abs((int) cord.real)>10000000000 ||
        abs((int) cord.imag)>10000000000 ||
        sqrt(norm(cord))    >10000000000 ){
      cord.root = 0;
      cord.iteration = jx;
      //printf("Inf\n");
      return cord;
    }
    for (size_t ix = 0; ix < numb_degree; ix++) {
      if (dist(cord,roots[ix]) < TOL){
        cord.root = ix+1;
        cord.iteration = jx;
        //printf("Yay\n");
        return cord;
      }
    }
    cord = newton_iteration(cord,numb_degree);
  }
  //printf("Doesn't converge\n");
  cord.root = 0;
  return cord;
}

void * newton_f(void * arg){
  int * block_n = (int*) arg;
  for (size_t ix = *block_n; ix < *(block_n+1); ix++) {
    for (size_t jx = 0; jx < numb_pixel; jx++) {
      pixel[ix][jx] = newton(pixel[ix][jx], roots, numb_degree);
    }
  }
  pthread_exit(0);
}



int main(int argc,char * argv[]){
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

  roots = (cplx*) malloc(sizeof(cplx) *numb_degree);
  for (size_t ix = 0; ix < numb_degree; ix++) {
    roots[ix].real = cos(ix*2*M_PI/(numb_degree));
    roots[ix].imag = sin(ix*2*M_PI/(numb_degree));
  }
  printf("%f,%f\n",roots[0].real,roots[0].imag);

  cplx* mem = (cplx*) malloc(sizeof(cplx)*numb_pixel*numb_pixel);
  pixel = (cplx **) malloc(sizeof(cplx*)*numb_pixel);

  for (size_t ix = 0, jx=0; ix < numb_pixel; ix++, jx+=numb_pixel) {
    pixel[ix] = mem+jx;
  }

  double pixel_size = 4.0/numb_pixel;
  for (size_t ix = 0; ix < numb_pixel; ix++) {
    for (size_t jx = 0; jx < numb_pixel; jx++) {
      pixel[ix][jx].real = -2.0 + jx*pixel_size;
      pixel[ix][jx].imag = 2.0 - ix*pixel_size;
    }
  }

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

  int h;
  for (size_t ix = 0; ix < numb_pixel; ix++) {
    for (size_t jx = 0; jx < numb_pixel; jx++) {
      if (pixel[ix][jx].iteration == 0){
        //printf("%f,%f\n", pixel[ix][jx].real,pixel[ix][jx].imag);
        h++;
      }
    }
  }

  printf("woops %d\n",h);

  colours = (char**) malloc(sizeof(char*)*10);
  for (size_t i = 0; i < 10; i++) {
    colours[i] = (char *) malloc(sizeof(char)*6);
  }
  colours[0] = "5 1 1 ";
  colours[1] = "1 5 1 ";
  colours[2] = "1 1 5 ";
  colours[3] = "1 1 5 ";
  colours[4] = "1 5 1 ";
  colours[5] = "1 5 5 ";
  colours[6] = "5 1 5 ";
  colours[7] = "5 5 1 ";
  colours[8] = "5 5 5 ";
  colours[9] = "9 5 1 ";


  pthread_t painter;

  pthread_create(&painter, &attr, ppm_writer, NULL);
  pthread_join(painter,NULL);





  //pthread_mutex_init(&mutex_sum, NULL);
  //
  // for (tx=0, ix=0; tx < n_threads; ++tx, ix+=block_size) {
  //   double ** arg = malloc(2*sizeof(double*));
  //   arg[0] = a+ix; arg[1] = b+ix;
  //   if (ret = pthread_create(threads+tx, NULL, newton, (void*)arg)) {
  //     printf("Error creating thread: %\n", ret);
  //     exit(1);
  // }
  //}





  // for (size_t ix = 0; ix < numb_pixel; ix++) {
  //   for (size_t jx = 0; jx < numb_pixel; jx++) {
  //     newton(pixel[ix][jx],roots,numb_degree);
  //   }
  // }

  return 0;


}

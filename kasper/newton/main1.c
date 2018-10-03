#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#define TOL 0.001

typedef struct cplx{
  double real;
  double imag;
  int root;
  int iteration;
}cplx;

struct cplx cplx_mul(struct cplx b, struct cplx c) {
      struct cplx a;
      a.real = -(b.imag * c.imag) + b.real * c.real;
      a.imag = b.real * c.imag + b.imag * c.real;
      return a;
    }

struct cplx cplx_pow(struct cplx a,long int d){
  if (d<1)
    return a;
  d = d-1;
  struct cplx c = a;
  for (size_t ix = 0; ix < d; ix++) {
    c = cplx_mul(c,a);
  }
  return c;
}

double norm(struct cplx a){
  double norm = 1.0;
  norm = norm*(a.real*a.real+a.imag*a.imag);
  return norm;
}

double dist(struct cplx a, struct cplx b){
  double dist;
  dist = (a.real-b.real)*(a.real-b.real)+(a.imag-b.imag)*(a.imag-b.imag);
  dist = sqrt(dist);
  return dist;
}

struct cplx cplx_conj(struct cplx a) {
  a.imag = -a.imag;
  return a;
}


struct cplx *newton_iteration(struct cplx *cord,long int numb_degree){
  struct cplx *a = (struct cplx *) malloc(sizeof(struct cplx));
  struct cplx *c = (struct cplx *) malloc(sizeof(struct cplx));
  struct cplx *d = (struct cplx *) malloc(sizeof(struct cplx));
  struct cplx* e = (struct cplx *) malloc(sizeof(struct cplx));
  const struct cplx *ptr = cord;


  memcpy(a,ptr,sizeof(struct cplx));
  memcpy(c,ptr,sizeof(struct cplx));
  memcpy(d,ptr,sizeof(struct cplx));


  *c = cplx_pow(*a,numb_degree);
  c->real = c->real - 1.0;

  *d = cplx_pow(*a,numb_degree-1);
  double norm1 = numb_degree*norm(*d);
  *d = cplx_conj(*d);

  *e = cplx_mul(*c,*d);
//  printf("%f,%f\n",e->real,e->imag);
//  printf("%f,%f\n",cord->real,cord->imag);
  e->real = (cord->real) - (e->real/norm1);
  e->imag = (cord->imag) - (e->imag/norm1);

  printf("%f,%f\n",e->real,e->imag);
  return e;

}

struct cplx * newton(struct cplx *cord, struct cplx * roots,long int numb_degree){
  int jx = 0;
  while(jx<100){
    jx++;
    if (sqrt(norm(*cord)) < TOL){
      cord->root = 0;
      cord->iteration = jx;
      printf("Zero\n");
      return cord;
    }
    if (abs((int) cord->real)>10000000000 ||
        abs((int) cord->imag)>10000000000 ||
        sqrt(norm(*cord))    >10000000000 ){
      cord->root = 100;
      cord->iteration = jx;
      printf("Inf\n");
      return cord;
    }
    for (size_t ix = 0; ix < numb_degree; ix++) {
      if (dist(*cord,roots[ix]) < TOL){
        cord->root = ix+1;
        cord->iteration = jx;
        printf("Yay\n");
        return cord;
      }

    }
    cord = newton_iteration(cord,numb_degree);
  }
  printf("Doesn't converge\n");
  exit(-1);
}





int main(int argc,char * argv[]){
  if (argc < 4) {
    printf("Wrong number of arguments\n");
    return 0;
  }

  long int numb_threads;
  if (argv[1][0] == '-' && argv[1][1] == 't') {
    numb_threads = strtol((argv[1]+2),NULL,10);
  }


  long int numb_pixel;
  if (argv[2][0] == '-' && argv[2][1] == 'l'){
    numb_pixel = strtol((argv[2]+2),NULL,10);
  }


  long int numb_degree;
  numb_degree = strtol(argv[3],NULL,10);

  struct cplx * roots = (struct cplx*) malloc(sizeof(struct cplx) *numb_degree);
  for (size_t ix = 0; ix < numb_degree; ix++) {
    roots[ix].real = cos(ix*2*M_PI/(numb_degree));
    roots[ix].imag = sin(ix*2*M_PI/(numb_degree));
  }
  printf("%f,%f\n",roots[0].real,roots[0].imag);

  cplx* mem = (cplx*) malloc(sizeof(cplx)*numb_pixel*numb_pixel);
  cplx ** pixel = (cplx **) malloc(sizeof(cplx*)*numb_pixel);

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


  struct cplx * cord = (struct cplx *) malloc(sizeof(struct cplx));
  cord->real = 2.0;
  cord->imag = 2.0;
  cord = newton(pixel[0]+999,roots,numb_degree);

  printf("real:%f imag:%f\n",cord->real,cord->imag);
  return 0;


}

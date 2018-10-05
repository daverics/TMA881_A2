#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<pthread.h>
#include<unistd.h>
#include<complex.h>

struct cplx {
  double re;
  double im;
  int root;
};

struct cplx mul_cpx( struct cplx b, struct cplx c) {
    struct cplx mul_ans;
    mul_ans.re = b.re*c.re - b.im*c.im;
    mul_ans.im = b.re*c.im + b.im*c.re;
    return mul_ans;
  }

struct cplx pow_cpx(struct cplx a, int d){
  struct cplx pow_ans = a;
  for (size_t i = 1; i < d; i++) {
    pow_ans = mul_cpx(pow_ans,a);
  }
  return pow_ans;
}

struct cplx conj_cpx(struct cplx a){
  a.im = -a.im;
  return a;
}

struct cplx cpx_add(struct cplx a,struct cplx b){
  struct cplx sum;
  sum.re=a.re+b.re;
  sum.im=a.im+b.re;
  return sum;
}

struct cplx cpx_min(struct cplx a,struct cplx b){
  struct cplx dif;
  dif.re=a.re-b.re;
  dif.im=a.im-b.re;
  return dif;
}



struct cplx newton_raphson_iteration(struct cplx startpoint,int d){
  struct cplx new_point, conj_startpoint_pow, pow_startpoint, help,help1,help2,help3,help4;
  conj_startpoint_pow = pow_cpx(conj_cpx(startpoint),d-1);
  pow_startpoint=pow_cpx(startpoint,d);
  pow_startpoint.re=pow_startpoint.re-1;
  help.re = d; help.im=0;
  help2 = mul_cpx(help,pow_startpoint);
  help3 = mul_cpx(help2,conj_startpoint_pow);
  help4 = cpx_min(startpoint,help3);
  return help4;
}













int main() {
 struct cplx cn1,cn_current;
 struct cplx cn2;
 struct cplx cn_mul;
 struct cplx cn_pow;

 cn1.re = 3;
 cn1.im = 10;
 cn2.re = 1;
 cn2.im = 2;
 cn_mul=mul_cpx(cn1,cn2);
 cn_pow=pow_cpx(cn2,2);
 cn_current = cn2;
 for (size_t i = 0; i < 9; i++) {
 cn_current = newton_raphson_iteration(cn_current,2);
}

 // printf("a is: %f %f\n", cn_pow.re,cn_pow.im);
 // cn_pow=conj_cpx(cn_pow);
 // printf("a is: %f %f\n", cn_pow.re,cn_pow.im);
 printf("After one iteration: %f %f\n", cn_current.re,cn_current.im);
  return 0;
}

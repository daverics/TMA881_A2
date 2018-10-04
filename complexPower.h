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

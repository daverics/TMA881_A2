
long int * benchmark(void (*f)(),char * str,long int n){
  struct timespec ts;
  struct timespec ts1;
  timespec_get(&ts,TIME_UTC);
  for (size_t ix = 0; ix < n; ix++) {
    (*f)();
  }
  timespec_get(&ts1,TIME_UTC);
  long int sec =+ ts1.tv_sec - ts.tv_sec;
  long int nsec = ts1.tv_nsec - ts.tv_nsec;

  if (nsec < 0){
    sec = sec - 1;
    nsec = 1000000000 + nsec;
  }

  printf("Calculation of %s took:\n%d seconds and %d nanoseconds\n",str,sec,nsec);
  long int *tid = (long int*) malloc(sizeof(long int)*2);
  (*tid) = sec;
  (*(tid+1)) = nsec;
  return tid;
}

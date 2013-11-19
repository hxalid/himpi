# include "config.h"
# include "hdnla_timer.h"
#ifdef HDNLA_SMPI
#include <mpi.h>
#endif
#include <stdio.h>

/* this could be specific for some processors
 * the default solution seems to be accurate enough
#define CLOCK_TIMER CLOCK_MONOTONIC_RAW
*/
#define CLOCK_TIMER CLOCK_MONOTONIC

inline double get_microsecond(struct timespec *res){
  return (res->tv_sec*1000000 + res->tv_nsec/1000);
}
inline double get_nanosecond(struct timespec *res){
  return (res->tv_sec*1000000000 + res->tv_nsec);
}
inline double get_second(struct timespec *res){
  return (res->tv_sec + res->tv_nsec/1000000000);
}

double get_res_real(){
  struct timespec start, end, dummy;
  int i;
  int iter = 100000;
  double eslape = 0;
  for(i =0; i < iter; i++){
    get_time(&start);
    get_time(&dummy);
    get_time(&end);
    eslape += get_timediff(&start, &end);
  }
  return eslape/iter;
}

void init_timer(){
  struct timespec res;
#ifdef HAVE_LIBRT
  if(get_res(&res)==0){
    //fprintf(stderr, "the clock resolution is equal to %f microseconds\n",
    //          get_res_real());
  }else{
    perror("the resolution of the clock lead to an error");
  }
#endif
}

int get_res(struct timespec *res){
#ifdef HAVE_LIBRT
  return clock_getres(CLOCK_TIMER, res);
#else
  return (int)-1;
#endif
}

inline int get_time(struct timespec *tp){
#ifdef HDNLA_SMPI
  double time = MPI_Wtime();
  time_t value = (time_t)floor(time);
  time -= (double) value;
  time = time * 1000000000;
  tp->tv_nsec = (long) time;
  tp->tv_sec = value ;
#else
#ifdef HAVE_LIBRT
  return clock_gettime(CLOCK_TIMER, tp);
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  tp->tv_sec = tv.tv_sec;
  tp->tv_nsec = tv.tv_usec*1000;
#endif
#endif
  return 0;
}

double get_timediff(struct timespec *start, struct timespec *end){
  return (double)(-start->tv_sec - ((double)start->tv_nsec)/1000000000 + end->tv_sec + ((double)end->tv_nsec)/1000000000);
}

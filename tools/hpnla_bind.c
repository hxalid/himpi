/*!
 * file to bind the process on the node
 *
 * Authors: Quintin Jean-Noël
 */

#include "hpnla_debug.h"
#define __USE_GNU

#include "config.h"
#ifndef _GNU_SOURCE
/* this should be done inside config.h */
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#endif

#include <sched.h>
#include <string.h>
#include <pthread.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
 
int hpnla_bind_process(char* str_mask){
#ifdef USE_SCHED_AFFINITY
  size_t row = 0, col = 0;
  if (str_mask == NULL) return -1;
  /*------------------------------binding on processors ---------------------*/
  cpu_set_t  mask;
  CPU_ZERO(&mask);
  char *range, *strtok_ptr, *strtok_ptr2;

  /* here we try to create the CPU Mask based on a string which could
   * be: 3;6-9;10 and so on.
   * like this it's suitable for a process with several threads
   */
  while((range=strtok_r(str_mask,",",&strtok_ptr))!=NULL){
    debug_print(2, row, col, "arg_p: %s arg : %s\n",range, str_mask);
    char * str_err = strdup(range);/* for the error message*/
    char * str_cpu = strtok_r(range,"-",&strtok_ptr2);
    int first_cpu = atoi(str_cpu), end_cpu;
    str_cpu = strtok_r(NULL,"-",&strtok_ptr2);
    if(str_cpu == NULL){
      end_cpu=first_cpu;
    }else{
      end_cpu = atoi(str_cpu);
      str_cpu = strtok_r(NULL,"-",&strtok_ptr2);
      if(str_cpu!=NULL || end_cpu < first_cpu){
        fprintf(stderr,"Error with the range of processors %s \n", str_err);
        //TODO delete this part
        return -1;
      }
    }
    int i;
    for(i=first_cpu;i<=end_cpu;i++){
      CPU_SET(i, &mask);
    }
    //arg = strstr(arg,",");
    str_mask=NULL;
  }
  /* the mask is created we can bind*/
  int pid = getpid();
  int result = sched_setaffinity(pid, sizeof(mask), &mask);
  debug_print(1, row, col, "I am on the processor : %d and "
      "the number of processors: %d\n",
      sched_getcpu(),
      CPU_COUNT(&mask));
  return result;
#else
  return -1;
#endif
}

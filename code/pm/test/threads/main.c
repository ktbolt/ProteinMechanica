//*============================================================*
//*                   pthread test                             *
//*============================================================*

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>

#define NUM_ITERATIONS  10000

typedef struct {
  int id;
  int num_proc;
  int dim;
  } ThreadData;

void *thread_func (void *arg);

//*============================================================*
//*                          main                              *
//*============================================================*

int 
main (int argc, char *argv[])
  {
  pthread_t *threads;
  pthread_attr_t  pthread_custom_attr;
  ThreadData *thread_data;
  int i;
  int num_threads;
  double startwtime, endwtime;

  if (argc != 2) {
    printf("Usage: %s n\n  where n is no. of thread\n", argv[0]);
    exit(1);
    }

  // spawn threads //

  num_threads = atoi(argv[1]);
  threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
  thread_data = (ThreadData*)malloc(num_threads * sizeof(ThreadData));
  pthread_attr_init(&pthread_custom_attr);

  startwtime = clock();

  for (i = 0; i < num_threads; i++) {
    thread_data[i].id = i;
    thread_data[i].num_proc = num_threads;
    pthread_create(&threads[i], &pthread_custom_attr, thread_func, 
                   (void*)(thread_data+i));
    }

  for (i = 0; i < num_threads; i++) {
    pthread_join(threads[i], NULL);
    }

  endwtime = clock();
  double dt = (endwtime - startwtime) / CLOCKS_PER_SEC;
  printf (">>> num it = %d \n", NUM_ITERATIONS); 
  printf (">>> time = %g \n", dt); 

  // do same calculation //

  startwtime = clock();
  double sum = 0.0;

  for (i = 0; i < NUM_ITERATIONS; i++) {
    sum += (float)i * (float)i; 
    }

  endwtime = clock();
  dt = (endwtime-startwtime) / CLOCKS_PER_SEC;
  printf (">>> sum = %g   time = %g \n", sum, dt); 

  exit(0);
  }

//*============================================================*
//*                          thread_func                       *
//*============================================================*

void * 
thread_func (void *arg)
  {
  ThreadData *tdata = (ThreadData *)arg;
  int id, num_threads, dn;
  int i, j, ibegin, iend;
  double sum;

  id = tdata->id;
  num_threads = tdata->num_proc;
  dn = NUM_ITERATIONS / num_threads;
  ibegin = id*dn;
  iend = ibegin + dn;
  sum = 0.0;

  for (i = ibegin; i < iend; i++) {
    sum += (float)i * (float)i; 
    }
  }


/*========================================================================*
 * test executing pm wihin a thread pthread)                              *
 *========================================================================*/

#include <stdio.h>
#include <sys/types.h>
#include <pthread.h>

#define MAX_THREAD 1000

typedef struct {
  int id;
  } parm;

typedef struct work {
  void (*func)();
  void *arg;
  struct work *next;
  } work_t;

typedef struct tpool {
  int num_threads;
  int max_size;
  int current_size;
  pthread_t *threads;
  work_t *head, *tail;
  pthread_mutex_t lock;
  pthread_cond_t not_empty, not_full;
  } tpool_t;

/*------------------------------------------------*
 *                    tpool_thread                *
 *------------------------------------------------*/

void *
tpool_thread (void *arg)
  {
  //fprintf(stderr, "\n\n>>> ***** tpool_thread ***** \n\n");
  //fflush (stderr);

  work_t *wl;
  tpool_t *tpl = (tpool_t*)arg;

  for (;;) {
    pthread_mutex_lock (&(tpl->lock));

    while (tpl->current_size == 0) {
      pthread_cond_wait (&(tpl->not_empty), &(tpl->lock));
      }

    wl = tpl->head;
    tpl->current_size -= 1;

    if (tpl->current_size == 0) {
      tpl->head = tpl->tail = NULL;
      }
    else {
      tpl->head = wl->next;
      }

    if (tpl->current_size == tpl->max_size-1) {
      pthread_cond_broadcast(&(tpl->not_full));
      }

    pthread_mutex_unlock (&(tpl->lock));

    (*(wl->func))(wl->arg);

    free (wl);
    }

  fprintf(stderr, "\n\n>>> ***** tpool_thread done ***** \n\n");
  fflush (stderr);
  }

/*------------------------------------------------*
 *                    tpool_insert                *
 *------------------------------------------------*/
void
tpool_insert (tpool_t *tpl, void (*func)(), void *arg)
  {
  parm *p = (parm*)arg;
  printf(">>> insert work item for node=%d\n", p->id);
  work_t *wl;

  pthread_mutex_lock (&(tpl->lock));
  
  while (tpl->current_size == tpl->max_size) {
    pthread_cond_wait (&(tpl->not_full), &(tpl->lock));
    }

  wl = (work_t*)malloc(sizeof(work_t));
  wl->func = func;
  wl->arg = arg;
  wl->next = NULL;

  if (tpl->current_size == 0) {
    tpl->head = tpl->tail = wl;
    pthread_cond_signal(&(tpl->not_empty)); 
    }
  else {
    tpl->tail->next = wl;
    tpl->tail = wl;
    }

  tpl->current_size += 1;

  pthread_mutex_unlock (&(tpl->lock));
  }

/*------------------------------------------------*
 *                    tpool_init                  *
 *------------------------------------------------*/

tpool_t*
tpool_init (int num_threads, int max_size)
  {
  int i;
  tpool_t *tpl;

  tpl = (tpool_t*)malloc(sizeof(tpool_t));
  tpl->num_threads = num_threads;
  tpl->max_size = max_size;
  tpl->current_size = 0;
  tpl->head = tpl->tail = NULL;

  pthread_mutex_init (&(tpl->lock), NULL);
  pthread_cond_init (&(tpl->not_empty), NULL);
  pthread_cond_init (&(tpl->not_full), NULL);

  tpl->threads = (pthread_t*)malloc(sizeof(pthread_t)*num_threads);

  for (i = 0; i < num_threads; i++) {
    pthread_create (&(tpl->threads[i]), NULL, tpool_thread, (void*)tpl);
    }

  return tpl;
  }

/*------------------------------------------------*
 *                    pmexec                      *
 *------------------------------------------------*/

void*
pmexec (void *arg)
  {
  char cmd[1000];

  parm *p=(parm *)arg;
  printf(">>> excute pm on node=%d\n", p->id);

  sprintf (cmd, "p sim_id=%d sim.pmc", p->id);
  system (cmd);

  printf(">>> done excute pm on node=%d\n", p->id);
  return (NULL);
  }

/*------------------------------------------------*
 *                    main                        *
 *------------------------------------------------*/

void 
main (int argc, char* argv[]) 
  {
  int n,i;
  pthread_t *threads;
  pthread_attr_t pthread_custom_attr;
  parm *p;

  if (argc != 2) {
    printf ("Usage: %s n\n  where n is no. of threads\n",argv[0]);
    exit(1);
    }

  n = atoi(argv[1]);

  if ((n < 1) || (n > MAX_THREAD)) {
    printf ("The no of thread should between 1 and %d.\n",MAX_THREAD);
    exit(1);
    }

  tpool_t* tpool = tpool_init (n, 1000);
  p = (parm*)malloc(sizeof(parm)*n);

  pthread_attr_init (&pthread_custom_attr);

  for (i = 0; i < n; i++) {
    p[i].id = i;
    tpool_insert (tpool, pmexec, (void *)(p+i));
    }

  /* synchronize the completion of each thread. */

  int err;

  for (i = 1; i < tpool->num_threads; i++) {
    err = pthread_join (tpool->threads[i], NULL);
    //printf (">>> err=%d \n", err);
    }

  }


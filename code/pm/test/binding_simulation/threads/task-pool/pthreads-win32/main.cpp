//*============================================================*
//* test master-worker task pool using c++           ==========*
//*============================================================*

#include <stdlib.h>
#include <string>
#include <sstream>
#include <list>
#include <iostream>

#include "pthread.h"
#include <signal.h>
//#include <sys/select.h>

using namespace std;

// the muticies, protectors of the shared resources //

pthread_mutex_t coutLock;
pthread_mutex_t task_queue_lock;
pthread_mutex_t finished_queue_lock;

// queue data //

list<string> task_queue;
list<string> finished_queue;

// information to pass to worker threads //

struct ThreadData { 
  unsigned long num;
  };

extern "C" {
  void *workerThread(void *threadarg);
  }

//*============================================================*
//*==========              workerThread              ==========*
//*============================================================*

void *
workerThread (void *threadarg) 
  {
  struct ThreadData *my_data;
  my_data = (ThreadData *) threadarg;
  int taskid = my_data->num;

  stringstream ss; 
  ss << taskid; 
  string taskString = ss.str();
  bool more_tasks = true;

  // keep on working until task_queue is empty  //

  while (more_tasks) {
    pthread_mutex_lock (&task_queue_lock);
    string workOnMe;

    if (task_queue.size() == 0) { 
      more_tasks = false; 
      }
    else {
      workOnMe = task_queue.front();
      task_queue.pop_front();
      }

    pthread_mutex_unlock (&task_queue_lock);

    if (!more_tasks) {
      break;
      }

    // execute pm cmd // 

    system (workOnMe.c_str());

    workOnMe = "thread " + taskString + " worked on " + workOnMe;

    // let's pretend this takes some time, add a delay to the computation

    /*
    struct timeval timeout;
    timeout.tv_sec = 0;
    timeout.tv_usec = 100000; // 0.1 second delay
    select (0, NULL, NULL, NULL, & timeout);
    */

    pthread_mutex_lock (&finished_queue_lock);
    finished_queue.push_back (workOnMe);
    pthread_mutex_unlock (&finished_queue_lock);
    }

  pthread_exit (NULL);
  }

//*============================================================*
//*==========                  main                  ==========*
//*============================================================*

int 
main (int argc, char *argv[])
  {
  unsigned long comp_DONE = 0; 
  unsigned long comp_START = 0;

  // set-up the mutexes //

  pthread_mutex_init (&coutLock,     NULL);
  pthread_mutex_init (&task_queue_lock,  NULL);
  pthread_mutex_init (&finished_queue_lock, NULL);

  if (argc != 3) { 
    cout << "Program requires two arguments: (1) number of threads to use,"
            " and (2) tasks to accomplish.\n"; 
    exit(1); 
    }

  unsigned long num_threads = atoi(argv[1]);
  unsigned long num_tasks = atoi(argv[2]);

  cout << ">>> number of threads=" << num_threads << endl;
  cout << ">>> number of tasks="   << num_tasks << endl;

  // fill task_queue with pm cmds // 

  for (unsigned long i = 0; i < num_tasks; i++) {
    stringstream ss; 
    ss << "p sim_id=" << i << " sim.pmc"; 
    task_queue.push_back(ss.str());
    }

  // start the worker threads //

  list<pthread_t*>    threadIdList; 
  list<ThreadData> thread_table; 

  // create threads //

  for (unsigned long i = 0; i < num_threads; i++) {
    pthread_t *tId = new pthread_t;   
    threadIdList.push_back(tId);
    ThreadData Y; 
    Y.num = i; 
    thread_table.push_back(Y);

    int rc = pthread_create (tId, NULL, workerThread, (void*)(&(thread_table.back())));

    if (rc) { 
      cout << "ERROR; return code from pthread_create() " << comp_START << "\n"; 
      cout.flush();
      exit(-1); 
      }
    }

  // now we wait for the threads to terminate //

  string stringOut;

  // poll the queue to get a status update on computation //

  while (comp_DONE != num_tasks) {
    pthread_mutex_lock (&task_queue_lock);
    comp_START = num_tasks - task_queue.size();
    pthread_mutex_unlock (&task_queue_lock);

    pthread_mutex_lock (&finished_queue_lock);
    comp_DONE = finished_queue.size();
    pthread_mutex_unlock (&finished_queue_lock);
    } 

  // call join to kill all worker threads //

  list<pthread_t*>::iterator i = threadIdList.begin();

  while (i != threadIdList.end()) {
    if (pthread_join(*(*i), NULL)!=0) { 
      cout << "Thread join error!\n"; 
      exit(1); 
      }

    delete (*i);
    threadIdList.erase(i++);  
    }

  cout << "\n";

  // let the user know what happened //

  for (list<string>::iterator i = finished_queue.begin(); i != finished_queue.end(); i++) {
    cout << (*i) << "\n";
    }

  // clean-up //

  pthread_mutex_destroy (&coutLock);
  pthread_mutex_destroy (&task_queue_lock);  
  pthread_mutex_destroy (&finished_queue_lock);  

  // pthread_exit(NULL);
  }


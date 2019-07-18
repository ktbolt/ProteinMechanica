
#include<stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

/*
main ()
  {
  system("D:/cygwin/home/dparker/bin/pnog.exe D:/cygwin/home/dparker/software/protmech/protmech/trunk/pm/test/exec/d");

  } 
*/

main ()
  {
  int fork_return;
  int count = 0;

  /* getpid() returns the process id of this process. */

  printf("Process %d about to fork a child.\n", getpid() );

  fork_return = fork();

  if (fork_return < 0) {
    printf("Unable to create child process, exiting.\n");
    exit(-1);
    }

  /* fork_return is the pid of the child process and I am
     the parent. Start printing a's. */

  if (fork_return > 0) {
    printf("Created child process %d.\n", fork_return);

    fork_return = fork();

    if (fork_return > 0) {
      while( count++ < 800) {
        putchar('a');

        if (count % 80 == 0){
          putchar('\n');
          sleep(2);
          }
        }
      }
    else {
      while( count++ < 800) {
        putchar('c');

        if (count % 80 == 0){
          putchar('\n');
          sleep(2);
          }
        }
      }
    }

  /* A 0 return tells me that I am the child. Print b's */ 

  else {

    while(count++ < 800) {
      putchar('b');

      if (count % 80 == 0) {
        putchar('\n');
        sleep(2);
        }
      }
    }
  }
  

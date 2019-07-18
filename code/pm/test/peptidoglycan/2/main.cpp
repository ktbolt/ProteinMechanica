
//*============================================================*
//* app:    p r o t e i n   m o d e l e r   a p p              *
//*============================================================*

#include "pm.h"

using namespace ProteinModeler;

int
main (int num_args, char **args) 
  {

  char *script;

  if (num_args == 2) {
    script = args[1];
    }
  else {
    script = NULL;
    }


  // initialize protein modeler

  pmSystem.init();


  // enter event processing loop

  pmSystem.procEvents (script);

  return (0);
  }


 
/* -------------------------------------------------------------------------- *
 *                              Protein Mechanica                             *
 * -------------------------------------------------------------------------- *
 * This is part of the Protein Mechanica coarse-grained molecular motor       *
 * modeling application originating from Simbios, the NIH National Center for *
 * Physics-Based Simulation of Biological Structures at Stanford, funded      *
 * under the NIH Roadmap for Medical Research, grant U54 GM072970.            *
 * See https://simtk.org.                                                     *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: David Parker                                                      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//*============================================================*
//* msg:                   m e s s a g e s                     *
//*============================================================*

#include "pm/pm.h"
#include <stdarg.h>

namespace ProteinMechanica {

///////////////////////////////////////////////////////////////
//        I n f o r m a t i o n a l   M e s s a g e s       //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========               pm_PrintMsg              ==========*
//*============================================================*

void
pm_PrintMsg (const char *prompt, char *format, ...)
  {

  va_list ap;
  int c, ia, i, n;
  float fa;
  char *sa;
  char des[3];
  char msg[1000], str[1000], tf[1000];

  if (!pmSystem.getCmdVerbose()) {
    return;
    }

  i = 0;
  n = 0;
  va_start (ap, format);
  strcpy (tf, format);
  //fn = va_arg (ap, char *);
  sprintf (msg, "%s", prompt);

  /*  parse the message format. when a 
      format (%) character is encountered
      extract an argument from the arglist
      of this function and print it into
      msg[] using the format parsed so far.  */
      
  while ((c = *format++) != 0) {
    if (i == 0) {
      if (c == '%') {
        des[i++] = c;
        }
      tf[n++] = c;
      }
    else if (i == 1) {
      des[i++] = c;
      des[i] = '\0';
      tf[n++] = c;
      tf[n] = '\0';


      /*  int format.  */

      if (!strcmp(des, "%d")) {
	ia = va_arg (ap, int);
        sprintf (str, tf, ia);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  float format.  */

      else if (!strcmp(des, "%f") || !strcmp(des, "%g")) {
	fa = va_arg (ap, double);
        sprintf (str, tf, fa);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  string format.  */

      else if (!strcmp(des, "%s")) {
	sa = va_arg (ap, char *);
        sprintf (str, tf, sa);
        sprintf (msg, "%s%s", msg, str);
	}

      i = 0;
      n = 0;
      }
    }

  va_end (ap);
  tf[n] = '\0';


  /*  finally, print the entrire message.  */

  fprintf (stderr, "%s%s \n", msg, tf); 
  }

///////////////////////////////////////////////////////////////
//               E r r o r   M e s s a g e s                //
/////////////////////////////////////////////////////////////

#define DM_ERROR_MAX  500
static int error_print;
static char error_msg[DM_ERROR_MAX];


/*------------------------------------------------------------*
 *                                                            *
 *               ****  pm_ErrorInit  ****                     *
 *                                                            *
 * initialize the error message queu.                         *
 *------------------------------------------------------------*/

void
pm_ErrorInit ()
  {

 /**************
  ***  body  ***
  **************/

  error_print = true;
  error_msg[0] = '\0';
  }

/*------------------------------------------------------------*
 *                                                            *
 *               ****  pm_ErrorReport  ****                   *
 *                                                            *
 * report an error message.                                   *
 *------------------------------------------------------------*/

void
pm_ErrorReport (char *mod, char *format, ...)
  {

  va_list ap;

  int c, ia, i, n;

  float fa;

  char *sa, *fn, ic;

  char des[3];

  char msg[1000], str[1000], tf[1000];

  i = 0;
  n = 0;
  va_start (ap, format);
  strcpy (tf, format);
  fn = va_arg (ap, char *);

  if (*fn == '*') {
    sprintf (msg, "\n **** %s ERROR: ", mod);
    }
  else {
    sprintf (msg, "\n **** %s ERROR [%s]: ", mod, fn);
    }


  /*  parse the message format. when a 
      format (%) character is encountered
      extract an argument from the arglist
      of this function and print it into
      msg[] using the format parsed so far.  */
      
  while ((c = *format++) != 0) {
    if (i == 0) {
      if (c == '%') {
        des[i++] = c;
        }
      tf[n++] = c;
      }
    else if (i == 1) {
      des[i++] = c;
      des[i] = '\0';
      tf[n++] = c;
      tf[n] = '\0';


      /*  int format.  */

      if (!strcmp(des, "%d")) {
	ia = va_arg (ap, int);
        sprintf (str, tf, ia);
        sprintf (msg, "%s%s", msg, str);
	}

      else if (!strcmp(des, "%c")) {
	ic = va_arg (ap, int);
        sprintf (str, tf, ic);
        sprintf (msg, "%s%s", msg, str);
	}

      /*  float format.  */

      else if (!strcmp(des, "%f")) {
	fa = va_arg (ap, double);
        sprintf (str, tf, fa);
        sprintf (msg, "%s%s", msg, str);
	}

      /*  string format.  */

      else if (!strcmp(des, "%s")) {
	sa = va_arg (ap, char *);
        sprintf (str, tf, sa);
        sprintf (msg, "%s%s", msg, str);
	}

      i = 0;
      n = 0;
      }
    }

  va_end (ap);
  tf[n] = '\0';


  /*  finally, print the entrire message.  */

  fprintf (stderr, "%s%s \n", msg, tf); 
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  pm_ErrorJmp  ****                      *
 *                                                            *
 * error long jump.                                           *
 *------------------------------------------------------------*/

void
pm_ErrorJmp (int error)
  {

 /**************
  ***  body  ***
  **************/

  /*
  longjmp (GrErrorEnv, error);
  */
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  pm_ErrorLog  ****                      *
 *                                                            *
 * log an error message.                                      *
 *------------------------------------------------------------*/

void
pm_ErrorLog (char *format, ...)
  {

  va_list ap;

  int c, ia, i;

  float fa;

  char *sa, *fn;

  char des[3], msg[1000], str[1000], tf[1000];

  int n;

 /**************
  ***  body  ***
  **************/

  i = 0;
  n = 0;
  va_start (ap, format);
  strcpy (tf, format);

  fn = va_arg (ap, char *);

  if (*fn == '*') {
    sprintf (msg, "\n **** DM ERROR: ");
    }
  else {
    sprintf (msg, "\n **** DM ERROR [%s]: ", fn);
    }


  /*  parse the message format. when a 
      format (%) character is encountered
      extract an argument from the arglist
      of this function and print it into
      msg[] using the format parsed so far.  */
      
  while ((c = *format++) != 0) {
    if (i == 0) {
      if (c == '%') {
        des[i++] = c;
        }
      tf[n++] = c;
      }
    else if (i == 1) {
      des[i++] = c;
      des[i] = '\0';
      tf[n++] = c;
      tf[n] = '\0';


      /*  int format.  */

      if (!strcmp(des, "%d")) {
	ia = va_arg (ap, int);
        sprintf (str, tf, ia);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  float format.  */

      else if (!strcmp(des, "%f")) {
	fa = va_arg (ap, double);
        sprintf (str, tf, fa);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  string format.  */

      else if (!strcmp(des, "%s")) {
	sa = va_arg (ap, char *);
        sprintf (str, tf, sa);
        sprintf (msg, "%s%s", msg, str);
	}

      i = 0;
      n = 0;
      }
    }

  va_end (ap);
  tf[n] = '\0';
  sprintf (error_msg, "%s%s \n", msg, tf); 
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  pm_ErrorLogPr  ****                    *
 *                                                            *
 * print a logged error message.                              *
 *------------------------------------------------------------*/

void
pm_ErrorLogPr ()
  {

  if (error_msg[0] == '\0') {
    return;
    }

  fprintf (stderr, "%s\n", error_msg); 
  error_msg[0] = '\0';
  }

/*------------------------------------------------------------*
 *                                                            *
 *               ****  pm_ErrorWarnReport  ****               *
 *                                                            *
 * report a warning message.                                  *
 *------------------------------------------------------------*/

void
pm_ErrorWarnReport (char *mod, char *format, ...)
  {

  va_list ap;

  int c, ia, i;

  float fa;

  char *sa, *fn;

  char des[3];

  char msg[1000], str[1000], tf[1000];

  int n;

  i = 0;
  n = 0;
  va_start (ap, format);
  strcpy (tf, format);

  fn = va_arg (ap, char *);

  if (*fn == '*') {
    sprintf (msg, "\n **** %s WARNING: ", mod);
    }
  else {
    sprintf (msg, "\n **** %s WARNING [%s]: ", mod, fn);
    }  


  /*  parse the message format. when a 
      format (%) character is encountered
      extract an argument from the arglist
      of this function and print it into
      msg[] using the format parsed so far.  */
      
  while ((c = *format++) != 0) {
    if (i == 0) {
      if (c == '%') {
        des[i++] = c;
        }
      tf[n++] = c;
      }
    else if (i == 1) {
      des[i++] = c;
      des[i] = '\0';
      tf[n++] = c;
      tf[n] = '\0';


      /*  int format.  */

      if (!strcmp(des, "%d")) {
	ia = va_arg (ap, int);
        sprintf (str, tf, ia);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  float format.  */

      else if (!strcmp(des, "%f")) {
	fa = va_arg (ap, double);
        sprintf (str, tf, fa);
        sprintf (msg, "%s%s", msg, str);
	}


      /*  string format.  */

      else if (!strcmp(des, "%s")) {
	sa = va_arg (ap, char *);
        sprintf (str, tf, sa);
        sprintf (msg, "%s%s", msg, str);
	}

      i = 0;
      n = 0;
      }
    }

  va_end (ap);
  tf[n] = '\0';


  /*  finally, print the entire message.  */

  fprintf (stderr, "%s%s \n", msg, tf); 
  }

}

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
//* csys        c o m p u t e r   s y s t e m                  *
//*============================================================*

#include "csys.h" 

namespace ProteinMechanica {


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_OsUserNameGet  ****              *
 *                                                            *
 * get user name.                                             *
 *------------------------------------------------------------*/

void
pm_getUserName (char *name)
  {
  char *s;
  s = getenv ("USER");
  strcpy (name, s);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_SysFileExists  ****              *
 *                                                            *
 *------------------------------------------------------------*/

int 
pm_FileExists (char *name)
  {

#if defined(_DM_UNIX_) || defined (_DM_CYGWIN_)

  int res;

  struct stat sbuf;

 /**************
  ***  body  ***
  **************/

  res = stat (name, &sbuf);

#else

  int res;

  res = 1;

#endif

  return (res != -1);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_SysDirRead  ****                 *
 *                                                            *
 *------------------------------------------------------------*/

void
pm_SysDirRead (char *dir_name, int *p_num_files, char ***p_file_names)
  {

#if defined(_DM_UNIX_) || defined (_DM_CYGWIN_)

  DIR *dirp;

#if defined(_DM_UNIX_) 
  struct direct *dp;
#else 
  struct dirent *dp;
#endif

  int num_files;

  char **file_names;

 /**************
  ***  body  ***
  **************/

  *p_num_files = 0;
  *p_file_names = NULL;
  dirp = opendir (dir_name);

  if (dirp == NULL) {
    return;
    }

  num_files = 0;

  while ((dp = readdir(dirp)) != NULL) {
    num_files += 1;
    }

  closedir (dirp);
  dirp = opendir (dir_name);
  file_names = (char**)malloc(sizeof(char*) * num_files);
  num_files = 0;

  while ((dp = readdir(dirp)) != NULL) {
    file_names[num_files++] = strdup (dp->d_name);
    }

  closedir (dirp);
  *p_num_files = num_files;
  *p_file_names = file_names;

#endif

  }

}

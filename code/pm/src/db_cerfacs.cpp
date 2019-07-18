
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

/////////////////////////////////////////////////////////////////
//        c e r f a c s    I n t e r f a c e                  //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmDbCerfacsInterface::PmDbCerfacsInterface(const string name) {
  this->name = name;
  }

//*============================================================*
//*==========                 open                   ==========*
//*============================================================*

void
PmDbCerfacsInterface::open (const string name, PmDbModeType otype, string dtype) {
  fprintf (stderr, ">>>>>> CERFACS: open [%s] \n", name.c_str());

  if (otype == PM_DB_MODE_READ) {
    this->fp = fopen (name.c_str(), "r");
    }
  else {
    this->fp = fopen (name.c_str(), "w");
    }
  }

//*============================================================*
//*==========                 getModes               ==========*
//*============================================================*

void 
PmDbCerfacsInterface::getModes(vector<PmMode>& modes)
  {
  char *s, line[1000];
  int n, vec_len, num;
  bool all_modes = false;
  int id;
  float val;
  PmMode mode;

  if (!fp) return; 
  fprintf (stderr, "    >>> reading cerfacs mode file \n");

  pm_DbLineGet (this->fp, line);
  s = (char *)strtok (line, " ");

  if (strcmp(s, "VECTOR")) {
    fprintf (stderr, "\n  ****  error: file not in cerfacs format. \n");
    return;
    }

  pm_DbLineGet (fp, line);
  vec_len = 0;

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    s = (char *)strtok (line, " ");

    if (*s == 'V') {
      break;
      }

    vec_len += 1;
    }

  fprintf (stderr, "    >>> eigenvector length = %d \n", vec_len);

  if (all_modes) {
    fprintf (stderr, "    >>> reading all modes. \n");
    n = 1;
    while (pm_DbLineGet(fp, line)) {
      n += 1;
      }

    num = n / (vec_len + 2)  + 1;
    }
  else {
    fprintf (stderr, "    >>> reading 1st 50 modes \n");
    num = 50;
    }

  for (int i = 0; i < num; i++) {
    if (!pm_DbLineGet (fp, line)) {
      num = i;
      break;
      }

    sscanf (line, "VECTOR %d VALUE %f\n", &id, &val);
    pm_DbLineGet (fp, line);

    mode.size = vec_len;
    mode.eigen_value = val;
    mode.eigen_vectors = new PmVector3[vec_len];

    for (int j = 0; j < vec_len; j++) {
      fscanf (fp, "%f %f %f \n", &mode.eigen_vectors[j][0], 
                                 &mode.eigen_vectors[j][1],
                                 &mode.eigen_vectors[j][2]);
      }
   
    modes.push_back(mode);
    }
  }

//*============================================================*
//*==========                 close                  ==========*
//*============================================================*

void
PmDbCerfacsInterface::close() {
  fclose (this->fp);
  }

//*============================================================*
//*==========                 read                   ==========*
//*============================================================*

void
PmDbCerfacsInterface::read(const string mol_name) {
  }



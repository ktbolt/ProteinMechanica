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
//    P r o t e i n   M o d e l e r   I n t e r f a c e       //
///////////////////////////////////////////////////////////////


// define pm file processing table

typedef void (PmDbPmInterface::*Pmfp)(char*);

typedef struct PmDbPmProc {
  char *name;
  Pmfp func;
  } PmDbPmProc;

static PmDbPmProc  pm_table[] = {
    {"grid",     &PmDbPmInterface::procGrid},
    {"particle", &PmDbPmInterface::procParticle},
    {"surface",  &PmDbPmInterface::procSurface},
    {"trace",    &PmDbPmInterface::procTrace},
    {NULL,       NULL}};

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmDbPmInterface::PmDbPmInterface(const string name) 
  {
  this->name = name;
  grid = NULL;
  mol = NULL;
  particle = NULL;
  surf = NULL;
  trace = NULL;
  }

//*============================================================*
//*==========              open                      ==========*
//*============================================================*
// open pm file.

void
PmDbPmInterface::open (const string name, PmDbModeType otype, string dtype) 
  {
  //fprintf (stderr, "    >>> PM: open file \"%s\" \n", name.c_str());

  if (otype == PM_DB_MODE_READ) {
    this->fp = fopen (name.c_str(), "rb");
    }
  else {
    this->fp = fopen (name.c_str(), "wb");
    }
  }

//*============================================================*
//*==========              close                     ==========*
//*============================================================*

void
PmDbPmInterface::close() {
  fclose (this->fp);
  }

//*============================================================*
//*==========              getMolecule               ==========*
//*============================================================*

PmMolecule *
PmDbPmInterface::getMolecule(int model, PmMoleculeType type) {
  return mol;
  }

PmMolecule *
PmDbPmInterface::getMolecule(string model, PmMoleculeType type) {
  return mol;
  }

//*============================================================*
//*==========              getGrid                   ==========*
//*============================================================*

PmGrid *
PmDbPmInterface::getGrid() {
  return grid;
  }

//*============================================================*
//*==========              getParticle               ==========*
//*============================================================*

PmParticle*
PmDbPmInterface::getParticle() {
  return particle;
  }

//*============================================================*
//*==========              getSurface                ==========*
//*============================================================*

PmSurface *
PmDbPmInterface::getSurface() {
  return surf;
  }

//*============================================================*
//*==========              getTrace                  ==========*
//*============================================================*

PmTrace *
PmDbPmInterface::getTrace() {
  return trace;
  }

//*============================================================*
//*==========              read                      ==========*
//*============================================================*
// read a pm file.

void
PmDbPmInterface::read(const string dbname) 
  {
  char line[DB_MAX_LINE];

  if (!fp) {
    return;
    }

  if (dbname != "") {
    name = dbname;
    }

  //fprintf (stderr, "    >>> read pm file\n");

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    for (int i = 0; pm_table[i].name != NULL; i++) {
      int n = strlen (pm_table[i].name);

      if (!strncasecmp(line, pm_table[i].name, n)) {
        (this->*pm_table[i].func) (line);
        }
      }
    }
  //fprintf (stderr, "    done.\n");
  }


//*============================================================*
//*==========              writeParticle             ==========*
//*============================================================*
// write a particle pm file.

void
PmDbPmInterface::writeParticle(string name, vector<PmVector3> verts, vector<float> data)
  {
  if (verts.size() == 0) { 
    return;
    }

  if (!fp) {
    return;
    }

  int num_verts; 
  PmVector3 v;

  num_verts = verts.size();

  fprintf (fp, "# pm particle file \n");
  fprintf (fp, "particles \"%s\" \n", name.c_str());
  fprintf (fp, "units nm \n");

  fprintf (fp, "binary\n");
  fwrite (&num_verts, sizeof(int), 1, fp);

  for (int i = 0; i < num_verts; i++) {
    v = verts[i];
    fwrite (&verts[i], sizeof(PmVector3), 1, fp);
    }

  
  if (data.size() == 0) { 
    return;
    }

  fwrite (&num_verts, sizeof(int), 1, fp);

  for (int i = 0; i < num_verts; i++) {
    fwrite (&data[i], sizeof(float), 1, fp);
    }
  }

//*============================================================*
//*==========              writeSurface              ==========*
//*============================================================*
// wrte a surface pm file.

void
PmDbPmInterface::writeSurface(string name, vector<PmVector3> verts, vector<int> polys)
  {
  if ((verts.size() == 0) || (polys.size() == 0)) {
    return;
    }

  if (!fp) {
    return;
    }

  int num_verts, num_tri;
  PmVector3 v;

  num_verts = verts.size();
  num_tri = polys.size();

  fprintf (fp, "# pm surface file \n");
  fprintf (fp, "# source \"%s\"\n", "isosurface");
  fprintf (fp, "surface \"%s\" \n", name.c_str());
  fprintf (fp, "units nm \n");

  fprintf (fp, "binary\n");
  fwrite (&num_verts, sizeof(int), 1, fp);
  fwrite (&num_tri, sizeof(int), 1, fp);

  for (int i = 0; i < num_verts; i++) {
    v = verts[i];
    fwrite (&verts[i], sizeof(PmVector3), 1, fp);
    }

  for (int i = 0; i < num_tri; i++) {
    fwrite (&polys[i], sizeof(int), 1, fp);
    }
  }

////////////////////////////////////////////////////////////////
//         functions to process pm file entries              //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pm_DbPmGridProc           ==========*
//*============================================================*
// process pm grid object.

void
PmDbPmInterface::procGrid(char *line)
  {
  /*
  fprintf (stderr, "\n>>>>>> PmDbPmInterface::procGrid\n");
  fprintf (stderr, ">>>>>> line = [%s] \n", line);
  */
  char *s;
  string grid_name;
  int nx, ny, nz;
  float dx, dy, dz, *data, vmin, vmax;
  PmExtent extent;

  s = strtok (line, "\"");
  s = strtok (NULL, "\"");

  if (name != "") {
    grid_name = name;
    }
  else {
    grid_name = s;
    }

  fprintf (stderr, "    >>> read grid name = %s \n", grid_name.c_str());

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    if (!strncmp(line, "binary", 6)) {
      fprintf (stderr, "    >>> binary file \n");
      fread  (&nx, sizeof(int), 1, fp);
      fread  (&ny, sizeof(int), 1, fp);
      fread  (&nz, sizeof(int), 1, fp);
      fprintf (stderr, "    >>> nx = %d  ny = %d  nz = %d \n", nx, ny, nz);

      fread  (&dx, sizeof(float), 1, fp);
      fread  (&dy, sizeof(float), 1, fp);
      fread  (&dz, sizeof(float), 1, fp);
      fprintf (stderr, "    >>> dx = %f  dy = %f  dz = %f \n", dx, dy, dz);

      data = new float[nx*ny*nz];
      fread  (data, sizeof(float), nx*ny*nz, fp);

      vmin = vmax = data[0];

      for (int i = 0; i < nx*ny*nz; i++) {
        if (data[i] > vmax) vmax = data[i];
        if (data[i] < vmin) vmin = data[i];
        }

      fprintf (stderr, "    >>> vmin = %f  vmax = %f  \n", vmin, vmax); 
      }
    }

  extent.min[0] = 0.0; extent.max[0] = dx*nx;
  extent.min[1] = 0.0; extent.max[1] = dy*ny;
  extent.min[2] = 0.0; extent.max[2] = dz*nz;

  grid = new PmGrid(grid_name, nx, ny, nz, dx, dy, dz, data);
  grid->setExtent (extent);
  }

//*============================================================*
//*==========              procParticle              ==========*
//*============================================================*
// process pm particle object.

void
PmDbPmInterface::procParticle(char *line)
  {
  char *s; 
  string part_name;
  float scale;

  s = strtok (line, "\"");
  s = strtok (NULL, "\"");

  if (name != "") {
    part_name = name;
    }
  else {
    part_name = s;
    }

  fprintf (stderr, "    >>> read particle name = %s \n", part_name.c_str());
  scale = 0.1;

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    if (!strncmp(line, "binary", 6)) {
      fprintf (stderr, "    >>> binary file \n");
      readBinaryParticle(part_name, scale); 
      return;
      }

    else if (!strncmp(line, "units", 5)) {
      if (!strncmp(line, "units nm", 8)) {
        scale = 1.0;
        }
      }
    }
  }

//*============================================================*
//*==========              readBinaryParticle        ==========*
//*============================================================*
// read binary particle file.  

void
PmDbPmInterface::readBinaryParticle(const string part_name, float scale)
  {
  int num_verts;
  PmVector3 *verts;
  float *data;
  int *conn, size;
  int i;
  PmExtent extent;

  fread (&num_verts, sizeof(int), 1, fp);
  verts = new PmVector3[num_verts];
  fread (verts, sizeof(PmVector3), num_verts, fp);

  for (i = 0; i < num_verts; i++) {
    verts[i][0] *= scale;
    verts[i][1] *= scale;
    verts[i][2] *= scale;

    if (i == 0) {
      extent.min = verts[i];
      extent.max = verts[i];
      }
    else {
      extent.update (verts[i]);
      }
    }

  fprintf (stderr, "    >>> read %d vertices. \n", num_verts);
  particle = new PmParticle(part_name, num_verts, verts);
  particle->setExtent (extent);

  // read mass data //

  fread (&num_verts, sizeof(int), 1, fp);

  if (num_verts) {
    fprintf (stderr, "    >>> read %d data. \n", num_verts);
    data = new float[num_verts];
    fread (data, sizeof(float), num_verts, fp);
    particle->setMasses(data);
    }
  else {
    fprintf (stderr, "    >>> no data. \n");
    }
  }

//*============================================================*
//*==========              pm_DbPmSurfProc           ==========*
//*============================================================*
// process pm surface object.

void
PmDbPmInterface::procSurface (char *line)
  {
  //fprintf (stderr, "\n>>>>>> PmDbPmInterface::procSurface \n");
  //fprintf (stderr, ">>>>>> line = [%s] \n", line);
  char *s; 
  string surf_name;
  float scale;
  bool has_charge = false;

  s = strtok (line, "\"");
  s = strtok (NULL, "\"");

  if (name != "") {
    surf_name = name;
    }
  else {
    surf_name = s;
    }

  //fprintf (stderr, "    >>> read surface name = %s \n", surf_name.c_str());
  scale = 0.1;

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    if (!strncmp(line, "binary", 6)) {
      //fprintf (stderr, "    >>> binary file \n");
      readBinarySurf(surf_name, scale, has_charge); 
      return;
      }

    else if (!strncmp(line, "charge", 6)) {
      has_charge = true;
      }

    else if (!strncmp(line, "units", 5)) {
      if (!strncmp(line, "units nm", 8)) {
        scale = 1.0;
        }
      }
    }
  }

//*============================================================*
//*==========              readBinarySurf            ==========*
//*============================================================*
// read binary surface file.  

void
PmDbPmInterface::readBinarySurf(const string surf_name, float scale, bool has_charge)
  {
  int num_verts;
  PmVector3 *verts;
  int num_tri;
  int *conn, size;
  int i;
  PmExtent extent;
  float *charge;

  fread (&num_verts, sizeof(int), 1, fp);
  fread (&num_tri, sizeof(int), 1, fp);
  verts = new PmVector3[num_verts];
  fread (verts, sizeof(PmVector3), num_verts, fp);

  size = 3*num_tri;
  conn = new int[size];
  fread (conn, sizeof(int), size, fp);

  for (i = 0; i < num_verts; i++) {
    verts[i][0] *= scale;
    verts[i][1] *= scale;
    verts[i][2] *= scale;

    if (i == 0) {
      extent.min = verts[i];
      extent.max = verts[i];
      }
    else {
      extent.update (verts[i]);
      }
    }

  //fprintf (stderr, "    >>> read %d vertices. \n", num_verts);
  //fprintf (stderr, "    >>> read %d triangles. \n", num_tri);
  surf = new PmSurface(surf_name, num_verts, verts, num_tri, conn, 3);
  surf->setExtent (extent);

  if (has_charge) {
    fprintf (stderr, ">>> read charge \n");
    charge = new float[num_verts];
    fread (charge, sizeof(float), num_verts, fp);
    surf->addVertexData("charge", charge);
    }
  }

//*============================================================*
//*==========              procTrace                 ==========*
//*============================================================*
// process pm trace object.

void
PmDbPmInterface::procTrace(char *line)
  {
  //fprintf (stderr, "\n>>>>>> PmDbPmInterface::procTrace\n");
  //fprintf (stderr, ">>>>>> line = [%s] \n", line);

  char *s;
  string trace_name;
  int num_pts;
  //PmVector3 *points; 
  //vector<PmVector3> pts; 
  PmExtent extent;
  float t; 
  PmVector3 pt;

  s = strtok (line, "\"");
  s = strtok (NULL, "\"");

  if (name != "") {
    trace_name = name;
    }
  else {
    trace_name = s;
    }

  fprintf (stderr, "    >>> read trace name = %s \n", trace_name.c_str());
  trace = new PmTrace(trace_name);
  num_pts = 0;

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    if (!strncmp(line, "binary", 6)) {
      fprintf (stderr, "    >>> binary file \n");
      }

    else { 
      sscanf (line, "%f %f %f %f\n", &t, &pt[0], &pt[1], &pt[2]);
      trace->addPoint(t, pt);

      if (num_pts) {
        extent.update (pt);
        }
      else {
        extent.min = pt;
        extent.max = pt; 
        }

      num_pts += 1;
      }
    }

  trace->setExtent(extent);
  }

//*============================================================*
//*==========              readEnergy                ==========*
//*============================================================*

void
PmDbPmInterface::readEnergy(const string fname, vector<string>& names,
                            vector<vector<float> >& values)
  {
  #ifdef dbg_PmDbPmInterface_readEnergy
  fprintf (stderr, ">>>>>> PmDbPmInterface::readEnergy \n");
  #endif
  char *s, line[DB_MAX_LINE];
  FILE *fp;
  bool read_data = false;
  int row = 0;
  int num_cols, n;
  names.clear();
  values.clear();

  num_cols = 0;
  fp = fopen (fname.c_str(), "r");

  while (1) {
    if (!pm_DbLineGet(fp, line)) {
      break;
      }

    #ifdef dbg_PmDbPmInterface_readEnergy
    fprintf (stderr, ">>> line=\"%s\" \n", line);
    #endif

    if (read_data && !isdigit(line[0])) {
      break;
      }

    else if (isdigit(line[0])) {
      read_data = true;
      }

    if (read_data) {
      values.push_back(vector<float>());

      for (int i = 0; i < num_cols; i++) {
        if (i == 0) {
          s = strtok (line, " ");
          }
        else {
          s = strtok (NULL, " ");
          }

        values[row].push_back(atof(s));
        }

      //fprintf (stderr, ">>> row=%d \n", row);
      row += 1;
      }

    if (line[0] == '#') {
      continue;
      }

    if (!strncmp(line, "interaction names", 17) ||
        !strncmp(line, "restraint names", 15)) {
      s = strtok (line, "=");

      while (1) {
        s = strtok (NULL, " ");

        if (s[0] == '}') {
          if (!strncmp(line, "restraint names", 15)) {
            num_cols += names.size();
            }
          else {
            num_cols += names.size()+2;
            }
          break; 
          }

        if (s[0] != '{') {
          #ifdef dbg_PmDbPmInterface_readEnergy
          fprintf (stderr, ">>> \"%s\" \n", s);
          #endif
          names.push_back(s);
          }
        }
      }
    }

  #ifdef dbg_PmDbPmInterface_readEnergy
  fprintf (stderr, "\n\n------ values ------- \n");
  fprintf (stderr, ">>> num cols=%d \n", num_cols);

  for (int i = 0; i < row; i++) {
    fprintf (stderr, "%d: ", i+1);
    for (int j = 0; j < num_cols; j++) {
      fprintf (stderr, "%f ", values[i][j]);
      }
    fprintf (stderr, "\n");
    }
  #endif
  }


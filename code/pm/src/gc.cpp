
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
//* gc:         g e o m e t r i c    c o m p u t a t i o n     *
//*============================================================*

#include "gc.h"

namespace ProteinMechanica {

typedef struct _coord_bin_ {
  int id;
  int count;
  float x, y, z, d;
  struct _coord_bin_ *link;
  } GcCoordBin;


//*============================================================*
//*==========              pm_GcBasisComp            ==========*
//*============================================================*
// compute an orthonormal basis from the given normal vector.

void 
pm_GcBasisComp (PmVector3& normal, PmVector3& u, PmVector3& v)
  {

  u[0] = -normal[1];
  u[1] =  normal[0];
  u[2] =  0.0;

  if ((u[0] == 0.0) && (u[1] == 0.0)) {
    u[0] = 1.0;
    }

  v = normal.cross(u);
  float mag = u.length();

  if (mag == 0.0) {
    return;
    }

  u = (1/mag) * u;
  mag = v.length();
  v = (1/mag) * v;
  }

//*============================================================*
//*==========              pm_GcCylLineInt             ==========*
//*============================================================*
// compute the intersection of a cylinder with a line.

void
pm_GcCylLineInt (float radius, float length, PmVector3& origin, PmVector3& axis, 
               PmVector3 line[2], float tol, bool& intersect, PmVector3& ipt)
  {

  #define ndbg_pm_GcCylLineInt 
  #ifdef dbg_pm_GcCylLineInt 
  fprintf (stderr, ">>>>>> pm_GcCylLineInt: \n");
  fprintf (stderr, "   >>> line[0] (%g %g %g) \n", line[0][0], line[0][1], line[0][2] );
  fprintf (stderr, "   >>> line[1] (%g %g %g) \n", line[1][0], line[1][1], line[1][2] );
  #endif

  PmVector3 h, c, v, d, o, oc, dh, f, g, vec;
  PmVector3 cpt, cpt1, cpt2; 
  float odh, vdh, cdh, c1, c2, c3, desc, k1, k2, t1, t2, t;
  float tc, tc1, tc2;
  bool cpt_in, spt_in;

  intersect = false;
  ipt.set(0, 0, 0);

  /*
  line[0].set(4.0, 0.0, 0.0);
  line[1].set(9.0, 0.0, 0.0);
  */

  h = axis;
  h.normalize();
  c = origin - 0.5*length*h;
  d = origin + 0.5*length*h;
  o = line[0];
  v = line[1] - line[0];
  oc = o - c;
  //fprintf (stderr, "   >>> d (%g %g %g) \n", d[0], d[1], d[2] );

  odh = o*h;
  vdh = v*h;
  cdh = c*h;

  // check intersection with caps //

  cpt_in = false; 
  tc1 = 1e6; 

  if (vdh != 0.0) {
    tc1 = h*(c - o) / vdh;
    cpt1 = o + tc1*v;
    vec = cpt1 - c;

    if (vec.length() <= radius) {
      cpt_in = true; 
      tc = tc1;
      cpt = cpt1;
      }
    
    tc2 = h*(d - o) / vdh;
    cpt2 = o + tc2*v;
    vec = cpt2 - d;

    if (vec.length() <= radius) {
      if (cpt_in) {
        if (tc2 < tc1) { 
          tc = tc2;
          cpt = cpt2;
          }
        }
      else {
        cpt_in = true; 
        tc = tc2;
        cpt = cpt2;
        }
      }
    }

  if (cpt_in) {
    ipt = cpt;
    intersect = true;
    }

  // check intersection with sides //

  spt_in = false;
  f = oc - (odh - cdh)*h;
  g = v - vdh*h;
  c3 = f*f - radius*radius;
  c2 = 2*f*g;
  c1 = g*g;
  desc = c2*c2 - 4.0*c1*c3;

  if (desc < 0.0) {
    return;
    }

  t1 = (-c2 + sqrt(desc)) / (2*c1);
  t2 = (-c2 - sqrt(desc)) / (2*c1);

  if ((t1 < 0.0) && (t2 < 0.0)) {
    return;
    }

  k1 = odh + t1*vdh - cdh;
  k2 = odh + t2*vdh - cdh;

  if ((k1 <= length) && (k1 >= 0.0)) {
    t = t1;
    spt_in = true;
    }
  else {
    t1 = -1.0;
    }

  if ((k2 <= length) && (k2 >= 0.0)) {
    t = t2;
    spt_in = true;
    }
  else {
    t2 = -1.0;
    }

  if ((t1 != -1.0) && (t2 != -1.0)) {
    if (t1 < t2)  {
      t = t1;
      }
    else {
      t = t2;
      }
    }

  if (spt_in) { 
    if (cpt_in && (tc < t)) {
      ipt = cpt;
      }
    else {
      ipt = o + t*v;
      intersect = true;
      }
    }
  }

//*============================================================*
//*==========          LineLineDist                  ==========*
//*============================================================*
// compute the distance between two lines.

void
pm_GcLineLineDist(PmVector3& p1, PmVector3& v1, PmVector3& p2, PmVector3& v2, 
                PmVector3& ipt, float& p_dist)
  {
  float u, t, dist; 
  float s1, s2;
  int i;

  float a =  v1 * v1;
  float b = -v1 * v2;
  float c = -b;
  float d = -v2 * v2;
  float e = v1*p2 - v1*p1;
  float f = v2*p2 - v2*p1;
  float det = a*d - c*b;

  if (det == 0.0) {
    p_dist = -1.0;
    return;
    }

  t = (e*d - f*b) / det;
  u = (a*f - e*c) / det;
  dist = 0.0;

  if ((t < 0.0) || (t > 1.0)) {
    p_dist = -1.0;
    return;
    }
  t = (e*d - f*b) / det;
  u = (a*f - e*c) / det;
  dist = 0.0;

  if ((t < 0.0) || (t > 1.0)) {
    p_dist = -1.0;
    return;
    }

  if ((u < 0.0) || (u > 1.0)) {
    p_dist = -1.0;
    return;
    }

  for (i = 0; i < 3; i++) {
    s1 = p1[i] + t*v1[i];
    s2 = p2[i] + u*v2[i];
    dist += (s1 - s2) * (s1 - s2);
    ipt[i] = s1;
    }

  p_dist = sqrt (dist);
  }

//*============================================================*
//*==========              pm_GcLinePointProj          ==========*
//*============================================================*
// compute the projection of a point onto a line. 

void
pm_GcLinePointProj (PmVector3 line[2], PmVector3& pt, float& dist, PmVector3& proj_pt)
  {

  float dp, mag, mag2;

  PmVector3 u, v;

  /*
  fprintf (stderr, "--------------------------\n");
  fprintf (stderr, "   >>> line[0] (%g %g %g) \n", line[0][0], line[0][1], line[0][2] );
  fprintf (stderr, "   >>> line[1] (%g %g %g) \n", line[1][0], line[1][1], line[1][2] );
  fprintf (stderr, "   >>> pt (%g %g %g) \n", pt[0], pt[1] , pt[2]); 
  */

  u = line[1] - line[0];
  mag = u.length();
  mag2 = mag * mag;
  v = pt - line[0];
  dp = u * v;

  proj_pt[0] = line[0][0] + (dp / mag2) * u[0];
  proj_pt[1] = line[0][1] + (dp / mag2) * u[1];
  proj_pt[2] = line[0][2] + (dp / mag2) * u[2];

  v = pt - proj_pt;
  dist = v.length();
  }

//*============================================================*
//*==========              pm_GcPolyNormComp           ==========*
//*============================================================*
// compute an average normal for a set of points.

void
pm_GcPolyNormComp (int n, PmVector3 *pts, PmVector3& norm)
  {

  int i, j;

  float nx, ny, nz, mag;

  nx = ny = nz = 0.0;

  for (i = 0; i < n; i++) {
    if (i == (n - 1)) {
      j = 0;
      }
    else {
      j = i + 1;
      }

    nx += pts[j][2]*pts[i][1] - pts[j][1]*pts[i][2];
    ny += pts[j][0]*pts[i][2] - pts[j][2]*pts[i][0];
    nz += pts[j][1]*pts[i][0] - pts[j][0]*pts[i][1];
    }

  mag = sqrt(nx*nx + ny*ny + nz*nz);

  if (mag != 0.0) {
    norm[0] = nx / mag;
    norm[1] = ny / mag;
    norm[2] = nz / mag;
    }
  else {
    norm[0] = 0.0;
    norm[1] = 0.0;
    norm[2] = 0.0;
    }
  }

//*============================================================*
//*==========              pm_GcPolyLineInt            ==========*
//*============================================================*
// compute the intersection of a polygon with a line.

void 
pm_GcPolyLineInt (int n, PmVector3 *pts, PmVector3& norm, PmVector3 line[2], float tol, 
                bool& intersect, PmVector3& ipt)
  {

  PmVector3 dir, v;

  float t, ndp;

  int i0, i1, i2;

  #ifdef dbg_pm_GcPolyLineInt 
  fprintf (stderr, ">>>>>> pm_GcPolyLineInt: \n");
  fprintf (stderr, "   >>> norm (%g %g %g) \n", norm[0], norm[1] , norm[2]); 
  fprintf (stderr, "   >>> line[0] (%g %g %g) \n", line[0][0], line[0][1], line[0][2] );
  fprintf (stderr, "   >>> line[1] (%g %g %g) \n", line[1][0], line[1][1], line[1][2] );
  #endif 

  intersect = false;
  ipt.set(0, 0, 0);
  dir = line[0] - line[1];
  ndp = dir * norm;

  #ifdef dbg_pm_GcPolyLineInt 
  fprintf (stderr, "   >>> ndp [%g]  \n", ndp); 
  #endif

  if (fabs(ndp) < tol) {
    return;
    }

  v = pts[0] - line[1];
  t = (v * norm) / ndp;

  if ((t + tol < 0) || (t-tol > 1.0)) {
    return;
    }

  ipt[0] = line[1][0]  +  t*dir[0];
  ipt[1] = line[1][1]  +  t*dir[1];
  ipt[2] = line[1][2]  +  t*dir[2];

  // check for triangle or a quad.

  if (n == 3) {
    pm_GcVecMaxProj (norm, i0, i1, i2);
    pm_GcPolyPtIn3 (pts, ipt, i0, i1, i2, intersect);
    }

  else if (n == 4) {
    PmVector3 tri[3];
    tri[0] = pts[0]; tri[1] = pts[1]; tri[2] = pts[2];
    pm_GcVecMaxProj (norm, i0, i1, i2);
    pm_GcPolyPtIn3 (tri, ipt, i0, i1, i2, intersect);

    if (!intersect) {
      tri[1] = pts[2]; tri[2] = pts[3];
      pm_GcPolyPtIn3 (tri, ipt, i0, i1, i2, intersect);
      }
    }

  // general polygon                   

  else {
    pm_GcPolyPointClassify (n, pts, norm, tol, ipt, intersect);
    }
  }

//*============================================================*
//*==========           pm_GcPolyFaceNormsComp         ==========*
//*============================================================*
// get face normals for a set of polygons.

void
pm_GcPolyFaceNormsComp (int num_poly, int *conn, int hgeom, int num_verts,
                       PmVector3 *verts, PmVector3 **p_fnorms)
  {
  PmVector3 pverts[40], v1, v2, norm, *fnorms;
  float nx, ny, nz, mag;
  int i, j, k;

  int n = hgeom;
  fnorms = new PmVector3[num_poly];

  for (i = 0; i < num_poly; i++) {
    if (!hgeom) {
      n = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      k = *conn;
      pverts[j][0] = verts[k][0];
      pverts[j][1] = verts[k][1];
      pverts[j][2] = verts[k][2];
      conn++;
      }

    v1 = pverts[1]   - pverts[0];
    v2 = pverts[n-1] - pverts[0];
    nx = v1[1]*v2[2] - v1[2]*v2[1];
    ny = v1[2]*v2[0] - v1[0]*v2[2];
    nz = v1[0]*v2[1] - v1[1]*v2[0];
    mag = sqrt(nx*nx + ny*ny + nz*nz);

    if (mag != 0.0) {
      //if (rev) mag = -mag;
      nx = nx / mag;
      ny = ny / mag;
      nz = nz / mag;
      }

    fnorms[i][0] = nx;
    fnorms[i][1] = ny;
    fnorms[i][2] = nz;
    }

  *p_fnorms = fnorms;
  }

//*============================================================*
//*==========          pm_GcPolyVertexNormsComp      ==========*
//*============================================================*
// get vertex normals.             

void
pm_GcPolyVertexNormsComp (int num_poly, int *conn, int hgeom, int num_verts, 
                        PmVector3 *verts, PmVector3 *fnorms, PmVector3 **p_vnorms)
  {

  PmVector3 norm, *vnorms;

  int i, j;

  int ptr, next_ptr; 

  int node, face;

  vector<int> edge_table(num_verts, -1);

  // compute face adjacency list

  int n = hgeom;
  int size = num_poly;
  vector<int> face_list(size);
  vector<int> link(size); 
  next_ptr = 0;

  for (i = 0; i < num_poly; i++) {
    if (!hgeom) {
      n = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      node = *conn;
      conn++;

      if (edge_table[node] == -1) {
        ptr = next_ptr++;

        if (ptr == size) {
          size += 1000;
          face_list.resize(size); 
          link.resize(size); 
          }

        edge_table[node] = ptr;
        face_list[ptr] = i;
        link[ptr] = -1;
        }
      else {
        ptr = next_ptr++;

        if (ptr == size) {
          size += 1000;
          face_list.resize(size); 
          link.resize(size); 
          }

        link[ptr] = edge_table[node];
        edge_table[node] = ptr;
        face_list[ptr] = i;
        }
      }
    }


  // average face normals //

  vnorms = new PmVector3[num_verts];

  for (i = 0; i < num_verts; i++) {
    ptr = edge_table[i];
    n = 0;
    norm.set(0.0, 0.0, 0.0);

    while (ptr != -1) {
      face = face_list[ptr];
      norm = norm + fnorms[face];
      ptr = link[ptr];
      n += 1;
      }

    vnorms[i][0] = norm[0] / n;
    vnorms[i][1] = norm[1] / n;
    vnorms[i][2] = norm[2] / n;
    }

  *p_vnorms = vnorms;
  }

//*============================================================*
//*==========              pm_GcSphereGen              ==========*
//*============================================================*
// generate geometry for a sphere.                          

void
pm_GcSphereGen (PmGcTri *f, int iterations, int *num)
  {

  int i,it;

  double a;

  PmVector3 pa, pb, pc;

  int nt = 0, ntold;

  PmVector3 p[6] = {PmVector3(0,0,1),  PmVector3(0,0,-1), PmVector3(-1,-1,0), 
                    PmVector3(1,-1,0), PmVector3(1,1,0),  PmVector3(-1,1,0)};

  a = 1 / sqrt(2.0);

  for (i = 0; i < 6; i++) {
    p[i][0] *= a;
    p[i][1] *= a;
    }

  f[0].p1 = p[0]; 
  f[0].p2 = p[3]; 
  f[0].p3 = p[4];

  f[1].p1 = p[0]; 
  f[1].p2 = p[4]; 
  f[1].p3 = p[5];

  f[2].p1 = p[0]; 
  f[2].p2 = p[5]; 
  f[2].p3 = p[2];

  f[3].p1 = p[0]; 
  f[3].p2 = p[2]; 
  f[3].p3 = p[3];

  f[4].p1 = p[1]; 
  f[4].p2 = p[4]; 
  f[4].p3 = p[3];

  f[5].p1 = p[1]; 
  f[5].p2 = p[5]; 
  f[5].p3 = p[4];

  f[6].p1 = p[1]; 
  f[6].p2 = p[2]; 
  f[6].p3 = p[5];

  f[7].p1 = p[1]; 
  f[7].p2 = p[3]; 
  f[7].p3 = p[2];

  nt = 8;

  if (iterations < 1) {
    *num = nt;
    return;
    }


  // bisect each edge and move to the surface of a unit sphere // 

  for (it = 0; it < iterations; it++) {
    ntold = nt;

    for (i = 0; i< ntold; i++) {
      pa[0] = (f[i].p1[0] + f[i].p2[0]) / 2;
      pa[1] = (f[i].p1[1] + f[i].p2[1]) / 2;
      pa[2] = (f[i].p1[2] + f[i].p2[2]) / 2;
      pb[0] = (f[i].p2[0] + f[i].p3[0]) / 2;
      pb[1] = (f[i].p2[1] + f[i].p3[1]) / 2;
      pb[2] = (f[i].p2[2] + f[i].p3[2]) / 2;
      pc[0] = (f[i].p3[0] + f[i].p1[0]) / 2;
      pc[1] = (f[i].p3[1] + f[i].p1[1]) / 2;
      pc[2] = (f[i].p3[2] + f[i].p1[2]) / 2;

      pa.normalize();
      pb.normalize();
      pc.normalize();

      f[nt].p1 = f[i].p1; 
      f[nt].p2 = pa; 
      f[nt].p3 = pc; 
      nt++;

      f[nt].p1 = pa; 
      f[nt].p2 = f[i].p2; 
      f[nt].p3 = pb; 
      nt++;

      f[nt].p1 = pb; 
      f[nt].p2 = f[i].p3; 
      f[nt].p3 = pc; 
      nt++;

      f[i].p1 = pa;
      f[i].p2 = pb;
      f[i].p3 = pc;
      }
    }

  *num = nt;
  }

//*============================================================*
//*==========            pm_GcVecMaxProj             ==========*
//*============================================================*
// compute the maximum projection for a vector.  

void 
pm_GcVecMaxProj (PmVector3& normal, int& i0, int& i1, int& i2)
  {

  float max_proj, nx, ny, nz;

  nx = normal[0];
  ny = normal[1];
  nz = normal[2];

  max_proj = (fabs(nx) > fabs(ny) ? fabs(nx) : fabs(ny));
  max_proj = (fabs(nz) > max_proj ? fabs(nz) : max_proj);

  if (fabs(nx) == max_proj) {
    i0 = 0;
    i1 = 1;
    i2 = 2;
    }
  else if (fabs(ny) == max_proj) {
    i0 = 1;
    i1 = 0;
    i2 = 2;
    }
  else if (fabs(nz) == max_proj) {
    i0 = 2;
    i1 = 0;
    i2 = 1;
    }
  }

//*============================================================*
//*==========            pm_GcTriGetIndexes          ==========*
//*============================================================*
// compute the indexes for a set of triangles.  

void
pm_GcTriGetIndexes (int num_tri, PmVector3 *verts, int **p_conn, int *p_num_iverts, 
                  PmVector3 **p_iverts) 
  {

  int num_verts; 
  int *conn;
  PmVector3 *iverts;
  PmVector3 max, min;
  float xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, max_dim;
  float x, y, z, tx, ty, tz;
  int size;
  GcCoordBin **hash_table, *ptr, *pptr;
  int i, j, num;
  int offset, id_count, index, found;
  float f = 0.00001;

  num_verts = num_tri * 3;

  // compute extent of vertices  //

  xmin = xmax = verts[0][0];
  ymin = ymax = verts[0][1];
  zmin = zmax = verts[0][2];

  PmExtent extent(verts[0], verts[0]);

  for (i = 1; i < num_verts; i++) {
    extent.update(verts[i]);
    }

  extent.get(min, max);
  xmin = min[0]; ymin = min[1]; zmin = min[2];
  xmax = max[0]; ymax = max[1]; zmax = max[2];
  dx = (xmax - xmin);
  dy = (ymax - ymin);
  dz = (zmax - zmin);
  max_dim = (dx > dy ? dx : dy);
  max_dim = (dz > max_dim ? dz : max_dim);

  size = num_verts / 10;
  hash_table = new GcCoordBin*[size+1];

  for (i = 0; i <= size; i++) {
    hash_table[i] = NULL;
    }

  offset = (int)sqrt(static_cast<double>(num_verts));
  id_count = 0;
  tx = dx * f;
  ty = dy * f;
  tz = dz * f;

  for (i = 0; i < num_verts; i++) {
    x = verts[i][0];
    y = verts[i][1];
    z = verts[i][2];
    index = int(offset*((x-xmin)/dx + offset*((y-ymin)/dy + offset*((z-zmin)/dz))));
    index = index % size;
    ptr = hash_table[index];
    found = 0;

    while (ptr) {
      if ((fabs(ptr->x - x) <= tx) &&
          (fabs(ptr->y - y) <= ty) &&
          (fabs(ptr->z - z) <= tz)) {
        found = 1;
        break;
        }

      ptr = ptr->link;
      }

    if (!found) {
      ptr = new GcCoordBin;
      ptr->x = x;
      ptr->y = y;
      ptr->z = z;
      ptr->id = id_count++;
      ptr->link = hash_table[index];
      hash_table[index] = ptr;
      }
    }

  //fprintf (stderr, "\n>>>>>> id_count [%d] \n", id_count);
  iverts = new PmVector3[id_count];
  conn = new int[4*num_tri];
  num = 0;
  num_verts = 0;

  for (j = 0; j < num_tri; j++) {
    conn[num++] = 3;

    for (i = 0; i < 3; i++) {
      x = verts[num_verts][0];
      y = verts[num_verts][1];
      z = verts[num_verts][2];
      num_verts += 1;
      index = int(offset*((x-xmin)/dx + offset*((y-ymin)/dy + offset*((z-zmin)/dz))));
      index = index % size;
      ptr = hash_table[index];
      found = 0;

      while (ptr) {
        if ((fabs(ptr->x - x) <= tx) &&
            (fabs(ptr->y - y) <= ty) &&
            (fabs(ptr->z - z) <= tz)) {
          iverts[ptr->id][0] = x;
          iverts[ptr->id][1] = y;
          iverts[ptr->id][2] = z;
          conn[num++] = ptr->id;
          break;
          }

        pptr = ptr;
        ptr = ptr->link;
        }
      }
    }

  // NOTE: delete ptrs //
  delete[] hash_table;
  *p_conn = conn;
  *p_num_iverts = id_count;
  *p_iverts = iverts;
  }

//*============================================================*
//*==========            pm_GcPolyPtIn3              ==========*
//*============================================================*
// determine if a point is within a triangle.   

void
pm_GcPolyPtIn3 (PmVector3 *tri, PmVector3& pt, int i0, int i1, int i2, bool& in)
  {

  float u0, u1, u2, v0, v1, v2, a, b;

  u0 = pt[i1] - tri[0][i1];
  v0 = pt[i2] - tri[0][i2];

  u1 = tri[1][i1] - tri[0][i1];
  v1 = tri[1][i2] - tri[0][i2];

  u2 = tri[2][i1] - tri[0][i1];
  v2 = tri[2][i2] - tri[0][i2];

  in = false;

  if (u1 == 0.0) {
    b = u0 / u2;

    if ((b >= 0.0) && (b <= 1.0)) {
      a = (v0 - b*v2) / v1;
      in = ((a >= 0.0) && ((a + b) <= 1.00001));
      }
    }

  else {
    b = (v0*u1 - u0*v1) / (v2*u1 - u2*v1);

    if ((b >= 0.0) && (b <= 1.00001)) {
      a = (u0 - b*u2) / u1;
      in = ((a >= 0.0) && ((a + b) <= 1.00001));
      }
    }
  }

//*============================================================*
//*==========            pm_GcPolyPointClassify      ==========*
//*============================================================*
// classify a point in a polygon.

void
pm_GcPolyPointClassify (int n, PmVector3 *pts, PmVector3& norm, float tol, PmVector3& pt, 
                      bool& res)
  {

  float wnum;

  pm_GcPolyWindingComp (n, pts, norm, pt, wnum);

  res = (wnum > 2.0*M_PI);
  }

//*============================================================*
//*==========            pm_GcPolyWindingComp        ==========*
//*============================================================*
// compute the winding number of a polygon.

void
pm_GcPolyWindingComp (int num_pts, PmVector3 *pts, PmVector3& normal, PmVector3& pt, 
                    float& wnum)
  {

  float vec1[2], vec2[2];

  float angle, len, dp;

  int j, i0, i1, i2;

  pm_GcVecMaxProj (normal, i0, i1, i2);
  vec1[0] = pts[num_pts-1][i1] - pt[i1];
  vec1[1] = pts[num_pts-1][i2] - pt[i2];
  len = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1]);

  if (len <= 0.0) {
    wnum = 1;
    return;
    }

  vec1[0] /= len;
  vec1[1] /= len;

  // sum the angles and see if answer mod 2*PI > PI
  angle = 0.0 ;

  for (j = 0; j < num_pts; j++) {
    vec2[0] = pts[j][i1] - pt[i1];
    vec2[1] = pts[j][i2] - pt[i2];
    len = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1]);
    vec2[0] /= len;
    vec2[1] /= len;
    dp = vec1[0]*vec2[0] + vec1[1]*vec2[1];

    if (dp < -1.0) {
      angle += M_PI;
      }

    else if (dp < 1.0) {
      if (vec1[0]*vec2[1] - vec2[0]*vec1[1] >= 0.0) {
        angle += acos (dp);
        }
      else {
        angle -= acos (dp);
        }
      }

    vec1[0] = vec2[0];
    vec1[1] = vec2[1];
    }

  wnum = fmod (fabs(angle) + M_PI, 4.0*M_PI);
  }


//*============================================================*
//*==========            pm_GcTriInscribedCircle     ==========*
//*============================================================*
// compute the center and radius of the incribed circle for a 
// triangle. 

void
pm_GcTriInscribedCircle(PmVector3 pts[3], PmVector3& center, float& radius, float& area) 
  {
  PmVector3 v1, v2, v3;
  float a, b, c, s, p;

  v1 = pts[1] - pts[0];
  v2 = pts[2] - pts[1];
  v3 = pts[0] - pts[2];

  a = v1.length();
  b = v2.length();
  c = v3.length();
  p = (a + b + c);
  s = p / 2.0;

  area = sqrt(s*(s-a)*(s-b)*(s-c));
  radius = area / s;

  /*
  for (int i = 0; i < 3; i++) {
    center[i] = (a*pts[0][i] + b*pts[1][i] + c*pts[2][i]) / p;
    }
  */

  center = (pts[0] + pts[1] + pts[2]) / 3.0;

  }


}

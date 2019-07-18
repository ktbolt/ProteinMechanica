
/*------------------------------------------------------------*
 *                                                            *
 *                      ****  poly  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

#include "gr.h"

typedef struct GeomPolygonEdge {
  int count;
  int poly;
  int node1, node2;
  struct GeomPolygonEdge *next;
  } GeomPolygonEdge;


int g_plist[50000];

void
PolyClosed (DmObj *geom);

void
GeomIndexPolygonsReOrient (int p, int num_polys, int *conn, int hgeom, int n, 
                           int *poly_table, GeomPolygonEdge **edge_table);

void 
GeomIndexPolygonsOrientTrav (int p, int num_polys, int *conn, 
                             GeomPolygonEdge **edge_table,
                             int hgeom, int n, int *poly_table, int poly_proc[]);


/*------------------------------------------------------------*
 *                                                            *
 *                      ****  PolyReorient  ****              *
 *                                                            *
 *------------------------------------------------------------*/

void 
PolyReorient (DmObj *geom)
  {

  int num_verts;

  DmPoint3 *verts;

  int num_polys;
 
  int *conn, *pconn; 

  int n, ne;

  int i, j;

  int n1, n2, first_node, min_n, nlist[10];

  int nn1, nn2;

  GeomPolygonEdge **edge_table, *ptr, *nptr;

  int size;

  int hgeom;

  int *poly_table;

  int *poly_proc;

  int p, pn, pcount;

  static char *fn = "GeomIndexPolygonsOrient";

 /**************
  ***  body  ***
  **************/

  PolyClosed (geom);

  fprintf (stderr, "\n---------- PolyReorient ---------- \n");

  gr_GeomIndexPolygonsGet (geom, &num_verts, &verts, &num_polys, &conn);

  hgeom = DM_TRUE;
  n = 3;

  edge_table = dm_MemAlloc(0, GeomPolygonEdge*, num_verts);

  for (i = 0; i < num_verts; i++) {
    edge_table[i] = DmNil(GeomPolygonEdge*);
    }

  for (i = 0; i < num_polys; i++) {
    g_plist[i] = -1;
    }

  poly_proc = dm_MemAlloc(0, int, num_polys);

  for (i = 0; i < num_polys; i++) {
    poly_proc[i] = 0;
    }


  /*  compute poly edge list.  */

  if (!hgeom) { 
    poly_table = dm_MemAlloc(0, int, num_polys);
    }

  pconn = conn;

  for (i = 0; i < num_polys; i++) {
    if (!hgeom) { 
      n = *pconn;
      pconn++;

      if (i == 0) {
        poly_table[i] = 0;
        }
      else {
        poly_table[i] = poly_table[i-1] + n + 1;
        }
      }

    first_node = *pconn;

    for (j = 0; j < n; j++) {
      n1 = *pconn;
      pconn++;

      if (j == n - 1) {
        n2 = first_node;
        }
      else {
        n2 = *pconn;
        }

      min_n = dm_MinVal (n1, n2); 
      ptr = dm_MemAlloc(0, GeomPolygonEdge, 1);
      ptr->poly = i;
      ptr->node1 = n1;
      ptr->node2 = n2;
      ptr->next = edge_table[min_n]; 
      edge_table[min_n] = ptr; 
      }
    }


  for (i = 0; i < num_verts; i++) {
    ptr = edge_table[i];
    ne = 0;

    while (ptr) { 
      ptr = ptr->next;
      ne += 1;
      }

    if (ne == 1) {
      fprintf (stderr, "\n\n ***** poly not closed ***** \n"); 
      break;
      } 
    }

  p = 4873;
  p = 11451;
  poly_proc[p] = 1;

  GeomIndexPolygonsOrientTrav (p, num_polys, conn, edge_table, hgeom, n, poly_table,
                               poly_proc);


#ifdef dbg_GeomIndexPolygonsOrient 
  fprintf (stderr, "\n---------- new conn -----\n", i);

  for (i = 0; i < num_polys; i++) {
    if (!hgeom) { 
      n = *conn;
      conn++;
      }

    fprintf (stderr, "%d: ", i);

    for (j = 0; j < n; j++) {
      fprintf (stderr, "%d ", *conn);
      conn++;
      }

    fprintf (stderr, "\n");
    }
#endif

  for (i = 0; i < num_polys; i++) {
    if (poly_proc[i] == 0) {
      /*
    if (g_plist[i] == -1) {
      fprintf (stderr, "%d: %d ", i, poly_proc[i]);
      */
      fprintf (stderr, "%d ", i);
      }
    }

  dm_MemFree (0, edge_table);
  }


/*------------------------------------------------------------*
 *                                                            *
 *          ****  GeomIndexPolygonsOrientTrav  ****           *
 *                                                            *
 *------------------------------------------------------------*/

void 
GeomIndexPolygonsOrientTrav (int p, int num_polys, int *conn, 
                             GeomPolygonEdge **edge_table,
                             int hgeom, int n, int *poly_table, int poly_proc[])
  {

  GeomPolygonEdge *ptr, *nptr;

  int *pconn;

  int i, j;

  int n1, n2, first_node, max_n, min_n, nlist[10];

  int nn1, nn2, max_nn;

  int pn, pcount, plist[10];

  int reorient, pr;

  int olist[20];

 /**************
  ***  body  ***
  **************/

  g_plist[p] = 1;
  /*
  poly_proc[p] = 1;
  */
  
  pr = 0;
  if (p == 4070) pr = 1;

  if (pr) fprintf (stderr, "------- GeomIndexPolygonsOrientTrav: poly [%d] ------- \n", p);

  if (!hgeom) {
    i = poly_table[p];
    n = conn[i];
    pconn = conn+i+1;
    }
  else {
    i = p*n;
    pconn = conn+i;
    }

  first_node = *pconn;
  pcount = 0;

  for (j = 0; j < n; j++) {
    n1 = *pconn;
    pconn++;

    if (j == n - 1) {
      n2 = first_node;
      }
    else {
      n2 = *pconn;
      }

    min_n = dm_MinVal (n1, n2); 
    max_n = dm_MaxVal (n1, n2); 

    if (pr) fprintf (stderr, ">>>>>> search edge  %d(%d %d): ", p, n1, n2);

    for (ptr = edge_table[min_n]; ptr; ptr = ptr->next) {
      pn = ptr->poly;
      /*
      if (pn == p) {
      */


      if ((pn == p) || poly_proc[pn]) {
        continue;
        }

      nn1 = ptr->node1;
      nn2 = ptr->node2;
      max_nn = dm_MaxVal (nn1, nn2); 

      if (max_n == max_nn) {
        reorient = 0;

        /*
        if (pn == 4070) {
          fprintf (stderr, "******************************************\n"); 
          fprintf (stderr, "p [%d]  n1 [%d]  n2 [%d] \n", p, n1, n2);
          fprintf (stderr, "pn [%d]  nn1 [%d]  nn2 [%d] \n", pn, nn1, nn2);
          fprintf (stderr, "******************************************\n"); 
          }
        */

        if (nn1 == n1) {
          /*
          GeomIndexPolygonsReOrient (pn, num_polys, conn, hgeom, n, poly_table, 
                                        edge_table);
          */
          reorient = 1;
          }

        if (pr) fprintf (stderr, " >>> match %d(%d %d) %d <<< ", pn, nn1, nn2, reorient);

        /*
        if (!poly_proc[pn]) { 
          plist[pcount++] = pn;
          }

        poly_proc[pn] = 1;
        */ 
        poly_proc[pn] = 1;
        olist[pcount] = reorient;
        plist[pcount++] = pn;
        }
      }

    if (pr) fprintf (stderr, "\n"); 
    }

  /*
  fprintf (stderr, " >>> pcount %d \n ", pcount); 
  */

  if (pcount > 3) {
    fprintf (stderr, " >>> p [%d]  pcount %d: ", p, pcount); 

    for (i = 0; i < pcount; i++) {
      pn = plist[i];
      fprintf (stderr, " %d ", pn); 
      }

    fprintf (stderr, "\n"); 
    }

  for (i = 0; i < pcount; i++) {
    p = plist[i];

    if (olist[i]) {
      /*
      GeomIndexPolygonsReOrient (p, num_polys, conn, hgeom, n, poly_table, edge_table);
      */
      GeomIndexPolygonsReOrient (p, num_polys, conn, hgeom, n, poly_table, edge_table);
      }

    GeomIndexPolygonsOrientTrav (p, num_polys, conn, edge_table, hgeom, n, poly_table,
                                    poly_proc);
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *          ****  GeomIndexPolygonsReOrient  ****             *
 *                                                            *
 *------------------------------------------------------------*/

void
GeomIndexPolygonsReOrient (int p, int num_polys, int *conn, int hgeom, int n, 
                           int *poly_table, GeomPolygonEdge **edge_table)
  {

  int *pconn;

  int i, j;

  int *ptr, nlist[10];

  GeomPolygonEdge *eptr;

  int first_node, n1, n2, min_n, nt;

 /**************
  ***  body  ***
  **************/

  if (!hgeom) {
    i = poly_table[p];
    n = conn[i];
    pconn = conn+i+1;
    }
  else {
    i = p*n;
    pconn = conn+i;
    }


  /*  change connectivity.  */

  ptr = pconn;

  for (i = 0; i < n; i++) {
    nlist[i] = *ptr;
    ptr++;
    }

  ptr = pconn;

  for (i = 0; i < n; i++) {
    *ptr = nlist[n - i - 1];
    ptr++;
    }


  /*  update edge table with new edge orientation.  */

  first_node = *pconn;

  for (j = 0; j < n; j++) {
    n1 = *pconn;
    pconn++;

    if (j == n - 1) {
      n2 = first_node;
      }
    else {
      n2 = *pconn;
      }

    min_n = dm_MinVal (n1, n2); 
    eptr = edge_table[min_n];

    while (eptr) {
      if (eptr->poly == p) {
        eptr->node1 = n2; 
        eptr->node2 = n1; 
        }

      eptr = eptr->next; 
      }
    }
  }



/*------------------------------------------------------------*
 *                                                            *
 *                      ****  PolyClosed  ****                *
 *                                                            *
 *------------------------------------------------------------*/

void 
PolyClosed (DmObj *geom)
  {

  int num_verts;

  DmPoint3 *verts;

  int num_polys;
 
  int *conn, *pconn; 

  int n, ne;

  int i, j;

  int n1, n2, first_node, min_n, nlist[10], max_n;

  int nn1, nn2;

  int min_nn, max_nn;

  GeomPolygonEdge **edge_table, *ptr, *nptr;

  int size;

  int hgeom;

  int p, pn, pcount;

  static char *fn = "GeomIndexPolygonsOrient";

 /**************
  ***  body  ***
  **************/

  fprintf (stderr, "\n---------- PolyClosed ---------- \n");

  gr_GeomIndexPolygonsGet (geom, &num_verts, &verts, &num_polys, &conn);

  hgeom = DM_TRUE;
  n = 3;

  edge_table = dm_MemAlloc(0, GeomPolygonEdge*, num_verts);

  for (i = 0; i < num_verts; i++) {
    edge_table[i] = DmNil(GeomPolygonEdge*);
    }


  /*  compute edge table  */

  pconn = conn;

  for (i = 0; i < num_polys; i++) {
    if (!hgeom) { 
      n = *pconn;
      pconn++;
      }

    first_node = *pconn;

    for (j = 0; j < n; j++) {
      n1 = *pconn;
      pconn++;

      if (j == n - 1) {
        n2 = first_node;
        }
      else {
        n2 = *pconn;
        }

      min_n = dm_MinVal (n1, n2); 
      max_n = dm_MaxVal (n1, n2); 
      ptr = edge_table[min_n];

      while (ptr) { 
        if ((ptr->node1 == min_n) && (ptr->node2 == max_n)) {
          ptr->count++;
          break;       
          } 
        ptr = ptr->next;
        }

      if (!ptr) { 
        ptr = dm_MemAlloc(0, GeomPolygonEdge, 1);
        ptr->count = 1;
        ptr->poly = i;
        ptr->node1 = min_n;
        ptr->node2 = max_n;
        ptr->next = edge_table[min_n]; 
        edge_table[min_n] = ptr; 
        }
      }
    }


  /*  check for closed poly  */

  for (i = 0; i < num_verts; i++) {
    ptr = edge_table[i];

    while (ptr) { 
      if (ptr->count != 2) {
        break;
        }

      ptr = ptr->next;
      }

    if (ptr) { 
      fprintf (stderr, "\n\n ***** poly not closed ***** \n"); 
      ptr = edge_table[i];
      fprintf (stderr, "\n >>>> pos [%d]: ", i); 

      while (ptr) { 
        p = ptr->poly;
        fprintf (stderr, "%d(%d %d)%d ", p, ptr->node1, ptr->node2, ptr->count); 
        ptr = ptr->next;
        }
      } 
    }

  }




/*------------------------------------------------------------*
 *                                                            *
 *                      ****  poly  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

#include "sfc2pm.h" 

int *g_plist;

/*------------------------------------------------------------*
 *                                                            *
 *                      ****  PolyReorient  ****              *
 *                                                            *
 *------------------------------------------------------------*/

void 
PolyReorient (int num_polys, int *conn, int num_verts, DmPoint3 *verts)
  {
  int *pconn; 
  int n, ne;
  int i, j;
  int n1, n2, first_node, min_n, nlist[100];
  int nn1, nn2;
  GeomPolygonEdge **edge_table, *ptr, *nptr;
  int size;
  int hgeom;
  int *poly_table;
  int *poly_proc;
  int p, pn, pcount;

  hgeom = 1;
  n = 3;


  /* build edge table */

  edge_table = (GeomPolygonEdge**)malloc(sizeof(GeomPolygonEdge*)* num_verts);

  for (i = 0; i < num_verts; i++) {
    edge_table[i] = NULL; 
    }

  poly_proc = (int*)malloc(sizeof(int)* num_polys);
  g_plist = (int*)malloc(sizeof(int)* num_polys);

  for (i = 0; i < num_polys; i++) {
    poly_proc[i] = 0;
    g_plist[i] = -1;
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
      ptr = (GeomPolygonEdge*)malloc(sizeof(GeomPolygonEdge));
      ptr->poly = i;
      ptr->node1 = n1;
      ptr->node2 = n2;
      ptr->next = edge_table[min_n]; 
      edge_table[min_n] = ptr; 
      }
    }


  /* check for closed surface */

  for (i = 0; i < num_verts; i++) {
    ptr = edge_table[i];
    ne = 0;

    while (ptr) { 
      ptr = ptr->next;
      ne += 1;
      }

    if (ne == 1) {
      fprintf (stderr, "\n\n ***** surface not closed ***** \n"); 
      break;
      } 
    }

  p = 1000;
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
  int n1, n2, first_node, max_n, min_n; 
  int nn1, nn2, max_nn;
  int pn, pcount, plist[10];
  int reorient;
  int olist[4];
  int pr = 0;
  static int level = 1;

  g_plist[p] = 1;
  /*
  poly_proc[p] = 1;
  */
  
  level += 1;
  //fprintf (stderr, " level %d \n", level);
  //fprintf (stderr, " %d ", p);

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

    for (ptr = edge_table[min_n]; ptr; ptr = ptr->next) {
      pn = ptr->poly;

      if ((pn == p) || poly_proc[pn]) {
        continue;
        }

      nn1 = ptr->node1;
      nn2 = ptr->node2;
      max_nn = dm_MaxVal (nn1, nn2); 

      if (max_n == max_nn) {
        reorient = 0;

        if (nn1 == n1) {
          reorient = 1;
          }

        // fprintf (stderr, " >>> match %d(%d %d) %d <<< ", pn, nn1, nn2, reorient);

        poly_proc[pn] = 1;
        olist[pcount] = reorient;
        plist[pcount++] = pn;
        //fprintf (stderr, "(%d) ", pcount); 
        }
      }
    }

  if (pcount > 3) {
    fprintf (stderr, " >>> pcount > 3:  p = %d  pcount = %d: ", p, pcount); 

    for (i = 0; i < pcount; i++) {
      pn = plist[i];
      fprintf (stderr, " %d ", pn); 
      }

    fprintf (stderr, "\n"); 
    }

  for (i = 0; i < pcount; i++) {
    p = plist[i];

    if (olist[i]) {
      GeomIndexPolygonsReOrient (p, num_polys, conn, hgeom, n, poly_table, edge_table);
      }

    GeomIndexPolygonsOrientTrav (p, num_polys, conn, edge_table, hgeom, n, poly_table,
                                    poly_proc);
    }

  level -= 1;
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
  int *ptr, nlist[100];
  GeomPolygonEdge *eptr;
  int first_node, n1, n2, min_n, nt;

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


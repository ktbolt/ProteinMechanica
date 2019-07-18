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
//* win_glut:           g l u t   w in d o w                   *
//*============================================================*

#include "win_glut.h"
#include "jpeg.h"

namespace PmGraphics {

static int glut_init = false;

static bool g_has_init_cmd = false;
static char g_init_cmd[1000];

//#define  dbg


//*============================================================*
//*==========              createSysWindow           ==========*
//*============================================================*
// create a system specific window.

void
GrWindow::sys_win_createWindow ()
  {
  #ifdef dbg
  fprintf (stderr, "\n>>>>>> GrWindow::create: glut window [%s]  \n", name.c_str());
  #endif
  WinPrvData* win_data = (WinPrvData*)prv_data;
  char *argv[2], str[80];

  GrGlutWin *glut_win_data = new GrGlutWin;
  win_data->sys_win_data = glut_win_data;
  memset (str, 0, 80);
  name.copy(str, 60);

  if (!glut_init) {
    int argc = 1;
    argv[0] = str;
    glutInit (&argc, argv);
    glut_init = true;
    }

  glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize (width, height);
  glutInitWindowPosition (x, y);

  // open a window //

  int wid = glutCreateWindow (str);
  //fprintf (stderr, "  >>>> glut window id [%d] \n", wid);
  glut_win_data->glut_win = wid;
  GrGlutWinMap *win_map = new GrGlutWinMap;
  win_map->wid = wid;
  win_map->win = this;
  win_table.push_back (win_map);

  // register the function to do all our opengl drawing //

  glutDisplayFunc (&glut_WinExpose);

  // set idle func //

  glutIdleFunc (&glut_WinIdle);

  // keyboard proc func //

  glutKeyboardFunc (&glut_WinKeyProc);

  // special proc func//

  glutSpecialFunc (&glut_WinSpecialProc);

  // mouse button proc func //

  glutMouseFunc (&glut_EventButp);

  // mouse motion proc func //

  glutMotionFunc (&glut_EventMotion);

  // set reshape func //

  glutReshapeFunc (&glut_WinReshape);

  glutSetWindow (wid);
  }

//*============================================================*
//*==========          sys_win_processEvents         ==========*
//*============================================================*
// process events.                   

void
GrWindow::sys_win_processEvents(char *cmd) 
  {
  //fprintf (stderr, ">>>>> GrWindow::sys_win_processEvents  cmd = %x \n", cmd); 

  if (cmd) {
    glut_WinExpose ();
    //grSystem.processCommand (cmd);
    g_has_init_cmd = true;
    strcpy (g_init_cmd, cmd);
    }

  glutPostRedisplay();
  glutMainLoop();
  }

void
GrWindow::sys_win_init_display()
  {
  glutPostRedisplay();
  glutMainLoop();
  }

//*============================================================*
//*==========       sys_win_setCurrentWindow         ==========*
//*============================================================*
// set the cutrrent window.

void
GrWindow::sys_win_setCurrentWindow() 
  {
  WinPrvData* win_data = (WinPrvData*)prv_data;
  GrGlutWin *glut_win_data = (GrGlutWin *)win_data->sys_win_data;
  //fprintf (stderr, ">>>>> GrWindow::setCurrentWindow [%d] \n", glut_win_data->glut_win);
  glutSetWindow (glut_win_data->glut_win);
  g_curr_win = glut_win_data->glut_win;
  }

//*============================================================*
//*==========        sys_win_getGlutWin              ==========*
//*============================================================*
// get the window object for the current glut window.

void
sys_win_getGlutWin (GrWindow **win) 
  {

  int wid = glutGetWindow();
  *win = NULL;

  for (int i = 0; i < win_table.size(); i++) {
    if (win_table[i]->wid  == wid) {
      *win = win_table[i]->win;
      return;
      }
    }
  }

//*============================================================*
//*==========          sys_win_swapBuffers           ==========*
//*============================================================*
// swap window buffers.                                     

void
GrWindow::sys_win_swapBuffers () {
  glutSwapBuffers();
  }

//*============================================================*
//*==========          sys_win_writeImage            ==========*
//*============================================================*
// write the window image to a file.

void 
GrWindow::sys_win_writeImage(const string name, GrImageFormat format, int *status)
  {
  char fname[100];

  FILE *fp;

  unsigned char rgb[3], *pixels;

  int width, height;

  GLenum type;

  float frgb[3];

  fprintf (stderr, "    >>> write window image \n");
  *status = 1;
  width = glutGet (GLUT_WINDOW_WIDTH);
  height = glutGet (GLUT_WINDOW_HEIGHT);
  fprintf (stderr, "    >>> width [%d] \n", width);
  fprintf (stderr, "    >>> height [%d] \n", height);

  if (format == GR_IMAGE_FORMAT_PPM) {
    sprintf (fname, "%s.ppm", name.c_str());

    if ((fp = fopen(fname, "w")) == NULL) {
      *status = 0;
      return;
      }

    fprintf (fp, "P6\n" );
    fprintf (fp, "# %s.ppm from protein modeler \n", name.c_str());
    fprintf (fp, "%d %d\n", width, height );
    fprintf (fp, "255\n" );

    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        glReadPixels (i, j, 1, 1, GL_RGB, GL_FLOAT, frgb);
        rgb[0] = (unsigned char)(255.0 * frgb[0]);
        rgb[1] = (unsigned char)(255 * frgb[1]);
        rgb[2] = (unsigned char)(255 * frgb[2]);
        fwrite (rgb, sizeof(unsigned char), 3, fp);
        }
      }

    fclose (fp);
    }

  else if (format == GR_IMAGE_FORMAT_JPEG) {
    sprintf (fname, "%s.jpg", name.c_str());
    pixels = new unsigned char[3*width];
    gr_JpegInit();
    gr_JpegStartCompress(fname, width, height);

    for (int i = height-1; i >= 0; i--) {
      glReadPixels(0, i, width, 1, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      gr_JpegWriteScanline(i, width, pixels);
      }

    gr_JpegFinishCompress();
    delete[] pixels;
    }
  }

//////////////////////////////////////////////////////////////
//                    c a l l b a c k s                    //
////////////////////////////////////////////////////////////

//*============================================================*
//*==========              glut_WinExpose            ==========*
//*============================================================*

void
glut_WinExpose ()
  {
  //fprintf (stderr, ">>>>> glut_WinExpose:   g_has_init_cmd = %d \n", g_has_init_cmd);
  GrWindow *win;

  // get current glut window //

  sys_win_getGlutWin (&win);

  /*
  fprintf (stderr, "\n >>>>>> GlutWinExpose:  win [%x]  g_curr_win [%d] \n", 
           win, g_curr_win);
  */

  if (!win) {
    return;
    }

  win->redraw();

  if (g_curr_win != -1) {
    glutSetWindow (g_curr_win);
    }

  if (g_has_init_cmd) {
    grSystem.processCommand (g_init_cmd);
    g_has_init_cmd = false;
    }
  }

//*============================================================*
//*==========              glut_WinIdle              ==========*
//*============================================================*

void
glut_WinIdle() 
  {
  GrWindow *win;
  sys_win_getGlutWin (&win);

  void (*func)(); 
  win->getIdleFunction(&func);

  if (func) {
    (*func)(); 
    }
  }

//*============================================================*
//*==========              glut_WinReshape           ==========*
//*============================================================*

void
glut_WinReshape (int width, int height)
  {

  GrWindow *win;
  sys_win_getGlutWin (&win);

  if (!win) {
    return;
    }

  win->setDimensions (width, height);
  }

//*============================================================*
//*==========              glut_WinKeyProc           ==========*
//*============================================================*

void
glut_WinKeyProc (unsigned char key, int x, int y)
  {

  int cmd_input; 
  GrCmdPrompt prompt;
  GrWindow *win; 
  static char cmd[1000];
  static int n = 0;

  //fprintf (stderr, "\n------ gr_GlutWinKeyProc ---------- \n"); 

  sys_win_getGlutWin (&win);

  if (!win) {
    fprintf (stderr, "**** WARNING [glut_WinKeyProc]: no window. \n"); 
    return;
    }

  if (key == KEY_ESCAPE) {
    cmd[0] = key;
    cmd[1] = '\0';

    //fprintf (stderr, "\n>>>>>> gr_GlutWinKeyProc: cmd [%s] \n", cmd); 
    grSystem.processCommand (cmd);

    win->getCmdInfo (&cmd_input, prompt);
    fprintf (stdout, "%s ", prompt); 
    fflush (stdout);

    //exit (0);
    return;
    }

#ifdef dbg_gr_GlutWinKeyProc 
  fprintf (stderr, "key [%d] \n", key);
#endif

  if (key == KEY_RETURN) {
    win->getCmdInfo (&cmd_input, prompt);
    fprintf (stdout, "\n"); 

    if (cmd_input) {
      cmd[n] = '\0';
      grSystem.processCommand (cmd);
      fprintf (stdout, "%s ", prompt); 
      }

    n = 0;
    }

  else if (key == KEY_BACKSPACE) {
    fprintf (stdout, "\b \b");
    n -= 1;
    if (n < 0) n = 0;
    }

  else { 
    fprintf (stdout, "%c", key);
    cmd[n++] = key;
    }

  fflush (stdout);
  }

//*============================================================*
//*==========              glut_WinSpecialProc       ==========*
//*============================================================*

void
glut_WinSpecialProc (int key, int x, int y)
  {

  GrWindow *win;

  int i;

  static char cmd[100];

  static char *key_name[] = {"left", "up", "right", "down"};

#define ndbg_gr_GlutWinSpecialProc
#ifdef dbg_gr_GlutWinSpecialProc
  fprintf (stderr, "\n------ glut_WinSpecialProc ------- \n");
#endif

  sys_win_getGlutWin (&win);

  //usleep (100);

#ifdef dbg_gr_GlutWinSpecialProc
  fprintf (stderr, "key [%d] \n", key);
#endif

  if ((key >= KEY_ARROW_LEFT) && (key <= KEY_ARROW_DOWN)) {
    i = key - KEY_ARROW_LEFT;
    sprintf (cmd, "keyboard arrow_%s",  key_name[i]);
    //fprintf (stderr, ">>> key=%d %s \n", i, key_name[i]);
    grSystem.processCommand (cmd);
    }
  }

//*============================================================*
//*==========              glut_EventButr            ==========*
//*============================================================*
// process mouse button release event.                      

static void
glut_EventButr (int mx, int my, GrEventButtonId button_id)
  {
  //fprintf (stderr, "---------- gr_GlutEventButr ---------- \n");
  GrWindow *win;
  sys_win_getGlutWin (&win);
  win->processButtonReleaseEvent (mx, my, button_id);
  }

//*============================================================*
//*==========              glut_EventButp            ==========*
//*============================================================*
// process mouse button press event.                      

void
glut_EventButp (int button, int state, int x, int y)
  {

  GrWindow *win;

  int mx, my, key_mod, shift, ctrl; 

  GrEventButtonId button_id;

#define ndbg_glut_EventButp 
#ifdef dbg_glut_EventButp 
  fprintf (stderr, "\n --------- glut_EventButP --------- \n"); 
  fprintf (stderr, " >>>>>> button [%d] \n", button); 
  fprintf (stderr, " >>>>>> state  [%d] \n", state); 
#endif

  sys_win_getGlutWin (&win);


  // determine key modifiers 

  key_mod = glutGetModifiers();
  shift = (key_mod == GLUT_ACTIVE_SHIFT); 
  ctrl = (key_mod == GLUT_ACTIVE_CTRL); 

  switch (button) {
    case GLUT_LEFT_BUTTON:  
      if (shift) 
        button_id = GR_BUTTON_SHIFT_LEFT; 
      else if (ctrl) 
        button_id = GR_BUTTON_CTRL_LEFT; 
      else 
        button_id = GR_BUTTON_LEFT; 
    break;

    case GLUT_MIDDLE_BUTTON: 
      if (shift) 
        button_id = GR_BUTTON_SHIFT_MIDDLE;
      else if (ctrl) 
        button_id = GR_BUTTON_CTRL_MIDDLE;
      else
        button_id = GR_BUTTON_MIDDLE; 
    break;

    case GLUT_RIGHT_BUTTON:
      if (shift)  
        button_id = GR_BUTTON_SHIFT_RIGHT; 
      else if (ctrl) 
        button_id = GR_BUTTON_CTRL_RIGHT; 
      else
        button_id = GR_BUTTON_RIGHT; 
    break;
    }

  mx = x;
  my = y;

  if (state == GLUT_UP) {
    glut_EventButr (x, y, button_id);
    }
  else {
    win->processButtonPressEvent (mx, my, button_id);
    }
  }

//*============================================================*
//*==========              glut_EventMotion          ==========*
//*============================================================*

void
glut_EventMotion (int x, int y)
  {
  //fprintf (stderr, "\n --------- gr_GlutEventMotion --------- \n");
  GrWindow *win;
  sys_win_getGlutWin (&win);
  win->processMotionEvent (x, y);
  }


}  // PmGraphics




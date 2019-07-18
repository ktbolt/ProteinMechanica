
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
//* win:                 w i n d o w                           *
//*============================================================*

#include "pm/gr/win.h"
#include "win_prv.h"
#include "pm/gr/scene.h"

namespace PmGraphics {

vector<class GrWindow*> GrWindow::win_objs(10);

static bool debug = 0;


//*============================================================*
//*==========      constructors / detructor          ==========*
//*============================================================*

GrWindow::GrWindow(const string name, int x, int y, int width, int height, int type) :
                   x(x), y(y), width(width), height(height), type(type)
  {
  this->name = name;
  init();
  if (debug) fprintf (stderr, ">>>>>> GrWindow full ctor: name [%s] \n", this->name.c_str());
  }


// initialize window data.

void
GrWindow::init() 
  {
  WinPrvData *win_data = new WinPrvData;
  prv_data = win_data; 

  GrVector3 min(0.0, 0.0, 0.0);
  GrVector3 max(1.0, 1.0, 1.0);
  
  extent.set(min, max);
  orig_extent.set(min, max); 
  scene = NULL;
  idle_func = NULL;
  color.set(0,0,0);

  // init private data

  win_data->do_pick = false;
  win_data->do_keyp = false;
  win_data->ren_init = false;

  win_data->record.enabled = false;
  win_data->record.params.format = GR_IMAGE_FORMAT_UNKNOWN;
  win_data->record.count = 0;

  win_data->event.mx = 0;
  win_data->event.my = 0;
  win_data->event.mouse_pos_set = false;

  win_data->event.type = GR_EVENT_NONE;
  win_data->event.button_state = GR_EVENT_BUTTON_STATE_NONE;
  win_data->event.button_id = GR_BUTTON_UNKNOWN;
  win_data->event.action = GR_WINDOW_ACTION_NONE;
  win_data->event.last = new GrWindowEvent;


  //  set mapping of device units to mouse transformations.  

  win_data->xform_map.rot = 0.5;
  win_data->xform_map.translate[0] = 0.1;
  win_data->xform_map.translate[1] = 0.1;
  win_data->xform_map.translate[2] = 0.1;
  win_data->xform_map.scale = 1.1;

  // button actions

  int n = 0;
  button_actions = new GrWindowActionfp[11];
  
  button_actions[n++] = &GrWindow::actionNoOp;
  button_actions[n++] = &GrWindow::actionRotxy;
  button_actions[n++] = &GrWindow::actionNoOp;
  button_actions[n++] = &GrWindow::actionTransxy;
  button_actions[n++] = &GrWindow::actionRotz;
  button_actions[n++] = &GrWindow::actionNoOp;
  button_actions[n++] = &GrWindow::actionTransz;
  button_actions[n++] = &GrWindow::actionPick;
  button_actions[n++] = &GrWindow::actionNoOp;
  button_actions[n++] = &GrWindow::actionScale;
  button_actions[n++] = &GrWindow::actionNoOp;


  //win_win_datas.push_back (this);

  grSystem.addWindow (this);

  } 

GrWindow::~GrWindow() { 
  //fprintf (stderr, ">>>>>> GrWindow dtor:  name [%s] \n", name.c_str());
  //fprintf (stderr, "   >>> scene [%x] \n", scene);
  delete scene;
  }

//////////////////////////////////////////////////////////////
//                    p u b l i c                          //
////////////////////////////////////////////////////////////

//*============================================================*
//*==========                 getColor               ==========*
//*============================================================*
// get widow background color.

void
GrWindow::getColor(GrColor& color) {
  color = this->color;
  }

//*============================================================*
//*==========                 setColor               ==========*
//*============================================================*
// set widow background color.

void
GrWindow::setColor(const GrColor& color) {
  this->color = color;
  }

//*============================================================*
//*==========                 open                   ==========*
//*============================================================*
// open a window.

void
GrWindow::open () {
  sys_win_createWindow();
  }

//*============================================================*
//*==========                 redraw                 ==========*
//*============================================================*
// redraw a window.

void
GrWindow::redraw () {

  if (!this->scene) { 
    return;
    }

  //fprintf (stderr, ">>>>>> GrWindow redraw \n");
  this->scene->render(); 
  }

//*============================================================*
//*==========           setCurrentWindow             ==========*
//*============================================================*
// set this window as the current drawing window.

void
GrWindow::setCurrentWindow () { 
  sys_win_setCurrentWindow();
  //scene->setViewport (width, height, extent);
  }

//*============================================================*
//*==========              getDimensions             ==========*
//*============================================================*
// get window dimensions.                         

void
GrWindow::getDimensions (int& width, int& height) {
  x = this->x;
  y = this->y;
  width = this->width;
  height = this->height;
  }

//*============================================================*
//*==========              setDimensions             ==========*
//*============================================================*
// set window dimensions.

void
GrWindow::setDimensions (const int& width, const int& height) {
  this->width = width;
  this->height = height;

  float xmin, xmax, ymin, ymax, zmin, zmax;
  float cx, cy, dx, dy, f;
  GrVector3 min, max;

  extent.get (min, max);
  xmin = min[0], ymin = min[1], zmin = min[2];
  xmax = max[0], ymax = max[1], zmax = max[2];

  cx = (xmax + xmin) / 2.0;
  cy = (ymax + ymin) / 2.0;
  dx = xmax - xmin;
  dy = ymax - ymin;

  if (width <= height) {
    f = (float)height / (float)width;
    dy = dx * f;
    ymin = cy - dy / 2.0;
    ymax = cy + dy / 2.0;
    }
  else {
    f = (float)width / (float)height;
    dx = dy * f;
    xmin = cx - dx / 2.0;
    xmax = cx + dx / 2.0;
    }

  min[0] = xmin, min[1] = ymin;
  max[0] = xmax, max[1] = ymax;
  extent.set (min, max);
  scene->setViewport (width, height, extent);
  }

//*============================================================*
//*==========              processEvents             ==========*
//*============================================================*
// process window events.

void
GrWindow::processEvents (const char *script_file, int cmd_input, const char *prompt,
                         char *cmd)
  {
  //fprintf (stdout, "\n>>>>>> GrWindow::processEvents \n");
  WinPrvData* data = (WinPrvData*)prv_data;

  if (cmd_input) {
    this->cmd_input = true;

    if (prompt) {
      strcpy (this->cmd_prompt, prompt);
      }

    fprintf (stdout, "\n");
    fprintf (stdout, "%s ", prompt);
    fflush (stdout);
    }

  sys_win_processEvents(cmd);
  }

//*============================================================*
//*==========              processMotionEvent        ==========*
//*============================================================*
// process window motion event.

void
GrWindow::processMotionEvent (int mx, int my)
  {
  int dx, dy;
  WinPrvData* prvd = (WinPrvData*)prv_data;
  prvd->event.type = GR_EVENT_MOUSE_MOTION; 
  if (debug) fprintf (stderr, ">>>>>> GrWindow::processMotionEvent:  \n");

  // save the last event //

  saveEvent();

  if (prvd->event.mouse_pos_set) { 
    dx = mx - prvd->event.mx;
    dy = my - prvd->event.my;

    if (prvd->event.action != GR_WINDOW_ACTION_NONE) {
      (this->*event_action)(mx, my, dx, dy);
      }
    }

  else {
    prvd->event.mx = mx;
    prvd->event.my = my;
    prvd->event.mouse_pos_set = true;
    }

  prvd->event.mx = mx;
  prvd->event.my = my;
  }

//*============================================================*
//*==========       processButtonReleaseEvent        ==========*
//*============================================================*
// process window button release event.

void
GrWindow::processButtonReleaseEvent (int mx, int my, int button_id) 
  {
  //fprintf (stderr, ">>>>>> GrWindow::processButtonReleaseEvent:  \n");
  WinPrvData* prvd = (WinPrvData*)prv_data;

  // save the last event
  saveEvent();

  prvd->event.mx = mx;
  prvd->event.my = my;

  prvd->event.mouse_pos_set = false;
  prvd->event.type = GR_EVENT_NONE;
  prvd->event.action = GR_WINDOW_ACTION_NONE;
  }

//*============================================================*
//*==========         processButtonPressEvent        ==========*
//*============================================================*
// process window button up event.

void
GrWindow::processButtonPressEvent (int mx, int my, int button_id)
  {
  //fprintf (stderr, ">>>>>> GrWindow::processButtonPressEvent:  \n");
  int dx, dy;
  WinPrvData* prvd = (WinPrvData*)prv_data;

  // save the last event
  saveEvent();

  prvd->event.type = GR_EVENT_BUTTON_PRESS; 
  prvd->event.button_id = (GrEventButtonId)button_id; 

  if (prvd->event.mouse_pos_set) {
    dx = mx - prvd->event.mx;
    dy = my - prvd->event.my;

    if ((button_id >= 0) && (button_id < GR_BUTTON_MAX)) {
      (this->*button_actions[button_id])(mx, my, dx, dy);
      }
    }

  else {
    prvd->event.mouse_pos_set = true;
    dx = 0;
    dy = 0;

    if ((button_id >= 0) && (button_id < GR_BUTTON_MAX)) {
      (this->*button_actions[button_id])(mx, my, dx, dy);
      }
    }

  prvd->event.mx = mx;
  prvd->event.my = my;
  }

//*============================================================*
//*==========                saveEvent               ==========*
//*============================================================*
// save last window event.

void
GrWindow::saveEvent ()
  {
  WinPrvData* prvd = (WinPrvData*)prv_data;
  prvd->event.last->mx = prvd->event.mx;
  prvd->event.last->my = prvd->event.my;
  prvd->event.last->type = prvd->event.type;
  prvd->event.last->button_state = prvd->event.button_state;
  }

//*============================================================*
//*==========              setExtent                 ==========*
//*============================================================*
// set window extent.     

void
GrWindow::setExtent (const GrExtent& extent) {
  this->extent = extent;
  }

//*============================================================*
//*==========              swapBuffers               ==========*
//*============================================================*
// swap window buffers.

void
GrWindow::swapBuffers () {
  sys_win_swapBuffers();
  }

//*============================================================*
//*==========                 writeImage             ==========*
//*============================================================*
// write the window image to a file.

void 
GrWindow::writeImage(const string name, GrImageFormat format, int *status)
  {
  sys_win_writeImage(name, format, status);
  }

//*============================================================*
//*==========              setRecordParams           ==========*
//*============================================================*
// set parameters to output a record of window images.

void
GrWindow::setRecordParams(const GrWindowRecordParams& params) 
  {
  WinPrvData* data = (WinPrvData*)prv_data;
  data->record.params = params; 
  }

//*============================================================*
//*==========              setEnableRecord           ==========*
//*============================================================*
// enable / disbale recording of window images.

void 
GrWindow::setEnableRecord(const bool flag)
  {
  //fprintf (stderr, "    >>> set record enabled [%d] \n", flag);
  WinPrvData* data = (WinPrvData*)prv_data;
  data->record.enabled = flag; 
  data->record.count = 1; 
  }

//*============================================================*
//*==========              getEnableRecord           ==========*
//*============================================================*
// enable / disbale recording of window images.

bool
GrWindow::getEnableRecord()
  {
  WinPrvData* data = (WinPrvData*)prv_data;
  return data->record.enabled;
  }

//*============================================================*
//*==========              record                    ==========*
//*============================================================*
// record window image.

void
GrWindow::record()
  {
  WinPrvData* data = (WinPrvData*)prv_data;
  int status;
  char count[100];
  string name;

  if (!data->record.enabled) {
    return;
    }

  if (data->record.params.prefix == "") {
    return;
    }

  name = data->record.params.prefix;
  sprintf (count, "_%d", data->record.count);
  name.append(count, strlen(count));
  writeImage (name, data->record.params.format, &status);
  fprintf (stderr, "    >>> record window into file \"%s\" \n", name.c_str()); 
  data->record.count += 1;
  }

////////////////////////////////////////////////////////////////
//           b u t t o n   e v e n t   f u n c s             //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              actionNoOp                ==========*
//*============================================================*

void
GrWindow::actionNoOp (int mx, int my, int dx, int dy) { }


//*============================================================*
//*==========              actionRotxy               ==========*
//*============================================================*
// process rotation x-y.

void
GrWindow::actionRotxy (int mx, int my, int dx, int dy)
  {
  GrVector3 rot;
  WinPrvData* prvd = (WinPrvData*)prv_data;

  prvd->event.action = GR_WINDOW_ACTION_ROTATE_XY;
  this->event_action = &GrWindow::actionRotxy;

  rot[0] = (float)dy * prvd->xform_map.rot;
  rot[1] = (float)dx * prvd->xform_map.rot;
  rot[2] = 0.0;

  if (prvd->event.type == GR_EVENT_BUTTON_PRESS) {
    return;
    }

 if ((fabs(rot[0]) < 1.0) && (fabs(rot[1]) < 1.0)) {
   return;
   }

  this->scene->rotate (true, rot); 
  }

//*============================================================*
//*==========              actionRotz                ==========*
//*============================================================*
// process rotation z.

void
GrWindow::actionRotz (int mx, int my, int dx, int dy)
  {
  GrVector3 rot;
  WinPrvData* prvd = (WinPrvData*)prv_data;

  prvd->event.action = GR_WINDOW_ACTION_ROTATE_Z;
  this->event_action = &GrWindow::actionRotz;
  rot[0] = 0.0; 
  rot[1] = 0.0;
  rot[2] = (float)(dx+dy)*prvd->xform_map.rot;

  if (prvd->event.type == GR_EVENT_BUTTON_PRESS) {
    return;
    }

 if (fabs(rot[2]) < 1.0) {
   return;
   }

  this->scene->rotate (true, rot);
  }

//*============================================================*
//*==========              actionPick                ==========*
//*============================================================*
// process pick.

void
GrWindow::actionPick (int mx, int my, int dx, int dy)
  {
  //fprintf (stderr, ">>>>>> GrWindow::actionPick:  \n");

  GrPickResult pick;
  int screen_pt[2];
  WinPrvData* prvd = (WinPrvData*)prv_data;

  prvd->event.action = GR_WINDOW_ACTION_PICK;
  prvd->event.action_func = &GrWindow::actionPick;

  prvd->do_pick = true;
  prvd->pick_x = mx;
  prvd->pick_y = this->height - my;

  screen_pt[0] = mx;
  screen_pt[1] = this->height - my;

  this->scene->performPick(screen_pt, pick);
  }

//*============================================================*
//*==========              actionTransxy             ==========*
//*============================================================*
// process translate x-y.

void
GrWindow::actionTransxy (int mx, int my, int dx, int dy)
  {
  GrVector3 trans;
  WinPrvData* prvd = (WinPrvData*)prv_data;

  prvd->event.action = GR_WINDOW_ACTION_TRANSLATE_XY;
  this->event_action = &GrWindow::actionTransxy;

  trans[0] = (float)dx  * prvd->xform_map.translate[0];
  trans[1] = (float)-dy * prvd->xform_map.translate[1];
  trans[2] = 0.0;


  // translate the scene
  this->scene->translate (true, trans);
  }

//*============================================================*
//*==========              actionTransz              ==========*
//*============================================================*
// process translate z.

void
GrWindow::actionTransz (int mx, int my, int dx, int dy)
  {

  }



//*============================================================*
//*==========              actionScale               ==========*
//*============================================================*
// process scale.

void
GrWindow::actionScale (int mx, int my, int dx, int dy)
  {

  float d, scale;
  WinPrvData* prvd = (WinPrvData*)prv_data;

  prvd->event.action = GR_WINDOW_ACTION_SCALE;
  this->event_action = &GrWindow::actionScale;
  d = (float)(dx + dy);

  if (d < 0.0) {
    scale = 1.0 / prvd->xform_map.scale;
    }
  else if (d == 0.0) {
    scale = 1.0;
    }
  else {
    scale = prvd->xform_map.scale;
    }

  this->scene->scale (true, scale);
  }


} // PmGraphics



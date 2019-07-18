
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
//* win:                  w i n d o w                          *
//*============================================================*

#ifndef _WINDOW_PRV_GR_H_
#define _WINDOW_PRV_GR_H_

#include "pm/gr/gr.h"

namespace PmGraphics {


#define GR_WINDOW_SYSTEM_NT


//  window key functions  

#define  GR_WIN_KEY_ZOOM  'z'


/*  button actions.  */

typedef enum {
  GR_WINDOW_ACTION_NONE,
  GR_WINDOW_ACTION_SCALE,
  GR_WINDOW_ACTION_SCALE_X,
  GR_WINDOW_ACTION_SCALE_Y,
  GR_WINDOW_ACTION_SCALE_Z,
  GR_WINDOW_ACTION_ROTATE_XY,
  GR_WINDOW_ACTION_ROTATE_X,
  GR_WINDOW_ACTION_ROTATE_Y,
  GR_WINDOW_ACTION_ROTATE_Z,
  GR_WINDOW_ACTION_TRANSLATE_XY,
  GR_WINDOW_ACTION_TRANSLATE_X,
  GR_WINDOW_ACTION_TRANSLATE_Y,
  GR_WINDOW_ACTION_TRANSLATE_Z,
  GR_WINDOW_ACTION_LOOKAT,
  GR_WINDOW_ACTION_PICK,
  GR_WINDOW_ACTION_KEY_PRESS,
  GR_WINDOW_ACTION_MAX
  } GrWindowAction;


/*  button ids.  */

typedef enum {
  GR_BUTTON_UNKNOWN,
  GR_BUTTON_LEFT,
  GR_BUTTON_MIDDLE,
  GR_BUTTON_RIGHT,
  GR_BUTTON_SHIFT_LEFT,
  GR_BUTTON_SHIFT_MIDDLE,
  GR_BUTTON_SHIFT_RIGHT,
  GR_BUTTON_CTRL_LEFT,
  GR_BUTTON_CTRL_MIDDLE,
  GR_BUTTON_CTRL_RIGHT,
  GR_BUTTON_MAX
  } GrEventButtonId;

typedef enum {
  GR_EVENT_BUTTON_STATE_NONE,
  GR_EVENT_BUTTON_STATE_SHIFT,
  GR_EVENT_BUTTON_STATE_CONTROL,
  GR_EVENT_BUTTON_STATE_CONTROL_SHIFT
  } GrEventButtonState;

typedef enum {
  GR_EVENT_NONE,
  GR_EVENT_BUTTON_PRESS,
  GR_EVENT_BUTTON_RELEASE,
  GR_EVENT_MOUSE_MOTION,
  GR_EVENT_WINDOW_EXPOSE,
  GR_EVENT_KEY_PRESS,
  GR_EVENT_KEY_RELEASE
  } GrEventType;


//  button action functions 

typedef void (GrWindow::*GrWinActionfp)(int, int, int, int);


//  stores window event info  

typedef struct GrWindowEvent {
  int mx, my;
  int mouse_pos_set;
  GrEventType type;
  GrEventButtonState button_state;
  GrEventButtonId button_id;
  GrWindowAction action;
  GrWinActionfp action_func;
  struct GrWindowEvent *last;
  } GrWindowEvent; 


//  stores window transform mapping. 

typedef struct GrWindowXformMap {
  float rot;
  GrVector3 translate;
  float scale;
  } GrWindowXformMap; 


// pick  object  

typedef struct WinPickObj {
  int pick_x, pick_y;
  GrVector3 front_pt, back_pt;
  GrVector3 world_front_pt, world_back_pt;
  } WinPickObj; 


typedef struct WinRecord {
  int enabled;
  int count;
  GrWindowRecordParams params;
  } WinRecord;


typedef struct WinPrvData {
  GrWindowEvent event;
  bool do_pick;
  int pick_x, pick_y;
  void (*pick_func)();

  bool do_keyp;
  int keyp_x, keyp_y;
  char keyp_string[80];

  int font_id;
  WinPickObj pick_info;
  int ren_init;
  GrColor color;
  WinRecord record;

  GrWindowXformMap xform_map;

  // windowing system data
  void *sys_win_data;
  } WinPrvData ;


} // PmGraphics 


#endif

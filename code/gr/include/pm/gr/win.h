
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

#ifndef _WINDOW_GR_H_
#define _WINDOW_GR_H_

#include "pm/gr/gr.h"
//#include "scene.h"

namespace PmGraphics {

class GrScene;


// GrWindow
// ------------

 class PM_EXPORT GrWindow {

   public:

     GrWindow(const string name, int x, int y, int width, int height, int type);
     ~GrWindow();

     void getColor(GrColor& color); 
     void setColor(const GrColor& color); 
     void setCurrentWindow(); 

     void getCmdInfo (int *cmd_input, GrCmdPrompt prompt) { 
            *cmd_input = this->cmd_input;
            strcpy (prompt, this->cmd_prompt);
            }

     void getDimensions (int& width, int& height);
     void setDimensions (const int& width, const int& height);
     void setIdleFunction ( void (*func)()) { idle_func = func; }
     void getIdleFunction ( void (**func)()) { *func = idle_func; }

     void setExtent (const GrExtent& extent);

     void processMotionEvent (int x, int y);
     void processButtonReleaseEvent (int mx, int my, int button_id);
     void processButtonPressEvent (int mx, int my, int button_id);

     void setScene (GrScene* scene) { this->scene = scene; }

     void open();
     void redraw();
     void processEvents (const char *script_file, int cmd_input, const char *prompt,
                         char *cmd);
     void swapBuffers();
     void writeImage(const string name, GrImageFormat format, int *status);
     void setRecordParams(const GrWindowRecordParams& params); 
     void setEnableRecord(const bool flag); 
     bool getEnableRecord();
     void record();

   private:

     string name;
     int type;
     int x, y;
     int width, height;
     int win_id;
     GrRenderMode render_mode;
     GrExtent extent;
     GrExtent orig_extent;
     GrScene *scene;
     GrColor color;

     int cmd_input;
     GrCmdPrompt cmd_prompt;
     void (*idle_func)();

     static vector<class GrWindow*> win_objs;

     void *prv_data;

     void init();
     void saveEvent();

     // button actions
     void actionNoOp (int, int, int, int);
     void actionRotxy (int, int, int, int);
     void actionRotz (int, int, int, int);
     void actionPick (int, int, int, int);
     void actionTransxy (int, int, int, int);
     void actionTransz (int, int, int, int);
     void actionScale (int, int, int, int);

     typedef void (GrWindow::*GrWindowActionfp)(int, int, int, int);
     GrWindowActionfp *button_actions, event_action;

     // window system interface //
     void sys_win_init_display();
     void sys_win_createWindow();
     void sys_win_processEvents(char *cmd);
     void sys_win_setCurrentWindow();
     void sys_win_swapBuffers();
     void sys_win_writeImage(const string name, GrImageFormat format, int *status);
   };

}


#endif





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

#ifndef _GR_WINDOW_GLUT_H_
#define _GR_WINDOW_GLUT_H_

#include "pm/gr/gr.h"  
#include "pm/gr/win.h"  
#include "win_prv.h"  
#include <stdio.h>

#ifdef USE_CYGWIN_SPECIAL_HEADERS
#include "glut.h"
#else
#include <GLUT/glut.h>
#endif

namespace PmGraphics {

#define min(a,b) ((a)<(b))?(a):(b)


//  ascii codes for keys    

#define KEY_ESCAPE     27
#define KEY_RETURN     13
#define KEY_BACKSPACE   8

#define KEY_ARROW_LEFT   100
#define KEY_ARROW_UP     101
#define KEY_ARROW_RIGHT  102
#define KEY_ARROW_DOWN   103


// window globals  

static void *display = NULL;
static int screen = 0;


// glut window data 

typedef struct GrGlutWin {
  int x, y;
  int type;
  int width, height;
  int font_id;
  int glut_win;
  int cmd_input;
  char cmd_prompt[100];
  } GrGlutWin;


// widow table

typedef struct GrGlutWinMap {
  int wid;
  GrWindow *win;
  } GrGlutWinMap;


static vector<GrGlutWinMap*>win_table;


/*------------------------------------------------------------*
 *               ****  colormaps ****                         *
 *------------------------------------------------------------*/

static float color_maps[] = {
    0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 
    0.000000, 1.000000, 0.333333, 0.776471, 0.443137, 0.556863, 
    0.443137, 0.556863, 0.219608, 0.666667, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.039216, 0.078431, 0.117647, 0.156863, 
    0.200000, 0.239216, 0.278431, 0.317647, 0.356863, 0.400000, 
    0.439216, 0.478431, 0.517647, 0.556863, 0.600000, 0.639216, 
    0.678431, 0.717647, 0.756863, 0.800000, 0.839216, 0.878431, 
    0.917647, 0.956863, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000, 
    1.000000, 1.000000, 0.000000, 0.000000, 1.000000, 1.000000, 
    0.333333, 0.443137, 0.776471, 0.556863, 0.443137, 0.219608, 
    0.556863, 0.666667, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.039216, 0.078431, 0.117647, 0.156863, 0.200000, 0.239216, 
    0.278431, 0.317647, 0.356863, 0.400000, 0.439216, 0.478431, 
    0.517647, 0.556863, 0.600000, 0.639216, 0.678431, 0.717647, 
    0.756863, 0.800000, 0.839216, 0.878431, 0.917647, 0.956863, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.333333, 0.443137, 
    0.443137, 0.219608, 0.776471, 0.556863, 0.556863, 0.666667, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.039216, 0.078431, 
    0.117647, 0.156863, 0.200000, 0.239216, 0.278431, 0.317647, 
    0.356863, 0.400000, 0.439216, 0.478431, 0.517647, 0.556863, 
    0.600000, 0.639216, 0.678431, 0.717647, 0.756863, 0.800000, 
    0.839216, 0.878431, 0.917647, 0.956863, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
};

//////////////////////////////////////////////////////////////
//                    c a l l b a c k s                    //
////////////////////////////////////////////////////////////

void glut_WinExpose ();
void glut_WinIdle ();
void glut_WinReshape (int width, int height);
void glut_WinKeyProc (unsigned char key, int x, int y);
void glut_WinSpecialProc (int key, int x, int y);

void glut_EventButp (int button, int state, int x, int y);
void glut_EventButr (int button, int state, int x, int y);
void glut_EventMotion (int x, int y);


/*------------------------------------------------------------*
 *               ****  event callbacks  ****                  *
 *------------------------------------------------------------*/

//static void (*win_event_func[10000])();
static int g_curr_win = -1; 

}

#endif


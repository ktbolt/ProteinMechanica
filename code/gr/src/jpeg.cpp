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
//* jpeg:        j p e g       i n e r f a c e                 *
//*============================================================*

#include "jpeg.h"
#ifdef GR_USE_JPEG
#include "jpeglib.h"
#include <setjmp.h>
#endif

namespace PmGraphics {

#ifdef GR_USE_JPEG
static struct jpeg_compress_struct cinfo;
static struct jpeg_error_mgr jerr;
static FILE *g_fp; 

struct my_error_mgr {
    struct jpeg_error_mgr pub;    /* "public" fields */
    jmp_buf setjmp_buffer;        /* for return to caller */
    };

METHODDEF(void)
my_error_exit (j_common_ptr cinfo) { 
  fprintf (stderr, "\n\n>>>>> jpeg read error. \n\n"); 
  }
#endif

//*============================================================*
//*==========              gr_JpegInit               ==========*
//*============================================================*
// initialize jpeg.  

void 
gr_JpegInit()
  {
#ifdef GR_USE_JPEG
  cinfo.err = jpeg_std_error (&jerr);
  jpeg_create_compress (&cinfo);
#endif
  }

//*============================================================*
//*==========         gr_JpegStartCompress           ==========*
//*============================================================*
// initialize jpeg.

void 
gr_JpegStartCompress(char *name, int width, int height) 
  {

#ifdef GR_USE_JPEG
  FILE *fp;

  static int quality = 95;

  static char *fn = "gr_JpgFileSetup";

  fp = fopen (name, "wb");

  if (!fp) {
    return;
    }

  g_fp = fp;
  jpeg_stdio_dest (&cinfo, fp);

  cinfo.image_width      = width;
  cinfo.image_height     = height;
  cinfo.input_components = 3;
  cinfo.in_color_space   = JCS_RGB;

  jpeg_set_defaults (&cinfo);
  jpeg_set_quality (&cinfo, quality, TRUE);
  jpeg_start_compress (&cinfo, TRUE);
#endif
  }

//*============================================================*
//*==========         gr_JpegWriteScanline           ==========*
//*============================================================*

void 
gr_JpegWriteScanline (int id, int width, void *rgb) 
  {
#ifdef GR_USE_JPEG

  JSAMPROW row_pointer[1];

  row_pointer[0] = (JSAMPROW)rgb; 

  (void) jpeg_write_scanlines (&cinfo, row_pointer, 1);
#endif
  }

//*============================================================*
//*==========         gr_JpegFinishCompress          ==========*
//*============================================================*

void
gr_JpegFinishCompress()
  {
#ifdef GR_USE_JPEG
  jpeg_finish_compress (&cinfo);
  fclose (g_fp);
  jpeg_destroy_compress (&cinfo);
#endif
  }

//*============================================================*
//*==========         gr_JpegFinishCompress          ==========*
//*============================================================*

void 
gr_JpegWriteImage(char *name, int width, int height, void *image) 
  {
#ifdef GR_USE_JPEG

  struct jpeg_compress_struct cinfo;

  struct jpeg_error_mgr jerr;

  FILE *fp;

  int row_stride;

  JSAMPROW row_pointer[1];

  unsigned char *idata;

  int i;

  int quality;

  cinfo.err = jpeg_std_error (&jerr);
  jpeg_create_compress (&cinfo);
  fp = fopen (name, "wb");

  if (!fp) {
    return;
    }

  jpeg_stdio_dest (&cinfo, fp);

  cinfo.image_width      = width;
  cinfo.image_height     = height;
  cinfo.input_components = 3;
  cinfo.in_color_space   = JCS_RGB;

  jpeg_set_defaults (&cinfo);

  quality = 95;
  jpeg_set_quality (&cinfo, quality, TRUE);

  jpeg_start_compress (&cinfo, TRUE);

  idata = (unsigned char*)image;
  row_stride = width * 3;
  i = height - 1;

  while (cinfo.next_scanline < cinfo.image_height) {
    /*
    fprintf (stderr, "scanline %d \n", i);
    */
    row_pointer[0] = & idata[i * row_stride];
    (void) jpeg_write_scanlines (&cinfo, row_pointer, 1);
    i -= 1;
    }

  jpeg_finish_compress (&cinfo);
  fclose (fp);
  jpeg_destroy_compress (&cinfo);
#endif
  }

//*============================================================*
//*==========         gr_JpegReadImage               ==========*
//*============================================================*

void 
gr_JpegReadImage(char *name, int *width, int *height, void **image)
  {
#ifdef GR_USE_JPEG

  struct jpeg_decompress_struct cinfo;

  struct my_error_mgr jerr;

  FILE *fp;

  int row_stride;

  JSAMPARRAY buffer;

  unsigned char *idata, *ptr;

  int i, j, iw, ih;

  *width = 0;
  *height = 0;
  *image = NULL;
  fp = fopen (name, "rb");

  if (!fp) {
    return;
    }

  cinfo.err = jpeg_std_error (&jerr.pub);
  jerr.pub.error_exit = my_error_exit;

  if (setjmp(jerr.setjmp_buffer)) {
    jpeg_destroy_decompress (&cinfo);
    fclose (fp);
    return;
    }

  jpeg_create_decompress (&cinfo);

  jpeg_stdio_src (&cinfo, fp);

  (void) jpeg_read_header (&cinfo, TRUE);

  (void) jpeg_start_decompress (&cinfo);

  iw = cinfo.output_width;
  ih = cinfo.output_height;
  fprintf (stderr, ">>>>>> iw [%d] \n", iw);
  fprintf (stderr, ">>>>>> ih [%d] \n", ih);
  fprintf (stderr, ">>>>>> sizeof(JSAMPARRAY) [%d] \n", sizeof(JSAMPARRAY));
  idata = new unsigned char[iw*ih*3]; 

  row_stride = cinfo.output_width * cinfo.output_components;
  buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
  i = ih - 1;

  while (cinfo.output_scanline < cinfo.output_height) {
    (void) jpeg_read_scanlines (&cinfo, buffer, 1);
    ptr = buffer[0];

    for (j = 0; j < row_stride; j++) {
      /*
      idata[i*row_stride + j] = buffer[j];
      */
      idata[i*row_stride + j] = ptr[j];
      }

    i -= 1;
    }

  *width  = iw;
  *height = ih;
  *image = idata;
#endif
  }

}


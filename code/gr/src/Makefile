
#SHELL = /usr/bin/sh

include gr_src.mk

# define include dirs
#--------------------
glut_inc = -framework GLUT 
gr_inc = -framework OpenGL 
jpeg_inc= -I../jpg/jpeg-8b/ 
inc = -I.  -I../../include \
      -I../include \
      -I../../common/include \
      $(jpeg_inc) $(gr_inc) $(glut_inc)

hdr = ../include/pm/gr/

common_dir = ../../common/src/
common_objs = $(common_dir)/common.o $(common_dir)/mth.o
common_lib = /Users/parkerda/software/ProteinMechanica/lib/libpm_common.dll

# define libs
#------------
glut_lib = -framework GLUT 
opengl_lib = -framework OpenGL 
x_lib = -L/usr/X11R6/lib  -lGlU -lXext -lX11
x_lib = 
jpeg_lib =  ../jpg/jpeg-8b/.libs/libjpeg.a
gr_libs = $(glut_lib) $(opengl_lib) $(x_lib) $(jpeg_lib)
#libs = $(gr_libs) $(common_objs) -lm
libs = $(gr_libs) $(common_lib) 

# special cywin flags for my pos computer
#----------------------------------------
#cwf = -DUSE_CYGWIN_SPECIAL_HEADERS
cwf =

# windows flags
#--------------
wflags = 

# use for debug: code crashes with -g -funroll-loops -march=pentium4 -mfpmath=sse
#copt = -g
#flags = $(inc) $(copt)  $(cwf)

# jpeg flags 
#-----------
jpeg = -DGR_USE_JPEG

# use for optimized execution 
copt = -O2 -pedantic -Wno-long-long
flags = -DFREEGLUT_STATIC $(inc) $(cwf) $(copt) $(wflags) $(jpeg)

sw_lib = /Users/parkerda/software/ProteinMechanica/lib/
lib = $(sw_lib)/libpm_gr.dll
link_flags = 

$(lib): $(obj)
	g++ -shared -o $(lib) $(obj) $(libs) $(link_flags) 

.cpp.o: 
	g++ $(flags) -c $<

clean:
	rm -f *.o;


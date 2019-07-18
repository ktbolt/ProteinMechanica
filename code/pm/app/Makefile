#----------------------------------------------#
# pm makefile                                  #
#----------------------------------------------#
cwd = `pwd`

src = \
      main.cpp
obj = $(src:.cpp=.o)

# define include dirs
#--------------------
pm_dir = ../src/
inc = -I$(pm_dir)  -I../include   -I../../common/include


# define libs
#------------
glut_lib = -framework GLUT 
opengl_lib =  -framework OpenGL
x_lib = 
lib_dir = /Users/parkerda/software/ProteinMechanica/lib/
sw_libs = $(lib_dir)/libpm_gr.dll $(lib_dir)/libpm_core.dll $(lib_dir)/libpm_common.dll
gr_libs = $(glut_lib) $(opengl_lib) $(x_lib)
libs = $(sw_libs) $(gr_libs) -lm

copt = -O2
copt = -g

# windows flags
flags = $(inc)  

exe = p

$(exe): $(obj) 
	g++ -o $(exe) $(obj) $(libs) 
	@cp $(exe) $(HOME)/bin

.cpp.o:
	g++ $(copt) $(flags) -c $<


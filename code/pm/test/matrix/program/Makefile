
SHELL = /usr/bin/sh

src = \
      main.cpp 
obj = $(src:.cpp=.o) 


# define include dirs

glut_inc = -I/usr/include/GL//include
opengl_inc = -I/usr/X11R6/include/
gr_inc = $(glut_inc) $(opengl_inc) 
inc =     \
      -I/home/davidp/research/software/molecular_simulation/pms++/pm++/src/      \
      -I/home/davidp/research/software/molecular_simulation/pms++/gr++/include   \
      -I/home/davidp/research/software/molecular_simulation/pms++/pm++/include/  \
      -I/home/davidp/research/software/molecular_simulation/pms++/common/include/


# define libs 

glut_lib = -L/lib -lglut32
opengl_lib = -lglu32 -lopengl32
x_lib = -L/usr/X11R6/lib  -lGL -lGlU -lXext -lX11
gr_libs = $(glut_lib) $(opengl_lib) $(x_lib)
sw_libs = -L/home/davidp/research/software/lib/ -lgr++ -lpm++
libs =  $(sw_libs) $(gr_libs) -lstdc++ -lComctl32 -lkernel32 -luser32 -lgdi32 -lm

copt = -O3
copt = -g
flags = $(inc) $(copt) 

exe = p 

$(exe): $(obj) 
	g++ -o $(exe) $(obj) $(libs)

.cpp.o: 
	g++ $(flags) -c $<

clean:
	rm -f *.o;


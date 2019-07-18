
# SHELL = /usr/bin/sh

# include sources
#----------------
include pm_src.mk


# define include dirs
#--------------------
glut_inc = -framework GLUT 
opengl_inc = -framework OpenGL 
gr_inc = $(glut_inc) $(opengl_inc)
inc = -I. \
      -I../../gr/include   \
      -I../include/        \
      -I../../common/include/ 

common_dir = ../../common/src/
common_objs = $(common_dir)/common.o \
              $(common_dir)/mth.o

# define libs
#------------
ode_path = /usr/local/
ode_lib = $(ode_path)/lib/
ode_libs = $(ode_lib)/libode.a 

glut_lib = -framework GLUT 
opengl_lib = -framework OpenGL 

x_lib = -L/usr/X11R6/lib  -lGL -lGlU -lXext -lX11
x_lib = 
gr_libs = $(glut_lib) $(opengl_lib) $(x_lib)
lib_dir = /Users/parkerda/software/ProteinMechanica/lib/
sw_libs = $(lib_dir)/libpm_common.dll $(lib_dir)/libpm_gr.dll 

libs =  $(sw_libs) $(gr_libs) $(ode_libs) $(cpp_libs) -lm

# special cywin flags for my pos computer
#----------------------------------------
copt = -g -O2
flags = -w $(inc) $(copt) -DUSE_GRAPHICS   

sw_lib = /Users/parkerda/software/ProteinMechanica/lib/
lib = $(sw_lib)/libpm_core.dll

$(lib): $(obj)
	g++ -shared -o $(lib) $(obj) $(common_objs) $(libs)

.cpp.o: 
	g++ $(flags) -c $<

atoms.o: atoms.cpp  atoms.h 
	g++ $(flags) -c $<

atom.o: atom.cpp atom.h atom_prv.h
	g++ $(flags) -c $<

body.o: body.cpp  body.h 
	g++ $(flags) -c $<

cmd.o: cmd.cpp cmd.h $(cmd_deps) 
	g++ $(flags) -c $<

cmpd.o: cmpd.cpp cmpd.h 
	g++ $(flags) -c $<

csys.o: csys.cpp csys.h 
	g++ $(flags) -c $<

db.o: db.cpp  db.h db_cerfacs.cpp db_gro.cpp db_pm.cpp db_pdb.cpp
	g++ $(flags) -c $<

force.o: force.cpp force.h 
	g++ $(flags) -c $<

gc.o: gc.cpp gc.h 
	g++ $(flags) -c $<

geom.o: geom.cpp geom.h 
	g++ $(flags) -c $<

graphics.o: graphics.cpp  graphics.h 
	g++ $(flags) -c $<

grid.o: grid.cpp grid.h 
	g++ $(flags) -c $<

joint.o: joint.cpp joint.h 
	g++ $(flags) -c $<

mesh.o: mesh.cpp mesh.h 
	g++ $(flags) -c $<

model.o: model.cpp model.h 
	g++ $(flags) -c $<

mol.o: mol.cpp mol.h 
	g++ $(flags) -c $<

motor.o: motor.cpp motor.h 
	g++ $(flags) -c $<

msg.o: msg.cpp 
	g++ $(flags) -c $<

msr.o: msr.cpp msr.h
	g++ $(flags) -c $<

ode_solv.o: ode_solv.cpp ode_solv.h  rbsolv.h
	g++ $(flags) -DdDOUBLE -c $<

part.o: part.cpp part.h 
	g++ $(flags) -c $<

pobj.o: pobj.cpp pobj.h
	g++ $(flags) -c $<

pot.o: pot.cpp pot.h
	g++ $(flags) -c $<

pm.o: pm.cpp ../include/pm/pm.h ../../common/include/pm/common.h
	g++ $(flags) -c $<

rbsim.o: rbsim.cpp rbsim.h sim.o
	g++ $(flags) -c $<

rbsolv.o: rbsolv.cpp rbsolv.h
	g++ $(flags) -c $<

res.o: res.cpp res.h
	g++ $(flags) -c $<

rest.o: rest.cpp rest.h
	g++ $(flags) -c $<

sd.o: sd.cpp sd.h 
	g++ $(flags) -c $<

sim.o: sim.cpp sim.h 
	g++ $(flags) -c $<

sobj.o: sobj.cpp sobj.h 
	g++ $(flags) -c $<

solid.o: solid.cpp solid.h 
	g++ $(flags) -c $<

state.o: state.cpp state.h
	g++ $(flags) -c $<

surf.o: surf.cpp surf.h
	g++ $(flags) -c $<

pmsys.o: pmsys.cpp pmsys.h ../include/pm/pm.h
	g++ $(flags) -c $<

trace.o: trace.cpp trace.h
	g++ $(flags) -c $<

clean:
	rm -f *.o;


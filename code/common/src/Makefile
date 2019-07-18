#SHELL = /usr/bin/sh

src = common.cpp mth.cpp
obj = $(src:.cpp=.o) 

# define include dirs

inc = -I. \
      -I../include  \
      -I../../pm++/include/ 

copt = -g
copt = -O2
flags = $(inc) $(copt) -DUSE_GRAPHICS  

#obj: $(obj)

#sw_lib = /Users/parkerda/software/protein-mechanica/protmech/protmech/trunk/lib/
sw_lib = /Users/parkerda/software/ProteinMechanica/lib/
lib = $(sw_lib)/libpm_common.dll
link_flags =

$(lib): $(obj)
	g++ -shared -o $(lib) $(obj) $(libs) $(link_flags)

.cpp.o: 
	g++ $(flags) -c $<

clean:
	rm -f *.o;


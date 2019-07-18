src = \
      main.cpp 
obj = $(src:.cpp=.o) 


# define include dirs
inc = -I../../../common/include/pm \
      -I../../../common/include/\
      -I../../src/  \
      -I../../include/pm/gr \
      -I../../include/pm/ \
      -I../../include/


# define libs 

sw_libs = ../../../lib/libcommon.dll  \
          ../../../lib/libgr++.dll

libs =  $(sw_libs) -lm

copt = -O3
copt = -g
flags = $(inc) $(copt) $(wflags)

exe = p 

$(exe): $(obj) 
	g++ -o $(exe) $(obj) ../../src/gr.o $(libs) 


.cpp.o: 
	g++ $(flags) -c $<

clean:
	rm -f *.o;


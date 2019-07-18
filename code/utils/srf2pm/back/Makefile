
include $(MAKE_INCLUDE)
SHELL = /usr/bin/sh

src =  main.c  
obj = $(src:.c=.o) 

inc = $(dm_inc)
copt = -g
copt = -O3

flags = $(inc) $(copt) -DWINDOWS
flags = $(inc) $(copt) 

libs = $(dm_slib)

exe = srf2pm

$(exe): main.o
	$(mach_cc) $(flags) -o $(exe) main.c $(libs)
	@cp srf2pm.exe $(sw_bin)

.c.o:
	$(mach_cc) $(flags) -c $< 

clean:
	rm -f *.o


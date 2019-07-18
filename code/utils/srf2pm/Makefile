
include $(MAKE_INCLUDE)
SHELL = /usr/bin/sh

src =  main.c  
obj = $(src:.c=.o) 

inc = $(dm_inc)
copt = -O3
copt = -g

flags = $(inc) $(copt) -DWINDOWS
flags = $(inc) $(copt) 

libs = 

exe = srf2pm

$(exe): main.o
	gcc $(flags) -o $(exe) main.c $(libs)
	@cp srf2pm.exe /home/dparker/bin

.c.o:
	gcc $(flags) -c $< 

clean:
	rm -f *.o


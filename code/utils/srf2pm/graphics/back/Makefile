
include $(MAKE_INCLUDE)
SHELL = /usr/bin/sh

src =  main.c read.c poly.c
obj = $(src:.c=.o) 

inc = $(dm_inc) $(gr_inc)  $(glut_inc)
opt = -g
flags = $(gr_inc) $(inc) $(opt)

libs = $(gr_lib)      \
       $(dm_slib)      \
       $(glut_lib)     \
       $(mach_gr_lib)  \
       $(mach_x_lib)   \
       -lm

libs = $(gr_slib)       \
       $(dm_slib)      \
       $(glut_lib)     \
       $(mach_gr_lib)  \
       $(mach_x_lib)   \
       -lm

exe = p    

$(exe): $(obj) 
	$(mach_cc) -o $(exe) $(obj) $(libs) 

.c.o: 
	$(mach_cc) $(flags) -c $<


CC = g++
FLAGS = -g -w -mno-cygwin  -DWIN32
LIBS = -lpthread
LIBS = pthreadVC1.dll
EXECS = pt

pt:main.cpp
	$(CC) $(FLAGS) -o $@ $? $(LIBS) -lm -mno-cygwin


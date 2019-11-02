#
# Makefile for Linux machines---the location of libraries and header files
# may be different depending on your flavor of linux and where the libraries
# have been installed.
#
# Rename to "Makefile" and type "make" to compile the simple example gtk
# programs.
#

CC = gcc

CFLAGS = -g

INCLUDE = -I. \
  -I/usr/include/ -I/usr/include/X11/ \
  -I/usr/include/gtk-1.2 \
  -I/usr/include/glib-1.2 -I/usr/lib64/glib/include

LDFLAGS =  -L. \
  -L/usr/lib64 \
  -L/usr/X11R6/lib

LDLIBS = \
  -lGL -lGLU \
  -lX11 -lXext -lXmu \
  -lgtk -lgdk -lglib -lgtkgl \
  -lm

SRCS = \
gtkglsimple.c

EX3_OBJS = \
gtkglsimple.o

.c.o:
	$(CC) $(DEFS) $(INCLUDE) $(CFLAGS) -c $<

all: gtkglsimple


gtkglsimple: $(EX3_OBJS) gtkglsimple.c
	$(CC) $(CFLAGS) $(INCLUDE) $(EX3_OBJS) $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f *.o core gtkglsimple

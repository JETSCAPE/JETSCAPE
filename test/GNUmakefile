# ===========================================================================
#  Makefile for MUSIC                                 Chun Shen Dec. 12, 2014
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++-mp-5
CFLAGS= -std=c++11 -Wall -O3 $(shell gsl-config --cflags) -fopenmp
#CFLAGS= -Wall -g $(shell gsl-config --cflags) -fopenmp

music_path="./music/src"

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	mpihydro
endif

#SRC		=	test_MUSIC_main.cpp music.cpp ../src/fluid_dynamics.cpp \
#			$music_path/grid.cpp $music_path/eos.cpp $music_path/init.cpp \
#            $music_path/util.cpp $music_path/evolve.cpp $music_path/advance.cpp \
#			$music_path/u_derivative.cpp \
#			$music_path/reconst.cpp $music_path/dissipative.cpp \
#			$music_path/minmod.cpp $music_path/grid_info.cpp \
#			$music_path/cornelius.cpp
#
#INC		= 	music.h $music_path/grid.h $music_path/eos.h $music_path/init.h \
#			$music_path/util.h $music_path/data.h \
#			$music_path/evolve.h $music_path/advance.h \
#			$music_path/u_derivative.h $music_path/reconst.h \
#			$music_path/dissipative.h \
#			$music_path/minmod.h $music_path/grid_info.h \
#			$music_path/cornelius.h
#
# -------------------------------------------------

OBJDIR		=	obj
TARGET		=	$(MAIN)
INSTPATH	=	../

OBJECTS		= $(OBJDIR)/test_MUSIC_main.o $(OBJDIR)/music.o $(OBJDIR)/fluid_dynamics.o

($OBJDIR)/test_MUSIC_main.o: test_MUSIC_main.cpp
	$(CC) $(CFLAGS) -c test_MUSIC_main.cpp -o $(OBJDIR)/test_MUSIC_main.o
($OBJDIR)/music.o: music.cpp
	$(CC) $(CFLAGS) -c music.cpp -o $(OBJDIR)/music.o
($OBJDIR)/fluid_dynamics.o: ../src/fluid_dynamics.cpp
	$(CC) $(CFLAGS) -c ../src/fluid_dynamics.cpp -o $(OBJDIR)/fluid_dynamics.o


.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

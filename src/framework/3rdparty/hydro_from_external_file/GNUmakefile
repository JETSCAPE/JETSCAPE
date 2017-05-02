# ===========================================================================
#  Makefile photon_emission                           Chun Shen May 5, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := h5c++
CFLAGS = -O3 -Wall

RM		=	rm -f
O               =       .o
HDF5LD          = 
HDF5FLAGS       = 
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	hydro_analysis.e
endif

SRC		=	main.cpp Hydroinfo_h5.cpp FluidcellStatistic.cpp \
            ParameterReader.cpp arsenal.cpp SurfaceFinder.cpp \
            cornelius.cpp Hydroinfo_MUSIC.cpp

INC		= 	Hydroinfo_h5.h Stopwatch.h FluidcellStatistic.h \
            ParameterReader.h arsenal.h SurfaceFinder.h cornelius.h \
			Hydroinfo_MUSIC.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(HDF5LD) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(LDFLAGS) $(HDF5FLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
main.cpp : Hydroinfo_h5.h Hydroinfo_MUSIC.h Stopwatch.h FluidcellStatistic.h ParameterReader.h
Hydroinfo_h5.cpp : Hydroinfo_h5.h
Hydroinfo_MUSIC.cpp : Hydroinfo_MUSIC.h
FluidcellStatistic.cpp : FluidcellStatistic.h Hydroinfo_h5.h Hydroinfo_MUSIC.h
ParameterReader.cpp : ParameterReader.h arsenal.h
arsenal.cpp : arsenal.h
SurfaceFinder.cpp : Hydroinfo_h5.h Hydroinfo_MUSIC.h ParameterReader.h SurfaceFinder.h cornelius.h
cornelius.cpp : cornelius.h

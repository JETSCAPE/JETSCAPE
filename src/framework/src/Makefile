#
# Simple Makefile for test only !!!
#
######################

CXX           = g++
LD            = g++
AR	      = ar

#############################

os = $(shell uname -s)

ifeq ($(os),Linux)
CXXFLAGS      = -g -O -fPIC -std=c++11
#-I$(ROOTSYS)/include
else
CXXFLAGS      = -O -fPIC -pipe -Wall -g -std=c++11 -stdlib=libc++ -pthread
# -Wno-deprecated-writable-strings -Wno-deprecated-declarations -Wno-unused-value -Wno-unused-variable -Wno-overloaded-virtual -Wno-return-type-c-linkage #-I$(ROOTSYS)/include
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -dynamiclib -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ROOTLIBS      = -lpthread -stdlib=libc++ -lm -ldl #$(shell root-config --libs)

LIBS          =  $(ROOTLIBS)


HDRS          = JetScape.h JetScapeModuleBase.h JetScapeTask.h JetScapeLogger.h JetScapeXML.h JetEnergyLoss.h JetEnergyLossManager.h JetScapeSignalManager.h FluidDynamics.h ElossModulesTest.h tinyxml2.h
SRCS          = JetScape.cc JetScapeModuleBase.cc JetScapeTask.cc JetScapeLogger.cc JetScapeXML.cc JetEnergyLoss.cc JetEnergyLossManager.cc JetScapeSignalManager.cc FluidDynamics.cc ElossModulesTest.cc tinyxml2.cc

OBJS          = $(SRCS:.cc=.o)


SLIB          = libJetScape.so

all: $(SLIB)

$(SLIB): $(OBJS)
	@echo "Linking $(SLIB) for "$(os)
#	@echo $(LD) $(LDFLAGSS) $(OBJS) $(LIBS) -o $(SLIB)
	$(LD) $(LDFLAGSS) $(OBJS) $(LIBS) -o ../lib/$(SLIB)
	@echo "done"


##########################################################

.cc.o:
	$(CXX) $(CXXFLAGS) -I./src/ -c -o $@ $<

.PHONY: clean cleanObj

clean :
	rm -f *.o ../lib/$(SLIB)
	@echo "All cleaned-up :-)"

cleanObj:
	rm -f *.o

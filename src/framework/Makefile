######################

CXX           = g++
LD            = g++

#############################

os = $(shell uname -s)

ifeq ($(os),Linux)
CXXFLAGS      = -O -Wall -I./lib
CXXFLAGSSTAT  = -g -O -std=c++11 -I./lib
else
CXXFLAGS      = -O -fPIC -pipe -Wall -std=c++11 -stdlib=libc++ -pthread -I./lib
CXXFLAGSSTAT  = -O -std=c++11 -I./lib
endif

ifeq ($(os),Linux)
LDFLAGS       = -O
LDFLAGSSTAT    = -O -static
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSSTAT   = -O -std=c++11 -Bstatic
endif

testJetScape :
	$(CXX) $(CXXFLAGS) -c src/testJetScape.cc
	$(LD) $(LDFLAGS) testJetScape.o -L./lib -lJetScape -lpthread -stdlib=libc++ -lm -ldl -o bin/testJetScape
	rm testJetScape.o

static :
	$(CXX) $(CXXFLAGSSTAT) -c src/testJetScape.cc
#	$(LD) $(LDFLAGSSTAT) testJetScape.o -L./lib -lJetScape.a -lpthread -stdlib=libc++ -lm -ldl -o bin/testJetScapeStatic
	$(LD) $(LDFLAGSSTAT) testJetScape.o ./lib/libJetScape.a -lpthread -stdlib=libc++ -lm -ldl -o bin/testJetScapeStatic
	rm testJetScape.o

clean :
	rm -f *.o ./bin/testJetScape
	@echo "All cleaned-up :-)"


# DO NOT DELETE

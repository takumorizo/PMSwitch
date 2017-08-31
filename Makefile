.PHONY: clean
BOOST=libs/boost_1_45_0_subset/
SRCDIR=src
OBJDIR=obj
BINDIR=bin
CXX=g++
#CXX=clang++

CPPFLAGS=-std=c++11 -DNDEBUG
#-DLOGDEBUG -DDEBUGREADS #-DNDEBUG -DMMTEST
# CXXFLAGS=-I$(SAMTOOLDIR) -I$(BOOST) -I$(INCLUDEDIR) -I$(SRCDIR) -Wno-deprecated -O3
CXXFLAGS=-I$(BOOST) -I$(INCLUDEDIR) -I$(SRCDIR) -O3 -g3
LDFLAGS=

SOURCES := $(shell find $(SRCDIR) -name "*.cpp")
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all: pre-build main-build

pre-build:
	mkdir -p obj bin

main-build: pmswitch

pmswitch: $(OBJECTS)
	$(CXX) -o $(BINDIR)/$@ $(CXXFLAGS) $(CPPFLAGS) $(OBJECTS) $(LDFLAGS)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

dependencies:
	cd libs; make

clean:
	rm -f $(OBJECTS)

test:
	cd tests;make test

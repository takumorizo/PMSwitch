.PHONY: clean
BOOST=../libs/boost_1_45_0_subset/
GTESTDIR=libs/gtest-1.7.0
GTEST=$(GTESTDIR)/include

SRCDIR=../src
INCLUDEDIR=../include
OBJDIR=obj
BINDIR=bin

CXX=g++
CPPFLAGS=-std=c++11 #-DLOGDEBUG -DDEBUGREADS #-DNDEBUG -DMMTEST
CXXFLAGS=-I$(BOOST) -I$(GTEST) -I$(INCLUDEDIR) -I$(SRCDIR) -O3 -Wno-deprecated
LDFLAGS=-L$(GTESTDIR) -lgtest

SOURCES := $(shell find $(SRCDIR) -name "*.cpp" ! -name "main.cpp")
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all: pre-build main-build

pre-build:
	mkdir -p $(OBJDIR) $(BINDIR)

main-build: pmswitch_test

dependencies:
	cd libs;make

test: clean all
	bin/pmswitch_test data

pmswitch_test: $(OBJECTS) obj/main.o
	$(CXX) -o bin/$@ $(CXXFLAGS) $(CPPFLAGS) $(OBJECTS) obj/main.o $(LDFLAGS)

obj/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*

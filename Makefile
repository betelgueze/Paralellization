##
# Makefile for compiling GAL project
# Author: Jiri Hon, xhonji01@stud.fit.vutbr.cz
# Date: 2014-12-09

CXX=g++
CXXFLAGS=-O2 -fopenmp
LDFLAGS=-O2  
TARGET=econn
OGDF_LIB=./lib/libOGDF.a
OGDF_INCLUDE=./include

all: $(TARGET)

$(TARGET): $(TARGET).o $(OGDF_LIB)
	$(CXX) $(LDFLAGS) $^ -o $@ $(OGDF_LIB) -lgomp -pthread

$(TARGET).o : $(TARGET).cpp $(OGDF_LIB)
	$(CXX) $(CXXFLAGS) -I$(OGDF_INCLUDE) -c $(TARGET).cpp

$(OGDF_LIB):
	wget http://www.stud.fit.vutbr.cz/~xhonji01/GAL/resources.tar.gz
	tar -xzf resources.tar.gz

resources:
	tar -czf resources.tar.gz experiments graphs include lib
clean:
	rm *.o $(TARGET)

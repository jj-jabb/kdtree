SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
DS := $(patsubst %.cpp,%.d,$(wildcard *.cpp))
EXOBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CXX = g++
CC = gcc
CXXFLAGS:= -std=c++0x -O3 -Wall -c -fmessage-length=0 -I./
LDLIBS =-L./ -L/usr/lib/ -L/usr/local/lib/ -L/usr/lib/x86_64-linux-gnu/
LIBS =-lstdc++ -lm
RM=rm -f
all : $(EXOBJS)
	$(CXX) $(LDLIBS) -o kdtree $(EXOBJS) $(LIBS)
		
clean :
	clear
	$(RM) $(OBJS) $(EXOBJS) $(DS)
	$(RM) kdtree
	
default : all

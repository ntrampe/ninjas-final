#
# Nicholas Trampe
# 02/17/2015
# CS 5201
# Assignment 5 - Makefile

# The purpose of this makefile is to build an executable for Assignment 5

# Generic GNU g++ compiler flags
CXX = /usr/bin/g++

CXXFLAGS = -std=c++11 -g -Wall -W -pedantic-errors -O3 -pthread
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS =  $(SOURCES:%.cpp=%.o)

.SUFFIXES: .cpp

%.o: %.cpp
	@echo "Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

default: all

all: driver

driver: ${OBJECTS}
	@echo "Building $@"
	@$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@
	@echo "Done"

clean:
	-@ rm -f core > /dev/null 2>&1
	-@ rm -f driver > /dev/null 2>&1
	-@ rm -f ${OBJECTS} > /dev/null 2>&1

# Dependency lists for the various source files

depend: $(SOURCES) $(HEADERS)
	@echo "Generating dependencies"
	@$(CXX) -MM *.cpp > $@

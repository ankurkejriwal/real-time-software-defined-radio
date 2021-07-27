# this makefile is intended for g++ on Linux

CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -g -O3 
LDFLAGS = -pthread
SOURCES = radio.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp support.cpp threadfunc.cpp processing.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = radio

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS)

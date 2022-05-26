CC=g++
LDFLAGS=-std=c++11 -O3 -lm -Wall
SOURCES=src/legalizer.cpp src/main.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=macroLegalizer
#INCLUDES=src/cell.h src/net.h src/partitioner.h
INCLUDES=src/legalizer.h src/macro.h

all: $(SOURCES) bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o:  %.c  ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o bin/$(EXECUTABLE) *.plt vertical* horizontal*

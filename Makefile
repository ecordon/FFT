# Makefile for FFT program

CC = g++
CFLAGS = -c -Wall
OBJECTS = FFT.o Complex.o
EXECUTABLE = FFT

all: $(OBJECTS) $(EXECUTABLE)
clean: 
	rm $(OBJECTS) $(EXECUTABLE)

FFT: $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE)

Complex.o: Complex.cpp Complex.h
	$(CC) $(CFLAGS) Complex.cpp

FFT.o: FFT.cpp Complex.h Input.h
	$(CC) $(CFLAGS) FFT.cpp

FFT1.o: FFT1.cpp Complex.h Input.h
	$(CC) $(CFLAGS) FFT1.cpp


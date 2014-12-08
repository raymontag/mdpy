# Makefile for molecular dynamics program
#
# Eric Jeckelmann -- Uni Mainz -- October 2003
#

CC=g++  -Wall -O3 -lm #-lg2c
#CC=g++  -Wall -g -DDEBUG -lm -lg2c

#HEADERS= particle.h setup.h state.h dynamics.h
SOURCES= main.cc setup.cc state.cc 

clean:
	rm -f hcmd fhcmd *.o

hcmd: 
	$(CC)  -o hcmd  $(SOURCES) dynamics.cc 

fhcmd: 
	$(CC)  -o fhcmd  $(SOURCES) fastdynamics.cc  

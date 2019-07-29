CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas -lconfig++
CC = g++


pmcmc: pmcmc.cc 
	${CC} ${CFLAGS} -o pmcmc pmcmc.cc ${LDFLAGS} 

sim: simulate.cc 
	${CC} ${CFLAGS} -o simulate simulate.cc ${LDFLAGS} 

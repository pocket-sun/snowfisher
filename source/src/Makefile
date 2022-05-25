prefix = ${GLB_DIR}
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CC = g++
CFLAGS = -Wall
local_CFLAGS = -g -O4
debug_CFLAGS = -DNODEBUG -DNOERROR -DNOPARADEBUG
OPT_FLAGS = -march=native -ffp-contract=fast -DARMA_NO_DEBUG -fopenmp

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
ARMA_INCLUDE_PATH = /usr/include
MCMC_HEADER_DIR = /home/cls/local/mcmc/include
LIBS= -L/home/cls/local/mcmc -lmcmc -lblas -llapack


BIN = supernova 
OBJ = detector.o pinched.o pysnow.o main.o 

all: $(BIN)

supernova:
	$(CC) $(debug_CFLAGS) $(CFLAGS) $(local_CFLAGS) $(OPT_FLAGS) -o $(BIN) *.cpp \
	-I$(ARMA_INCLUDE_PATH) -I$(MCMC_HEADER_DIR) $(INCFLAGS) $(local_LDFLAGS) $(LIBS)

install:
	cp $(BIN) ../

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) $(SOBJ)

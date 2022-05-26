prefix = ${GLB_DIR}
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CC = g++
CFLAGS = -Wall
local_CFLAGS = -g -O4
OPT_FLAGS = -march=native -ffp-contract=fast -DARMA_NO_DEBUG -fopenmp


INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
ARMA_INCLUDE_PATH = /usr/include
MCMC_HEADER_DIR = /home/cls/local/mcmc/include
SNOWFLUX_HEADER = ${SNOWGLOBEPATH}src
LIBS= -L/home/cls/local/mcmc -lmcmc -lblas -llapack -L${SNOWFLUX_HEADER} -lsnowflux


BIN = fisher
OBJ = fisher.o mcmc.o

all: $(BIN)

#fisher:
#	$(CC) $(CFLAGS) $(local_CFLAGS) $(OPT_FLAGS) -o $(BIN) *.cpp \
#	-I$(ARMA_INCLUDE_PATH) -I$(MCMC_HEADER_DIR) -I$(SNOWFLUX_HEADER) \
#	-Wl,-rpath $(SNOWFLUX_HEADER) $(LIBS)
fisher: fisher.o mcmc.o
	$(CC) -Wl,-rpath,$(SNOWFLUX_HEADER) mcmc.o fisher.o -o $(BIN) $(LIBS) $(local_LDFLAGS) $(OPT_FLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(local_CFLAGS) $(OPT_FLAGS) -c \
	$< $(INCFLAGS) -I$(ARMA_INCLUDE_PATH) -I$(MCMC_HEADER_DIR) -I$(SNOWFLUX_HEADER) \
	 -Wl,-rpath,$(SNOWFLUX_HEADER) $(local_LDFLAGS) $(LIBS)

install:

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ) $(SOBJ)

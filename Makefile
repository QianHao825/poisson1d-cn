##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
#
# -- option to dedicated machine
#
HOSTNAME?=$(shell hostname)
include $(HOSTNAME).mk

#
# -- Compiler Option
OPTC=${OPTCLOCAL}

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=${LIBSLOCAL}

# -- Include directories
INCLATLAS=${INCLUDEBLASLOCAL}
INCL= -I $(TPDIR)/include $(INCLATLAS)

#
#################################################################
# makefile
############
#
SOL?=

OBJENV= tp_env.o

# ---- Libraries (add sparse for EX10) ----
OBJLIBPOISSON = lib_poisson1D$(SOL).o \
                lib_poisson1D_writers.o \
                lib_poisson1D_richardson$(SOL).o \
                lib_poisson1D_sparse.o

# ---- Programs ----
OBJTP2ITER   = $(OBJLIBPOISSON) tp_poisson1D_iter.o
OBJTP2DIRECT = $(OBJLIBPOISSON) tp_poisson1D_direct.o
OBJTP2SPARSE = $(OBJLIBPOISSON) tp_poisson1D_sparse.o

.PHONY: all run testenv tpPoisson1D_iter tpPoisson1D_direct tpPoisson1D_sparse \
        run_testenv run_tpPoisson1D_iter run_tpPoisson1D_direct run_tpPoisson1D_sparse clean

all: bin/tp_testenv bin/tpPoisson1D_iter bin/tpPoisson1D_direct bin/tpPoisson1D_sparse

run: run_testenv run_tpPoisson1D_iter run_tpPoisson1D_direct run_tpPoisson1D_sparse

testenv: bin/tp_testenv
tpPoisson1D_iter: bin/tpPoisson1D_iter
tpPoisson1D_direct: bin/tpPoisson1D_direct
tpPoisson1D_sparse: bin/tpPoisson1D_sparse

%.o : $(TPDIRSRC)/%.c
	$(CC) $(OPTC) -c $(INCL) $<

bin/tp_testenv: $(OBJENV)
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tpPoisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tpPoisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tpPoisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tpPoisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

# ---- EX10 ----
bin/tpPoisson1D_sparse: $(OBJTP2SPARSE)
	$(CC) -o bin/tpPoisson1D_sparse $(OPTC) $(OBJTP2SPARSE) $(LIBS)

run_testenv:
	bin/tp_testenv

run_tpPoisson1D_iter:
	bin/tpPoisson1D_iter
	bin/tpPoisson1D_iter 1
	bin/tpPoisson1D_iter 2

run_tpPoisson1D_direct:
	bin/tpPoisson1D_direct
	bin/tpPoisson1D_direct 1
	bin/tpPoisson1D_direct 2

# ---- EX10 ----
run_tpPoisson1D_sparse:
	bin/tpPoisson1D_sparse

clean:
	rm -f *.o bin/*

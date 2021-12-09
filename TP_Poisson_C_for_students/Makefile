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
CC=gcc

# 
# -- Compiler Option
#
OPTC=-O3 -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=-llapacke -lblas -lm

# -- Include directories
INCLBLASLAPACK= -I /usr/include

INCL= -I $(TPDIR)/include $(INCLBLASLAPACK) 
#
#################################################################
# makefile
############
#
OBJENV= tp_env.o
OBJTP2ITER= lib_poisson1D.o tp2_poisson1D_iter.o
OBJTP2DIRECT= lib_poisson1D.o tp2_poisson1D_direct.o
#

all: bin/tp_testenv bin/tp2poisson1D_iter bin/tp2poisson1D_direct

testenv: bin/tp_testenv

tp2poisson1D_iter: bin/tp2poisson1D_iter

tp2poisson1D_direct: bin/tp2poisson1D_direct

tp_env.o: $(TPDIRSRC)/tp_env.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c 

lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

tp2_poisson1D_iter.o: $(TPDIRSRC)/tp2_poisson1D_iter.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_iter.c  

tp2_poisson1D_direct.o: $(TPDIRSRC)/tp2_poisson1D_direct.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_direct.c  

bin/tp_testenv: $(OBJENV) 
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tp2poisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tp2poisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tp2poisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tp2poisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

run_testenv:
	bin/tp_testenv

run_tp2poisson1D_iter:
	bin/tp2poisson1D_iter

run_tp2poisson1D_direct:
	bin/tp2poisson1D_direct

clean:
	rm *.o bin/*

PYTHON_INCLUDES:=$(shell python-config --includes)
PYTHON_LIBS:=$(shell python-config --ldflags)
NUMPY_INCLUDES:=-I$(shell python -c "import numpy; print numpy.get_include()")
FFTW_LIB_PATH=/usr/lib/
CC=gcc
CFLAGS=-Wall -std=c99 -O3 -funroll-loops -fPIC

UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  EXTRA_LDFLAGS=
else
  EXTRA_LDFLAGS=
endif

LDFLAGS=-L$(FFTW_LIB_PATH) -lfftw3 -lfftw3_threads -lm -lgsl -lblas $(PYTHON_LIBS) $(EXTRA_LDFLAGS)

all: _grid_model.so

_grid_model.so: src/grid_model.c src/grid_model.h src/regional.c src/regional.h
	mkdir -p bin
	swig -python src/grid_model.i
	$(CC) -fPIC -O3 -c -fno-strict-aliasing -DDSFMT_MEXP=19937 -o bin/dsfmt.o src/dSFMT/dSFMT.c
	$(CC) $(CFLAGS) -c src/regional.c -o bin/regional.o
	$(CC) $(CFLAGS) -c src/grid_model_wrap.c -o bin/grid_model_wrap.o $(PYTHON_INCLUDES) $(NUMPY_INCLUDES)
	$(CC) $(CFLAGS) -c -DDSFMT_MEXP=19937 src/grid_model.c -o bin/grid_model.o
	$(CC) $(CFLAGS) -shared bin/regional.o bin/dsfmt.o bin/grid_model.o bin/grid_model_wrap.o -o _grid_model.so $(LDFLAGS)
	mv src/grid_model.py .

clean:
	rm -f *.pyc
	rm -rf bin/
	rm src/grid_model_wrap.c
	rm _grid_model.so grid_model.py

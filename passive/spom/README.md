A stochastic patch occupancy model
==================================

 This is an modified implementation of the stochastic patch occupancy model used in 

  - Rybicki and Hanski. Species–area relationships and extinctions caused by habitat loss and fragmentation.
    [Ecology Letters (2013) 16: 27–38.][1]

  This package provides a C implementation of the computational model together with a minimal Python interface for running simulations. Two example scripts illustrate how to run a simulation on a randomly generated auto-correlated habitat. For details of the model, see the article and its supplementary material. 

  If you use the model in your work, please cite the above article.

  [1]: http://dx.doi.org/10.1111/ele.12065

Authors
-------

  - Joel Rybicki joel.rybicki@aalto.fi
  - Ilkka Hanski ilkka.hanski@helsinki.fi

Requirements
------------

  The simulation software depends on the following libraries:
  - numpy http://www.numpy.org/
  - SWIG 2.0 http://www.swig.org/
  - FFTW 3.3 http://www.fftw.org/
  - GNU Scientific Library http://www.gnu.org/software/gsl/
  - dSFMT random number generator (included) http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
  - matplotlib (for plotting, optional) http://matplotlib.org/

  If you are running OS X with homebrew and pip installed, you can install the dependencies by running:
    $ brew install swig fftw gsl && pip install numpy matplotlib

  The software should compile on GNU/Linux and OS X.

Installation
------------

  To compile the C module, run 'make'. If this does not work, you may need to modify 'model/Makefile'. In this case, see 'model/README' for further instructions.

Usage
-----

  When you have compiled the module, you can run the example scripts with the following commands:
    
    python example.py 

    python steps.py

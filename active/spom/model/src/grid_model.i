/* SWIG interface file for convolution.h */
%module grid_model

%{
#define SWIG_FILE_WITH_INIT
#include "grid_model.h"
%}

%include "numpy.i"


%init %{
    import_array();
%}

%apply (double* IN_ARRAY2, int DIM1, int DIM2 ) { (double* habitat_m, int h, int w) }
%apply (unsigned char* IN_ARRAY2, int DIM1, int DIM2 ) { (unsigned char* arr, int h, int w) }
%apply (unsigned char** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2) { (unsigned char **view_array, int *h, int *w) }
%apply (double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2) { (double **view_array, int *h, int *w)}
%apply (unsigned char** ARGOUTVIEW_ARRAY1, int *DIM1) { (unsigned char **view_array, int *n) }

%apply (int DIM1, int* ARGOUT_ARRAY1) { (int steps, int *return_vals)}
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int *out, int n)}

%include "grid_model.h"
#clear  (double* habitat_m, int h, int w);
#clear (unsigned char *arr, int h, int w);
#clear  (unsigned char **view_array, int *h, int *w) 
#clear  (double **view_array, int *h, int *w);
#clear  (int steps, int *return_vals);
#clear  (int *out, int n);
#clear (unsigned char **view_array, int *n);

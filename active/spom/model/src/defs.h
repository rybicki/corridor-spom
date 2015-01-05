#ifndef __DEFS_H_
#define __DEFS_H_

#define ZERO_THRESHOLD 1e-15 

/* Uncomment to use the Gaussian kernel instead of exponential kernel? */
#define USE_GAUSSIAN_KERNEL

/* Are we using long doubles for computation? SWIG does not support these. */
typedef long double t_float;
#define __USING_LONG_DOUBLES

#endif

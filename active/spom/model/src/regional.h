#ifndef __REGIONAL_H_
#define __REGIONAL_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "defs.h"

typedef struct stochasticity_grid {
    t_float *arr;
    int height, width;
    double weight_factor;
    double stddev;
    double mean;
    int depth;
    gsl_rng *rng;
    double *weights;
} t_stochasticity_grid;

t_stochasticity_grid *create_stochasticity_grid(t_float *arr, int height, int width, double weight_factor, double mean, double stddev, long int seed);
void free_stochasticity_grid(t_stochasticity_grid *g);
void generate_stochasticity(t_stochasticity_grid *g);

#endif /* __REGIONAL_H_ */

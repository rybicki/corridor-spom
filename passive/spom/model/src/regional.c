#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "regional.h"

void init_weights(t_stochasticity_grid *g) {
    /* The weights will be normalized s.t. the sum of squared weights is one */
    g->weights[0] = 1.0;
    double sum = 1.0;
    for (int i = 1; i<g->depth; i++) {
        double val = g->weights[i-1] / g->weight_factor;
        g->weights[i] = val;
        sum += val*val;
    }
    double N = sqrt(sum);
    for (int i = 0; i < g->depth; i++) {
        g->weights[i] /= N;
    }
}

/*
 'arr' should point to a double array of size at least height*width.
*/
t_stochasticity_grid *create_stochasticity_grid(double *arr, int height, int width, double weight_factor, double mean, double stddev, long int seed) {
    t_stochasticity_grid *g = malloc( sizeof(t_stochasticity_grid) );
    g->arr = arr;
    g->height = height;
    g->width = width;

    /* It is best to use array dimensions of power of 2 */

    int d = height > width ? height : width; // Take the max 
    g->depth = ceil(log2(d)) + 1;
    g->stddev = stddev;
    g->mean = mean;
    g->weight_factor = weight_factor;

    g->rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(g->rng, seed);

    g->weights = malloc( sizeof(double) * g->depth );
    init_weights(g);

    return g;
}

void free_stochasticity_grid(t_stochasticity_grid *g) {
    gsl_rng_free(g->rng);
    free(g->weights);
    g->weights = NULL;
    g->arr = NULL; 
    free(g);
}

static inline double gaussian(const gsl_rng *rng, double stddev) {
    return gsl_ran_gaussian_ziggurat(rng, stddev);
}

void fill_quadrant(const t_stochasticity_grid *g, int y0, int x0, int y1, int x1, int level) {
    double weight = g->weights[level];
    double rval = gaussian(g->rng, g->stddev);
    for(int i = y0; i < y1; i++) {
        for(int j = x0; j < x1; j++) {
            int offset = i * g->width + j;
            double val = rval*weight;
            g->arr[offset] += val;
        }
    }
}

void recursive_division(const t_stochasticity_grid *g, int level, int y0, int x0, int y1, int x1) {
    if (level >= g->depth) return;
    if (x1-x0 == 1 && y1-y0 == 1) return; // Skip 1x1 grids

    if (y0 >= y1) return;
    if (x0 >= x1) return; 
    //printf("Depth: %i. (%i, %i, %i, %i)\n", level, y0, x0, y1, x1);

    fill_quadrant(g, y0, x0, y1, x1, level);

    if (level + 1 == g->depth) return; // Avoid unnecessary recursion; skip level with 1x1 grids

    int ydist = (y1-y0)/2.0 + 0.5;
    int xdist = (x1-x0)/2.0 + 0.5;
    int ymid = y0+ydist;
    int xmid = x0+xdist;

    /* Upper left */
    recursive_division(g, level+1, y0, x0, ymid, xmid);

    /* Upper right */
    recursive_division(g, level+1, ymid, x0, y1, xmid);

    /* Bottom left */
    recursive_division(g, level+1, y0, xmid, ymid, x1);

    /* Bottom right */
    recursive_division(g, level+1, ymid, xmid, y1, x1);
}

/**
 Generate spatially autocorrelated coefficients using a recursive quadtree subdivision. 
 The pointer 'arr' is a 2-dimensional output array of size 'height' * 'width'. 
 The variable 'variance' gives the variance of the log-norm distribution with mean and cut-off at 1.0.
*/
void generate_stochasticity(t_stochasticity_grid *g) {
    /* Initialize the array by generating the last level (i.e., 1x1 grids) */
    for(int i = 0; i < g->height; i++) {
        for(int j = 0; j < g->width; j++) {
            double weight = g->weights[g->depth-1];
            g->arr[i*g->width + j] = weight*gaussian(g->rng, g->stddev);
        }
    }

    recursive_division(g, 0, 0, 0, g->height, g->width);

    for(int i = 0; i < g->height; i++) {
        for(int j = 0; j < g->width; j++) {
            double y = g->arr[i*g->width + j];
            y = exp(y + g->mean);
            if (y > 1.0) y = 1.0;
            g->arr[i*g->width + j] = y;
        }
    }

}


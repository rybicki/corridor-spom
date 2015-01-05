/**
 Module containing C structures for the grid patch occupancy model.
 The module implements the basic dynamics of the model and a method 
 for computing connectivity in the grid via Fourier transforms.
*/

/* Use the native complex type */
#include <complex.h>
#include "fftw3.h"

#include "defs.h"
#include "regional.h"

#ifndef _GRID_MODEL_H_
#define _GRID_MODEL_H_

#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )

typedef struct {
    /* Model parameters */
    double alpha;       /* Inverse of the average dispersion distance. */
    double phenotype; /* The phenotype parameter (phi) */
    double colonization_rate; /* The colonization parameter (c) */
    double extinction_rate; /* Extinction rate (e) */
    double fitness_bell_width; /* Width of the gaussian bell (gamma) */
    int connectivity_iterations; /* NEW: Directed active dispersal */
} t_species_parameters;

typedef struct {
    /* Regional stochasticity parameters */
    double mean;
    double stddev;
    double weight_factor;
} t_stochasticity_parameters;

/*
 The grid object is used to manipulate the landscape grid.
 The object contains the h*w landscape matrix and necessary
 structures for computing connectivity fast.

 The physical size of the occupancy, kernel and connectivity
 matrices are 2h * 2w. These dimensions are denoted as n0 and n1. 
*/
typedef struct {
    int w, h;         /* The logical grid dimensions */
    int n0, n1;       /* The grid dimensions with padding */
    unsigned char *occ_m; /* Occupancy (bit) matrix */
   
    t_species_parameters species; /* The struct containig species parameters. */ 
    t_float *fit_m;    /* The fitness matrix */

    t_float *dir_m; /* Normalization factors for directed dispersal */

    /* Parameters related to regional stochasticity affecting fitness */
    int use_regional_stochasticity; /* If true compute regional stochasticity */
    t_stochasticity_parameters regional_stochasticity; /* Parameters for regional stochasticity */
    t_stochasticity_grid *stochasticity_grid; 

    /* Internal stuff */

    t_float *work_m;    /* Some working space for computing connectivity */
    t_float *conn_m;   /* Connectivity matrix */

    double *rand_buf; /* A buffer for random number generation */
    unsigned int rand_buf_len; /* Length of the above buffer */

    int fft_n;        /* The length of the FFTd complex arrays */
    fftwl_complex *occ_fft;  /* FFTd occupancy and fitness matrix */
    fftwl_complex *ker_fft;  /* FFTd kernel matrix */
    
    /* FFTW plans */
    fftwl_plan p_occ_fwd, /* Plan for FFTing the occupancy matrix. */
             p_occ_bwd,  /* Plan for inverse FFT. */
             p_ker_fwd;  /* Plan for FFTing the kernel matrix. */

} t_grid;

t_species_parameters create_species_parameters(double phenotype, double colonization_rate, double extinction_rate, double alpha, double fitness_bell_width, int connectivity_iterations);

/* Initialize this module. The function will initialize the FFTW and dSFMT libraries. 
   - 'nthreads' denotes the number of maximum threads available to FFTW
   - 'rng_seed' is the seed for dSFMT random number generator
*/
void init_fg(int nthreads, long rng_seed);
void init_fg_with_wisdom(int nthreads, long rng_seed, const char *fftw_wisdom_file);
/* This will free up any resources used by the module */
void finalize_fg(void);

/* Constructors and destructors for grid structures */
t_grid *create_grid(int g0, int g1, double phenotype, double colonization_rate, double extinction_rate, double alpha, double fitness_bell_width, int connectivity_iterations); 
t_grid *create_grid_from_struct(int g0, int g1, t_species_parameters species_params);
void free_grid(t_grid *g); 

/* Setup regional stochasticity */
void enable_regional_stochasticity(t_grid *g, t_stochasticity_parameters parameters);
void disable_regional_stochasticity(t_grid *g);

/* Calculate fitness from the given habitat_m matrix. 
   The size should be g->h x g->w. */
void calculate_fitness(t_grid *g, double *habitat_m, int h, int w);
void init_direction_factors(t_grid *g);

void copy_occupancy_from(t_grid *g, unsigned char *arr, int h, int w);

/* Give access to the internal matrices. Be careful with these! */
void give_occupancy_matrix(t_grid *g, unsigned char **view_array, int *h, int *w);

void give_fitness_matrix(t_grid *g, double **view_array, int *h, int *w);
void give_connectivity_matrix(t_grid *g, double **view_array, int *h, int *w);
void give_work_matrix(t_grid *g, double **view_array, int *h, int *w);

void give_conn_double_matrix(t_grid *g, double **view_array, int *h, int *w);

/* Accessories for the species fitness values */ 
double get_fitness(t_grid *g, int y, int x);
void set_fitness(t_grid *g, int y, int x, double val);

/* Accessories for the species occupancy matrix */
int is_occupied(t_grid *g, int y, int x);
void set_occupied(t_grid *g, int y, int x, unsigned char val);

/* Call this *before* any calls to simulate_step(). */
void initialize_simulation(t_grid *g);

/* This advance the simulation by one iteration. Call this after setting up the fitness matrix and initial occupancy values */ 
int simulate_step(t_grid *g);

void many_steps(t_grid *g, int steps, int *return_vals);

#endif

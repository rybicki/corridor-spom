#include "grid_model.h"
#include "defs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//#define NDEBUG

#include <assert.h>

#include "dSFMT/dSFMT.h"

/* This is the state object from the dSFMT random number generator */
dsfmt_t rng_state;

void init_fg(int nthreads, long rng_seed) {
    init_fg_with_wisdom(nthreads, rng_seed, NULL);
}

/*
 Initialize the FFTW library. The nthreads parameters
 gives an upper bound for the threads to be used by FFTW.
*/
void init_fg_with_wisdom(int nthreads, long rng_seed, const char *fftwl_wisdom_file) {
    /* Initialize FFTW */
    fftwl_init_threads();
    fftwl_plan_with_nthreads(nthreads);

    fftwl_import_system_wisdom();

    if (fftwl_wisdom_file != NULL) { 
        /* FIXME: Does importing wisdom ever work? */
       fftwl_import_wisdom_from_filename(fftwl_wisdom_file);
    }
    /* Initialize the random number generator */
    dsfmt_init_gen_rand(&rng_state, rng_seed);

}

/*
 Free resources allocated by the FFTW library.
*/
void finalize_fg(void) {
    fftwl_cleanup_threads();
    fftwl_cleanup();
}

/*
 Calls the FFTW library functions for initializing 
 the Fourier transform plans. The parameter
 flags is passed to fftwl_plan_dft_xxx() calls.
*/
void init_grid_plans(t_grid *g, int flags) {
    assert(g != NULL);
//    printf("init_grid_plans\n");
    g->p_occ_fwd = fftwl_plan_dft_r2c_2d(
                   g->n0, g->n1, g->work_m, g->occ_fft, flags
    );
    g->p_occ_bwd = fftwl_plan_dft_c2r_2d(
                   g->n0, g->n1, g->occ_fft, g->conn_m, flags
    );
    g->p_ker_fwd = fftwl_plan_dft_r2c_2d(
                   g->n0, g->n1, g->work_m, g->ker_fft, flags
    );
}

void free_grid_plans(t_grid *g) {
    fftwl_destroy_plan(g->p_occ_fwd);
    fftwl_destroy_plan(g->p_occ_bwd);
    fftwl_destroy_plan(g->p_ker_fwd);
}

t_species_parameters create_species_parameters(double phenotype, double colonization_rate, double extinction_rate, double alpha, double fitness_bell_width, int connectivity_iterations) {
    t_species_parameters params;
    params.phenotype = phenotype;
    params.colonization_rate = colonization_rate;
    params.extinction_rate = extinction_rate;
    params.alpha = alpha;
    params.fitness_bell_width = fitness_bell_width;
    params.connectivity_iterations = connectivity_iterations;
    return params;
}

t_grid *create_grid(int g0, int g1, double phenotype, double colonization_rate, double extinction_rate, double alpha, double fitness_bell_width, int connectivity_iterations) {
    t_species_parameters params = create_species_parameters(phenotype, colonization_rate, extinction_rate, alpha, fitness_bell_width, connectivity_iterations);
    return create_grid_from_struct(g0, g1, params);
}

/*
 Creates a new grid object for species with parameters given in the struct
 species_params.  The parameters g0 and g1 respectively 
 correspond to the height and width of the created grid. 
*/
t_grid *create_grid_from_struct(int g0, int g1, t_species_parameters species_params) {
    assert(g0 > 0);
    assert(g1 > 0);
    /* We use a matrix of size n0 x n1 where 
       n0 = 2*width, n1 = 2*height to calculate
       connectivity with fast fourier transiform. */
    const int n1 = 2*g1, n0 = 2*g0;
    const int N = n0*n1;
    const int FFT_N = n0 * (n1/2 + 1);
    
    t_grid *g = malloc( sizeof(t_grid) );
    g->w = g1;
    g->h = g0;
    g->n0 = n0;
    g->n1 = n1;
    g->fft_n = FFT_N;

    g->species = species_params;

    g->occ_m = malloc( sizeof(unsigned char) * g0 * g1 );
    g->fit_m = malloc( sizeof(t_float) * g0 * g1);
    g->dir_m = malloc( sizeof(t_float) * g0 * g1);

    /* Allocate memory for the matrices. We use FFTW's malloc and free
       to ensure proper SIMD alignment. */

    g->work_m = (t_float*) fftwl_malloc( sizeof(t_float)*N );
    g->conn_m = (t_float*) fftwl_malloc( sizeof(t_float)*N );

    /* Allocate a buffer for fast random number generation.
       It should always be at least of size g0*g1. */
    int buf_len = dsfmt_get_min_array_size();
    if (buf_len < g0*g1) {
        buf_len = g0*g1;
        if (buf_len % 2 != 0) buf_len++; // It seems to buffer needs to be of even size
    } 
    g->rand_buf = (double*) fftwl_malloc( sizeof(double) * buf_len);
    g->rand_buf_len = buf_len; 

    g->occ_fft = (fftwl_complex*) fftwl_malloc( sizeof(fftwl_complex)*FFT_N);
    g->ker_fft = (fftwl_complex*) fftwl_malloc( sizeof(fftwl_complex)*FFT_N);

    int fail = (g->occ_fft == 0 || g->ker_fft == 0 || g->work_m == 0 || g->occ_m == 0 || g->fit_m == 0);
    assert(!fail);

    /* Initialize FFTW plans for this grid */
    init_grid_plans(g, FFTW_MEASURE);

    /* Initialize the arrays */

    for(int i = 0; i<g0*g1; i++) {
        g->occ_m[i] = 0;
        g->fit_m[i] = 1;
	g->dir_m[i] = 0;
    }

    for(int i = 0; i<n0*n1; i++) {
        g->work_m[i] = 0;
        g->conn_m[i] = 0;
    }

    for(int i = 0; i<FFT_N; i++) {
        g->occ_fft[i] = 0;
        g->ker_fft[i] = 0;
    }

    /* By default, regional stochasticity is not enabled. */
    g->use_regional_stochasticity = 0;
    g->stochasticity_grid = NULL;

    return g; 
}

/*
  Free all the memory of a t_grid allocated by create_grid().
*/
void free_grid(t_grid *g) {
    disable_regional_stochasticity(g);

    //free_grid_plans(g); /* Calling this after finalize will cause a segmentation fault */

    free(g->occ_m);
    free(g->fit_m);
    free(g->dir_m);
    
    fftwl_free(g->rand_buf);
    fftwl_free(g->work_m);
    fftwl_free(g->occ_fft);
    fftwl_free(g->conn_m);
    fftwl_free(g->ker_fft);
    
    g->fit_m = 0;
    g->occ_m = 0;
    g->work_m = 0;
    g->conn_m = 0;
    g->occ_fft = 0;
    g->ker_fft = 0;

    g->w = 0;
    g->h = 0;
    g->n0 = 0;
    g->n1 = 0;

    free(g);
}

void enable_regional_stochasticity(t_grid *g, t_stochasticity_parameters parameters) {
    if (g->use_regional_stochasticity) { 
        /* Free previous resources if necessary */
        disable_regional_stochasticity(g);
    }
    // NOTE: We're using work_m for storing the generated stochasticity
    g->use_regional_stochasticity = 1;
    int seed = dsfmt_genrand_uint32(&rng_state);
    t_float *arr = g->work_m;
    
    t_stochasticity_grid *s = create_stochasticity_grid(arr, g->h, g->w, parameters.weight_factor, parameters.mean, parameters.stddev, seed);
    g->stochasticity_grid = s; 
}

void disable_regional_stochasticity(t_grid *g) {
    g->use_regional_stochasticity = 0; 
    if (g->stochasticity_grid != NULL) {
        free_stochasticity_grid(g->stochasticity_grid);
        g->stochasticity_grid = NULL;
    }
}

/*
 A helper for assigning val into g->arr[y][x].  
*/
void grid_set(t_grid *g, t_float *arr, int y, int x, t_float val) {
    assert(g != NULL);
    assert(g->conn_m == arr || g->work_m == arr);
    assert(x >= 0 && y >= 0);
    arr[y*g->n1 + x] = val;
}

/*
 A helper for retrieving the value of g->arr[y][x].
*/
t_float grid_get(t_grid *g, const t_float *arr, int y, int x) {
    assert(g != NULL);
    assert(g->conn_m == arr || g->work_m == arr);
    assert(x >= 0 && y >= 0);

    return arr[y*g->n1 + x];
}

int is_occupied(t_grid *g, int y, int x) {
    assert( y >= 0 && y < g->h );
    assert( x >= 0 && x < g->w );

    return g->occ_m[y * g->w + x];
}

void set_occupied(t_grid *g, int y, int x, unsigned char val) {
    assert( y >= 0 && y < g->h );
    assert( x >= 0 && x < g->w );
    
    g->occ_m[y * g->w + x] = val;
}

double fitness_function(double phenotype, double habitat_quality, double bell_width) {
    double val = (phenotype - habitat_quality);
    double c = bell_width;
	double result = 0;
    //return exp(- (val*val) / (2*c*c) );
	if(val == 0 ){
		result = 1;
	}
	else{
		result = 0;
	}
	if(phenotype == 0.5){
		result = 0.5;
	}
	return result;
}

void calculate_fitness(t_grid *g, double *habitat_m, int h, int w) {
    assert(h == g->h && w == g->w);
    for(int y = 0; y < g->h; y++) {
        for(int x = 0; x < g->w; x++) {
            double F = 0; /* Fitness value */
            double h = habitat_m[y * g->w + x];
            /* Negative habitat value implies zero fitness! */
            if (h >= 0) {
                double z = g->species.phenotype;
                double b = g->species.fitness_bell_width;
                F = fitness_function(z,h,b);
            }
            assert( F >= 0 );
            set_fitness(g, y, x, F);
        }
    }

    /* Update: we need to update dir_m now */
    init_direction_factors(g);
}

void copy_occupancy_from(t_grid *g, unsigned char *arr, int h, int w) {
    assert(h == g->h && w == g->w);
    for(int y = 0; y < g->h; y++) {
        for(int x = 0; x < g->w; x++) {
            int index = y * g->w + x;
            g->occ_m[index] = arr[index];
        }
    }
}

/* Since the fitness matrix has different dimensions 
   we use this to access the matrix */
double get_fitness(t_grid *g, int y, int x) {
    assert(y >= 0 && y < g->h);
    assert(x >= 0 && x < g->w);
    return g->fit_m[y * g->w + x];
}

void set_fitness(t_grid *g, int y, int x, double val) {
    assert(y >= 0 && y < g->h);
    assert(x >= 0 && x < g->w);
    g->fit_m[y * g->w + x] = val;
}

double get_dir_factor(t_grid *g, int y, int x) {
    assert(y >= 0 && y < g->h);
    assert(x >= 0 && x < g->w);
    double f = g->dir_m[y * g->w + x];
    return f;
}

void give_occupancy_matrix(t_grid *g, unsigned char **view_array, int *h, int *w) {
    *view_array = g->occ_m;
    *h = g->h;
    *w = g->w;
}

void give_fitness_matrix(t_grid *g, double **view_array, int *h, int *w) {
#ifdef __USING_LONG_DOUBLES
    assert(0);
#endif

   *view_array = g->fit_m;
    *h = g->h;
    *w = g->w;
}

void give_connectivity_matrix(t_grid *g, double **view_array, int *h, int *w) {
 #ifdef __USING_LONG_DOUBLES
    assert(0);
#endif

   *view_array = g->conn_m;
    *h = g->n0;
    *w = g->n1;
}

void give_work_matrix(t_grid *g, double **view_array, int *h, int *w) {
 #ifdef __USING_LONG_DOUBLES
    assert(0);
#endif

   *view_array = g->work_m;
    *h = g->n0;
    *w = g->n1;
}


void give_conn_double_matrix(t_grid *g, double **view_array, int *h, int *w) {
    int N = g->n0 * g->n1;
    *h = g->n0;
    *w = g->n1;

    double *tmp = (double*) fftwl_malloc( sizeof(t_float)*N );
    for (int i = 0; i<N; i++) {
	tmp[i] = (double) g->conn_m[i];
    }   
    *view_array = tmp; // YES THIS LEAKS
}


double exp_kernel(double dist, double alpha) {
    /* The normalization factor for the kernel */
    double PI = acos(-1.0);
    double N = (alpha * alpha) / (2*PI);
    return N*exp(-dist*alpha);
}

void create_exp_dispersal_kernel(t_grid *g, double alpha) {
   assert(alpha <= 1);
   assert(0 <= alpha);
   for(int i = 0; i<g->n0; i++) {
        for(int j = 0; j<g->n1; j++) {
            int dy = MIN(i, g->n0 - i);
            int dx = MIN(j, g->n1 - j);
            double dist = sqrt(dx*dx + dy*dy);
            double val = exp_kernel(dist, alpha);
	    assert(val <= 1);
            assert(val >= 0);
            grid_set(g, g->work_m, i, j, val);
        }
    }
}

void create_gaussian_dispersal_kernel(t_grid *g, double alpha) {
    /* The normalization factor for the kernel */
    double PI = acos(-1.0);
    assert(alpha > 0 );

    double N = 1 / (2*PI*alpha);
    for(int i = 0; i<g->n0; i++) {
        for(int j = 0; j<g->n1; j++) {
            int dy = MIN(i, g->n0 - i);
            int dx = MIN(j, g->n1 - j);
			double val = N*exp(-(dx*dx + dy*dy)/(2*alpha));
			assert(val <= 1);
            assert(val >= 0);
            grid_set(g, g->work_m, i, j, val);
        }
    }
}

/*
    This function will perform the Fourier transformation 
    for the dispersal kernel. After building the kernel into work_m
    call this function (once!) *before* any calculating convolutions
*/
void setup_kernel(t_grid *g) {
    /* Write the kernel into work matrix */
#ifdef USE_GAUSSIAN_KERNEL
    create_gaussian_dispersal_kernel(g, g->species.alpha);
#else
    create_exp_dispersal_kernel(g, g->species.alpha);
#endif
    /* FFT it into ker_fft */
    fftwl_execute(g->p_ker_fwd);
    /* Cleanup the work matrix to use for computing occupancy and fitness values */
    for(int i = 0; i< g->n0 * g->n1; i++) {
        g->work_m[i] = 0;
    }
}

/*
 Initialize the dir_m matrix with normalization factors for each patch.
 */
void init_direction_factors(t_grid *g) {
    assert(g != NULL);
    const int FFT_N = g->fft_n;

    /* Update the work matrix with fitness values */
    for(int i = 0; i < g->n0; i++) {
        for(int j = 0; j < g->n1; j++) {
            double f = 0.0;
            if (i < g->h && j < g->w) {
                f = get_fitness(g, i, j);
		assert(f >= 0);
            } /* Otherwise just clear rest of the work matrix */
            grid_set(g, g->work_m, i, j, f); 
        }
    }

    /* Compute the convolution via FFT */
    fftwl_execute(g->p_occ_fwd);
    for(int i = 0; i < FFT_N; i++) {
       g->occ_fft[i] *= g->ker_fft[i]; 
    }
    fftwl_execute(g->p_occ_bwd);

    /* Normalize the results. */
    for(int i = 0; i < g->n0*g->n1; i++) {
        g->conn_m[i] /= g->n0*g->n1;
   }

    /* Save to dir_m */
    for(int y = 0; y < g->h; y++) {
        for(int x = 0; x < g->w; x++) {
	    double s = grid_get(g, g->conn_m, y, x);
	    /* DEBUGGING 
	    if (s > 1) {
		printf("%i,%i ",y,x);
		printf("s=%f\n",s);
		} 
	    */
	    assert(s <= 2);
	    if (s < 0) { // NEW: Fixme, we get negative sign for some reason from FFTW
		s = 0.0;
	    }
	    int index = y * g->w + x;
	    g->dir_m[index] = s;
	}
    }
#if 0
    // ------ THE SLOW WAY; for debugging ---------

    /* For each patch j, compute N_j = sum_i (f_ji * q_i) */
    double alpha = g->species.alpha;
    for(int y = 0; y < g->h; y++) {
        for(int x = 0; x < g->w; x++) {
	    int index = y * g->w + x;
	    double p = 0;
	    /* Sum over other patches */
	    for(int j = 0; j < g->h; j++) {
		for(int k = 0; k < g->w; k++) {
		    int dy = y-j;
		    int dx = x-k;
		    //if (dy == 0 && dx == 0) continue;
		    double dist = sqrt(dx*dx + dy*dy);
		    double val = exp_kernel(dist, alpha);
		    double F = get_fitness(g, j, k);
		    p += val*F;
		}
	    }
	    //double S = grid_get(g, g->conn_m, y, x);
	    assert(abs(p - g->dir_m[index]) == 0);
	    //printf("%f = %f\n", p, g->dir_m[index]);
	}
    }
#endif
}

void initialize_simulation(t_grid *g) {
    //printf("Setup kernel\n");
    setup_kernel(g);
    //printf("Done");
}

/*
    Updates the g->conn_m matrix with new connectivity values.
*/
void calculate_connectivity(t_grid *g) {
    assert(g != NULL);
    const int FFT_N = g->fft_n;

    /* Update the work matrix with new occupancy and fitness values */
    for(int i = 0; i < g->n0; i++) {
        for(int j = 0; j < g->n1; j++) {
            double f = 0.0;
            if (i < g->h && j < g->w) {
                if (is_occupied(g, i,j)) {
                    f = get_fitness(g, i,j);
 
		    /* Updated: directed dispersal coefficient goes here */
		    double n = get_dir_factor(g, i, j);
		    if (n <= 0) {
			f = 0;
		    } else { 
			f /= n;
		    }
		}
            } /* Otherwise just clear rest of the work matrix */
            grid_set(g, g->work_m, i, j, f); 
        }
    }

    /* The convolution corresponds to point-wise multiplication
       of the FFTd elements. */
    fftwl_execute(g->p_occ_fwd);
    for(int i = 0; i < FFT_N; i++) {
       g->occ_fft[i] *= g->ker_fft[i]; 
    }
    fftwl_execute(g->p_occ_bwd);

    /* Normalize the results. */
    for(int i = 0; i < g->n0*g->n1; i++) {
        g->conn_m[i] /= g->n0*g->n1;
    }
}

void iterated_connectivity(t_grid *g, int iterations) {
    assert(iterations > 0);
    calculate_connectivity(g); // First iteration
    // g->conn_m now contains S^i

    const int FFT_N = g->fft_n;

    for(int its = 0; its<iterations-1; its++) {
	/* Copy current connectivity to work matrix */
	for(int i = 0; i < g->n0; i++) {
	    for(int j = 0; j < g->n1; j++) {
		double S = 0.0;
		if (i < g->h && j < g->w) {
		    S = grid_get(g, g->conn_m, i, j);
		    S *= get_fitness(g, i, j);
		    // We may get some numerical problems here. 
			if (S >= 2) {
				printf("S=%f\n",S);
				} 
		    assert(S < 2);
		    assert(S > -1); 
		    
		    /* HACK to deal with rounding errors */
		    if (S < 0) {
			//printf("S=%f\n", S);
			S = 0.0; // Treat negative numbers as zeros
		    }
		    if (S > 1) S = 1.0;

		    /* Updated: directed dispersal coefficient goes here */
		    double n = get_dir_factor(g, i, j);
		    if (n <= 0) {
			S = 0;
		    } else { 
			S /= n;
		    }
		} /* Otherwise just clear rest of the work matrix */

		assert(S >= 0);
		grid_set(g, g->work_m, i, j, S); 
	    }
	}

	/* The convolution corresponds to point-wise multiplication of the FFTd elements. */
	fftwl_execute(g->p_occ_fwd);
	for(int i = 0; i < FFT_N; i++) {
	    g->occ_fft[i] *= g->ker_fft[i]; 
	}
	fftwl_execute(g->p_occ_bwd);

	/* Normalize the results. */
	for(int i = 0; i < g->n0*g->n1; i++) {
	    g->conn_m[i] /= g->n0*g->n1;
	    if (g->conn_m[i] < ZERO_THRESHOLD) g->conn_m[i] = 0;
	}
    }	
}

/* Return the new number of occupied cells. */
int simulate_step(t_grid *g) {
    assert(g != NULL);

    iterated_connectivity(g, g->species.connectivity_iterations); // NEW

    int occupied_cells = 0;

    /* Get enough random numbers to complete this iteration */
//    dsfmt_fill_array_open_close(&rng_state, g->rand_buf, g->rand_buf_len); 
    /* Note: This will probably slow down the simulation with small grid sizes (
             less than 500 cells */

    /* If we're using regional stochasticity, calculate the stochasticity matrix */
    if (g->use_regional_stochasticity) {
        /* This assumes work_m is free to be used */
        generate_stochasticity(g->stochasticity_grid);
    }

    /* Update the status of the population at each patch (ui, uj) */
    for(int ui = 0; ui < g->h; ui++) {
        for(int uj = 0; uj < g->w; uj++) {
            double F = get_fitness(g, ui, uj);
            if (F <= 0) { 
                /* Non-positive fitness implies extinction */
                set_occupied(g, ui, uj, 0);
                continue;
            }

            /* Get a uniform random number from the range [0, 1). */
            //double random = g->rand_buf[ui*g->w + uj];
            double random = dsfmt_genrand_open_close(&rng_state);
            /* IEEE574 floating point comparison with small integers is safe */
            char occupied = is_occupied(g, ui, uj);
            if (!occupied) {
                /* Check for colonization */
                double S = grid_get(g, g->conn_m, ui, uj);

		/* Update: Scale colonization with the quality of patch to get the directed dispersal effect (with dir_m factors) */
		S *= F;

                if (S < 0) { /* Floating point arithmetic may yield (small) negative numbers.. */ 
                    S = 0.0;
                }

		assert(0 <= S);
                double c = g->species.colonization_rate;
                assert(c >= 0);
                double pr = 1 - exp(-c*S);
            
                assert(0 <= pr);
                assert(pr <= 1);

                if (random < pr) {
                    set_occupied(g, ui, uj, 1);
                    occupied_cells++; // This cell got a new population
                }
            } else {
                /* Check for extinction */
                assert(F >= 0);

                /* Are we using regional stochasticity for fitness? */
                if (g->use_regional_stochasticity) {
                    t_stochasticity_grid *sg = g->stochasticity_grid;
                    t_float *reg_arr = sg->arr;
                    t_float y = reg_arr[ui*sg->width + uj];
		    assert(y >= 0);
                    if (y < 1.0) { /* case y > 1 won't affect F */
                        F *= y; 
                    }
                }

                double e = g->species.extinction_rate;
                double pr = 1 - exp(-e/F);
		assert(0 <= pr);
                assert(pr <= 1);
                
                if (random < pr) {
                    set_occupied(g, ui, uj, 0);
                } else {
                    occupied_cells++; // This cell did not die
                }
            }
        } 
    }
    return occupied_cells;

}

void many_steps(t_grid* grid, int steps, int *return_vals) {
    for(int i = 0; i<steps; i++) {
        return_vals[i] = 0;
    }

    for(int i = 0; i<steps; i++) {
        return_vals[i] = simulate_step(grid);
        if (return_vals[i] == 0) { /* Extinction, no need to simulate further */
            break;
        }
    }
}

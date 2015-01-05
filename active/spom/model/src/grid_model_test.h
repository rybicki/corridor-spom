#include "grid_model.h"

void init_grid_plans(t_grid *g, int flags);
void grid_set(t_grid *g, double *arr, int y, int x, double val);
double grid_get(t_grid *g, const double *arr, int y, int x);
void create_exp_dispersal_kernel(t_grid *g, double alpha);
void setup_kernel(t_grid *g);
void initialize_simulation(t_grid *g);
void calculate_connectivity(t_grid *g);
int simulate_step(t_grid *g);


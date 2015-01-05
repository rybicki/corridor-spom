/*
 Some simple tests to ensure that the FFT computations yield correct results. 
*/
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "grid_model_test.h"

const int MAX_H = 25;
const int MAX_W = 25;

const int MIN_CON_H = 1;
const int MIN_CON_W = 1;
const int MAX_CON_H = 15;
const int MAX_CON_W = 15;

double MAGIC_Z = 0.512345;
double MAGIC_C = 0.598765;
double MAGIC_ALPHA = 0.5;
double MAGIC_E = 1.0;
double MAGIC_B = 0.2;

const double EPSILON = 0.00000000001;

int MAGIC_SEED = 100;
const int REPEATS = 10;

const int LARGE_H = 1024;
const int LARGE_W = 1024;
const int LARGE_REPEATS = 5000;

void test_create() {
    printf("test_create()..");
    init_fg(1, 1);

    for(int i = 1; i<MAX_H; i++) {
        for(int j = 1; j<MAX_W; j++) {
            t_grid *g = create_grid(i,j, MAGIC_Z, MAGIC_C, MAGIC_E, MAGIC_ALPHA, MAGIC_B);

            assert(g->w == j && g->h == i);
            assert(g->species.phenotype == MAGIC_Z);
            assert(g->species.colonization_rate == MAGIC_C);
            assert(g->species.alpha = MAGIC_ALPHA);

            assert(g->conn_m != NULL);
            assert(g->work_m != NULL);

            free_grid(g);
        }
    }

    finalize_fg();
    printf("...OK.\n");
}

int test_connectivity(int H, int W) {
    //printf("test_connectivity(%i, %i)..", H, W);
    char occupancy[H][W];
    double fitness[H][W];
    double connectivity[H][W];

    init_fg(1,1);
    t_grid *g = create_grid(H, W, MAGIC_Z, MAGIC_C, MAGIC_E, MAGIC_ALPHA, MAGIC_B);

    int occupied = 0;

    for(int i = 0; i<H; i++) {
        for(int j = 0; j<W; j++) {
            occupancy[i][j] = (rand() % 2 == 0);
            fitness[i][j] = 1.0; // 1.0 / (2 + rand() % 10);
            connectivity[i][j] = 0.0;
            occupied += occupancy[i][j];

            set_occupied(g, i, j, occupancy[i][j]);
            set_fitness(g, i, j, fitness[i][j]);
        }
    }
    //printf("\tGenerated %i occupied patches. ", occupied);

    /* Brute force connectivity */
    for(int y = 0; y<H; y++) {
        for(int x = 0; x<W; x++) {
            double S = 0.0;
            for(int i = 0; i<H; i++) {
                for(int j = 0; j<W; j++) {
                    int dy = abs(y-i);
                    int dx = abs(x-j);
                    double dist = sqrt(dy*dy + dx*dx);
                    double e = exp(-MAGIC_ALPHA*dist);
                    S += e * occupancy[i][j] * fitness[i][j];
                }
            }
            double PI = acos(-1.0);
            double scale_factor = (MAGIC_ALPHA * MAGIC_ALPHA) / (2*PI);
            connectivity[y][x] = scale_factor* S;
        }
    }

    //printf("\tSetting up kernel\n");
    setup_kernel(g); 
    //printf("\tCalculating connectivity\n");
    calculate_connectivity(g);
    //printf("\tChecking values..\n");

    int fails = 0;
    for(int i = 0; i<H; i++) {
        for(int j =0; j<W; j++) {
            double res = grid_get(g, g->conn_m, i, j);
            if (fabs(connectivity[i][j] - res) > EPSILON) {
                printf("Wrong result (bruteforce vs conv): %f vs %f", connectivity[i][j], res);
                fails++;
            }
            //assert( connectivity[i][j] - res < EPSILON );
        }
    }

    if (fails == 0) {
        //printf("\tConnectivity computation OK. Testing step()...");
    } else {
        printf("\n\t%i failures in connectivity!\n", fails);
    }


    simulate_step(g);

    free_grid(g);
    finalize_fg();
    if (fails == 0) {
        //printf("...OK.\n");
    } else printf("..FAILED.\n");
    return fails;
}

void test_large_input(int h, int w) {
    printf("test_large_input(%i, %i)..", h, w);
    init_fg(8,1);
    t_grid *g = create_grid(h, w, MAGIC_Z, MAGIC_C, MAGIC_E, MAGIC_ALPHA, MAGIC_B);

    for(int i = 0; i<MAX_H; i++) {
        for(int j = 0; j<MAX_W; j++) {
            set_occupied(g, i, j, (rand() % 2 == 0));
            set_fitness(g, i, j, 1.0 / (2 + rand() % 10));
        }
    }
    initialize_simulation(g);
    printf("Input of size %i x %i generated and simulator initialized\n", h, w);
    
    clock_t start = clock();

    for(int i = 0; i<LARGE_REPEATS; i++) {
        simulate_step(g);
    }
    clock_t end = clock();

    printf("Elapsed: %f seconds on average\n", (double)(end - start) / CLOCKS_PER_SEC / LARGE_REPEATS);
    
    free_grid(g);
    finalize_fg();
    printf("..OK\n");
}

int main(int argc, const char **argv) {
    if (argc > 1) {
        printf("Trying to use given argument as a seed..\n");
        MAGIC_SEED = atoi(argv[1]);
    }
    printf("Running tests with\n\tMAGIC_SEED=%i\n", MAGIC_SEED);
    srand( MAGIC_SEED );
    MAGIC_Z = 1.0 / (1 + rand() % 20);
    MAGIC_C = 1.0 / (1 + rand() % 20);
    MAGIC_ALPHA = 1.0 / (1 + rand() % 20);

    printf("\tREPEATS=%i\n\tMAX_H=%i\n\tMAX_W=%i", REPEATS, MAX_H, MAX_W);
    printf("\n\tMAGIC_Z=%f\n\tMAGIC_C=%f\n\tMAGIC_E=%f\n\tMAGIC_ALPHA=%f\n\tMAGIC_B=%f", MAGIC_Z, MAGIC_C, MAGIC_E, MAGIC_ALPHA, MAGIC_B);
    printf("\n\tEPSILON=%.15f (Maximum difference in results)\n", EPSILON);
    printf("This might take a while.\n");
    test_create();

    /*printf("Testing connectivity");
    test_connectivity(3,3);
    test_connectivity(4,3);*/

    printf("Testing connectivity (%ix%i..%ix%i): ", MIN_CON_H, MAX_CON_H, MAX_CON_H, MAX_CON_W);
    int iterations = 0, failures = 0;
    for(int h = MIN_CON_H; h<MAX_CON_H; h++) {
        for(int w = MIN_CON_W; w<MAX_CON_W; w++) {
            for(int i = 0; i<REPEATS; i++) {
                srand( MAGIC_SEED + i);
                failures += test_connectivity(h,w);
                iterations++;
            }
        }
    }
    printf("\nRan %i test iterations. %i floating point errors encountered\n", iterations, failures);

    test_large_input(128, 128);
    test_large_input(512, 512);
    test_large_input(LARGE_H, LARGE_W);
    return 0;
}

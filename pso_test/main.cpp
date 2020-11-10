//https://github.com/kkentzo/pso

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> // for printf
#include "pso.h"
#include "pso.c"
//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================

double pso_sphere(double *vec, int dim, void *params) {

    double sum = 0;
    int i;
    for (i=0; i<dim; i++)
        sum += pow(vec[i], 2);
    return sum;
}
double pso_rosenbrock(double *vec, int dim, void *params)
{
    double sum = 0;
    int i;
    for (i=0; i<dim-1; i++)
        sum += 100 * pow((vec[i+1] - pow(vec[i], 2)), 2) +	\
            pow((1 - vec[i]), 2);
    return - sum;
}

double pso_griewank(double *vec, int dim, void *params)
{
    double sum = 0.;
    double prod = 1.;
    for (int i=0; i<dim;i++) {
        sum += pow(vec[i], 2);
        prod *= cos(vec[i] / sqrt(i+1));
    }
    return - (sum / 4000 - prod + 1);
}
/********** Griewank function N Dimensions ***********/
//            size_t dim_in = 20;
//            auto xx = x;
//            //transfer interval from [0, 1] to [-10, 10] [-5, 5] [-6,6]
//            for (int i = 0; i < dim_in; i++)
//              xx[i] = 10. * x[i] - 5.;
//            double sum = 0.0;
//            double f = 1.0;
//            for (size_t i = 0; i < dim_in; ++i)
//            {
//              sum += xx[i] * xx[i]/4000;
//              f = f * std::cos(xx[i]/std::sqrt(i+1));
//            }
//            double obj = sum - f + 1;
//            return tools::make_vector(-obj); //maximum = 0 with (0, 0, 0, 0);

/********** Rastrigin function N Dimensions ***********/
double pso_rastrigin(double *vec, int dim, void *params)
{
//  auto xx = x;
////transfer interval from [0, 1] to [-3, 3] [-2, 2]
//  for (int i = 0; i < dim; i++)
//    xx[i] = 4. * x[i] - 2.;
  double f = 10. * dim;
  for (int i = 0; i < dim; ++i)
    f += vec[i] * vec[i] - 10. * std::cos(2 * M_PI * vec[i]);
  return -f; //maximum = 0 with (0, 0, 0, 0);
}

///********** 2# SCHWEFEL function N Dimensions **************/
double pso_schwefel(double *vec, int dim, void *params)
{
//  size_t dim_in = 20;
//  auto xx = x;
//  // transfer interval from [0, 1] to [-500, 500] [250, 500] [100,500]
//  for (int i = 0; i < dim_in; i++)
//  {
//  xx[i] = 300. * x[i] + 200.;
//  }
  double sum = 0.;
  for (size_t i = 0; i < dim; i++)
  {
    sum = sum + vec[i] * sin(sqrt(abs(vec[i])));
  }
  double obj = 418.9829 * dim - sum;
  return -obj; //maximum = 0 with (420.9687, ...,420.9687)
}

//==============================================================

int main(int argc, char **argv)
{
  pso_settings_t *settings = NULL;
  pso_obj_fun_t obj_fun = pso_rastrigin; //schwefel(200,500),griewank(-5,5),rastrigin(-2,2)

  settings = pso_settings_new(20, -2, 2);
//  printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",settings->dim, settings->size);

  // set some general PSO settings
  settings->goal = 1e-5;
  settings->size = 10;
  settings->nhood_strategy = PSO_NHOOD_RING;
  settings->nhood_size = 5;
  settings->w_strategy = PSO_W_LIN_DEC;

  // initialize GBEST solution
  pso_result_t solution;
  // allocate memory for the best position buffer
  solution.gbest = (double *)malloc(settings->dim * sizeof(double));

  int iters = 10;
  for (int run = 0; run < iters; run++) {
    // run optimization algorithm
    pso_solve(obj_fun, NULL, &solution, settings);
  }

  // free the gbest buffer
  free(solution.gbest);

  // free the settings
  pso_settings_free(settings);

  return 0;
}
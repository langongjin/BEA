// please see the explanation in the documentation
// http://www.resibots.eu/limbo

#include <iostream>

// you can also include <limbo/limbo.hpp> but it will slow down the compilation
#include "limbo/bayes_opt/boptimizer.hpp"

#include <sys/time.h>

using namespace limbo;

using vec_t = Eigen::VectorXd;
using mat_t = Eigen::MatrixXd;

#define ASSIGNMENT_OP <<

#define USE_NLOPT;

// Here we define the parameters
struct Params {
    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer {
    };
//    struct opt_nloptnograd : public defaults::opt_nloptnograd { };
// depending on which internal optimizer we use, we need to import different parameters
#ifdef USE_NLOPT
    struct opt_nloptnograd : public defaults::opt_nloptnograd {
    };
#elif defined(USE_LIBCMAES)
    struct opt_cmaes : public defaults::opt_cmaes {
    };
#else
    struct opt_gridsearch : public defaults::opt_gridsearch {
    };
#endif

  struct kernel
      : public defaults::kernel
  {
    BO_PARAM(double, noise, 0.000001);
    BO_PARAM(bool, optimize_noise, false);
  };

  struct bayes_opt_bobase
      : public defaults::bayes_opt_bobase
  {
    BO_PARAM(bool, stats_enabled, true);
    BO_PARAM(bool, bounded, true); //false
  };

//  struct kernel_maternthreehalves
//      : public defaults::kernel_maternthreehalves
//  {
//    BO_PARAM(double, sigma_sq, 0.001);
//    BO_PARAM(double, l, 0.1);
//  };

  struct kernel_maternfivehalves : public defaults::kernel_maternfivehalves
  {
    BO_PARAM(double, sigma_sq, 0.01);
    BO_PARAM(double, l, 0.5); //rastrigin:0.1,tang:0.2,schew:0.5,griewank:0.1
  };

  struct init_randomsampling : public defaults::init_randomsampling
  {
    BO_PARAM(int, samples, 10);
  };

  struct stop_maxiterations : public defaults::stop_maxiterations
  {
    //1s -> 19*10 (180+10), 10s -> 25*10 (240+10), 100s -> 28 (270+10)
    BO_PARAM(int, iterations, 241); //stop after iterations (10+10)
  };

  // we use the default parameters for acqui_ucb
  struct acqui_ucb : public defaults::acqui_ucb
  {
    //UCB(x) = \mu(x) + \alpha \sigma(x). high alpha have high exploration
    //iterations is high, alpha can be low for high accuracy in enoughiterations.
    // In contrast, the low iterations should have high alpha for high
    // searching in limited iterations, which guarantee to optimal.
    BO_PARAM(double, alpha, 5); // it depends on the mu(x) and sigma(x)
  };

  // we use the default parameters for acqui_gpucb
  struct acqui_gpucb : public limbo::defaults::acqui_gpucb
  {
    //UCB(x) = \mu(x) + \kappa \sigma(x). high alpha have high exploration
    BO_PARAM(double, delta, 5); // default alpha = 0.1
  };
};


// Here we define the evaluation function
struct eval_func {
    // number of input dimension (x.size())
    BO_PARAM(size_t, dim_in, 20); //change i at row 156 chromosome.hpp
    // number of dimenions of the result (res.size())
    BO_PARAM(size_t, dim_out, 1);

    // the function to be optimized
    Eigen::VectorXd operator()(const Eigen::VectorXd& x) const
    {
      /********** 1# ACKLEY function N Dimensions **************/
//      size_t dim_in = 2;
//      auto xx = x;
//      // transfer interval from [0, 1] to [-32.768, 32.768]
//      for (int i = 0; i < dim_in; i++)
//      {
//        xx[i] = 65.536 * x[i] - 32.768;
//      }
//      const double a = 20.;
//      const double b = 0.2;
//      const double c = 2 * M_PI;
//      double sum1 = 0.;
//      double sum2 = 0.;
//      for (size_t i = 0; i < dim_in; i++)
//      {
//        sum1 = sum1 + xx[i] * xx[i];
//        sum2 = sum2 + std::cos(c * xx[i]);
//      }
//      double term1 = -a * std::exp(-b * std::sqrt(sum1 / dim_in));
//      double term2 = -std::exp(sum2 / dim_in);
//      double obj = term1 + term2 + a + std::exp(1);
//      return tools::make_vector(-obj); //max = 0, at (0,...,0)

      /********** 2# SCHWEFEL function N Dimensions **************/
      size_t dim_in = 20;
      auto xx = x;
      // transfer interval from [0, 1] to [-500, 500] [250, 500] [100,500]
      for (int i = 0; i < dim_in; i++)
      {
        xx[i] = 300. * x[i] + 200.;
      }
      double sum = 0.;
      for (size_t i = 0; i < dim_in; i++)
      {
        sum = sum + xx[i] * sin(sqrt(abs(xx[i])));
      }
      double obj = 418.9829 * dim_in - sum;
      return tools::make_vector(-obj); //maximum = 0 with (420.9687, ...,420.9687)

      /********** 3# Ellipsoid function N Dimensions **************/
//            size_t dim_in = 5;
//            double inner = 0., outer = 0.;
//            for (size_t i = 0; i < dim_in; ++i)
//            {
//              for(size_t j = 0; j < i; j++)
//              {
//                inner = inner + std::pow((131.072 * x[j] - 65.536), 2); //(-65.536, 65.536)
//              }
//              outer = outer + inner;
//            }
//            return tools::make_vector(-outer); //maximum = 0 at (0, ..., 0)

      /********** 4# Sphere function N Dimensions **************/
//            size_t dim_in = 5;
//            double inner = 0.;
//            for (size_t i = 0; i < dim_in; ++i)
//            {
//              inner = inner + std::pow((10. * x[i] - 5.), 2);
//            }
//            return tools::make_vector(-inner); //max = 0 with (0, 0)

      /********** 5# Rosenbrock  function N Dimensions **************/
//            size_t dim_in = 5;
//            auto xx = x;
//            // transfer interval from [0, 1] to [-2.048, 2.048]
//            for (int i = 0; i < dim_in; i++)
//              xx[i] = 4.096 * x[i] - 2.048;
//            double sum = 0.;
//            double term = 0.;
//            double xnext = 0.;
//            for(size_t i = 0; i < (dim_in - 1); i++)
//            {
//              xnext = xx[i + 1];
//              term = 100. * std::pow((xnext - xx[i] * xx[i]), 2.0) + std::pow((xx[i] - 1), 2.0);
//              sum = sum + term;
//            }
//            double obj = 0.001 * sum;
//      //      double obj = (sum - 382700)/375500.; //rescale
//            return tools::make_vector(-obj); //maximum = 0 with (1,...,1)

      /********* 6# Michalewicz function N = 2/5/10 Dimensions **********/
//            size_t dim_in = 5; //todo not good results
//            auto xx = x;
//            // transfer interval from [0, 1] to [0, pi]
//            for (int i = 0; i < dim_in; i++)
//              xx[i] = M_PI * x[i];
//            double sum = 0.;
//            double term = 0.;
//            double m = 10.;
//            for(size_t i = 0; i < dim_in; i++)
//            {
//              term = std::sin(xx[i]) * std::pow(std::sin(i * xx[i] * xx[i]/M_PI), 2 * m);
//              sum = sum + term;
//            }
//            double obj = sum;
//            return tools::make_vector(obj); //max= -1.8013(2D) at (2.20,1.57)/-4.687658(5D)/-9.66015(10D)

      /********** 7# StyblinskiTang function N Dimensions ****************/
//            size_t dim_in = 20;
//            auto xx = x;
//            // transfer interval from [0, 1] to [-5, 5] [-4,4]
//            for (int i = 0; i < dim_in; i++)
//              xx[i] = 10. * x[i] - 5.;
//            double sum = 0.;
//            double term;
//            for(size_t i = 0; i < dim_in; i++)
//            {
//              term = std::pow(xx[i], 4.0) - 16 * xx[i] * xx[i] + 5 * xx[i];
//              sum = sum + term;
//            }
//            double obj = sum/2.0;
//            //max= 39.16599 * d, (5D:195.82995) at (-2.903534,...,-2.903534)
//            return tools::make_vector(-obj);
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
//            size_t dim_in = 20;
//            auto xx = x;
//            //transfer interval from [0, 1] to [-3, 3] [-2, 2]
//            for (int i = 0; i < dim_in; i++)
//              xx[i] = 4. * x[i] - 2.;
//            double f = 10. * dim_in;
//            for (size_t i = 0; i < dim_in; ++i)
//              f += xx[i] * xx[i] - 10. * std::cos(2 * M_PI * xx[i]);
//            return tools::make_vector(-f); //maximum = 0 with (0, 0, 0, 0);
    }
};


//s.member等价于(&s)->member。反之：(*p).member则是(&(*p))->member，即是p->member

int main()
{
  struct timeval timeStart, timeEnd;
  double timeDiff;
  gettimeofday(&timeStart, NULL);

  /*********** configuration ************/
  //  using Kernel_t = kernel::SquaredExpARD<Params>;
  using Kernel_t = kernel::MaternFiveHalves <Params>;
  using Mean_t = mean::Data<Params>;
  using GP_t = model::GP<Params, Kernel_t, Mean_t>;
  using Acqui_t = acqui::GP_UCB<Params, GP_t>;

  int iters = 1;
  for (int run = 0; run < iters; run++)
  {
    // we use the default acquisition function / model / stat / etc.
//    bayes_opt::BOptimizer< Params > boptimizer;
    bayes_opt::BOptimizer<Params, modelfun<GP_t>, acquifun<Acqui_t>> boptimizer;
    // run the evaluation
    boptimizer.optimize(eval_func());
  }
  // the best sample found
  //    std::cout << "Best sample: " << boptimizer.best_sample().transpose() << " - Best observation: " << boptimizer.best_observation()(0) << std::endl;
  gettimeofday(&timeEnd, NULL);

  timeDiff = 1000000 * (timeEnd.tv_sec - timeStart.tv_sec) + timeEnd.tv_usec - timeStart.tv_usec; //tv_sec: value of second, tv_usec: value of microsecond
  timeDiff /= 1000;
//  std::cout << "timeDiff : " << timeDiff << std::endl;

  return 0;
}

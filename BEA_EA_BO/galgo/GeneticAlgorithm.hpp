//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

#include <sys/time.h>

namespace galgo {

enum class MutationType {
    MutationSPM,
    MutationBDM,
    MutationUNM,
    MutationGAM_UncorrelatedOneStepSizeFixed,
    MutationGAM_UncorrelatedOneStepSizeBoundary,
    MutationGAM_UncorrelatedNStepSize,
    MutationGAM_UncorrelatedNStepSizeBoundary,
    MutationGAM_sigma_adapting_per_generation,
    MutationGAM_sigma_adapting_per_mutation,
};

template <typename T>
struct MutationInfo
{
    MutationInfo<T>() :
        _type(MutationType::MutationSPM),
        _sigma(1.0),
        _ratio_boundary(1.0 / 10.0), //default: 1.0/6.0
        _sigma_lowest(0.001) // default: 0.0001
    {
    }

    MutationType _type;
    double _sigma;
    double _ratio_boundary;
    double _sigma_lowest;
};

template <typename ParamTYPE>
struct ConfigInfo
{
    ConfigInfo() : mutinfo()
    {
        // DEFAULTS
        covrate = 0.1;  // cross-over rate
        mutrate = 0.8;  // mutation rate usually is 1.0 for real-valued
        SP = 1.5;
        tolerance = 0.0;
        recombination_ratio = 0.60; //Real Valued crossover ratio, can't be 0.5 because 0.5 will generate two same offsprings after Xover

        elitpop = 1;
        tntsize = 2; // k-tournament size k=2/4, higher value higher pressure
        precision = 10;

        Objective = nullptr;
        Selection = TNT; // RWS, TNT, RNK, RSP, TRS
        CrossOver = RealValuedSimpleArithmeticRecombination; //
        //Mutation = SPM; // derived from by mutinfo._type
        Adaptation = nullptr;
        Constraint = nullptr;
        FixedValue = nullptr;
        StopCondition = nullptr;

        nbgen = 125; //The num of gens for EA, 131+19, 125+25, 122+28
        popsize = 10;
        output = false;
    }

    MutationInfo<ParamTYPE> mutinfo;

    double covrate;
    double mutrate;
    double SP;
    double tolerance;
    double recombination_ratio;

    int elitpop;
    //int matsize; // set to popsize when ga is constructed, maybe change by ga.matsize = ... after constructor and before ga.run()
    int tntsize;
    int precision;

    std::vector<double> (*Objective)(const std::vector<ParamTYPE>&);
    void (*Selection)(Population<ParamTYPE>&);
    void (*CrossOver)(const Population<ParamTYPE>&, CHR<ParamTYPE>&, CHR<ParamTYPE>&);
    void (*Mutation)(CHR<ParamTYPE>&);
    void (*Adaptation)(Population<ParamTYPE>&) = nullptr;
    std::vector<double> (*Constraint)(const std::vector<ParamTYPE>&);
    void (*FixedValue)(Population<ParamTYPE>&, int k);
    bool (*StopCondition)(galgo::GeneticAlgorithm<ParamTYPE>&);

    std::vector<bool>       force_value_flag;
    std::vector<ParamTYPE>  force_value;

    int nbgen;
    int popsize;
    bool output;
};


/*-------------------------------------------------------------------------------------------------*/

template <typename T>
class GeneticAlgorithm
{
   template <typename K> friend class Population;
   template <typename K> friend class Chromosome;
   template <typename Z>  using FuncKT = std::vector<double>(*)(const std::vector<Z>&);

protected:
   Population<T> pop;             // population of chromosomes
   std::vector<PAR<T>> param;     // parameter(s)
   std::vector<T> lowerBound;     // parameter(s) lower bound
   std::vector<T> upperBound;     // parameter(s) upper bound
   std::vector<T> initialSet;     // initial set of parameter(s)
   std::vector<int> idx;          // indexes for chromosome breakdown

public:
   // objective function pointer
   FuncKT<T> Objective;

   // selection method initialized to roulette wheel selection
   void (*Selection)(Population<T>&) = TNT;

   // cross-over method initialized to 1-point cross-over
   void(*CrossOver)(const Population<T>&, CHR<T>&, CHR<T>&) = RealValuedSimpleArithmeticRecombination;

   // mutation method initialized to single-point mutation
   void (*Mutation)(CHR<T>&) = GAM_UncorrelatedNStepSizeBoundary;

   // adaptation to constraint(s) method
   void (*Adaptation)(Population<T>&) = nullptr;

   // constraint(s)
   std::vector<double> (*Constraint)(const std::vector<T>&) = nullptr;

   MutationInfo<T> mutinfo;

   double covrate = 0.5;   // cross-over rate
   double mutrate = 0.8;   // mutation rate
   double SP = 1.5;        // selective pressure for RSP selection method
   double tolerance = 0.0; // terminal condition (inactive if equal to zero)
   double recombination_ratio = 0.60; // Real Valued crossover ratio number = popsize * recombination_ratio

   int elitpop = 1;   // elit population size
   int matsize;       // mating pool size, set to popsize by default
   int tntsize = 2;  // k-tournament size k=2/4, higher value higher pressure
   int precision = 10; // precision for outputting results
   bool output;   // control if results must be outputted

   // Prototype to set fixed value of parameters while evolving
   void (*FixedValue)(Population<T>&, int k) = nullptr;
   std::vector<bool> force_value_flag;
   std::vector<T> force_value;

   bool (*StopCondition)(GeneticAlgorithm<T>&) = nullptr;

   // constructor
   template <int...N> GeneticAlgorithm(FuncKT<T> objective, int popsize, int nbgen, bool output, MutationInfo<T> mutinfo, const Parameter<T,N>&...args);
   template <int...N> GeneticAlgorithm(const ConfigInfo<T>& config, const Parameter<T, N>&...args);

   // run genetic algorithm
   void run();

   // return best chromosome
   const CHR<T>& result() const;

   // print results for each new generation
   void print(bool force  =  false) const;

   Population<T>& get_pop() { return pop; }

protected:
   void setMutation(const MutationInfo<T>& mt)
   {
       mutinfo = mt;
       if (mt._type == MutationType::MutationSPM)       { Mutation = SPM; }
       else if (mt._type == MutationType::MutationBDM)  { Mutation = BDM; }
       else if (mt._type == MutationType::MutationUNM)  { Mutation = UNM; }
       else if (mt._type == MutationType::MutationGAM_UncorrelatedOneStepSizeFixed)     { Mutation = GAM_UncorrelatedOneStepSizeFixed; }
       else if (mt._type == MutationType::MutationGAM_UncorrelatedOneStepSizeBoundary)  { Mutation = GAM_UncorrelatedOneStepSizeBoundary; }
       else if (mt._type == MutationType::MutationGAM_UncorrelatedNStepSize)            { Mutation = GAM_UncorrelatedNStepSize; }
       else if (mt._type == MutationType::MutationGAM_UncorrelatedNStepSizeBoundary)    { Mutation = GAM_UncorrelatedNStepSizeBoundary; }
       else if (mt._type == MutationType::MutationGAM_sigma_adapting_per_generation)    { Mutation = GAM_sigma_adapting_per_generation; }
       else if (mt._type == MutationType::MutationGAM_sigma_adapting_per_mutation)      { Mutation = GAM_sigma_adapting_per_mutation; }
       else Mutation = SPM;
   }

   GeneticAlgorithm(const ConfigInfo<T>& config); // No parameters set

   int nbbit;     // total number of bits per chromosome
   int nbgen;     // number of generations
   int nogen = 0; // numero of generation
   int nbparam;   // number of parameters to be estimated
   int popsize;   // population size

   // end of recursion for initializing parameter(s) data
   template <int I = 0, int...N>
   typename std::enable_if<I == sizeof...(N), void>::type init(const TUP<T,N...>&);

   // recursion for initializing parameter(s) data
   template <int I = 0, int...N>
   typename std::enable_if<I < sizeof...(N), void>::type init(const TUP<T,N...>&);

   // check inputs validity
   void check() const ;

   void init_from_config(const ConfigInfo<T>& config);
};

/*-------------------------------------------------------------------------------------------------*/
//template <typename T, int PARAM_NBIT>
//class GeneticAlgorithmN : public GeneticAlgorithm<T>
//{
//public:
//    // constructor
//    GeneticAlgorithmN(const ConfigInfo<T>& config,
//        std::vector<T>& _lowerBound,
//        std::vector<T>& _upperBound,
//        std::vector<T>& _initialSet);
//};

/*-------------------------------------------------------------------------------------------------*/
template <typename T>
void GeneticAlgorithm<T>::init_from_config(const ConfigInfo<T>& config)
{
    setMutation(config.mutinfo); // Mutation is set here

    Objective = config.Objective;
    Selection = config.Selection;
    CrossOver = config.CrossOver;
    Adaptation = config.Adaptation;
    Constraint = config.Constraint;
    FixedValue = config.FixedValue;
    StopCondition = config.StopCondition;

    covrate = config.covrate;
    mutrate = config.mutrate;
    SP = config.SP;
    tolerance = config.tolerance;
    recombination_ratio = config.recombination_ratio;

    elitpop = config.elitpop;
    tntsize = config.tntsize;
    precision = config.precision;

    force_value_flag = config.force_value_flag;
    force_value = config.force_value;

    nbgen = config.nbgen;
    popsize = config.popsize;
    matsize = popsize;      // matsize default to popsize
    output = config.output;

    nogen = 0;
}

template <typename T> template <int...N>
GeneticAlgorithm<T>::GeneticAlgorithm(const ConfigInfo<T>& config, const Parameter<T, N>&...args)
{
    init_from_config(config);

    nbbit = sum(N...);
    nbparam = sizeof...(N);
    TUP<T, N...> tp(args...);
    init(tp);
}

template <typename T>
GeneticAlgorithm<T>::GeneticAlgorithm(const ConfigInfo<T>& config)
{
    init_from_config(config);
}

// constructor
template <typename T> template <int...N>
GeneticAlgorithm<T>::GeneticAlgorithm(FuncKT<T> objective, int popsize, int nbgen, bool output, MutationInfo<T> mutinfo, const Parameter<T,N>&...args)
{
   setMutation(mutinfo);
   Objective = objective;

   // getting total number of bits per chromosome
   nbbit = sum(N...);
   nbgen = nbgen;

   // getting number of parameters in the pack
   nbparam = sizeof...(N);

   popsize = popsize;
   matsize = popsize;
   output = output;

   // unpacking parameter pack in tuple
   TUP<T,N...> tp(args...);

   // initializing parameter(s) data
   init(tp);
}

//// constructor
//template <typename T, int PARAM_NBIT>
//GeneticAlgorithmN<T, PARAM_NBIT>::GeneticAlgorithmN(const ConfigInfo<T>& _config, std::vector<T>& _lowerBound, std::vector<T>& _upperBound, std::vector<T>& _initialSet)
//    : GeneticAlgorithm<T>(_config)
//{
//    for (int i = 0; i<_lowerBound.size(); i++)
//    {
//        std::vector<T> w;
//        w.push_back(_lowerBound[i]); w.push_back(_upperBound[i]); w.push_back(_initialSet[i]);
//        Parameter<T, PARAM_NBIT> p(w);
//        param.emplace_back(new decltype(p)(p));
//
//        if (i == 0) idx.push_back(0);
//        else idx.push_back(idx[i - 1] + PARAM_NBIT);
//    }
//    lowerBound = _lowerBound;
//    upperBound = _upperBound;
//    initialSet = _initialSet;
//
//    nbbit = (int)_lowerBound.size()*PARAM_NBIT;
//    nbparam = (int)_lowerBound.size();
//};

/*-------------------------------------------------------------------------------------------------*/

// end of recursion for initializing parameter(s) data
template <typename T> template <int I, int...N>
inline typename std::enable_if<I == sizeof...(N), void>::type
GeneticAlgorithm<T>::init(const TUP<T,N...>& tp) {}

// recursion for initializing parameter(s) data
template <typename T> template <int I, int...N>
inline typename std::enable_if<I < sizeof...(N), void>::type
GeneticAlgorithm<T>::init(const TUP<T,N...>& tp)
{
   // getting Ith parameter in tuple
   auto par = std::get<I>(tp);

   // getting Ith parameter initial data
   const std::vector<T>& data = par.getData();

   // copying parameter data
   param.emplace_back(new decltype(par)(par));

   lowerBound.push_back(data[0]);
   upperBound.push_back(data[1]);

   // if parameter has initial value
   if (data.size() > 2) {
      initialSet.push_back(data[2]);
   }
   // setting indexes for chromosome breakdown
   if (I == 0) {
      idx.push_back(0);
   } else {
      idx.push_back(idx[I - 1] + par.size());
   }
   // recursing
   init<I + 1>(tp);
}

/*-------------------------------------------------------------------------------------------------*/

// check inputs validity
template <typename T>
void GeneticAlgorithm<T>::check() const
{
   if (!initialSet.empty()) {
      for (int i = 0; i < nbparam; ++i) {
         if (initialSet[i] < lowerBound[i] || initialSet[i] > upperBound[i]) {
            throw std::invalid_argument("Error: in class galgo::Parameter<T,N>, initial parameter value cannot be outside the parameter boundaries, please choose a value between its lower and upper bounds.");
         }
      }
      if (initialSet.size() != (unsigned)nbparam) {
         throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, initial set of parameters does not have the same dimension than the number of parameters, please adjust.");
      }
   }
   if (SP < 1.0 || SP > 2.0) {
      throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, selective pressure (SP) cannot be outside [1.0,2.0], please choose a real value within this interval.");
   }
   if (elitpop > popsize || elitpop < 0) {
      throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, elit population (elitpop) cannot outside [0,popsize], please choose an integral value within this interval.");
   }
   if (covrate < 0.0 || covrate > 1.0) {
      throw std::invalid_argument("Error: in class galgo::GeneticAlgorithm<T>, cross-over rate (covrate) cannot outside [0.0,1.0], please choose a real value within this interval.");
   }
}

/*-------------------------------------------------------------------------------------------------*/

// run genetic algorithm
template <typename T>
void GeneticAlgorithm<T>::run()
{
   // checking inputs validity
   check();

   // setting adaptation method to default if needed
   if (Constraint != nullptr && Adaptation == nullptr) {
      Adaptation = DAC;
   }

   // initializing population
   pop = Population<T>(*this);

   if (output) {
      std::cout << "\n Running Genetic Algorithm...\n";
      std::cout << " ----------------------------\n";
   }

   // creating population
   pop.creation();

   // initializing best result and previous best result
   double bestResult = pop(0)->getTotal();
   double prevBestResult = bestResult;

   // outputting results
   if (output) print();

  struct timeval timeStartEA, timeEndEA;
  double timeDiffEA;

  std::ofstream ctime;
  ctime.open("../ctime.txt", std::ios::app);

   // starting population evolution
   for (nogen = 1; nogen <= nbgen; ++nogen)
   {
     gettimeofday(&timeStartEA,NULL);

     // evolving population
     pop.evolution(); //nogen = 0

     // initializing best result and previous best result
     bestResult = pop(0)->getTotal();

     // outputting results
     if (output)
     {
       print();
     }

     // checking convergence
     if (tolerance != 0.0)
     {
       if (fabs(bestResult - prevBestResult) < fabs(tolerance))
       {
         break;
       }
       prevBestResult = bestResult;
     }

     if (StopCondition != nullptr)
     {
       if (StopCondition(*this) == true)
       {
         break;
       }
     }
     gettimeofday(&timeEndEA,NULL);
     //tv_sec: value of second, tv_usec: value of microsecond
     timeDiffEA = 1000000 * (timeEndEA.tv_sec - timeStartEA.tv_sec)
                  + timeEndEA.tv_usec - timeStartEA.tv_usec;
     timeDiffEA/=1000;
     //std::cout << "timeDiffEA: " << timeDiffEA << std::endl;

     ctime << std::fixed << timeDiffEA << std::endl;
   }

   // outputting contraint value
   if (Constraint != nullptr)
   {
      // getting best parameter(s) constraint value(s)
      std::vector<double> cst = pop(0)->getConstraint();
      if (output) {
         std::cout << "\n Constraint(s)\n";
         std::cout << " -------------\n";
         for (unsigned i = 0; i < cst.size(); ++i) {
            std::cout << " C";

            if (nbparam > 1) {
               std::cout << std::to_string(i + 1);
            }
            std::cout << "(x) = " << std::setw(6) << std::fixed << std::setprecision(precision) << cst[i] << "\n";
         }
         std::cout << "\n";
      }
   }
}

/*-------------------------------------------------------------------------------------------------*/

// return best chromosome
template <typename T>
inline const CHR<T>& GeneticAlgorithm<T>::result() const
{
   return pop(0);
}

/*-------------------------------------------------------------------------------------------------*/

// print results for each new generation
template <typename T>
void GeneticAlgorithm<T>::print(bool force) const
{
  // getting best parameter(s) from best chromosome
  std::vector< T > bestParam = pop(0)->getParam();
  std::vector< double > bestResult = pop(0)->getResult();

  std::cout << " Generation = " << std::setw(std::to_string(nbgen).size())
            << nogen << " |";
  for (int i = 0; i < nbparam; ++i)
  {
    std::cout << " X";
    if (nbparam > 1)
    {
      std::cout << std::to_string(i + 1);
    }
    std::cout << " = " << std::setw(2 + precision) << std::fixed
              << std::setprecision(precision) << bestParam[i] << " |";
  }
  for (unsigned i = 0; i < bestResult.size(); ++i)
  {
    std::cout << " F";
    if (bestResult.size() > 1)
    {
      std::cout << std::to_string(i + 1);
    }
    std::cout << "(x) = " << std::setw(12) << std::fixed
              << std::setprecision(precision) << bestResult[i];

    if(nogen not_eq 0)
    {
      //save the best value
      std::ofstream value;
      value.open("../value.txt", std::ios::app);
      value << std::fixed << bestResult[i] << std::endl;
    }

    if (i < bestResult.size() - 1)
    {
      std::cout << " |";
    }
    else
    {
      std::cout << "\n";
    }
  }
  std::cout << std::endl;
}

//=================================================================================================

}

#endif

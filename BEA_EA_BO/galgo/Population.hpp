//=================================================================================================
//                  Copyright (C) 2018 Alain Lanthier - All Rights Reserved  
//                  License: MIT License    See LICENSE.md for the full license.
//                  Original code 2017 Olivier Mallet (MIT License)              
//=================================================================================================

#ifndef POPULATION_HPP
#define POPULATION_HPP

#include <iostream>
#include <fstream>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/core.hpp"
#include <sys/time.h>
#include <math.h>

namespace galgo
{

  template <typename T> class Population
  {
    public:
    // nullary constructor
    Population()
    {
    }

    // constructor
    Population(const GeneticAlgorithm <T> &ga);

    // create a population of chromosomes
    void creation();

    // evolve population, get next generation
    void evolution();

    // access element in current population at position pos
    const CHR <T> &operator()(int pos) const;

    // access element in mating population at position pos
    const CHR <T> &operator[](int pos) const;

    // return iterator to current population beginning
    typename std::vector< CHR < T>>::

    iterator begin();

    // return const iterator to current population beginning
    typename std::vector< CHR < T>>::

    const_iterator cbegin() const;

    // return iterator to current population ending
    typename std::vector< CHR < T>>::

    iterator end();

    // return const iterator to current population ending
    typename std::vector< CHR < T>>::

    const_iterator cend() const;

    // select element at position pos in current population and copy it into mating population
    void select(int pos);

    // set all fitness to positive values
    void adjustFitness();

    // compute fitness sum of current population
    double getSumFitness() const;

    // get worst objective function total result from current population
    double getWorstTotal() const;

    // return population size
    int popsize() const;

    // return mating population size
    int matsize() const;

    // return tournament size
    int tntsize() const;

    // return numero of generation
    int nogen() const;

    // return number of generations
    int nbgen() const;

    // return selection pressure
    double SP() const;

    const GeneticAlgorithm <T> *ga_algo()
    {
      return ptr;
    }

    std::vector< CHR < T>>&

    get_newpop()
    {
      return newpop;
    }

    private:
    std::vector< CHR < T>> curpop;               // current population
    std::vector< CHR < T>> matpop;               // mating population
    std::vector< CHR < T>> newpop;               // new population

    const GeneticAlgorithm <T>
        *ptr = nullptr; // pointer to genetic algorithm
    int nbrcrov;                              // number of cross-over
    int matidx;                               // mating population index

    // elitism => saving best chromosomes in new population
    void elitism();

    // create new population from recombination of the old one
    void recombination();

    // complete new population randomly
    void completion();

    // update population (adapting, sorting)
    void updating();
  };

  /*-------------------------------------------------------------------------------------------------*/

  // constructor
  template <typename T>
  Population< T >::Population(const GeneticAlgorithm <T> &ga)
  {
    ptr = &ga;
    nbrcrov = (int)ceil(ga.covrate * (ga.popsize - ga.elitpop));

    // adjusting nbrcrov (must be an even number)
    if (nbrcrov % 2 != 0)
    {
      nbrcrov += 1;
    }

    // for convenience, we add elitpop to nbrcrov
    nbrcrov += ga.elitpop;

    // allocating memory
    curpop.resize(ga.popsize);
    matpop.resize(ga.matsize);
  }

  /*-------------------------------------------------------------------------------------------------*/

  // create a population of chromosomes
  template <typename T> void Population< T >::creation()
  {
    std::cout << "population creation in the first generation" << std::endl;
    int start = 0;
    // initializing first chromosome
    if (!ptr->initialSet.empty())
    {
      curpop[0] = std::make_shared< Chromosome< T>>(*ptr);
      curpop[0]->initialize();
      curpop[0]->evaluate();
      start++;
    }
    // getting the rest
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif

    /*** loading the populations from B.O. txt ***/
    std::vector< float > indiv;
    std::ifstream infile;
    int popsize = 10;
    int clusterCount = popsize;
    int dimens = 20; //dimensions
    std::cout << "reading txt ... " << std::endl;
    infile.open("../individuals.txt");
    double tmp;
    while (not infile.eof()) /////does the end of the file?
    {
      infile >> tmp;
      indiv.push_back(tmp);
    }
    //transmit indiv into array indivs[][]
    float indivs[indiv.size()/(dimens+1)][dimens];
    float fitness[indiv.size()/(dimens+1)]; //nrow
    for(int i = 0; i < indiv.size(); i++)
    {
      int m = floor(i/(dimens+1));
      int n = i%(dimens+1);
      if(n < dimens)
        indivs[m][n] = indiv[i];
      else
        fitness[m] = indiv[i];
    }

    float MaxFitPerCluster[clusterCount]; //clusterCount
    float MaxIndivPerCluster[clusterCount][dimens]; //clusterCount, dimens

    ////************** loading the lastest individuals from BO ************///
//    std::cout << " ------ loading the latest individuals ------ " << std::endl;
//    for(int i = 0; i < popsize; i++)
//    {
//      int k = indiv.size() -2 - i * (dimens + 1);
//      std::cout << "pop_" << i << "#  k=" << k << " ";
//      MaxFitPerCluster[i] = indiv[k];
//      for(int j = 0; j < dimens; j++)
//      {
//        MaxIndivPerCluster[i][j] = indiv[indiv.size() -1 - (i + 1) * (dimens + 1) + j];
//        std::cout << " [" << i << "][" << j << "]: " << MaxIndivPerCluster[i][j] << "  ";
//      }
//      std::cout << "MaxFitPerCluster[" << i << "]: " << MaxFitPerCluster[i] << std::endl;
//    }

    ////************ load the best popsize individuals from BO ***********////
//    std::vector< double > fit;
////    int MaxBoGen = (indiv.size()-1) / (dimens + 1); // not including
//    // initial samples
//    for (int i = 0; i < (indiv.size()-1)/(dimens + 1); i++)
//    {
//      fit.push_back(indiv[i * (dimens + 1) + dimens]);
//    }
//    std::sort(fit.begin(), fit.end()); // ranking from min to max
//    //matching top popsize fitness for each individual
//    //   std::cout << std::endl << "tmpnumgen: " << tmpnumgen << std::endl;
//    for (int i = 0; i < fit.size(); i++) // MaxBoGen = fit.size()
//    {
////      std::cout << "fit[" << (MaxBoGen - tmpnumgen -1) << "]: " << fit[MaxBoGen - tmpnumgen -1] << "  ";
//      for(int j = 0; j < popsize; j++)
//      {
//        //matching all to top j (in fit)
//        if (indiv[i * (dimens + 1) + dimens] == fit[fit.size() - 1 - j])
//        {
//          MaxFitPerCluster[j] = fit[fit.size() - 1 - j];
//          for(int k = 0; k < dimens; k++)
//          {
//            MaxIndivPerCluster[j][k] = indiv[i * (dimens + 1)+ k];
//          }
//          std::cout << std::endl;
//        }
//      }
//    }

    ///*********************** k-means individuals ************************///
    //transmit array into Mat
//    cv::Mat points(indiv.size() / (dimens + 1), //rows
//                   dimens, //colums
//                   CV_32F, //float, CV_64F: double
//                   indivs); //indivsTop
//    //    std::cout << "Data: " << points << std::endl;
//
//    cv::Mat centers(clusterCount, 1, points.type());
//
//    // values of 1st half of data set is set to 10
//    //change the values of 2nd half of the data set; i.e. set it to 20
//    std::cout << "----------- k-means clustering ---------- "
//              << "points.rows: " << points.rows << "  points.cols: "
//              << points.cols << std::endl;
//
//    cv::Mat labels;
//    kmeans(points,  //data
//           clusterCount,  //clusterCount: searched clusters in k-means
//           labels, cv::TermCriteria(
//            //type of termination criteria : It has 3 flags as below:
//            //cv2.TERM_CRITERIA_EPS - stop the algorithm iteration if specified accuracy, epsilon, is reached.
//            //cv2.TERM_CRITERIA_MAX_ITER - stop the algorithm after the specified number of iterations, max_iter.
//            //cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER - stop the iteration when any of the above condition is met.
//            CV_TERMCRIT_EPS + CV_TERMCRIT_ITER,
//            10, //max_iter - An integer specifying maximum number of iterations
//            1.0), //epsilon - Required accuracy
//           3, //聚类3次，取结果最好的那次
//           cv::KMEANS_PP_CENTERS, //聚类的初始化采用PP特定的随机算法
//           centers);
//    // we can print the matrix directly.
////        std::cout << "Data: " << points << std::endl;
////        std::cout << "Center: " << centers.at<float>(0) << std::endl;
//    /***** pick the center individual of each cluster *****/
////    memset(MaxFitPerCluster, -9999.99, sizeof(MaxFitPerCluster)); //initial
//    std::cout << "Center: " << std::endl;
//    std::cout << centers << std::endl;
//    for(int i = 0; i < clusterCount; i++)
//    {
//      for(int j = 0; j < dimens; j++)
//      {
//        MaxIndivPerCluster[i][j] = centers.at<float>(i*dimens + j);
//      }
//    }
    /***** pick the best individual of each cluster *****/
//    memset(MaxFitPerCluster, -9999.99, sizeof(MaxFitPerCluster)); //initial
//    //std::cout << "labels.rows: " << labels.rows << std::endl;
//    for (int i = 0; i < labels.rows; i++) //labels.rows is the number of all indi
//    {
//      std::cout << "labels[" << i << "]: " << labels.at< int >(i) << "  ";
//      for (int j = 0; j < clusterCount; j++) //clusterCount equals popsize
//      {
//        if (labels.at< int >(i) == j) //labels.at<int>(i) is the label of the cluster
//        {
//          std::cout << std::fixed << "fitness[" << i << "]: " << fitness[i]
//                    << "  " << std::endl;
//          if (fitness[i] > MaxFitPerCluster[j])
//          {
//            MaxFitPerCluster[j] = fitness[i];
//            for (int k = 0; k < dimens; k++)
//            {
//              MaxIndivPerCluster[j][k] = indivs[i][k]; //save best indiv
//              std::cout << "MaxIndivPerCluster[" << j << "][" << k << "]: "
//                        << indivs[i][k] << "  ";
//            }
//            std::cout << " MaxFitPerCluster[" << j << "]: "
//                      << MaxFitPerCluster[j] << std::endl;
//          }
//        }
//      }
//    }

    /********************* k-means on top p% individuals *********************/
    /****** top 50% individuals ******/
    std::cout << std::endl;
    std::cout << "--------- k-means clustering on top p% --------- "
              << std::endl;
    std::vector<float> fitbuf;
    for(int i = 0; i < indiv.size()/(dimens+1); i++)
    {
      fitbuf.push_back(indiv[i * (dimens + 1) + dimens]);
    }
    std::sort(fitbuf.begin(), fitbuf.end());

    float percent = 0.4; //0.2, 0.4
    float indivsTop[(int)(fitbuf.size() * percent)][dimens]; //individuals with top 50%
    float fitnessTop[(int)(percent*indiv.size()/(dimens+1))];
    for(int i = 0; i < fitbuf.size(); i++)
    {
      for(int j = 0; j < fitbuf.size() * percent; j++)
      {
        //matching from fitbuf.size() ~ fitbuf.size() * percent
        if (indiv[i * (dimens + 1) + dimens] == fitbuf[fitbuf.size() - 1 - j])
        {
          for(int k = 0; k < dimens; k++)
          {
            //indivsTop: individuals with top 50%, pass to points
            indivsTop[j][k] = indiv[i * (dimens + 1)+ k];
            std::cout << "  indivsTop[" << j << "][" << k << "]: " << indivsTop[j][k];
          }
          fitnessTop[j] = fitbuf[fitbuf.size() - 1 - j];
          std::cout << " fitnessTop[" << j << "]: " << fitnessTop[j] << std::endl;
        }
      }
    }

    /**** k-means on top p% individuals *****/
    //transmit array into Mat
    cv::Mat pointsTop(percent*indiv.size()/(dimens+1), //rows
                      dimens, //colums
                      CV_32F, //float, CV_64F: double
                      indivsTop); //indivsTop
    //std::cout << "Data: " << pointsTop << std::endl;

    cv::Mat centersTop(clusterCount, 1, pointsTop.type());
    cv::Mat labelsTop;
    // values of 1st half of data set is set to 10
    //change the values of 2nd half of the data set; i.e. set it to 20
    std::cout << "pointsTop.rows: " << pointsTop.rows
              << "  pointsTop.cols: " << pointsTop.cols << std::endl;

    kmeans(pointsTop,  //data
           clusterCount,  //clusterCount: searched clusters in k-means
           labelsTop,
           cv::TermCriteria(
             //type of termination criteria : It has 3 flags as below:
             //cv2.TERM_CRITERIA_EPS - stop the algorithm iteration if specified accuracy, epsilon, is reached.
             //cv2.TERM_CRITERIA_MAX_ITER - stop the algorithm after the specified number of iterations, max_iter.
             //cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER - stop the iteration when any of the above condition is met.
             CV_TERMCRIT_EPS+CV_TERMCRIT_ITER,
             10, //max_iter - An integer specifying maximum number of iterations
             1.0), //epsilon - Required accuracy
           3, //聚类3次，取结果最好的那次
           cv::KMEANS_PP_CENTERS, //聚类的初始化采用PP特定的随机算法
           centersTop);
    // we can print the matrix directly.
    //    std::cout << "Data: " << points << std::endl;
    //    std::cout << "Center: " << centers.at<float>(0) << std::endl;
    /***** pick the center individual of each cluster ****/
    //std::cout << "CentersTop: " << centersTop << std::endl;

    /***** pick the best individual of each cluster ****/
//    float MaxTopPerCluster[clusterCount];
    memset(MaxFitPerCluster, -9999.99, sizeof(MaxFitPerCluster)); //initial
//    float MaxTopInPerCluster[clusterCount][dimens];
    //std::cout << "labels.rows: " << labels.rows << std::endl;
    for (int i =0; i < labelsTop.rows; i++) //the number of Top indivs
    {
      std::cout << "labelsTop[" << i << "]: " << labelsTop.at<int>(i) << "  ";
      for(int j=0; j < clusterCount; j++) //clusterCount equals popsize
      {
        if(labelsTop.at<int>(i) == j) //labels.at<int>(i) is the label of the cluster
        {
          if(fitnessTop[i] > MaxFitPerCluster[j])
          {
            MaxFitPerCluster[j] = fitnessTop[i];
            for(int k = 0; k < dimens; k++)
            {
              MaxIndivPerCluster[j][k] = indivsTop[i][k]; //save best top indivs
              std::cout << "MaxIndivPerCluster[" << j << "][" << k << "]: "
                        << indivsTop[i][k] << "  ";
            }
            std::cout << " MaxFitPerCluster[" << j << "]: "
                      << MaxFitPerCluster[j];
          }
        }
      }
      std::cout << std::endl;
    }


    for (int i = start; i < ptr->popsize; ++i)
    {
      curpop[i] = std::make_shared< Chromosome< T>>(*ptr);
      curpop[i]->create(i, MaxFitPerCluster[i], MaxIndivPerCluster[i]);
      curpop[i]->evaluate0gen(); //row 316 in Chromosome.hpp
    }
    this->updating(); // updating population
  }

  /*-------------------------------------------------------------------------------------------------*/

  // population evolution (selection, recombination, completion, mutation), get next generation
  template <typename T> void Population< T >::evolution()
  {
    // initializing mating population index
    matidx = 0;

    // selecting mating population
    // curpop[] -> matpop[]
    ptr->Selection(*this);

    // applying elitism if required
    // curpop[] -> newpop[0...elitpop-1]
    this->elitism();

    // crossing-over mating population
    // matpop[] -> newpop[elitpop...nbrcrov-1]
    this->recombination();

    // completing new population
    // matpop[] -> newpop[nbrcrov...popsize]
    this->completion();

    // moving new population into current population for next generation
    curpop = std::move(newpop);

    // updating population
    this->updating();
  }

  /*-------------------------------------------------------------------------------------------------*/

  // elitism => saving best chromosomes in new population, making a copy of each elit chromosome
  template <typename T> void Population< T >::elitism()
  {
    // (re)allocating new population
    newpop.resize(ptr->popsize);

    if (ptr->elitpop > 0)
    {
      // copying elit chromosomes into new population
      std::transform(curpop.cbegin(),
                     curpop.cend(),
                     newpop.begin(),
                     [](const CHR< T > &chr) -> CHR< T >
                     { return std::make_shared< Chromosome< T>>(*chr); });

      for (size_t i = 0; i < curpop.size(); i++)
      {
        transmit_sigma< T >(*curpop[i], *newpop[i]);
      }
    }
  }

  /*-------------------------------------------------------------------------------------------------*/

  // create new population from recombination of the old one
  template <typename T> void Population< T >::recombination()
  {
    // creating a new population by cross-over
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif

    for (int i = ptr->elitpop; i < nbrcrov; i = i + 2)
    {
      // initializing 2 new chromosome
      newpop[i] = std::make_shared< Chromosome< T>>(*ptr);
      newpop[i + 1] = std::make_shared< Chromosome< T>>(*ptr);

      // crossing-over mating population to create 2 new chromosomes
      ptr->CrossOver(*this, newpop[i], newpop[i + 1]);

      // mutating new chromosomes
      ptr->Mutation(newpop[i]);
      ptr->Mutation(newpop[i + 1]);

      if (ptr->FixedValue != nullptr)
      {
        ptr->FixedValue(*this, i);
        ptr->FixedValue(*this, i + 1);
      }

      // evaluating new chromosomes
      //std::cout << "Nu_evaluation_re: " << i;
      newpop[i]->evaluate();
      //std::cout << "Nu_evaluation_re- " << (i + 1);
      newpop[i + 1]->evaluate();
    }
  }

  /*-------------------------------------------------------------------------------------------------*/

  // complete new population
  template <typename T> void Population< T >::completion()
  {
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif
    for (int i = nbrcrov; i < ptr->popsize; ++i)
    {
      // selecting chromosome randomly from mating population
      int pos = uniform< int >(0, ptr->matsize);
      newpop[i] = std::make_shared< Chromosome< T>>(*matpop[pos]);
      transmit_sigma< T >(*matpop[pos], *newpop[i]);

      // mutating chromosome
      ptr->Mutation(newpop[i]);

      if (ptr->FixedValue != nullptr)
      {
        ptr->FixedValue(*this, i);
      }
      //std::cout << "Nu_evaluation_com: " << i;
      // evaluating chromosome
      newpop[i]->evaluate();
    }
  }

  /*-------------------------------------------------------------------------------------------------*/

  // update population (adapting, sorting)
  template <typename T> void Population< T >::updating()
  {
    // adapting population to constraints
    if (ptr->Constraint != nullptr)
    {
      ptr->Adaptation(*this);
    }
    // sorting chromosomes from best to worst fitness
    std::sort(curpop.begin(), curpop.end(), [](
        const CHR< T > &chr1,
        const CHR< T > &chr2) -> bool
    { return chr1->fitness > chr2->fitness; });
  }

  /*-------------------------------------------------------------------------------------------------*/

  // access element in current population at position pos
  template <typename T>
  const CHR <T> &Population< T >::operator()(int pos) const
  {
#ifndef NDEBUG
    if (pos > ptr->popsize - 1)
    {
      throw std::invalid_argument(
          "Error: in galgo::Population<T>::operator()(int), exceeding current population memory.");
    }
#endif

    return curpop[pos];
  }

  /*-------------------------------------------------------------------------------------------------*/

  // access element in mating population at position pos
  template <typename T>
  const CHR <T> &Population< T >::operator[](int pos) const
  {
#ifndef NDEBUG
    if (pos > ptr->matsize - 1)
    {
      throw std::invalid_argument(
          "Error: in galgo::Population<T>::operator[](int), exceeding mating population memory.");
    }
#endif

    return matpop[pos];
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return iterator to current population beginning
  template <typename T> inline typename std::vector< CHR < T>>

  ::iterator Population< T >::begin()
  {
    return curpop.begin();
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return const iterator to current population beginning
  template <typename T> inline typename std::vector< CHR < T>>

  ::const_iterator Population< T >::cbegin() const
  {
    return curpop.cbegin();
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return iterator to current population ending
  template <typename T> inline typename std::vector< CHR < T>>

  ::iterator Population< T >::end()
  {
    return curpop.end();
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return const iterator to current population ending
  template <typename T> inline typename std::vector< CHR < T>>

  ::const_iterator Population< T >::cend() const
  {
    return curpop.cend();
  }

  /*-------------------------------------------------------------------------------------------------*/

  // select element at position pos in current population and copy it into mating population
  template <typename T> inline void Population< T >::select(int pos)
  {
#ifndef NDEBUG
    if (pos > ptr->popsize - 1)
    {
      throw std::invalid_argument(
          "Error: in galgo::Population<T>::select(int), exceeding current population memory.");
    }
    if (matidx == ptr->matsize)
    {
      throw std::invalid_argument(
          "Error: in galgo::Population<T>::select(int), exceeding mating population memory.");
    }
#endif

    matpop[matidx] = curpop[pos];
    matidx++;
  }

  /*-------------------------------------------------------------------------------------------------*/

  // set all fitness to positive values (used in RWS and SUS selection methods)
  template <typename T> void Population< T >::adjustFitness()
  {
    // getting worst population fitness
    double worstFitness = curpop.back()->fitness;

    if (worstFitness < 0)
    {
      // getting best fitness
      double bestFitness = curpop.front()->fitness;
      // case where all fitness are equal and negative
      if (worstFitness == bestFitness)
      {
        std::for_each(curpop.begin(), curpop.end(), [](CHR< T > &chr) -> void
        { chr->fitness *= -1; });
      }
      else
      {
        std::for_each(curpop.begin(),
                      curpop.end(),
                      [worstFitness](CHR< T > &chr) -> void
                      { chr->fitness -= worstFitness; });
      }
    }
  }

  /*-------------------------------------------------------------------------------------------------*/

  // compute population fitness sum (used in TRS, RWS and SUS selection methods)
  template <typename T> inline double Population< T >::getSumFitness() const
  {
    return std::accumulate(curpop.cbegin(), curpop.cend(), 0.0, [](
        double sum,
        const CHR< T > &chr) -> double
    { return sum + chr->fitness; });
  }

  /*-------------------------------------------------------------------------------------------------*/

  // get worst objective function total result from current population (used in constraint(s) adaptation)
  template <typename T> inline double Population< T >::getWorstTotal() const
  {
    auto it = std::min_element(curpop.begin(), curpop.end(), [](
        const CHR< T > &chr1,
        const CHR< T > &chr2) -> bool
    { return chr1->getTotal() < chr2->getTotal(); });
    return (*it)->getTotal();
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return population size
  template <typename T> inline int Population< T >::popsize() const
  {
    return ptr->popsize;
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return mating population size
  template <typename T> inline int Population< T >::matsize() const
  {
    return ptr->matsize;
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return tournament size
  template <typename T> inline int Population< T >::tntsize() const
  {
    return ptr->tntsize;
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return numero of generation
  template <typename T> inline int Population< T >::nogen() const
  {
    return ptr->nogen;
  }


  /*-------------------------------------------------------------------------------------------------*/

  // return number of generations
  template <typename T> inline int Population< T >::nbgen() const
  {
    return ptr->nbgen;
  }

  /*-------------------------------------------------------------------------------------------------*/

  // return selection pressure
  template <typename T> inline double Population< T >::SP() const
  {
    return ptr->SP;
  }

  //=================================================================================================

}

#endif



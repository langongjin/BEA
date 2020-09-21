//| Copyright Inria May 2015
//| This project has received funding from the European Research Council (ERC) under
//| the European Union's Horizon 2020 research and innovation programme (grant
//| agreement No 637972) - see http://www.resibots.eu
//|
//| Contributor(s):
//|   - Jean-Baptiste Mouret (jean-baptiste.mouret@inria.fr)
//|   - Antoine Cully (antoinecully@gmail.com)
//|   - Konstantinos Chatzilygeroudis (konstantinos.chatzilygeroudis@inria.fr)
//|   - Federico Allocati (fede.allocati@gmail.com)
//|   - Vaios Papaspyros (b.papaspyros@gmail.com)
//|   - Roberto Rama (bertoski@gmail.com)
//|
//| This software is a computer library whose purpose is to optimize continuous,
//| black-box functions. It mainly implements Gaussian processes and Bayesian
//| optimization.
//| Main repository: http://github.com/resibots/limbo
//| Documentation: http://www.resibots.eu/limbo
//|
//| This software is governed by the CeCILL-C license under French law and
//| abiding by the rules of distribution of free software.  You can  use,
//| modify and/ or redistribute the software under the terms of the CeCILL-C
//| license as circulated by CEA, CNRS and INRIA at the following URL
//| "http://www.cecill.info".
//|
//| As a counterpart to the access to the source code and  rights to copy,
//| modify and redistribute granted by the license, users are provided only
//| with a limited warranty  and the software's author,  the holder of the
//| economic rights,  and the successive licensors  have only  limited
//| liability.
//|
//| In this respect, the user's attention is drawn to the risks associated
//| with loading,  using,  modifying and/or developing or reproducing the
//| software by the user in light of its specific status of free software,
//| that may mean  that it is complicated to manipulate,  and  that  also
//| therefore means  that it is reserved for developers  and  experienced
//| professionals having in-depth computer knowledge. Users are therefore
//| encouraged to load and test the software's suitability as regards their
//| requirements in conditions enabling the security of their systems and/or
//| data to be ensured and,  more generally, to use and operate it in the
//| same conditions as regards security.
//|
//| The fact that you are presently reading this means that you have had
//| knowledge of the CeCILL-C license and that you accept its terms.
//|
#ifndef LIMBO_STAT_CONSOLE_SUMMARY_HPP
#define LIMBO_STAT_CONSOLE_SUMMARY_HPP

#include "stat/stat_base.hpp"

namespace limbo {
    namespace stat {
        ///@ingroup stat
        ///write the status of the algorithm on the terminal
        template <typename Params>
        struct ConsoleSummary : public StatBase<Params> {
            template <typename BO, typename AggregatorFunction>
            void operator()(const BO& bo, const AggregatorFunction& afun)
            {
                if (!bo.stats_enabled() || bo.observations().empty())
                    return;

                std::cout << bo.total_iterations() << " new point: "
                          << bo.samples().back().transpose() //back(): returnthe last element.
                          << " value: " << afun(bo.observations().back()) //same with bo.observations().back()
                          << " best:" << afun(bo.best_observation(afun)) << std::endl;// same with bo.best_observation()
                //save fitness of initial samples
//                std::ofstream value;
//                value.open("../value.txt", std::ios::app);
//                value << std::fixed << afun(bo.observations().back()) << " -2-"<< std::endl;
                /******* save individuals for lastest **********/
                std::ofstream individuals;
                individuals.open("../individuals1.txt", std::ios::app);
                individuals << std::fixed << bo.samples().back().transpose() << " "<< afun(bo.observations().back()) << std::endl;//scientific notation

//                /******* save top pop-size individuals **********/
//                int eval = 0;
//                int Nomin = 0;
//                int Rankeval[] = 0;
//                double fitness[] = {0};
//                //std::vector <double> individuals;
//                double individuals[20][]; //pop-size of EA
//
//                if(eval < 20) //pop-size of EA
//                {
//                    individuals[eval] << bo.samples().back().transpose();
//                    std::cout << "individuals[" << eval << "] "
//                              << individuals[eval] << std::endl;
//                    if (eval == 1)
//                    {
//                        if(fitness[1] > fitness[0])
//                        {
//                            Nomin = 0;
//                            Rankeval[0] = 1; // ranking
//                            Rankeval[1] = 0; // the best is ranking 0
//                        }
//                        else
//                        {
//                            Nomin = 1;
//                            Rankeval[0] = 0;
//                            Rankeval[1] = 1; // the best is ranking 0
//                        }
//                    }
//                    if (eval > 1) //0,1, < 2,3
//                    {
//                        for{int i = 0; i < eval; i++}
//                        {
//                            if (fitness[eval] < fitness[Nomin])
//                            {
//                                Nomin = eval;
//                                Rankeval[eval] = eval;
//                            }
//                            if (fitness[eval] > fitness[i])
//                            {
//                                Rankeval[i]++;
//
//                            }
//                            else if (fitness[eval] < fitness[i])
//                            {
//                                Rankeval[eval] = Rankeval[i] + 1;
//                                if ()
//                            }
//
//                        }
//                        if(fitness[eval] < fitness[Nomin])
//                        {
//                            Nomin = eval;
//
//                        }
//                    }
//                }
//                if (i > 20)
//                {
//                    if(fitness[eval] > fitness[Nomin])
//                    {
//                        individuals[Nomin] = afun(bo.observations().back());
//
//                    }
//                }
//
//                eval ++;
            }
        };
    }
}

#endif

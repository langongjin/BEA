//#include "DifferentialEvolution.h"
//#include <ctime>
//
//class SimpleQuadriatic : public de::IOptimizable
//{
//public:
//  double EvaluteCost(std::vector<double> inputs) const override
//  {
//    assert(inputs.size() == 2);
//    double x = inputs[0];
//    double y = inputs[1];
//    return x * x + 2 * x * y + 3 * y * y;
//  }
//
//  unsigned int NumberOfParameters() const override
//  {
//    return 2;
//  }
//
//  std::vector<Constraints> GetConstraints() const override
//  {
//    std::vector<Constraints> constr(NumberOfParameters());
//    for (auto& c : constr)
//    {
//      c = Constraints(-100.0, 100.0, true);
//    }
//    return constr;
//  }
//};
//
//int main()
//{
//  SimpleQuadriatic cost;
//  de::DifferentialEvolution de(cost, 100, std::time(nullptr));
//  de.Optimize(1000, true);
//  return 0;
//}

#include "DifferentialEvolution.h"
#include "TestFunctions.h"

int main()
{
  // Creat Rastring's function in 5 dimensions
  de::Rastrigin cost(20);
  int iters = 10;
  for (int run = 0; run < iters; run++) {
    // Create Differential Evolution optimizer with population size of 50
    de::DifferentialEvolution de(cost, 10);

    // Optimize for 1000 iterations with verbose output.
    de.Optimize(150, true);
  }
  return 0;
}
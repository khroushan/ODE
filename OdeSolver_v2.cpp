// (setq scroll-step 1); keyboard scroll one line at a time 
#include <iostream>

class AbstractOdeSolver
{
protected:
  double stepSize;
  double initialTime;
  double finalTime;
  double initialValue;
public:
  void SetStepSize(double h);
  void SetTimeInterval(double t0, double t1);
  void SetInitialValue(double y0);
  double (*RightHandSide)(double y, double t);
  virtual double
  SolveEquation(double (*RightHandSide)(double y, double t)) = 0;
  virtual double
  RKsolver(double (*RightHandSide)(double y, double t)) = 0;
};

void AbstractOdeSolver::SetStepSize(double h)
{
  stepSize = h;
}

void AbstractOdeSolver::SetTimeInterval(double t0, double t1)
{
  initialTime = t0;
  finalTime = t1;
}

void AbstractOdeSolver::SetInitialValue(double y0)
{
  initialValue = y0;
}

// ************************************************************
// derived class
class ForwardEulerSolver: public AbstractOdeSolver
{
public:
  // double RightHandSide(double y, double t);
  double
  SolveEquation(double (*RightHandSide)(double y, double t));
  double
  RKsolver(double (*RightHandSide)(double y, double t));
  
};
// ================================================
// double ForwardEulerSolver::(*RightHandSide)(double y, double t)
// {
//   return 1 + y + 3.*t*t;
// }
// ***************************************
double ForwardEulerSolver::
SolveEquation(double (*RightHandSide)(double y, double t))
{
  int Nstep = (finalTime - initialTime)/stepSize;
  double y = initialValue;
  // implementing y_{i+1} = y_i + f(y_i,t_i)* dt
  
  for(int i=0; i<Nstep; i++)
  {
    double t = initialTime + i*stepSize;
    y += (*RightHandSide)(y,t)*stepSize;
    std::cout << t << "  " << y << "  \n";
  }
  return 0;
}
// ****************************************
double ForwardEulerSolver::
RKsolver(double (*RightHandSide)(double y, double t))
{
  //  Runge-Kuttal method to solve ODE
  int Nstep = (finalTime - initialTime)/stepSize;
  double y = initialValue;

  for(int i=0; i<Nstep; i++)
  {
    double t = initialTime + i*stepSize;
    double k1 = (*RightHandSide)(y      , t);
    double k2 = (*RightHandSide)(y+k1/2., t+stepSize/2.);
    double k3 = (*RightHandSide)(y+k2/2., t+stepSize/2.);
    double k4 = (*RightHandSide)(y+k3   , t+stepSize);
    y += (k1+2*k2+2*k3+k4)*(stepSize/6.);

    std::cout << t << "  " << y << "\n";
  }
  return 0;
}
// ***************************************
double func1(double y, double t)
{
  return y + t;
}



// ==============
//  Main Function
// ==============

int main(int argc, char* argv[])
{
  ForwardEulerSolver example;
  example.SetStepSize(0.1);
  example.SetTimeInterval(1, 5.0);
  example.SetInitialValue(1.0);
  // example.SolveEquation();
  example.RKsolver(func1);
  return 0;
}
  
  

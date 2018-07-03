#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 20;
double dt = 0.1;
/*
vars[0],..vars[9]->x1-x10
vars[10],..vars[19]->y1-y10
vars[20],..vars[29]->phi1-phi10
vars[30],..vars[39]->v1-v10
vars[40],..vars[49]->cte11-cte10
vars[50],..vars[59]->e_phi1-e_phi10
// you need 9 actautors for 10 time intervals
vars[60],..vars[68]->del1-del9
vars[69],..vars[67]->a1-a9
*/
int x_start=0;
int y_start=x_start+N;
int phi_start=y_start+N;
int v_start=phi_start+N;
int cte_start=v_start+N;
int ephi_start=cte_start+N;
int del_start=ephi_start+N;
int a_start=del_start+N-1;

//Reference velocity(try to be this fast)
//OBJECTIVES
double ref_v=90;
double ref_cte=0;
double ref_ephi=0;
// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    //The cost is stored in the first element of the vector fg;
    // Set it to zero at the start.
    fg[0]=0;

    //Cost function
    //TODO:Define the cost related reference state

    for (int t=0; t <N;t++){
      fg[0]+= 1600 * CppAD::pow(vars[cte_start+t]-ref_cte,2);
      fg[0]+= 1100 * CppAD::pow(vars[ephi_start+t]-ref_ephi,2);
      fg[0]+= CppAD::pow(vars[v_start+t]-ref_v,2);
    }

    //Minimize the use of actuators
    for (int t=0; t<N-1;t++){
      fg[0]+= 5 * CppAD::pow(vars[del_start+t],2);
      fg[0]+= 5 * CppAD::pow(vars[a_start+t],2);
      fg[0]+= 70 * CppAD::pow(vars[cte_start + t] * vars[v_start + t], 2);
    }

    //Minimize the value gap between sequential actuatations
    for (int t=0; t<N-2;t++){
      fg[0]+= 10 * CppAD::pow(vars[del_start+t+1]-vars[del_start+t],2);
      fg[0]+= 10 * CppAD::pow(vars[a_start+t+1] - vars[a_start+t],2);
    }

  //Setting up constraints

  fg[1+x_start] = vars[x_start];
  fg[1+y_start] = vars[y_start];
  fg[1+phi_start] = vars[phi_start];
  fg[1+v_start] = vars[v_start];
  fg[1+cte_start] = vars[cte_start];
  fg[1+ephi_start] = vars[ephi_start];

  for (int t = 1; t < N; t++) {
    // The state at time t+1 .
    AD<double> x1 = vars[x_start + t];
    AD<double> y1 = vars[y_start + t];
    AD<double> phi1 = vars[phi_start + t];
    AD<double> v1 = vars[v_start + t];
    AD<double> cte1 = vars[cte_start + t];
    AD<double> ephi1 = vars[ephi_start + t];

    // The state at time t.
    AD<double> x0 = vars[x_start + t - 1];
    AD<double> y0 = vars[y_start + t - 1];
    AD<double> phi0 = vars[phi_start + t - 1];
    AD<double> v0 = vars[v_start + t - 1];
    AD<double> cte0 = vars[cte_start + t - 1];
    AD<double> ephi0 = vars[ephi_start + t - 1];


    AD<double> del0 = vars[del_start + t - 1];
    AD<double> a0 = vars[a_start + t - 1];

    AD<double> f0 = coeffs[0] + coeffs[1] * x0;
    AD<double> phides0 = CppAD::atan(coeffs[1]);

    // Here's `x` to get you started.
    // The idea here is to constraint this value to be 0.
    //
    // Recall the equations for the model:
    // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
    // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
    // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
    // v_[t] = v[t-1] + a[t-1] * dt
    // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
    // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
    fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(phi0) * dt);
    fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(phi0) * dt);
    fg[1 + phi_start + t] = phi1 - (phi0 - (v0 * del0 / Lf )* dt);
    fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
    fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(ephi0) * dt));
    fg[1 + ephi_start + t] = ephi1 - ((phi0 - phides0) - (v0 * del0 / Lf) * dt);
  }
}
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  // we have 10 timesteps and 6 state variables(x,y,phi,v,cte,ephi) and 2 actuators
  size_t n_vars = 6 * N + 2 * (N-1);
  // TODO: Set the number of constraints
  // 6 equations for N times steps
  size_t n_constraints = N*6;



  double x = state[0];
  double y = state[1];
  double phi = state[2];
  double v = state[3];
  double cte = state[4];
  double ephi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }


  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[phi_start] = phi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[ephi_start] = ephi;


  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    for (int i = 0; i < del_start; i++) {
      vars_lowerbound[i] = -1.0e19;
      vars_upperbound[i] = 1.0e19;
    }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = del_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332*2.67;
    vars_upperbound[i] = 0.436332*2.67;
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -0.8;
    vars_upperbound[i] = 0.8;
  }


  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  // For the initial state ,we need to set special constraint(1 element of the N )
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[phi_start] = phi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[ephi_start] = ephi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[phi_start] = phi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[ephi_start] = ephi;
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  vector<double> result;
  result.push_back(solution.x[del_start]);
  result.push_back(solution.x[a_start]);

  for (int i = 0; i < N-1; i++) {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }
  return result;
}

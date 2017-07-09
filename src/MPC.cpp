#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
 public:
  FG_eval(Eigen::VectorXd coeffs,
          double Lf,
          size_t N,
          double dt,
          double ref_v,
          double cte_weight,
          double epsi_weight,
          double v_weight,
          double steer_weight,
          double a_weight,
          double delta_steer_weight,
          double delta_a_weight) :
      coeffs_{coeffs},
      Lf_{Lf},
      N_{N},
      dt_{dt},
      ref_v_{ref_v},
      cte_weight_{cte_weight},
      epsi_weight_{epsi_weight},
      v_weight_{v_weight},
      steer_weight_{steer_weight},
      a_weight_{a_weight},
      delta_steer_weight_{delta_steer_weight},
      delta_a_weight_{delta_a_weight},
      x_start_{0},
      y_start_{x_start_ + N_},
      psi_start_{y_start_ + N_},
      v_start_{psi_start_ + N_},
      cte_start_{v_start_ + N_},
      epsi_start_{cte_start_ + N_},
      delta_start_{epsi_start_ + N_},
      a_start_{delta_start_ + N_ - 1}
  {}

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < N_; t++) {
      fg[0] += cte_weight_ * CppAD::pow(vars[cte_start_ + t], 2);
      fg[0] += epsi_weight_ * CppAD::pow(vars[epsi_start_ + t], 2);
      fg[0] += v_weight_ * CppAD::pow(vars[v_start_ + t] - ref_v_, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < N_ - 1; t++) {
      fg[0] += steer_weight_ * CppAD::pow(vars[delta_start_ + t], 2);
      fg[0] += a_weight_ * CppAD::pow(vars[a_start_ + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N_ - 2; t++) {
      fg[0] += delta_steer_weight_ * CppAD::pow(vars[delta_start_ + t + 1] - vars[delta_start_ + t], 2);
      fg[0] += delta_a_weight_ * CppAD::pow(vars[a_start_ + t + 1] - vars[a_start_ + t], 2);
    }
    fg[1 + x_start_] = vars[x_start_];
    fg[1 + y_start_] = vars[y_start_];
    fg[1 + psi_start_] = vars[psi_start_];
    fg[1 + v_start_] = vars[v_start_];
    fg[1 + cte_start_] = vars[cte_start_];
    fg[1 + epsi_start_] = vars[epsi_start_];

    for (int t = 1; t < N_; t++) {
      AD<double> x1 = vars[x_start_ + t];
      AD<double> y1 = vars[y_start_ + t];
      AD<double> psi1 = vars[psi_start_ + t];
      AD<double> v1 = vars[v_start_ + t];
      AD<double> cte1 = vars[cte_start_ + t];
      AD<double> epsi1 = vars[epsi_start_ + t];

      // The state at time t.
      AD<double> x0 = vars[x_start_ + t - 1];
      AD<double> y0 = vars[y_start_ + t - 1];
      AD<double> psi0 = vars[psi_start_ + t - 1];
      AD<double> v0 = vars[v_start_ + t - 1];
      AD<double> cte0 = vars[cte_start_ + t - 1];
      AD<double> epsi0 = vars[epsi_start_ + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start_ + t - 1];
      AD<double> a0 = vars[a_start_ + t - 1];

      AD<double> f0 = coeffs_[0] + coeffs_[1] * x0 + coeffs_[2] * x0 * x0 + coeffs_[3] * x0 * x0 * x0;
      AD<double> psides0 = CppAD::atan(coeffs_[1] + (2 * coeffs_[2] * x0) + (3 * coeffs_[3] * (x0*x0)));

      fg[1 + x_start_ + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt_);
      fg[1 + y_start_ + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt_);
      fg[1 + psi_start_ + t] = psi1 - (psi0 + v0 * delta0 / Lf_ * dt_);
      fg[1 + v_start_ + t] = v1 - (v0 + a0 * dt_);
      fg[1 + cte_start_ + t] =  cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt_));
      fg[1 + epsi_start_ + t] = epsi1 - (psi0 - psides0 + v0 * delta0 / Lf_ * dt_);
    }
  }

 private:
  Eigen::VectorXd coeffs_;
  double Lf_;
  size_t N_;
  double dt_;
  double ref_v_;
  double cte_weight_;
  double epsi_weight_;
  double v_weight_;
  double steer_weight_;
  double a_weight_;
  double delta_steer_weight_;
  double delta_a_weight_;

  size_t x_start_;
  size_t y_start_;
  size_t psi_start_;
  size_t v_start_;
  size_t cte_start_;
  size_t epsi_start_;
  size_t delta_start_;
  size_t a_start_;

};

//
// MPC class definition implementation.
//
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;

  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = N_ * state.size() + (N_ - 1) * 2;
  // Number of constraints
  size_t n_constraints = N_ * state.size();

  // Initial value of the independent variables.
  // Should be 0 except for the initial values.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  // Set the initial variable values
  vars[x_start_] = x;
  vars[y_start_] = y;
  vars[psi_start_] = psi;
  vars[v_start_] = v;
  vars[cte_start_] = cte;
  vars[epsi_start_] = epsi;

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start_; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start_; i < a_start_; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start_; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start_] = x;
  constraints_lowerbound[y_start_] = y;
  constraints_lowerbound[psi_start_] = psi;
  constraints_lowerbound[v_start_] = v;
  constraints_lowerbound[cte_start_] = cte;
  constraints_lowerbound[epsi_start_] = epsi;

  constraints_upperbound[x_start_] = x;
  constraints_upperbound[y_start_] = y;
  constraints_upperbound[psi_start_] = psi;
  constraints_upperbound[v_start_] = v;
  constraints_upperbound[cte_start_] = cte;
  constraints_upperbound[epsi_start_] = epsi;

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs,
                  Lf_,
                  N_,
                  dt_,
                  ref_v_,
                  cte_weight_,
                  epsi_weight_,
                  v_weight_,
                  steer_weight_,
                  a_weight_,
                  delta_steer_weight_,
                  delta_a_weight_);

  // options
  std::string options;
  options += "Integer print_level  0\n";
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

  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Predicted X and Y path
  predicted_x_.clear();
  predicted_y_.clear();
  for (int i = 0; i < N_ - 1; i++) {
    predicted_x_.push_back(solution.x[x_start_ + i + 1]);
    predicted_y_.push_back(solution.x[y_start_ + i + 1]);
  }

  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {solution.x[delta_start_], solution.x[a_start_]};
}

#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

constexpr double miles_per_hour_to_meters_per_second_ = 1600.0/3600.0;

class MPC {
 public:
  MPC(double Lf,
      size_t N,
      double dt,
      double ref_v_mph,
      double cte_weight,
      double epsi_weight,
      double v_weight,
      double steer_weight,
      double a_weight,
      double delta_steer_weight,
      double delta_a_weight) :
      Lf_{Lf},
      N_{N},
      dt_{dt},
      ref_v_{ref_v_mph * miles_per_hour_to_meters_per_second_},
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
  {
  }

  ~MPC();

  const vector<double> &predicted_x() const { return predicted_x_;}
  const vector<double> &predicted_y() const { return predicted_y_;}

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

 private:

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

  vector<double> predicted_x_;
  vector<double> predicted_y_;

  size_t x_start_;
  size_t y_start_;
  size_t psi_start_;
  size_t v_start_;
  size_t cte_start_;
  size_t epsi_start_;
  size_t delta_start_;
  size_t a_start_;

};

#endif /* MPC_H */

#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  const size_t N = 10;
  const double latency = 0.1;
  const double Lf = 2.67;
  const double dt = 0.09;
  const double ref_v_mph = 70;
  const double cte_weight = 8;
  const double epsi_weight = 8;
  const double v_weight = 1;
  const double steer_weight = 1;
  const double a_weight = 1;
  const double delta_steer_weight = 1000;
  const double delta_a_weight = 1;
  
  constexpr double miles_per_hour_to_meters_per_second = 1600.0/3600.0;

  // MPC is initialized here!
  MPC mpc{Lf,
          N,
          dt,
          ref_v_mph,
          cte_weight,
          epsi_weight,
          v_weight,
          steer_weight,
          a_weight,
          delta_steer_weight,
          delta_a_weight};

  h.onMessage([&mpc,&Lf, &latency](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          const vector<double> ptsx = j[1]["ptsx"];
          const vector<double> ptsy = j[1]["ptsy"];
          const double px = j[1]["x"];
          const double py = j[1]["y"];
          const double psi = j[1]["psi"];
          const double v = static_cast<double>(j[1]["speed"]) * miles_per_hour_to_meters_per_second;
          const double steer = static_cast<double>(j[1]["steering_angle"]) * -1;
          const double throttle = j[1]["throttle"];

          Eigen::VectorXd xvals(ptsx.size());
          Eigen::VectorXd yvals(ptsy.size());
          for (size_t i = 0; i < ptsx.size(); i++) {
            const double dx = ptsx[i] - px;
            const double dy = ptsy[i] - py;
            xvals[i] = (dx * cos(psi)) + (dy * sin(psi));
            yvals[i] = (dy * cos(psi)) - (dx * sin(psi));
          }

          Eigen::VectorXd coeffs = polyfit(xvals, yvals, 3);

          double predicted_psi = (v / Lf) * steer * latency;
          double predicted_x = v * cos(predicted_psi) * latency;
          double predicted_y = v * sin(predicted_psi) * latency;
          double predicted_v = v + (throttle * latency);
          double epsi = -atan(coeffs[1]);
          double cte = polyeval(coeffs, 0) - 0 + (v * sin(epsi) * latency);
          epsi += v/Lf * steer * latency;

          Eigen::VectorXd state(6);
          state << predicted_x, predicted_y, predicted_psi, predicted_v, cte, epsi;

          auto steer_throttle = mpc.Solve(state, coeffs);

          double steer_value = -steer_throttle[0] / deg2rad(25.0);
          double throttle_value = steer_throttle[1];

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;
          msgJson["mpc_x"] = mpc.predicted_x();
          msgJson["mpc_y"] = mpc.predicted_y();

          vector<double> next_x_vals;
          vector<double> next_y_vals;
          double space = 3.0;
          for (int i = 1; i < 15 ; ++i) {
            next_x_vals.push_back(space * i);
            next_y_vals.push_back(polyeval(coeffs, space * i));
          }
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

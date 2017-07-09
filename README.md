# CarND-Controls-MPC

This project implements a Model Predictive Control (MPC) to guide a car in the simulator using 
2 actuators: steering angle and throttle/brake. The goal is to drive the car around the track safely.

The key idea of MPC is to use a model for the trajectory of the car that includes the variable 
actuators and then optimize the trajectory so that it matches a reference trajectory.
 
## The Model
I used a Kinematic Model in this project, extended with crosstrack error and orientation error:
 * x: car x-coordinate.
 * y: car y-coordinate.
 * psi: car orientation.
 * v: car velocity.
 * cte: crosstrack error. Difference between car position and expected car position.
 * epsi: orientation error. Difference between car orientation and expected car orientation.
  

The model is implemented with the following equations, where delta[t] is the steering angle and
a[t] is the throttle/break factor:

      x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      v_[t+1] = v[t] + a[t] * dt
      cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt

The model is implemented in FG_eval::operator() in the file MPC.cpp

## Cost Function
The objective for the MPC optimization uses 7 weighted components: 
* Crosstrack Error.
* Orientation Error.
* Target velocity.
* Steering angle.
* Acceleration/Brake.
* Delta steering angle.
* Delta acceleration.

The target velocity is set to 70 MPH, meaning the target velocity cost is minimal when velocity is 
70 MPH.

After trying several different weights I found the following values for the cost components provided 
a good result:
* Crosstrack Error = 8
* Orientation Error =  8
* Target velocity = 1
* Steering angle =  1
* Acceleration/Brake = 1
* Delta steering angle = 1000
* Delta acceleration = 1.

In order to have a smooth driving the delta steering angle, i.e. the difference between the current
and previous steering angle should have a big weight compared to the other metrics.

The cost function is implemented in FG_eval::operator().

## Length and Duration
The parameters N and dt control the prediction horizon T. 
N is the number of timesteps in the horizon, dt is how much time elapses between actuations. 
For example, if N were 10 and dt were 0.05, then T would be 0.5 seconds.

After trying with several different values for N and dt I decided to use N = 10 and dt = 0.09, 
which gives T = 0.9s and considering an average velocity of 65mph it results in an horizon of 26 meters.

## Polynomial Fitting and MPC Preprocessing
One of the inputs for the method implement in this project is the expected X and Y values for the 
track. Those values are fit to a 3rd degree polynomial and used in the MPC optimization.

During preprocessing I transform the X and Y values from track coordinates to car coordinates. 

Preprocessing and polynomial fitting are implemented in main.cpp. 

## Taking Latency into account
The simulator has a 100ms latency between sending commands and actually having them applied to the
car. This is taken into account in this project by predicting what would be the state of the car
after 100ms and using that as input for the MPC optimization. The latency prediction is implemented
in main.cpp

## Conclusion
This project successfully implements a Model Predict Control (MPC). The implementation has a target
velocity of 70 MPH, and the car drives around the track in the simulator at around 65 MPH.



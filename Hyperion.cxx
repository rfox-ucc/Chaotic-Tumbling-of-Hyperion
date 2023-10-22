
#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "gnuplot.cxx"
using namespace std;
//The function compute takes the inital conditions (velocities, positions andangular displacements) and calculates the position and angular displacement ofHyperion
// T: the final time of the simulation
// N: the number of times to be plotted, including t=0.
// name: an identifying string that is attached to graphs and file names
// difference: indicates whether the model should calculate the real difference between values of theta (true) or else format the differences to be between pi
// and -pi (false)
//
// Results are output as .jpg files
int compute( double vx_0, double vy_0, double x_0, double y_0, double T, double theta_0, double theta1_0, int N, string name, bool difference)
{
//the neccessary arrays are defined
double x[N];
double y[N];
double vx[N];
double vy[N];
double theta[N];
double omega[N];
double t[N];
double diff[N];
double theta1[N];
double omega1[N];
//the timestep size delta_t is determined, inital time is assumed to be 0
double delta_t = T / (N - 1);
//the initial conditions are inputted as the 0th element of respective arrays
//omega is set to be zero
x[0] = x_0;
y[0] = y_0;
vx[0] = vx_0;
vy[0] = vy_0;
theta[0] = theta_0;
omega[0] = 0;
omega1[0] = 0;
theta1[0] = theta1_0;
t[0] = 0;
//The difference between the two values of theta is the absolute value of their difference, it must be positive as we will take
//the log of these values to determine the lyapunov coefficient later
diff[0]=abs(theta1[0]-theta[0]);
//the variable r here corresponds to rc, the distance from the centre of mass to Saturn. The values of r are not stored in an array as they are not required
//later for graphing
double r;
int j = 0;
while(j < N - 1)
{
//add a timestep to the current time to find the time for the next step
t[j+1] = t[j] + delta_t;
//r is the distance between the point (xc,yc) and the origin, therefore the value of r is SQRT(x^2+y^2), the next line computes this
r = sqrt(x[j]*x[j] + y[j]*y[j]);
//The difference equations (11) from the report are implemented here
//note that theta and omega correspond to one pair of initial conditions and theta1 omega1 to another
//the trigonometry functions of the math library are used here
omega[j+1] = omega[j] - ((3*4*M_PI*M_PI)/(r*r*r*r*r))*(x[j]*sin(theta[j])- y[j]*cos(theta[j]))*(x[j]*cos(theta[j])+y[j]*sin(theta[j]))*delta_t;
omega1[j+1] = omega1[j] - ((3*4*M_PI*M_PI)/(r*r*r*r*r))*(x[j]*sin(theta1[j])-y[j]*cos(theta1[j]))*(x[j]*cos(theta1[j])+y[j]*sin(theta1[j]))*delta_t;
theta[j+1] = theta[j] + omega[j+1]*delta_t;
theta1[j+1] = theta1[j] + omega1[j+1]*delta_t;
vx[j+1] = vx[j] - ((4*M_PI*M_PI*x[j])/(r*r*r))*delta_t;
vy[j+1] = vy[j] - ((4*M_PI*M_PI*y[j])/(r*r*r))*delta_t;
x[j + 1] = x[j] + vx[j+1]*delta_t;
y[j + 1] = y[j] + vy[j+1]*delta_t;
//the next segment of the loop checks if the difference boolean is true or false
//if it is false, the values of theta are bounded between pi and -pi so as to make the graphs clearer and closer to those to be recreated from the book
//if the value of theta exceeds pi, then 2pi is subtracted from it, this keeps the values between pi and -pi
if (difference == false)
{
if (theta[j+1] > M_PI)
{
theta[j+1] = theta[j+1] - 2*M_PI;
}
if (theta1[j+1] > M_PI)
{
theta1[j+1] = theta1[j+1] - 2*M_PI;
}
diff[j+1] = abs(theta1[j+1] - theta[j+1]);
if (diff[j+1] > M_PI)
{
diff[j+1] = abs(diff[j+1] - 2*M_PI);
}
}
//if the difference boolean is true, than the values of theta are left unbounded so as to allow the true difference in the angles to be found
//the log of the difference is also taken to prepare the values for plotting against time to find the lyapunov exponent
else if (difference == true)
{
diff[j+1] = abs(theta1[j+1] - theta[j+1]);
diff[j+1] = log(diff[j+1]);
}
j = j + 1;
}
//these strings correspond to the titles of the various plots genrated
string tstr = "Trajectory of " + name;
string thstr = "Theta against time for " + name;
string ostr = "Omega against time for " + name;
string dstr = "Difference of theta for " + name;
//if the difference boolean is set to true, then only the graph of difference against time is generated to reduce clutter
if (difference == true)
{
gnuplot_one_function_general(dstr, "lines", "t (hyperion years)", "log(difference theta) ", t, diff, N,true,(dstr+".jpg"),false,false,0,10,0,3.5);
}
else
{
gnuplot_one_function_general(dstr, "lines", "t (hyperion years)", "difference theta", t, diff, N,true,(dstr+".jpg"),false,false,0,10,0,3.5);
//notice that the trajectory jpg is made square so as to more accurately represent the true trajectory
gnuplot_one_function_general(tstr, "lines", "x (HU)", "y (HU)", x, y, N,true,(tstr+".jpg"),true,true,-2,2,-2,2);
gnuplot_one_function_general(thstr, "lines", "t (Hyperion Years)", "theta", t, theta, N,true,(thstr+".jpg"),false,false,-30,30,-30,30);
gnuplot_one_function_general(ostr, "lines", "t (hyperion years)", "omega (1/second)", t, omega, N,true,(ostr+".jpg"),false,false,-30,30,-30,30);
}
return 0;
}
int main()
{
//the main program generates the trajectory, theta against time, omega against time and difference of theta for a circular and elliptical orbit
compute(0,2*M_PI,1,0,10,0,0.01,25000,"Circular Orbit",false);
compute(0,5,1,0,10,0,0.01,25000,"Elliptical Orbit",false);
//next the program computes the log of the difference in theta for two different inital values with orbits of increasing eccentricity
double e=0;
while (e < 0.9)
{
//the inital conditions are determined based on equations (4.11) from the Book Computational Physics by Giordano and Nakanishi
//x corresponds to the minimum distance from Saturn of Hyperion's orbit, at this time velocity is at its maximum
//here (1.77e19/5.683e26) is the ratio between the masses of Hyperion and Saturn Respectively
double vy_0 = 2*M_PI*sqrt(((1+e)/(1-e))*(1+(1.77e19/5.683e26)));
double x = (1-e);
//the compute function generates a graph based on the inital conditions provided and a name corresponding to its eccentricity
compute(0,vy_0,x,0,15,0,0.01,25000, ("Orbit e = " + to_string(e)),true);
e += 0.1;
}
return 0;
}



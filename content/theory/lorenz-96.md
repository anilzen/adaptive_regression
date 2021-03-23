---
title: "Lorenz-96 Model"
weight: -20
---

The [Lorenz-96 model](https://en.wikipedia.org/wiki/Lorenz_96_model) is a chaotic, continuous-in-time, discrete-in-space, dynamical system that was proposed by Lorenz in 1996 as a toy model for weather dynamics.[^fn1].  
$$ \frac{dx_j}{dt} = ( x_{j+1} - x_{j-2} ) x_{j-1} - x_j + F(t). $$

The system is implemented in `helpers.py` in class `Lorenz`. The arguments of the class are
- `N`: Number of grid points.
- `dt`: Time step, defaults to 0.05. This is hard-coded for forecasts, so be careful when you change it. The error growth time scale is assumed to be such that 0.05 corresponds to 6 hours in an operational weather forecast system.
- `forcing`: External forcing function, defaults to a constant of 8 which is the chaotic regime.

The interior equations are implemented through the `rhs` method using vectorized numpy expressions.

```py
dotx[2:-1] = (x[3:]-x[0:-3])*x[1:-2] - x[2:-1] + forcing[2:-1]
```
Periodic boundary conditions are imposed on the outer grid points.

The Lorenz-96 system is solved using the `solve` method, which takes the number of days and initial data as input. The method populates the `sol` attribute by integrating the system using a 4th order Runge-Kutta method   
```py
self.sol[0] = init_data
for i in range(1, self.nt+1):
    self.sol[i] = self.rk4(self.sol[i-1])
    self.t += self.dt
```
The total number of time steps `nt` is determined from the days and the time step via `self.nt = int(4*0.05/self.dt*days)` as `dt=0.05` corresponds to 6 hours.

An `animate` method is provided to view the evolution of the solution in time.


[^fn1]: Lorenz, E.N., 1996. [Predictability: A problem partly solved](https://www.ecmwf.int/en/elibrary/10829-predictability-problem-partly-solved). *ECMWF Proc. Seminar on predictability* (Vol. 1, No. 1).
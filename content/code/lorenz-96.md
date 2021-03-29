---
title: "Lorenz-96 Model"
weight: -20
---

The [Lorenz-96 model](https://en.wikipedia.org/wiki/Lorenz_96_model) is a chaotic, continuous-in-time, discrete-in-space, dynamical system that was proposed by Lorenz in 1996 as a toy model for weather dynamics[^fn1].  
$$ \frac{dx_j}{dt} = ( x_{j+1} - x_{j-2} ) x_{j-1} - x_j + F(t). $$
Even though the equations have been used and studied for a long time, their chaotic behavior has been mathematically proven only recently[^fn2].

The system is implemented in `helpers.py` in class `Lorenz`. 

#### Parameters
- `N`: Number of grid points, defaults to 40.
- `dt`: Time step, defaults to 0.05. This is hard-coded for forecasts, so be careful when you change it. The error growth time scale is assumed to be such that 0.05 corresponds to 6 hours in an operational weather forecast system.
- `forcing`: External forcing function, defaults to the traditional choice of 8, where the system behaves chaotically.

#### Methods
- `rk4` is a 4th order Runge-Kutta integrator for solving the Lorenz system of ODEs.
- `rhs` implements the right hand side using vectorized numpy expressions.
```py
dotx[2:-1] = (x[3:]-x[0:-3])*x[1:-2] - x[2:-1] + forcing[2:-1]
```
Periodic boundary conditions are imposed on the outer grid points.
- `solve` solves the system of ODEs for a number of days with given initial data. The method populates the `sol` attribute by integrating the system using a 4th order Runge-Kutta method. The core portion of the method is
```py
self.sol[0] = init_data
for i in range(1, self.nt+1):
    self.sol[i] = self.rk4(self.sol[i-1])
    self.t += self.dt
```
The total number of time steps `nt` is determined from the days and the time step via `self.nt = int(4*0.05/self.dt*days)` as `dt=0.05` corresponds to 6 hours.
- `animate` is provided to view the dynamic solution in time.

#### Usage
The default system can be initiated with random data, solved for 90 days, and visualized by
```py
N=40; F=8;
initial_data=np.random.normal(0.25*F,0.5*F,N)
Default=Lorenz()
Default.solve(days=90,init_data=initial_data)
anim=Default.animate()
```

[^fn1]: Lorenz, E.N., 1996. [Predictability: A problem partly solved](https://www.ecmwf.int/en/elibrary/10829-predictability-problem-partly-solved). *ECMWF Proc. Seminar on predictability* (Vol. 1, No. 1).
[^fn2]: Bedrossian, J., Blumenthal, A. and Punshon-Smith, S., 2020. [A regularity method for lower bounds on the Lyapunov exponent for stochastic differential equations](https://arxiv.org/abs/2007.15827). arXiv preprint arXiv:2007.15827.
---
title: "Lorenz-96 Model"
weight: -20
---

The [Lorenz-96 model](https://en.wikipedia.org/wiki/Lorenz_96_model) is a chaotic, continuous-in-time, discrete-in-space, dynamical system that was proposed by Lorenz in 1996 as a toy model for weather dynamics.[^fn1].  
$$ \frac{dx_j}{dt} = ( x_{j+1} - x_{j-2} ) x_{j-1} - x_j + F. $$

The system is implemented in `helpers.py` in class `Lorenz`. The arguments of the class are
- `N`: Number of grid points.
- `F`: Constant external forcing term, defaults to 8 which is the chaotic regime.
- `days`, `dt`: The total number of time steps, `nt`, is determined from total days and the time step via `self.nt = 4*0.05/dt*self.days`. This setting is because the error growth time scale is assumed to be such that a `dt` of 0.05 corresponds to 6 hours in an operational weather forecast system.
- `bias`: Amplitude of bias.
- `noise`: Amplitude of noise.
- `pert_type`: Type of perturbation. Perturbation types are implemented through the method `update_perturbation` in the class. The implemented types are
  - `"None"`: No bias, corresponds to the Nature solution.
  - `"Bias"`: A constant (in time and space) bias applied to all grid points with amplitude `bias`. Corresponds to a constant shift of the forcing term. 
  - `"Noise"`: Normally distributed bias that varies randomly at each time step but is constant across the grid, implemented via `np.random.normal(0., self.noise, self.N)`. 
  - `"GridBias"`: 
  - `"TimeBias"`:

The interior equations are implemented through the `rhs` method using vectorized numpy expressions.

```py
dotx[2:-1] = (x[3:]-x[0:-3])*x[1:-2] - x[2:-1] + self.F + \
    + self.pert[2:-1]
```
Periodic boundary conditions are imposed on the outer grid points.

The Lorenz-96 system is solved using the `solve` method, which has an additional input `init_data`. The method populates the `sol` attribute by integrating the system with updated perturbations using a 4th order Runge-Kutta method   
```py
self.sol[0] = init_data
for i in range(1, self.nt+1):
    self.update_perturbation()
    self.sol[i] = self.rk4(self.sol[i-1])
```

An `animate` method is provided to view the evolution of the solution in time.

[^fn1]: Lorenz, E.N., 1996. [Predictability: A problem partly solved](https://www.ecmwf.int/en/elibrary/10829-predictability-problem-partly-solved). *ECMWF Proc. Seminar on predictability* (Vol. 1, No. 1).
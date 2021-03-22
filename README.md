# Model Output Statistics and Adaptive Regression

This repository contains Python code, Jupyter notebook, LaTeX source, and documentation for comparing Model Output Statistics (MOS) and Adaptive Regression (AR) for the Lorenz-96 model. 

The Jupyter notebook sets up an Observing Systems Simulation Experiment (OSSE) where Nature is given by the Lorenz-96 system. Biased versions of this system are used as models with which forecasts are made. Bias-correction is performed both by MOS and AR, based on an idea by Eugenia Kalnay using Kalman filters. 

Bias can be constant or stochastic. It can vary in space across the grid with a normal distribution or vary in time mimicking seasonal patterns. One can also use blended models by combining forecasts using different types of biases (each specific bias choice with a type and amplitude correponds to a mode). 

Statistics are collected for a 5-year nature run and compared to both MOS and AR forecasts. Forecast quality significantly depends on the type of bias. AR is most accurate for time-varying biases, which is the most realistic scenario.

Fore more information on usage, see the [online documentation](https://anilzen.github.io/adaptive_kalman/). 
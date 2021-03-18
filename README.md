# Model Output Statistics and Adaptive Kalman Filter for the Lorenz-96 model

An Observing Systems Simulation Experiment for forecasts using Model Output Statistics (MOS) and Adaptive Regression (AR). AR forecasts use Kalman filters based on an idea by Eugenia Kalnay. 

The forecasts are made for the Lorenz-40 model with external forcing and various choices of bias. The bias can be constant or stochastic. It can vary in space across the grid with a normal distribution or vary in time mimicking seasonal patterns. 

Statistics are collected for a 5-year nature run and compared to both MOS and AR forecasts. Forecast quality significantly depends on the type of bias. AR is most accurate for time-varying biases, which is the most realistic scenario.

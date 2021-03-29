---
title: "Forecast"
weight: -10
---

Forecasts are made using the `Forecast` class, which inherits from `Lorenz`. The time step and the number of grid points in the related Lorenz system are hard-coded to 40 and 0.05 respectively. 

#### Parameters
Forecasts are made from initial data provided by Nature at certain stations for a certain number of forecast days using a biased model. Therefore, the parameters are

- `Nature`: Nature reference array (`Lorenz.sol`). Forecast quality is compared to Nature. Nature also provides initial data for the forecasts.
- `stations`: Stations are represented as grid points of the model. Forecasts are made at those grid points. Since the model is periodic, the specific distribution of stations turns out to be not so important.
- `forecast_days`: The range of each forecast.
- `forcing`: The forcing for the forecast model should be different from the Nature forcing to introduce deliberate bias into the forecasts. 

### Making a deliberately biased forecast

The method `make_forecast` solves the biased Lorenz model up to the forecast day for each day in the Nature run. It populates the attributes `forecast` and `station_forecast`. For each day within the Nature run, the biased model is solved for the number of forecast days and the full solution is stored in `forecast`. So `Forecast.forecast` is an array shaped like `(total_days+1, 4*forecast_days+1, N)`. A subset of this array is stored under `station_forecast` corresponding to the location of the stations.

### Making corrected forecasts using MOS

The method `regression_pars` performs a linear, two-parameter regression. The predictor is the station forecast and the predictand is the Nature solution evaluated at these stations. The comparison between the Nature solution and the forecast provides the regression parameters. 


- `make_mos_forecast`

### Making corrected forecasts using AR


### Computing the quality of a forecast

- `error_norms`
- `process` is a shorthand method to trigger each method of the forecast so all forecast attributes are populated. 

#### Usage


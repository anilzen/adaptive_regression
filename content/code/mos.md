---
title: "MOS Forecast"
weight: 0
---

A MOS forecast starts with the computation of the regression parameters using the method `regression_pars`, which performs a linear, two-parameter regression. The predictor is the station forecast and the predictand is the Nature solution evaluated at these stations. The comparison between the Nature solution and the forecast provides the regression parameters. 

Writing $x$ as the biased forecasts for each forecast hour (one predictor), $\hat{y}$ as the corrected forecast (predictand), and $y$ as the Nature solution, we have 
$$ \hat{y} = \sum_k b_k x_k, \qquad y = \hat{y} + \epsilon, $$  
where $\epsilon$ is the forecast error and $k$ runs through 0,1 with $x_0=1$ and $x_1=x$. The forecast parameters $b_k$ are computed by requiring that the least-square error is minimized. Doing this for each forecast hour and representing the variables through matrices, we obtain for the matrix of regression parameters
$$ B = (X^T X)^{-1} X^T Y,$$ 
where $Y$ is the vector of Nature solutions for each forecast hour, and $X$ is the matrix of biased forecasts (predictor) with 1's at the first row. This is roughly implemented as follows.
```py
nt = 4*forecast_days+1
predictor = station_forecast
for k in range(len(predictor)):
    predictand[k] = Nature[4*k:4*k+nt, station]
for i in range(nt):
    X = vstack((ones(len(predictor)), predictor[:, i])).T
    b[i] = dot(dot(linalg.inv(dot(X.T, X)), X.T), predictand[:, i])
```
The main difficulty in the code is just accounting the rows and columns of the matrices.

The method `make_mos_forecast` uses the regression parameters to compute the forecast vector $\hat{Y}$ for each station via $\hat{Y}=X B$ which corresponds to the code
```py
mos_forecast = b[:, :, 0].T[np.newaxis, :]+b[:, :, 1].T*forecast[:, :, stations]
```

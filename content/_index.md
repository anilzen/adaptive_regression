# Adaptive Regression with Kalman Filters

This documentation describes an adaptive regression code that postprocesses numerical model output to correct systematic biases. The method is demonstrated on the [Lorenz-96 model](code/lorenz-96). Forecasts based on [Model Output Statistics](code/mos) and [Adaptive Regression](code/adaptive-regression) are compared in a 5-year run. 


The idea goes back to [Eugenia Kalnay](https://www2.atmos.umd.edu/~ekalnay/) who first described it in her [lecture notes](https://www2.atmos.umd.edu/~ekalnay/syllabi/AOSC630/METO630ClassNotes9MOS-KF.pdf) on statistical forecasting. Experiments with the Lorenz model support her conclusions at the end of her lecture notes: 

> Kalman Filtering provides a simple algorithm for adaptive regression. It requires little training so that it is able to adapt rather quickly to changes in the model, and to long-lasting weather regimes. It is particularly good in correcting model biases. However, in general it is not as good as regression based on long dependent samples.

Model Output Statistics performs better under ideal conditions with long dependent samples. However, there are a few reasons to prefer adaptive regression in an operational setting:
- Climate change is causing long-lasting weather regimes that are different from the usual weather patterns of the past.
- Numerical weather forecasting is in rapid development. Models are upgraded frequently to incorporate more detailed physics. Long-term training data for a given model may not be available and/or may be costly to obtain.
- Adaptive regression incorporates "errors of the day".

We explore these conclusions quantitatively using a simple Python implementation of adaptive regression. The documentation describes the parameters and methods implemented in the code. For usage, please refer to the Jupyter notebook [Lorenz.ipynb](https://github.com/anilzen/adaptive_regression/blob/main/code/Lorenz.ipynb).
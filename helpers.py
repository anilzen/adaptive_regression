import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np


class Lorenz:
    """
    Implements perturbed Lorenz model.
    """

    def __init__(self, days, N=40, dt=0.05, F=8., bias=0., noise=0., pert_type="None"):
        self.days = days
        self.N = N
        self.dt = dt
        self.nt = 4*self.days
        self.F = F
        self.noise = noise
        self.bias = bias
        self.gridbias = np.random.normal(0., noise, self.N)
        self.pert_type = pert_type
        self.pert = np.zeros(self.N)
        self.sol = np.zeros((self.nt+1, self.N))

    # Fourth order Runge-Kutta integrator for method of lines
    def rk4(self, cur):
        dt = self.dt
        # Runge-Kutta RK4
        k1 = dt * self.rhs(cur)
        k2 = dt * self.rhs(cur+0.5*k1)
        k3 = dt * self.rhs(cur+0.5*k2)
        k4 = dt * self.rhs(cur+k3)
        return cur+(k1+2.*(k2+k3)+k4)/6.

    # Right hand side for Lorenz 40 model
    def rhs(self, x):
        dotx = np.zeros(self.N)
        # Boundary equations (periodic boundary conditions)
        dotx[0] = (x[1]-x[-2])*x[-1] - x[0] + self.F + self.pert[0]
        dotx[1] = (x[2]-x[-1])*x[0] - x[1] + self.F + self.pert[1]
        dotx[-1] = (x[0]-x[-3])*x[-2] - x[-1] + self.F + self.pert[2]
        # Interior equations
        dotx[2:-1] = (x[3:]-x[0:-3])*x[1:-2] - x[2:-1] + self.F + \
            + self.pert[2:-1]
        return dotx

    def update_perturbation(self):
        if self.pert_type == "None":
            pass
        elif self.pert_type == "Bias":
            self.pert = self.bias*np.ones(self.N)
        elif self.pert_type == "Noise":
            self.pert = np.random.normal(0., self.noise, self.N)
        elif self.pert_type == "GridBias":
            self.pert = self.gridbias
        elif self.pert_type == "TimeBias":
            print("TimeBias not yet implemented!")

    def solve(self, init_data):
        if len(init_data) != self.N:
            print("STOP! Inconsistent data!")
        self.sol[0] = init_data
        for i in range(1, self.nt+1):
            self.update_perturbation()
            self.sol[i] = self.rk4(self.sol[i-1])

    def animate(self):
        fig = plt.figure()
        ax = plt.axes(xlim=(0, self.N), ylim=(self.sol.min(), self.sol.max()))
        line, = ax.plot(np.arange(self.N), self.sol[0], color='k', linewidth=1, marker='o', markersize=2)
        plt.xlabel('j')
        plt.ylabel('x')
        plt.legend(["Lorenz"])  # loc=3, frameon=False)

        def animate(i):
            line.set_ydata(self.sol[i])
        # call the animator.  blit=True means only re-draw the parts that have changed.
        return animation.FuncAnimation(fig, animate, frames=self.nt, interval=100, blit=False)


class Forecast(Lorenz):

    # Initialize forecast
    def __init__(self, total_days, forecast_days, Truth, stations, pert_type, bias=0., noise=0.):
        Lorenz.__init__(self, days=forecast_days, noise=noise, bias=bias, pert_type=pert_type)
        if Truth.shape[0]-4*forecast_days-1 != 4*total_days:
            print("Inconsistent Truth with total days.")
        self.total_days = total_days
        self.forecast_days = forecast_days
        self.Truth = Truth
        self.stations = stations
        self.forecast = np.zeros((total_days+1, 4*forecast_days+1, self.N))
        self.hours = 6*np.arange(4*self.forecast_days+1)

    # Make daily forecasts
    def make_forecast(self):
        for i in range(self.total_days+1):  # A forecast for each day
            self.solve(init_data=self.Truth[4*i])  # Perfect initial data (each day corresponds to 4 timesteps)
            self.forecast[i] = self.sol.copy()  # Biased evolution
        self.station_forecast = self.forecast[:, :, self.stations]

    def regression_pars(self):
        n_stat = len(self.stations)
        nt = 4*self.forecast_days+1
        self.b = np.zeros((n_stat, nt, 2))
        for j in np.arange(n_stat):
            predictor = self.station_forecast[:, :, j]
            predictand = np.zeros_like(predictor)
            for k in range(len(predictor)):
                predictand[k] = self.Truth[4*k:4*k+nt, self.stations[j]]
            for i in range(nt):
                H = np.vstack((np.ones(len(predictor)), predictor[:, i])).T
                self.b[j, i] = np.dot(np.dot(np.linalg.inv(np.dot(H.T, H)), H.T), predictand[:, i])

    def make_mos_forecast(self):
        self.mos_forecast = self.b[:, :, 0].T[np.newaxis, :]+self.b[:, :, 1].T*self.forecast[:, :, self.stations]

    def error_norms(self):
        errors = np.zeros_like(self.station_forecast)
        mos_errors = np.zeros_like(errors)
        for i in range(self.total_days+1):
            StationTruth = self.Truth[4*i:4*i+4*self.forecast_days+1, self.stations]
            errors[i] = self.station_forecast[i] - StationTruth
            mos_errors[i] = self.mos_forecast[i] - StationTruth
        self.errors = np.linalg.norm(errors, axis=2).mean(axis=0)
        self.mos_errors = np.linalg.norm(mos_errors, axis=2).mean(axis=0)

    def process(self):
        self.make_forecast()
        self.regression_pars()
        self.make_mos_forecast()
        self.error_norms()


def kalman(h, bak, Pak, yok, Qpar, rk):
    # Parameters of Kalman Filter
    # Q is a diagonal matrix for uncorrelated errors with the variance
    # of each predictor: KFNature[:,stations].var(axis=0)
    # Try 0 or identity matrix
    Q = Qpar*np.identity(2)
  #  Q  = np.diag((Qpar, Qpar))
    # rk is observational error covariance. Chosen 0 or very small.
    # Kalman Filter Algorithm
    yfk = np.dot(h.T, bak)
    ek = yok - yfk
    Pfk = Pak + Q
    wk = np.dot(np.dot(h.T, Pfk), h) + rk
    kk = np.dot(Pfk, h)/wk
    bakup = bak + np.dot(kk, ek)
    Pakup = Pfk - np.dot(np.dot(kk, h.T), Pfk)
    return ek, bakup, Pakup


def adaptive_regression(Forecast, Truth):
    period, nt, nstat = Forecast.shape
    bak = np.zeros((period+1, nstat, nt, 2))
    bak[0] = regression(models[0][:101], Nature[:419, stations])
    Pak = np.zeros((nstat, nt, 2, 2))
    for i in range(nstat):
        for j in range(nt):
            Pak[i, j] = np.identity(2)
    ek = np.zeros((period+1, nstat, nt))

    for j in range(n_stat):
        predor = StationForecast[:, :, j]
        predant = np.zeros_like(predor)
        for k in range(len(predor)):
            predant[k] = Truth[4*k:4*k+nt, j]
        for i in range(nt):
            H = np.vstack((np.ones(len(predor)), predor[:, i])).T
            b[j, i] = np.dot(np.dot(np.linalg.inv(np.dot(H.T, H)), H.T), predant[:, i])
    return b

# The observe function adds random perturbations
# with a Gaussian distribution to the solution


def observe(Xt, sigma):
    Xo = Xt+np.random.normal(0, sigma, np.size(Xt))
    return Xo

# The H_observe function provides a partial observation
# with linear dependence on true values.
# The resulting vector y will have zero values
# H is a matrix
# def H_observe(Xt, H):
#    Yo=np.dot(H,Xt)
#    return Yo

# Construction of background covariance matrix from observation samples
# (equivalent to numpy function cov)


def bg_cov(samples):
    nx = np.size(samples[0])
    bg_covar = np.zeros((nx, nx))
    sample_mean = np.mean(samples, axis=0)
    for i in range(nx):
        for j in range(i, nx):
            a = 0
            for k in range(len(samples)):
                a += (samples[k, i]-sample_mean[i])*(samples[k, j]-sample_mean[j])
#                a+=(samples[k,i]-init[i])*(samples[k,j]-init[j])
            bg_covar[i, j] = a/(len(samples)-1)
    # Symmetrize
    return bg_covar + bg_covar.T - np.diag(bg_covar.diagonal())

# Check that the covariance matrix is symmetric, posdef, and full rank


def check_Pb(Pb):
    K = len(Pb[0])
    if ((Pb.T == Pb).all() and (np.linalg.eigvals(Pb) > 0).all()
            and np.linalg.matrix_rank(Pb) == K):
        return True
    else:
        return False

# Optimal Interpolation algorithm


def OI(Xb, Pb, Yo, Ro, H):
    K_OI = np.dot(np.dot(Pb, H.T), np.linalg.inv(np.dot(np.dot(H, Pb), H.T)+Ro))
    dob = Yo-np.dot(H, Xb)
    Xa = Xb+np.dot(K_OI, dob)
    # doa = Yo-np.dot(H,Xa)
    # Pa = np.dot((np.eye(np.size(Xb))-np.dot(K_OI,H)),Pb)
    return Xa  # , 0.5*(Pa.T+Pa), dob, doa

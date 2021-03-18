# Construct MOS for biased models
from helpers import l40_step, regression, \
    model_forecasts, compute_error_norms, kalman
import numpy as np
import matplotlib.pylab as plt
##############
# Parameters
##############
# set general parameters
dt = 0.05   # time step
K  = 40     # number of grixyd points 
F  = 8.     # F... forcing
alpha = 0.  # alpha... coefficient of noise
pars = [dt, K, F, alpha, np.zeros(K)]
# Set model parameters
stations = np.array([0,10,20,30])
bias=2; noise = 1000
gridbias = np.random.normal(0., 4., K)
# Construct parameter array for each model
parsB=pars[:]; parsN = pars[:]; parsNB=pars[:]; parsGB=pars[:]
parsB[2] += bias # bias
parsN[3] += noise # noise
parsNB[2] += bias # bias
parsNB[3] += noise # noise
parsGB[4] = gridbias
parsTest = parsB[:]
##############
# Initialization
##############
# Random initial data
X0=np.random.normal(0.25*F,0.5*F,K)
# Reinitialize by removing transient effects (90 day run)
nt=4*90+1 # Note the +1 for nt for the use of the range function.
for i in range(1,nt):
    X0 = l40_step(X0, pars)
##############
# Nature
##############
# Generate a 2 year nature solution and add 3 more days for forecasts
long = 2 * 360 * 4
nt = 5*4+1 # Note the +1 for nt for the use of the range function.
hours = 6*np.arange(nt)
Nature = np.zeros((long+nt, len(X0)))
Nature[0] = X0
for i in range(1,long+nt):
    Nature[i] = l40_step(Nature[i-1], pars)
##############
# Compute regression parameters
##############                                                                  
model = model_forecasts(long, nt, Nature, parsTest)[:,:,stations]
b = regression(model, Nature[:,stations])
error_norms = compute_error_norms(model, Nature[:,stations])
mos_forecast = b[:,:,0].T[np.newaxis,:]+b[:,:,1].T*model
mos_error_norm = compute_error_norms(mos_forecast, Nature[:,stations])
##############
# Construct Kalman regressionOOO
##############                                                                  
# Initialize and set Kalman parameters and matrices
short = 5*4
model = model_forecasts(short, nt, Nature, parsTest)[:,:,stations]
bak = np.zeros((int(long/4)+2, 4, nt, 2))
bak[0] = regression(model, Nature[:,stations])
#bak[0] = b
Pak = np.zeros((len(stations), nt, 2, 2))
ek  = np.zeros((int(long/4.)+1, len(stations), nt))
emodel = np.zeros_like(ek)
emos = np.zeros_like(ek)
for i in range(len(stations)):
    for j in range(nt):
        Pak[i,j] = np.identity(2)
KFmodel = np.zeros((int((long+nt)/4)+1, nt, pars[1]))
# Construct independent nature run for comparing MOS to adaptive Kalman
Nature[0] = Nature[-1]
for i in range(1,nt-1):
    Nature[i] = l40_step(Nature[i-1], pars)
    if i%4==0:
        k = int(i/4)
        for j in range(nt):
            if j==0:    KFmodel[k, j] = Nature[i]
            else:       KFmodel[k, j] = l40_step(KFmodel[k,j-1], parsTest)
for i in range(nt-1, long+nt):
    Nature[i] = l40_step(Nature[i-1], pars)
    if i%4==0:
        k = int(i/4)
        # For future prediction, run the model into the future
        for j in range(nt):
#            parsTest[2]+=.3*np.sin(2*np.pi*(i+j)/31)
            if j==0:    KFmodel[k, j] = Nature[i]
            else:       KFmodel[k, j] = l40_step(KFmodel[k,j-1], parsTest)
            for i_st,station in enumerate(stations):
                # Perfectly observed nature
                yok = Nature[4*(k-1)+j%4,station]
                y_model = KFmodel[k-1-int(j/4), j, station]
                # Model forecast error
                emodel[k-5,i_st,j] = yok - y_model
                # MOS forecast error
                mos_fcs = b[i_st, j, 0]+b[i_st, j, 1]*y_model
                emos[k-5,i_st,j] = yok - mos_fcs
                # KF update of regression parameters and forecast error
                h = np.array([1,y_model])
                ek[k-5,i_st,j], bak[k-5+1,i_st,j], Pak[i_st,j] = \
                    kalman(h, bak[k-5, i_st,j], Pak[i_st, j], yok, 0.0001*1.2**j, 1.)
    
hours = 6*np.arange(nt)
plt.plot(hours,np.linalg.norm(emodel[60:],axis=0).mean(axis=0), label='Model')
plt.plot(hours,np.linalg.norm(emos[60:],axis=0).mean(axis=0), label='MOS')
plt.plot(hours,np.linalg.norm(ek[60:],axis=0).mean(axis=0), label='Adaptive')
#plt.plot(hours,np.linalg.norm(ek_check[60:],axis=0).mean(axis=0), label='Check')
plt.title("Comparison of errors: Bias")
plt.ylabel("Forecast Error Norm")
plt.xlabel("Forecast Hour")
plt.legend(loc='upper left')
plt.xticks(hours[::2])
#plt.savefig("ComparisonBias.pdf")
plt.show()

#f, (ax1, ax2) = plt.subplots(1, 2, sharex = True, figsize=(12,4))
#for i in range(0,nt,2):
#    ind1=(0,i,0)
#    ind2=(0,i,1)
#    ax1.plot(bak.T[ind1[::-1]]-b[ind1])
#    ax2.plot(bak.T[ind2[::-1]]-b[ind2])
#ax1.set_xticks(np.arange(0,int(long/4+1),180))
#ax1.set_ylabel("Predictor 1: $b_0$")
#ax2.set_ylabel("Predictor 2: $b_1$")
#f.suptitle('Difference of KF predictors to MOS predictors at one station', fontsize=12, va='center')
#plt.tight_layout()
##plt.savefig("EvolutionBias.pdf")

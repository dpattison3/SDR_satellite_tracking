#%%
# import tools
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy
import plotly.graph_objects as go


#%%
# Define some useful constants
# how many samples to cut out of the frequency transitions
gate = 200
# sample rate of the measurements
sampleRate = 2e6


#%%
# import the binary file
rawData_benGS = np.fromfile(
    "/Users/benjaminpattison/Documents/Projects/satNav/SDR_satellite_tracking/data/test5_benGS.dat",
    dtype="uint8",
)
rawData_domGS = np.fromfile(
    "/Users/benjaminpattison/Documents/Projects/satNav/SDR_satellite_tracking/data/test5_domGS.dat",
    dtype="uint8",
)

# subtract 255/2 to go from unsigned to signed values
benData = rawData_benGS - 127.5
domData = rawData_domGS - 127.5

# split the data into its real and imaginary parts
benRealParts = benData[0::2]
benImagParts = benData[1::2]
domRealParts = domData[0::2]
domImagParts = domData[1::2]

# create a variable that has the real and imaginary parts together
benCplx = benRealParts + 1j * benImagParts
domCplx = domRealParts + 1j * domImagParts

#%%
# number of samples for each frequency
ben_nS = np.size(benCplx, 0)
dom_nS = np.size(domCplx, 0)
if ben_nS != dom_nS:
    print("You don't have the same number of samples from each ground station")
    exit()
else:
    nS_total = dom_nS
    nS_each = int(nS_total / 3)

#%%
# Correlation of the two reference signals using the differential phase method discussed here: http://www.panoradio-sdr.de/correlation-for-time-delay-analysis/
# and referencing some of the code from the matlab scripts in his library
"""
Transmitter location according to wikipedia
https://geohack.toolforge.org/geohack.php?pagename=KUFX&params=37.205_N_121.950_W_type:landmark_region:US-CA_source:FCC
repeaters here
    37.151583    -121.609917 (SE of SJ)
    37.659389    -121.933028 (NE of Fremont)
"""
# The first reference
# diff the signal
dPhase1_benGS = np.diff(np.unwrap(np.angle(benCplx[: nS_each - gate])))
dPhase1_domGS = np.diff(np.unwrap(np.angle(domCplx[: nS_each - gate])))

# remove the mean
dPhase1_benGS = dPhase1_benGS - np.mean(dPhase1_benGS)
dPhase1_domGS = dPhase1_domGS - np.mean(dPhase1_domGS)

# correlate the signals
dPhaseXcorr1 = signal.correlate(dPhase1_benGS, dPhase1_domGS)
dPhaseXcorr1_lags = signal.correlation_lags(len(dPhase1_benGS), len(dPhase1_domGS))
dPhaseXcorr1_max = np.max(dPhaseXcorr1)

refSig1Lag_samples = np.abs(dPhaseXcorr1_lags[np.argmax(dPhaseXcorr1)])
refSig1Lag_time = refSig1Lag_samples / sampleRate

# same process, for the second signal
dPhase2_benGS = np.diff(np.unwrap(np.angle(benCplx[-nS_each + gate :])))
dPhase2_domGS = np.diff(np.unwrap(np.angle(domCplx[-nS_each + gate :])))

dPhase_benGS = dPhase2_benGS - np.mean(dPhase2_benGS)
dPhase_domGS = dPhase2_domGS - np.mean(dPhase2_domGS)

dPhaseXcorr2 = signal.correlate(dPhase2_benGS, dPhase2_domGS)
dPhaseXcorr2_lags = signal.correlation_lags(len(dPhase2_benGS), len(dPhase2_domGS))
dPhaseXcorr_max2 = np.max(dPhaseXcorr1)

refSig2Lag_samples = np.abs(dPhaseXcorr2_lags[np.argmax(dPhaseXcorr2)])
refSig2Lag_time = refSig2Lag_samples / sampleRate

# print out the two lags (in units of seconds)
print(
    f"For the first reference: the lag in time is {refSig1Lag_time} and the lag in samples is {refSig1Lag_samples}"
)
print(
    f"For the second reference: the lag in time is {refSig2Lag_time} and the lag in samples is {refSig2Lag_samples}"
)

#%%
# plot the differential phase correlation
"""
xVals1 = dPhaseXcorr1_lags / 2e6
xVals2 = dPhaseXcorr1_lags / 2e6

fig1 = go.Figure()
fig1.add_trace(
    go.Scatter(x=xVals1, y=dPhaseXcorr1, mode="lines", name="Reference 1 correlation")
)
fig1.add_trace(
    go.Scatter(x=xVals2, y=dPhaseXcorr2, mode="lines", name="Reference 2 correlation")
)

fig1.show()

plt.plot(xVals, dPhaseXcorr1)
plt.plot(xVals, dPhaseXcorr2)
plt.xlabel("Lag [seconds]")
plt.ylabel("Not sure... correlation maybe?")
plt.xlim((-2, 2))
plt.show()
"""

#%%
"""#%%
# define center freq (Hz)
centerFreq = 98.5e6

# +/- bandwidth for plotting (Hz)
bwdth = 400e3

# define sample frequency
sampFreq = 2.0e6
# %%
plt.specgram(cplx, NFFT=4096 * 8, Fs=sampFreq, Fc=centerFreq)
plt.title("PSD of signal")
plt.xlabel("Time")
plt.ylabel("Frequency")
plt.ylim([centerFreq - bwdth, centerFreq + bwdth])
plt.show()

# %%
plt.psd(cplx, NFFT=4096 * 8, Fs=sampFreq, Fc=centerFreq)
plt.ylim([-42, 0])
plt.yticks(ticks=[-45, -40, -35, -30, -25, -20, -15, -10, -5, 0])
plt.title("PSD of signal")

# %%

# %%
"""
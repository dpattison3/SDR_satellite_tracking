#%%
# import statements
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy

#%%
# Define some useful constants
# how many samples to cut out of the frequency transitions
gate = 200

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
"""
Transmitter location according to wikipedia
https://geohack.toolforge.org/geohack.php?pagename=KUFX&params=37.205_N_121.950_W_type:landmark_region:US-CA_source:FCC
repeaters here
    37.151583    -121.609917 (SE of SJ)
    37.659389    -121.933028 (NE of Fremont)
"""
# correlation of the two signals for just the first reference sample. this method uses the standard complex correlation talked about here:
# http://www.panoradio-sdr.de/correlation-for-time-delay-analysis/
xCorr = signal.correlate(benCplx[:nS_each], domCplx[:nS_each])
lags = signal.correlation_lags(len(benCplx[:nS_each]), len(domCplx[:nS_each]))

#%%
# this method uses the differential phase talked about at the same link on the first reference frequency sample
dPhase_benGS = np.diff(np.unwrap(np.angle(benCplx[: nS_each - gate])))
dPhase_domGS = np.diff(np.unwrap(np.angle(domCplx[: nS_each - gate])))

# remove the mean
dPhase_benGS = dPhase_benGS - np.mean(dPhase_benGS)
dPhase_domGS = dPhase_domGS - np.mean(dPhase_domGS)

# append a zero to the beginning to make the length of the abs correlation (size fit)
dPhase_benGS = np.append(0, dPhase_benGS)
dPhase_domGS = np.append(0, dPhase_domGS)

dPhaseXcorr = signal.correlate(dPhase_benGS, dPhase_domGS)
dPhaseXcorr_lags = signal.correlation_lags(len(dPhase_benGS), len(dPhase_domGS))
ref1 = max(signal.correlate(dPhase_benGS, dPhase_benGS))
ref2 = max(signal.correlate(dPhase_domGS, dPhase_domGS))
dPhaseXcorr_max = np.max(dPhaseXcorr)

refSig1Lag = np.abs(dPhaseXcorr_lags[np.argmax(dPhaseXcorr)]) / 2e6

print(refSig1Lag)

#%%
# this method uses the differential phase talked about at the same link on the second reference frequency sample
dPhase2_benGS = np.diff(np.unwrap(np.angle(benCplx[-nS_each + gate :])))
dPhase2_domGS = np.diff(np.unwrap(np.angle(domCplx[-nS_each + gate :])))

# remove the mean
dPhase_benGS = dPhase2_benGS - np.mean(dPhase2_benGS)
dPhase_domGS = dPhase2_domGS - np.mean(dPhase2_domGS)

# append a zero to the beginning to make the length of the abs correlation (size fit)
dPhase2_benGS = np.append(0, dPhase2_benGS)
dPhase2_domGS = np.append(0, dPhase2_domGS)

dPhaseXcorr2 = signal.correlate(dPhase2_benGS, dPhase2_domGS)
dPhaseXcorr2_lags = signal.correlation_lags(len(dPhase2_benGS), len(dPhase2_domGS))
# ref1 = max(signal.correlate(dPhase_benGS, dPhase_benGS))
# ref2 = max(signal.correlate(dPhase_domGS, dPhase_domGS))
dPhaseXcorr_max2 = np.max(dPhaseXcorr)

refSig2Lag = np.abs(dPhaseXcorr2_lags[np.argmax(dPhaseXcorr2)]) / 2e6

print(refSig2Lag)
#%%
# plot the differential phase correlation
xVals = dPhaseXcorr_lags / 2e6 + refSig1Lag
plt.plot(xVals, dPhaseXcorr)
plt.plot(xVals, dPhaseXcorr2)
plt.xlabel("Lag [seconds]")
plt.ylabel("Not sure... correlation maybe?")
plt.xlim((-2, 2))
plt.show()


#%%
# plot the complex correlation.
xVals = np.linspace(0, np.size(xCorr), np.size(xCorr))
plt.plot(np.abs(xCorr))
plt.show()
plt.plot(lags, np.abs(xCorr))
plt.show()

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
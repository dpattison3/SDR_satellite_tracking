#%%
# import statements
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

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
sampleRate = 2e6
xCorr = signal.correlate(benCplx[:nS_each], domCplx[:nS_each])
lags = signal.correlation_lags(len(benCplx[:nS_each]), len(domCplx[:nS_each]))

#%%
# plot the correlation.
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
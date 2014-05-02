import matplotlib.pyplot as plt
import numpy as np
# import scipy as sp
from scipy import interpolate
import ReadData
import libMATERIALS


Ge_table, Ge_start, dU_Ge, N_Ge1, N_Ge2 = libMATERIALS.initialize_ge()
TGe = interpolate.interp1d(Ge_table[2, :], Ge_table[0, :], 'slinear')

Tx = np.zeros((100))
Ex = np.zeros((100))
ct2 = 0
for iter in range(10000, 173000, 1000):
    ct, L, dL, Ne, qfluxL, qfluxC, qfluxR = ReadData.readdata1(iter)
    qL, qC, qR, Tavg, Eavg = ReadData.readdata2(iter, Ne[0], Ne[1], Ne[2])
    ct2 = ct2 + 1
    for i in range(Ne[0]):
        Ex[i] = Ex[i] + np.sum(Eavg[i, :, :]) / (Ne[1] * Ne[2])


Ex = Ex / ct2 / (dL[0] * dL[1] * dL[2])
Tx = TGe(Ex)


x = np.linspace(0.5 * dL[0], L[0] - 0.5 * dL[0], Ne[0])


plt.figure()
plt.plot(x, Tx)
plt.show()

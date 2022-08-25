import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

pixels = 400
sigma_l = 10.0 / (2 * np.pi)

def dist(coord, norm=1./float(pixels), center=np.array((pixels/2 - 0.5, pixels/2 - 0.5))):
    spot = coord - center
    dist = np.hypot(spot[0], spot[1])
    norm_dist = dist * norm
    return norm_dist

def H(omega, sigma=sigma_l):
    nuem = 2.0 * np.sqrt(2.0 * np.pi) * omega * sigma
    denom = erf(np.sqrt(2.0) * np.pi * omega * sigma)
    return nuem / denom

filter = np.zeros((pixels, pixels))
for i in range(len(filter[:,0])):
    for j in range(len(filter[0,:])):
        filter[i,j] = H(dist(np.array((i,j))))

filter = np.fft.fftshift(filter)

fig, ax = plt.subplots()
ax.imshow(filter)
plt.show()
np.save("lor_rendering\\4_pi_filter", filter)

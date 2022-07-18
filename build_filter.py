import numpy as np
import matplotlib.pyplot as plt


# image = np.loadtxt("lor_rendering\\filter_dot.data")
image = np.loadtxt("lor_rendering\\220_filter.data")


# image = image.reshape((400,400))

image = image.reshape((220,220,220))
image = image[:,:,110]


transformed = np.fft.fftn(image)

inverse = 1 / transformed

np.save("lor_rendering\dot_inverse220", inverse)

fig, ax = plt.subplots()
fig, ay = plt.subplots()
fig, az = plt.subplots()
ax.imshow(np.real(np.fft.ifftn(1/inverse)))
ay.imshow(np.real(inverse))
az.imshow(np.imag(inverse))
plt.show()

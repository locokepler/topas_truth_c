from distutils.command.clean import clean
from inspect import cleandoc
import numpy as np
import matplotlib.pyplot as plt

# source_file = "lor_rendering\\500000_3cyl"
# source_file = "lor_rendering\\dot_9_y_stat"
# source_file = "lor_rendering\\Derenzo10000"
# source_file = "lor_rendering\\Derenzo1000000bg"
# source_file = "lor_rendering\\xcat_run1"
# source_file = "lor_rendering\\xcat_run2_sigma"
# source_file = "lor_rendering\\xcat_run3_2_sigma"
# source_file = "lor_rendering\\derenzo_600"
# source_file = "lor_rendering\\derenzo_600_thin"
# source_file = "lor_rendering\\derenzo_12"
# source_file = "lor_rendering\\derenzo_no_bg"
# source_file = "lor_rendering\\derenzo_2_500"
source_file = "lor_rendering\\derenzo_2_5000"



pixels = 400

image = np.load(source_file + ".npy")

print(np.shape(image))

slice = 100
if (len(np.shape(image)) > 2):
    image = image[:,:,slice]
else:
    image = image.reshape((400,400))



fig, pre = plt.subplots()
pre_ = pre.imshow(image)
pre_bar = plt.colorbar(pre_)

image_F = np.fft.fftn(image)

filter = np.load("lor_rendering\dot_inverse.npy")
# filter = np.load("lor_rendering\dot_inverse220.npy")

symmetry = 2

filter = 1 * filter

cutoff_low = symmetry
cutoff_high = pixels - symmetry
# filter[cutoff_low:cutoff_high,cutoff_low:cutoff_high] = 1
# filter[cutoff_high:,:] = 1
# filter[:,cutoff_high:] = 1
# filter[:cutoff_low, :] = 1
# filter[:, :cutoff_low] = 1

filter[cutoff_low:cutoff_high, :] = 1
filter[:, cutoff_low:cutoff_high] = 1

# now filter out the edge behavior

fig, ax = plt.subplots()
# fig, ay = plt.subplots()
ax.imshow(np.real(filter))
# ay.imshow(np.real(image_F))
# plt.show()

cleaned_F = np.multiply(image_F, filter)

symmetry = 30

cutoff_low = symmetry
cutoff_high = pixels - symmetry

cleaned_F[cutoff_low:cutoff_high, :] = 1
cleaned_F[:, cutoff_low:cutoff_high] = 1

cleaned = np.fft.ifftn(cleaned_F)

fig, az = plt.subplots()

cleaned = np.real(cleaned)
# cleaned = np.where(cleaned > 0,cleaned, 0.0)

az_ = az.imshow(cleaned)
az_bar = plt.colorbar(az_)
plt.show()

np.save(source_file + "_clean", cleaned)
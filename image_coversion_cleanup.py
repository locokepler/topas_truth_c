import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as image



# source_file = "lor_rendering\\xcat_run2_sigma"
# source_file = "lor_rendering\\xcat_run3_2_sigma"
# source_file = "lor_rendering\\derenzo_600"
# source_file = "lor_rendering\\derenzo_600_thin"
# source_file = "lor_rendering\\derenzo_12"
# source_file = "lor_rendering\\derenzo_no_bg"
# source_file = "lor_rendering\\derenzo_2_500"
source_file = "lor_rendering\\derenzo_2_5000"


pixels = 400

source_image = np.loadtxt(source_file + ".data")

print(np.shape(source_image))
source_image = source_image.reshape((pixels,pixels,np.shape(source_image)[1]))
source_image = image.gaussian_filter(source_image, 1)

# source_image = source_image.reshape((400,400))
# source_image = image.gaussian_filter(source_image, 1)

print(np.shape(source_image))


np.save(source_file, source_image)
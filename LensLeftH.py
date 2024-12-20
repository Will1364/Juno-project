import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter

img1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/jupRef/LeftH/IMD_580411682.672409.bmp",0)
img2=cv2.imread("C:/Users/Bruger/Desktop/30330 project/jupRef/LeftH/IMD_580411713.422424.bmp",0)
img3=cv2.imread("C:/Users/Bruger/Desktop/30330 project/jupRef/LeftH/IMD_667200558.422409.bmp",0)

scale1 = np.median(img3)/np.median(img1)
scale2 = np.median(img3)/np.median(img2)

Img2 = cv2.convertScaleAbs(img2, alpha=scale2)
Img1 = cv2.convertScaleAbs(img1, alpha=scale1)

stacked_areas = np.stack([Img1,Img2,img3], axis=0)
median_distortionMIN = np.min(stacked_areas, axis=0).astype(np.uint8)
median_distortionMED = np.median(stacked_areas, axis=0).astype(np.uint8)

plt.imshow(median_distortionMIN, cmap='gray')
plt.colorbar()
plt.show()


plt.imshow(median_distortionMED, cmap='gray')
plt.colorbar()
plt.show()

#cv2.imwrite('C:/Users/Bruger/Desktop/30330 project/LensLeftHMin.bmp', median_distortionMIN)

#cv2.imwrite('C:/Users/Bruger/Desktop/30330 project/LensLeftH.bmp', median_distortionMED)
import cv2
import numpy as np
import matplotlib.pyplot as plt

#Load images
PJ22N1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530842.422424.bmp",0)
PJ22N2=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530842.922424.bmp",0)
PJ22N3=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530872.922424.bmp",0)
PJ22N4=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530873.422409.bmp",0)

PJ23N1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ23 - JupNorth/IMD_626090721.172424.bmp",0)
PJ23N2=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ23 - JupNorth/IMD_626090721.672424.bmp",0)
PJ23N3=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ23 - JupNorth/IMD_626090722.172424.bmp",0)


J221=PJ22N1#[:, 0:50]
J222=PJ22N2#[:, 0:50]
J223=PJ22N3#[:, 0:50]
J224=PJ22N4#[:, 0:50]

J231=PJ23N1#[:, 0:50]
J232=PJ23N2#[:, 0:50]
J233=PJ23N3#[:, 0:50]

#Stack images to 3D array and collapse them into a median image
stacked_areas = np.stack([J221,J222,J223,J224,J231,J232,J233], axis=0)
median_distortion = np.median(stacked_areas, axis=0).astype(np.uint8)
plt.imshow(median_distortion, cmap='gray')
plt.colorbar()
plt.show()

J221m=PJ22N1[:, 400:579]
J222m=PJ22N2[:, 400:579]
J223m=PJ22N3[:, 400:579]
J224m=PJ22N4[:, 400:579]

J231m=PJ23N1[:, 400:579]
J232m=PJ23N2[:, 400:579]
J233m=PJ23N3[:, 400:579]

#Find the median of an area uncurrupted by the column
stacked_areas_m = np.stack([J221m,J222m,J223m,J224m,J231m,J232m,J233m], axis=0)
median_m = np.median(stacked_areas_m, axis=0).astype(np.uint8)
MEDIAN=np.median(median_m)

#Subtract the median from the column image
medSub=median_distortion-MEDIAN
plt.imshow(medSub, cmap='gray')
#plt.savefig('C:/Users/Bruger/Desktop/30330 project/søjle profiler/MedianToSubtract01.png', format='png')#, dpi=dpi, bbox_inches='tight')
plt.colorbar()
plt.show()
#Save image
#cv2.imwrite('C:/Users/Bruger/Desktop/30330 project/MedianToSubtract1.bmp', medSub)

#PJ22N1[:,0:50]=PJ22N1[:,0:50]-medSub
#plt.imshow(PJ22N1,cmap='gray')

PJ22N12=PJ22N1-medSub
plt.imshow(PJ22N12,cmap='gray')


PJ22N13=PJ22N1-median_distortion
PJ22N13[PJ22N1-median_distortion<1]=0
#plt.imshow(PJ22N13,cmap='gray')


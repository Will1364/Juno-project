import cv2
import numpy as np
import matplotlib.pyplot as plt

img1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/temp/PJ08 - Adrastea And Rings/IMD_557575033.922424.bmp",0)#Adresta and Rings
img2 = cv2.imread("C:/Users/Bruger/Desktop/30330 project/temp/PJ20 - Amalthea/IMD_612388582.422424.bmp",0)#Amalthea
img3 = cv2.imread("C:/Users/Bruger/Desktop/30330 project/temp/PJ21 - Ganymedes and JupNorth and JupSouth/IMD_616956708.672424.bmp",0)#Ganymede and Jup
img1 = cv2.flip(img1, 0)
img1=img1/5.4
img2=img2/1.4


#COL=cv2.imread("C:/Users/Bruger/Desktop/30330 project/MedianToSubtract1.bmp",0)*5

stacked_areas = np.stack([img1,img2,img3], axis=0)
median_distortionMIN = np.min(stacked_areas, axis=0).astype(np.uint8)
median_distortionMED = np.median(stacked_areas, axis=0).astype(np.uint8)
#median_distortion=np.median([median_distortionMIN,median_distortionMED], axis=0).astype(np.uint8)
#median_distortion[400:-1,:]=np.min([COL[400:-1,:],median_distortion[400:-1,:]], axis=0).astype(np.uint8)
Jup=median_distortionMED-median_distortionMIN
median_distortion=median_distortionMED
median_distortion[400:-1,:]=median_distortionMED[400:-1,:]-Jup[400:-1,:]
plt.imshow(median_distortion, cmap='gray')
plt.show()

#cv2.imwrite('C:/Users/Bruger/Desktop/30330 project/LensToSubtract.bmp', median_distortion)

PJ22N1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530902.922424.bmp",0) #J nederst venstre
plt.imshow(PJ22N1-median_distortion*1.5,cmap='gray')
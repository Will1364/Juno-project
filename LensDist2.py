import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter

#Load the images
img1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/temp/PJ08 - Adrastea And Rings/IMD_557575033.922424.bmp",0)#Adresta and Rings
img2 = cv2.imread("C:/Users/Bruger/Desktop/30330 project/temp/PJ20 - Amalthea/IMD_612388582.422424.bmp",0)#Amalthea
img3 = cv2.imread("C:/Users/Bruger/Desktop/30330 project/temp/PJ21 - Ganymedes and JupNorth and JupSouth/IMD_616956708.672424.bmp",0)#Ganymede and Jup
img1 = cv2.flip(img1, 0)


#Scale the images
scale1 = np.median(img1)/np.median(img2)
scale2 = np.median(img1)/np.median(img3)

Img2 = cv2.convertScaleAbs(img2, alpha=scale1)
Img3 = cv2.convertScaleAbs(img3, alpha=scale2)


#Stack images and collapse into min and med
stacked_areas = np.stack([img1,Img2,Img3], axis=0)
median_distortionMIN = np.min(stacked_areas, axis=0).astype(np.uint8)
median_distortionMED = np.median(stacked_areas, axis=0).astype(np.uint8)
#Subtract images to isolate jupiter
Jup=median_distortionMED-median_distortionMIN

#Scale jupiter and set small values to 0
Jup[Jup<90]=0
Jup[Jup>2]=Jup[Jup>2]-90
scaleJup=np.mean(Jup[400:-1,:])/np.mean(median_distortionMED[400:-1,:])


MED=cv2.convertScaleAbs(median_distortionMED, alpha=scaleJup)
Jup2=cv2.convertScaleAbs(Jup, alpha=scaleJup)

#Subtract Jupiter from median collapsed image
median_distortion = median_distortionMED - Jup
MED = MED - Jup2

plt.imshow(Jup, cmap='gray')
plt.colorbar()
plt.show()

plt.imshow(MED, cmap='gray')
plt.colorbar()
plt.show()
plt.imshow(median_distortion, cmap='gray')
plt.colorbar()
plt.show()

#Select column and apply a gaussian filter
y1=gaussian_filter(median_distortion[:,300], 10)
#Set threshold for clamping 
threshold = 15

# Find the index where the curve reaches close to the minimum
min_val = np.min(y1)
clamp_start_idx = np.where(y1 < min_val + threshold)[0][0]

# Blend towards the minimum smoothly
y_smooth_clamped1 = np.copy(y1)
blend_factor = np.linspace(1, 0, len(y1) - clamp_start_idx)
y_smooth_clamped1[clamp_start_idx:] = (
    y1[clamp_start_idx:] * blend_factor + min_val * (1 - blend_factor)
)

#Apply flat clamping
min_idx = np.argmin(y1)
y_clamped = np.copy(y1)
y_clamped[min_idx:] = y1[min_idx]

plt.plot(median_distortion[:,300],label="Original Curve")
plt.plot(y1, label="Gauss smoothed Curve")
plt.plot( y_smooth_clamped1, label="Smoothly Clamped Curve")
plt.plot(y_clamped, label="Flat Clamped Curve")
plt.legend()
plt.xlabel("Pixels in column")
plt.ylabel("Intensity")
plt.title("Intensity curve of column 300")
plt.show()


#Make ampty image and apply the same procedure on larger scale
#This is done several times to make it look smoother and to compare images
MEDNEW=np.full((580,752),0)

for x in range(len(MED[0,:])):
    y = savgol_filter(median_distortion[:,x], window_length=34, polyorder=2)  # Example curve with a minimum
    y=gaussian_filter(median_distortion[:,x], 10)


    threshold = 1  # Adjust this value as needed

    # Find the index where the curve reaches close to the minimum
    min_val = np.min(y)
    clamp_start_idx = np.where(y < min_val + threshold)[0][0]

    # Blend towards the minimum smoothly
    y_smooth_clamped = np.copy(y)
    blend_factor = np.linspace(1, 0, len(y) - clamp_start_idx)
    y_smooth_clamped[clamp_start_idx:] = (
        y[clamp_start_idx:] * blend_factor + min_val * (1 - blend_factor)
    )
    #MEDNEW=np.append(MEDNEW, y_smooth_clamped,axis=1)
    MEDNEW[:,x]=y_smooth_clamped

for x in range(len(MED[0,:])):
    y = savgol_filter(MEDNEW[:,x], window_length=34, polyorder=2)  # Example curve with a minimum

    y=gaussian_filter(MEDNEW[:,x], 10)

    threshold = 15  # Adjust this value as needed

    # Find the index where the curve reaches close to the minimum
    min_val = np.min(y)
    clamp_start_idx = np.where(y < min_val + threshold)[0][0]

    # Blend towards the minimum smoothly
    y_smooth_clamped = np.copy(y)
    blend_factor = np.linspace(1, 0, len(y) - clamp_start_idx)
    y_smooth_clamped[clamp_start_idx:] = (
        y[clamp_start_idx:] * blend_factor + min_val * (1 - blend_factor)
    )
    #MEDNEW=np.append(MEDNEW, y_smooth_clamped,axis=1)
    MEDNEW[:,x]=y_smooth_clamped

plt.imshow(MEDNEW, cmap='gray')
plt.title("Modelling Lens Reflection Using Smoothed Clamping")
plt.colorbar()
plt.show()

#Repeated with flat clamping
MEDNEW2=np.full((580,752),0)

for x in range(len(MED[0,:])):
    #y = savgol_filter(median_distortion[:,x], window_length=34, polyorder=2)  # Example curve with a minimum

    y=gaussian_filter(median_distortion[:,x], 10)

    min_idx = np.argmin(y)

    # Clamp the values after the minimum
    y_clamped = np.copy(y)
    y_clamped[min_idx:] = y[min_idx]
    MEDNEW2[:,x]=y_clamped
    # Plot to visualize
    
plt.imshow(MEDNEW2, cmap='gray')
plt.colorbar()
plt.show()

for x in range(len(MED[0,:])):
    #y = savgol_filter(MEDNEW2[:,x], window_length=34, polyorder=2)  # Example curve with a minimum

    y=gaussian_filter(MEDNEW2[:,x], 10)

    threshold = 15  # Adjust this value as needed

    # Find the index where the curve reaches close to the minimum
    min_val = np.min(y)
    clamp_start_idx = np.where(y < min_val + threshold)[0][0]

    # Blend towards the minimum smoothly
    y_smooth_clamped = np.copy(y)
    blend_factor = np.linspace(1, 0, len(y) - clamp_start_idx)
    y_smooth_clamped[clamp_start_idx:] = (
        y[clamp_start_idx:] * blend_factor + min_val * (1 - blend_factor)
    )
    #MEDNEW=np.append(MEDNEW, y_smooth_clamped,axis=1)
    MEDNEW2[:,x]=y_smooth_clamped

plt.imshow(MEDNEW2, cmap='gray')
plt.colorbar()
plt.title("Modelling Lens Reflection Using Flat Clamping")
plt.show()

#cv2.imwrite('C:/Users/Bruger/Desktop/30330 project/Lense2.bmp', MEDNEW2)

PJ22N1=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530902.922424.bmp",0)
A1=PJ22N1-MEDNEW2*0.3197969543147208
A1[A1<0]=0
plt.imshow(A1, cmap='gray')
plt.colorbar()
plt.show()

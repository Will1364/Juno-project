import cv2
import matplotlib.pyplot as plt 

#Requires photutils - can be acquired by pip install
from photutils.aperture import CircularAperture
from photutils.aperture import aperture_photometry
#from astropy import units as u


#Example code to count pixels in an apperture

Img=cv2.imread("C:/Users/Bruger/Desktop/30330 project/PJ22/PJ22 - JupNorth/IMD_621530842.422424.bmp",0)
#Img=Img*u.adu #Can in theory give the values the unit of adu, but then it makes a weird format

aperture_ref = CircularAperture([170, 400], 13) #x and y coord. of center and radius of aperture
#angle = Angle(0,'deg') #create degree of tilt if needed
#apertureB = EllipticalAperture([857, 673], 290, 190, angle) #x and y coord. of ellipse, semimajor axis, semiminor axis, tilt

phot_table = aperture_photometry(Img, aperture_ref) 
print(phot_table)
F_ref=phot_table['aperture_sum']


plt.imshow(Img, cmap='gray')
plt.title("Example Image")
ap = aperture_ref.plot(color='white', lw=2, label='Photometry aperture')

#aparent magnitude:
#m_object=-2.5*np.log10(flux_object/F_ref)+m_ref

#absolute magnitude (if relevant?)
#D=some distance from celestia
#M=m_object - 5*np.log10(D/10)


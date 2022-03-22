from plotting import plot2d
import matplotlib.pyplot as plt

ncfile = r'c:\Users\weste_bt\OneDrive - Stichting Deltares\Documents\GitHub\aeolis-python\examples\barchan\aeolis.nc' 

fig, ax = plot2d(ncfile, itimes=0, ifrac=0, param='zb', 
                 cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)

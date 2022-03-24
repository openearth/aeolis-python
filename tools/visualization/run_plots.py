from plotting import plot1d, plot2d
import matplotlib.pyplot as plt

ncfile = r'c:\Users\weste_bt\OneDrive - Stichting Deltares\Documents\GitHub\aeolis-python\examples\barchan\aeolis.nc' 


fig, ax = plot2d(ncfile, itimes=310, ifrac=0, param='zb', 
                  cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)

# fig, ax = plot2d(ncfile, itimes=-1, ifrac=0, param='ustars', 
                  # cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)
                  
# fig, ax = plot1d(ncfile, itransects=75, itimes=[0, 1], ifrac=0, params=['zb', 'zsep'])
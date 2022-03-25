from plotting import plot1d, plot2d, plot3d
import matplotlib.pyplot as plt

ncfile = r'c:\Users\weste_bt\OneDrive - Stichting Deltares\Documents\GitHub\aeolis-python\examples\parabolic\aeolis2.nc' 


# fig, ax = plot2d(ncfile, itimes=310, ifrac=0, param='zsep', 
                  # cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)

fig, ax = plot2d(ncfile, itimes=-1, ifrac=0, param='zb', 
                    cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)

# fig, ax = plot3d(ncfile, itimes=-1, scalez=1., vangle=60., hangle=30.)
                   
# fig, ax = plot1d(ncfile, itransects=75, itimes=-1, ifrac=0, params=['zb', 'zsep'])
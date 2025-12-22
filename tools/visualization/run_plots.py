from plotting import plot1d, plot2d, plot3d
import matplotlib.pyplot as plt

# ncfile = r'c:\Users\weste_bt\aeolis\Tests\RotatingWind\Barchan_Grid270\aeolis.nc'
# ncfile = r'c:\Users\weste_bt\aeolis\case6_AeolianTransport_2Bart\aeolis_base_ne.nc'
ncfile = r'c:\Users\weste_bt\aeolis\ZM\aeolis_small_1yr.nc'

# fig, ax = plot2d(ncfile, itimes=310, ifrac=0, param='zsep', 
                  # cmap=plt.cm.jet, clim='auto', delta=False, itimeDelta=0)

fig, ax = plot2d(ncfile, itimes=0, ifrac=0, param='zb', cmap=plt.cm.viridis, clim='auto', delta=True, itimeDelta=-1)
plt.show()

# fig, ax = plot3d(ncfile, itimes=-1, scalez=1., vangle=60., hangle=30.)
                   
# fig, ax = plot1d(ncfile, itransects=75, itimes=-1, ifrac=0, params=['zb', 'zsep'])
# fig, ax = plot1d(ncfile, itransects=75, itimes=[0, -1], ifrac=0, params='zb')
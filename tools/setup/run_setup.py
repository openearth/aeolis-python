from setup_tools import create_grd, create_z_cone

# path = r'c:\Users\weste_bt\OneDrive - Stichting Deltares\Documents\01-Projects\AeoLiS\Tests\RotatingWind\Barchan_Grid355' 
path = r'c:\Users\weste_bt\OneDrive - Stichting Deltares\Documents\GitHub\aeolis-python\examples\barchan'


create_grd(path, nx=300, ny=150, x0=0., y0=0., dx=1., dy=1., alfa=355., plotBool=True)

# create_z_cone(path, nx=300, ny=150, height=6., diameter=60., midx=75, midy=75, dx=1., dy=1.)



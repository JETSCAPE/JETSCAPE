import pandas as pd

#open the freezeout surface file (MUSIC format)
# tau, x , y, eta, ds0, ds1, ds2, ds3, ....
surf_tau = pd.read_csv( 'surface.dat', header=None, delimiter=' ', usecols=[0] )
surf_x = pd.read_csv( 'surface.dat', header=None, delimiter=' ', usecols=[1] , names=['surf_x'])
surf_y = pd.read_csv( 'surface.dat', header=None, delimiter=' ', usecols=[2] , names=['surf_y'])
surf_eta = pd.read_csv( 'surface.dat', header=None, delimiter=' ', usecols=[3] )

surf_x_list = surf_x['surf_x'].tolist()
surf_y_list = surf_y['surf_y'].tolist()

#default minimum value of r_max in fm
#every hydro grid will be at least (2 r_max)fm x (2 r_max)fm 
r_max = 5.0

x_max = 0.0
y_max = 0.0 

for i in range( 0, len(surf_x_list) ):
    x_max = max( x_max, abs( surf_x_list[i]) )

for i in range( 0, len(surf_y_list) ):
    y_max = max( y_max, abs( surf_y_list[i]) )


x_y_max = max(x_max, y_max)
r_max = max(r_max, x_y_max)

print( "x_max = " + str(x_max) + ", y_max = " + str(y_max) + ", r_max = " + str(r_max) )


import numpy as np

def calc_ecc():
    file = open('freestream_input', mode = 'r', encoding = 'utf-8-sig')
    counter=0
    for line in file:
        line = line.split(' ')
        if(counter==13): #get dx from file
           dx=line[1]
        if(counter==14): #get dy from file
           dy=line[1]
        counter+=1
    file.close()

    dx=np.float64(dx)
    dy=np.float64(dy)

    #e=np.loadtxt('test.dat')
    #load T^00 block file
    e=np.loadtxt('output/initial_e_projection.dat')

    N_x=np.size(e,1) # number of points in x
    N_y=np.size(e,0) # number of points in y

    zmax=6

    eps=np.zeros((zmax,zmax),dtype=complex)
    eps_r3=0
    eps_r5=0

    for iy in range(0,N_y):
     y=dy*(2*iy-N_y)/2
     for ix in range(0,N_x):
      x=dx*(2*ix-N_x)/2

      z=np.complex(x,y)
      zbar=np.conjugate(z)

      for i in range(0,zmax):
       for j in range(0,zmax):
        powz=np.power(z,i)
        powzbar=np.power(zbar,j)
        powz_powzbar=powz*powzbar
        eps[i][j]+=e[iy][ix]*powz*powzbar

    for iy in range(0,N_y):
     y=dy*(2*iy-N_y)/2
     for ix in range(0,N_x):
      x=dx*(2*ix-N_x)/2

      z=np.complex(x,y)
      eps_r3+=e[iy][ix]*np.power(np.absolute(z - (    eps[1][0]/ eps[0][0])),3)
      eps_r5+=e[iy][ix]*np.power(np.absolute(z - (    eps[1][0]/ eps[0][0])),5)


    eps*=dx*dy
    eps_r3*=dx*dy
    eps_r5*=dx*dy

    for i in range(0,zmax):
     for j in range(0,zmax):
      if(not(i==0 and j==0)):  eps[i][j] /= eps[0][0]
    
    f=open("ecc_new.dat","w")
    for i in range(0,zmax):
     for j in range(0,zmax):
      f.write(str(i)+"	"+str(j)+"	"+"("+str(np.real(eps[i][j]))+","+str(np.imag(eps[i][j]))+")"+"\n")
    f.close()

    f=open("ecc_r3_new.dat","w")
    f.write(str(eps_r3/eps[0][0]))
    f.close()

    f=open("ecc_r5_new.dat","w")
    f.write(str(eps_r5/eps[0][0]))
    f.close()
    
    return eps

if __name__ == '__main__':
    calc_ecc()

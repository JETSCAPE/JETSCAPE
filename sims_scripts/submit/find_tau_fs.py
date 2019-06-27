import numpy as np

def calc_tau_fs(tau_0, e_0, alpha):
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

    #load T^00 block file
    e=np.loadtxt('output/initial_e_projection.dat')

    #integrate over the transverse plane to find dE/d\eta_s
    dE_dn = np.sum(e) * dx * dy

    tau_fs = tau_0 * ( (dE_dn / e_0) ** alpha )

    return tau_fs

if __name__ == '__main__':
    calc_tau_fs()

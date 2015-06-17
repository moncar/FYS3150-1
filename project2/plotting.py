from pylab import*
import sys

def reader(filename):
    readfile = open(filename, 'r')
    
    values = []
    for line in readfile:
        if line.startswith("n"):
            n_line = line.strip()
            n = float(n_line[3:])

        elif line.startswith("rho_max"):
            rho_max_line = line.strip()
            rho_max = float(rho_max_line[9:])

        elif line.startswith("rho_min"):
            rho_min_line = line.strip()
            rho_min = float(rho_min_line[9:])

        else:
            values.append(float(line))

    readfile.close()

    return n, rho_max, rho_min, array(values)

def normalise(eigenvec, rho_max, rho_min, n):
    h = float((rho_max - rho_min))/n
    
    eig_square = eigenvec**2; nv = 0
    for i in range(len(eigenvec)-1):
        nv += (h/2)*(eig_square[i+1] + eig_square[i])
    return eig_square/sqrt(nv)

def set_values(eigval, eigvec, rho_max, rho_min, n):
    norm_eigvec = normalise(eigvec, rho_max, rho_min, n)

    return norm_eigvec

mode = 'd'

if mode == 'd':
    try:
        name = sys.argv[1]
        titlename = sys.argv[2]
    except IndexError:
        raise IndexError("specify name!")

    n, rho_max, rho_min, values = reader(name)

    eigenvalues = values[:n]
    eigenvector1 = values[n:2*n]
    normalised_eigenvec1 = set_values(eigenvalues, eigenvector1, 
            rho_max, rho_min, n)

    eigenvector2 = values[2*n:3*n]
    normalised_eigenvec2 = set_values(eigenvalues, eigenvector2, 
            rho_max, rho_min, n)

    eigenvector3 = values[3*n:4*n]
    normalised_eigenvec3 = set_values(eigenvalues, eigenvector3, 
            rho_max, rho_min, n)

    rho = values[4*n:]

    plot(rho, normalised_eigenvec1)
    hold('on')
    plot(rho, normalised_eigenvec2)
    plot(rho, normalised_eigenvec3)
    title(titlename)
    legend(["ground state","state2","state3"])
    show()

elif mode == 'c':
    try:
        namec = sys.argv[1:]
    except IndexError:
        raise IndexError("c failed")
        
    #omega_r=0.01
    n, rho_max, rho_min, values1 = reader(namec[0])
    eigenvalue = values1[:n]
    eigenvector1 = values1[n:2*n]
    norm_eig1 = set_values(eigenvalue, eigenvector1,
            rho_max, rho_min, n)
    
    rho = values1[4*n:]

    #omega_r=0.5
    a,b,c, values2 = reader(namec[1])
    eigenvector2 = values2[n:2*n]
    norm_eig2 = set_values(eigenvalue, eigenvector2,
            rho_max, rho_min, n)

    #omega_r=1.0
    a,b,c, values3 = reader(namec[2])
    eigenvector3 = values3[n:2*n]
    norm_eig3 = set_values(eigenvalue, eigenvector3,
            rho_max, rho_min, n)
    
    #omega_r=5.0
    a,b,c, values4 = reader(namec[3])
    eigenvector4 = values4[n:2*n]
    norm_eig4 = set_values(eigenvalue, eigenvector4,
            rho_max, rho_min, n)

    plot(rho, norm_eig1)
    hold('on')
    plot(rho, norm_eig2)
    plot(rho, norm_eig3) 
    plot(rho, norm_eig4)
    title("Varying omega_r")
    legend(["0.01","0.5","1.0","5.0"])
    show()

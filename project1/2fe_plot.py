from pylab import *
import sys

def opener(fname, test):
    openfile = open(fname, 'r')

    c = 0; 
    for line in openfile: 
        if line.startswith('n'): 
            continue
        c += 1
            
    openfile = open(fname, 'r')
            
    values = zeros(c)
    
    i = 0
    for line in openfile:
        if line.startswith('n'):
            n = int(line[4:])
            continue
        values[i] = float(line)
        i += 1
    
    openfile.close()
    
    if test == False:
        v = values[:n]; x = values[n:]
        return v, x
    elif test == True:
        return values

def comm(line):
    from subprocess import PIPE, Popen
    p = Popen(args=line, stdout=PIPE, shell=True)
    return p.communicate()[0]

try:
    filename = sys.argv[1:]
except IndexError:
    raise IndexError("specify results")

if filename[0] == "tri_res.dat":
    val, x10 = opener(filename[0], False)
    x_ec = linspace(0,1,10E6)
    u = 1 - (1 - exp(-10))*x_ec - exp(-10*x_ec)
    plot(x10,val)
    hold('on')
    plot(x_ec, u)
    xlabel("x=ih"); ylabel("u(x)")
    legend(["Approx","Exact"])

if filename[0] == "tri_error.dat":
    err_val = opener(filename[0], True)
    N = linspace(0,1000000,6)
    plot(N,err_val)
    hold('on')
    title("Error-estimate")
    xlabel("N"); ylabel("eps(N)")
show()

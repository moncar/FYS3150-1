from pylab import *
import sys

def analytic(x,t,N):
    """ Formula for closed-form solution """
    S = 1-x
    for i in range(1,N):
        S += -(2/(pi*i))*sin(i*pi*x)*exp(-i**2*pi**2*t) #closed-form solution

    return S

namefile = sys.argv[1]

filename = ["t1","t2","t3","t4","t5"]
for name in filename:
    """Fancy read from file and assign to array"""
    exec "{0} = fromfile(\"{0}.txt\",sep=\" \"); {0} = {0}.reshape((len({0})/1.,1))".format(name)

n = len(t1)
x = linspace(0,1,n) #x-values

plot(x,t1, 'r', x,t2,'r', x,t3,'r', x,t4,'r', x,t5, 'r')
#plot(x,analytic(x,0.002,100), x,analytic(x,0.05,100), 
 #       x,analytic(x,0.1,100), x,analytic(x,0.2,100), x,analytic(x,0.5,100))

xlabel('x'); ylabel('u(x,T)')
legend(["t=0.002","t=0.05","t=0.1","t=0.2","t=0.5", "te=0.002","te=0.05","te=0.1","te=0.2","te=0.5"])
savefig(namefile)
show()

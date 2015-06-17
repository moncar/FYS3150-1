from pylab import *
import scipy, sys

name = sys.argv[1]
exec "{0} = fromfile(\"{0}.txt\",sep=\" \"); {0} = {0}.reshape((len({0})/2.,2))".format(name)

def analytical_solution(x,t,N):
    S = 1-x
    for n in xrange(1,N):
        S += -(2./(pi*n))*sin(n*pi*x)*exp(-n**2*pi**2*t)
                            
    return S


x = test_1D[:,0]; xProbPos = test_1D[:,1];
probx = xProbPos/max(xProbPos)

time = float(sys.argv[2])
x2 = linspace(0,1,len(xProbPos))

probx_analytical = analytical_solution(x2,time,100)

ylab = "u(x,T)" + " " + "(" + "scaled by max:" + str(max(xProbPos)) + ")"
plot(x2,probx)
hold('on')
plot(x2,probx_analytical)
xlabel('position'); ylabel(ylab)
legend(["simulated", "analytical"])
axis([0,1,0,1])
#savefig(sys.argv[3])
show()

err_test = True
if err_test == True:
    plot(x2,abs(probx-probx_analytical))
    xlabel("position"); ylabel("difference from analytical")
    savefig(sys.argv[3])
    show()


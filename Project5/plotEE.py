from pylab import *
import sys

def solution(x,y,t):
    return (1-y)*exp(x+t)

method = sys.argv[1]

if method == "a":
    n = 100
    x = linspace(0,1,n)
    xi,yi = meshgrid(x,x);
    sol = solution(xi,yi,0.5)
    pcolor(x,x,sol)
    cb = colorbar()
    cb.set_label("u(x,y,t)")
    xlabel("x:"); ylabel("y")
    savefig("analytic2D")
    show()
elif method == "e":
    err_test = True
    filename = "expEuler2D"
    exec "{0} = fromfile(\"{0}.txt\",sep=\" \"); {0} = {0}.reshape((len({0})/1.,1))".format(filename)
    
    n = int(expEuler2D[0]); 
    
    expEuler2D = expEuler2D[1:]
    expEuler2D = expEuler2D.reshape(len(expEuler2D))
    
    expResults = zeros((n,n))
    for i in range(n-1):
        expResults[:,i] = expEuler2D[n*i:n*(i+1)]

    x = linspace(0,1,n)
    xi,yi = meshgrid(x,x);
    pcolor(xi,yi,expResults)
    cb = colorbar()
    cb.set_label("u(x,y,T)")
    savefig("explicit2D")
    show()
    if err_test == True:
        pcolor(xi,yi,abs(expResults-solution(xi,yi,0.5))) #plot difference
        cb = colorbar()
        savefig("explicit2D_error_05")
        show()
elif method == "i":
    err_test = True
    filename = "impJacobi2D"
    exec "{0} = fromfile(\"{0}.txt\",sep=\" \"); {0} = {0}.reshape((len({0})/1.,1))".format(filename)
    n = int(impJacobi2D[0]); 
    impJacobi2D = impJacobi2D[1:]
    impJacobi2D = impJacobi2D.reshape(len(impJacobi2D))
    
    impResults = zeros((n,n))
    for i in range(n-1):
        impResults[:,i] = impJacobi2D[n*i:n*(i+1)]
    
    x = linspace(0,1,n)
    xi,yi = meshgrid(x,x);
    pcolor(xi,yi,impResults)
    cb = colorbar()
    cb.set_label("u(x,y,T)")
    xlabel("x"); ylabel("y")
    savefig("implicit2D")
    show()
    if err_test == True:
        pcolor(xi,yi,abs(impResults-solution(xi,yi,0.5))) #plot difference
        cb = colorbar()
        savefig("limplicit2D_error_05")
        show()
elif method == "acid":
    x = linspace(0,1,exp_n)
    for i in range(exp_n):
        plot(x,expResults[:,i])
        hold('on')
    show()
elif method == "MC":
    filename = "mc2D"
    exec "{0} = fromfile(\"{0}.txt\",sep=\" \"); {0} = {0}.reshape((len({0})/1.,1))".format(filename)
    n =  int(mc2D[0])
    mc2D = mc2D[1:]
    mc2D = mc2D.reshape(len(mc2D))
    
    tResults = zeros((n,n))
    for i in range(n-1):
        tResults[:,i] = mc2D[n*i:n*(i+1)]
    
    for j in range(n-1):
        m = max(tResults[:,j])
        if m == 0:
            continue
        else:
            tResults[:,j] /= m

    x = linspace(0,1,n)
    xi,yi = meshgrid(x,x)
    pcolor(xi,yi,tResults)
    cb = colorbar()
    cb.set_label("u(x,y,T)")
    savefig("lMC2Dl_02")
    show()
else:
    print "specify method"

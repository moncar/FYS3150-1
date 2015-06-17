from pylab import *
import sys

filename = sys.argv[1]
en = fromfile(filename,sep=" ")

t = linspace(0,1,len(en))

hold('on')
xlabel('t (years)')
ylabel('Total energy')
plot(t,en)
show()

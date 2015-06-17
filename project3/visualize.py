from pylab import *

textfiles = ["sun","mercury","venus","earth","mars","jupiter","saturn","uranus","neptune","pluto"]

for t in textfiles:
    exec "{0} = fromfile(\"{0}.txt\",sep=\" \"); {0} = {0}.reshape((len({0})/3.,3))".format(t)

#xkcd(); #awesomeness
subplots_adjust(right=0.77)

plot(sun[:,0], sun[:,1], 'black')
plot(mercury[:,0], mercury[:,1], 'orange')
plot(venus[:,0], venus[:,1], 'green')
plot(earth[:,0], earth[:,1], 'blue')
plot(mars[:,0], mars[:,1], 'red')
plot(jupiter[:,0], jupiter[:,1], 'brown')
plot(saturn[:,0], saturn[:,1], 'yellow')
plot(uranus[:,0], uranus[:,1], 'grey')
plot(neptune[:,0], neptune[:,1], 'cyan')
plot(pluto[:,0], pluto[:,1], 'magenta')
xlabel('rx (AU)'); ylabel("ry (AU)")
axis('equal')

legend(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Neptune","Uranus","Pluto"],
        labelspacing=0.3, bbox_to_anchor=[1,1.05],
        loc=2, numpoints=1)

show()

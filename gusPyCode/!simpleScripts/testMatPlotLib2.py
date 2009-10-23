from pylab import *
t = arange(0.0, 5.2, 0.2)

# red dashes, blue squares and green triangles
plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
show()

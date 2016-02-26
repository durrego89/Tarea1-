import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------Runge Kutta


d = np.loadtxt("rk4.dat", skiprows=2,  unpack=True)


dato2 = d[2]

dato0 = d[3]

dato4 = d[4]


pos= abs(dato0)<0.1

dato2a = dato2[pos]
dato2b = dato4[pos]


#grafica 1
s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(-3.2, 3.2)
plt.ylim(-2.5, 2.5)
plt.scatter(dato2a,dato2b,s) # One plot

plt.title(u"Metodo Runge Kutta 4th")
plt.xlabel("$q_{3}$") 
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("rk4.png")
plt.close()

#----------------------------------------------------------grafica 2




s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(-2.9, -0.75)
plt.ylim(-1, 1)
plt.scatter(dato2a,dato2b,s) # One plot

plt.title(u"Zoom 1 Runge Kutta")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("Zoom_1_rk4.png")
plt.close()

#----------------------------------------------------------grafica 3

s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(-0.2, 0.75)
plt.ylim(-0.3, 0.3)
plt.scatter(dato2a,dato2b,s) # One plot

plt.title(u"Zoom 2 Runge Kutta")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("Zoom_2_rk4.png")
plt.close()


#--------------------------------------------------------- Metodo Simpletico

g = np.loadtxt("simpletico.dat", skiprows=2,  unpack=True)

dat0 = g[3]
dat2 = g[2]
dat4 = g[4]


pos= abs(dat0)<0.1

dat2a = dat2[pos]
dat2b = dat4[pos]



s=0.01
fig=plt.figure()

plt.xlim(-3.2, 3.2)
plt.ylim(-2.5, 2.5)

plt.scatter(dat2a,dat2b,s)

plt.title(u"Metodo LeapFrog")
plt.xlabel("$q_{3}$") 
plt.ylabel("$p_{3}$")

#plt.show()
fig.savefig("simpletico.png")
plt.close()

#----------------------------------------------------------grafica 2




s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(-2.9, -0.75)
plt.ylim(-1, 1)
plt.scatter(dat2a,dat2b,s) # One plot

plt.title(u"Zoom 1 Simpletico")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("Zoom_1_Simpletico.png")
plt.close()

#----------------------------------------------------------grafica 3

s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(-0.2, 0.75)
plt.ylim(-0.3, 0.3)
plt.scatter(dat2a,dat2b,s) # One plot

plt.title(u"Zoom 2 Simpletico")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("Zoom_2_Simpletico.png")
plt.close()

#--------------------------------------------------------- Metodo Simpletico a=0.4325

g = np.loadtxt("simpletico2.dat", skiprows=2,  unpack=True)

dat0 = g[3]
dat2 = g[2]
dat4 = g[4]


pos= abs(dat0)<0.1

dat21a = dat2[pos]
dat21b = dat4[pos]



s=0.01
fig=plt.figure()

plt.xlim(-2, 2)
plt.ylim(-2.5, 2.5)

plt.scatter(dat21a,dat21b,s)

plt.title(u"Metodo LeapFrog a = 0.4325")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")

#plt.show()
fig.savefig("simpletico_04325.png")
plt.close()

#----------------------------------------------------------grafica 2




s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(0.2, 0.7)
plt.ylim(-0.7,0.7)
plt.scatter(dat2a,dat2b,s) # One plot

plt.title(u"Zoom Simpletico (0.45,0) y a = 0.4325")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("Zoom_1_Simpletico_04325.png")
plt.close()

#----------------------------------------------------------grafica 3

s=0.03
fig=plt.figure() # Create and empty figure
plt.xlim(-0.08, 0.08)
plt.ylim(-0.1, 0.1)
plt.scatter(dat2a,dat2b,s) # One plot

plt.title(u"Zoom 2 Simpletico a = 0.4325")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")


#plt.show()
fig.savefig("Zoom_2_Simpletico_04325.png")
plt.close()

#--------------------------------------------------------- Metodo Simpletico a=0.425

g = np.loadtxt("simpletico3.dat", skiprows=2,  unpack=True)

dat0 = g[3]
dat2 = g[2]
dat4 = g[4]


pos= abs(dat0)<0.1

dat21a = dat2[pos]
dat21b = dat4[pos]



s=0.01
fig=plt.figure()

plt.xlim(-3, 3)
plt.ylim(-2.5, 2.5)

plt.scatter(dat21a,dat21b,s)

plt.title(u"Metodo LeapFrog a = 0.425")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")

#plt.show()
fig.savefig("simpletico_0425.png")
plt.close()

#--------------------------------------------------------- Metodo Simpletico a=0.46

g = np.loadtxt("simpletico4.dat", skiprows=2,  unpack=True)

dat0 = g[3]
dat2 = g[2]
dat4 = g[4]


pos= abs(dat0)<0.1

dat21a = dat2[pos]
dat21b = dat4[pos]



s=0.01
fig=plt.figure()

plt.xlim(-3, 3)
plt.ylim(-2.5, 2.5)

plt.scatter(dat21a,dat21b,s)

plt.title(u"Metodo LeapFrog a = 0.46")
plt.xlabel("$q_{3}$")
plt.ylabel("$p_{3}$")

#plt.show()
fig.savefig("simpletico_046.png")
plt.close()

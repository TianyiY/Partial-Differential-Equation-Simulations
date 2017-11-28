import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

x=[i for i in range(0, 1000)]

m=9.1e-31    # particle mass - kg
planck_constant=6.6e-34    # J*s
h_bar=planck_constant/2/math.pi

L=2.*math.pi*0.053e-9    # orbit length - m
B=math.sqrt(2./L)   # amplitude - m
n=1.  # energy level

E=(math.pi**2*h_bar**2*n**2/(2.*m*L**2))   # definite energy
k=math.sqrt(2.*m*E/h_bar**2)
p=k*h_bar   # definite momentum


wave_init=[]
for i in range(len(x)):
    res1=B*math.sin(math.pi*n*x[i]/L)
    wave_init.append(res1)

time=5000
I=complex(0, 1)
E_time=cmath.exp(-I*E*time/h_bar)

phi_xt=[]
for i in range(len(wave_init)):
    res2=E_time*wave_init[i]
    phi_xt.append(res2)

phi_xt_real=[]
for i in range(len(phi_xt)):
    res3=phi_xt[i].real
    phi_xt_real.append(res3)

phi_xt_imag=[]
for i in range(len(phi_xt)):
    res4=phi_xt[i].imag
    phi_xt_imag.append(res4)

phi_xt_prob=[]
total=0
for i in range(len(phi_xt)):
    res5=phi_xt_real[i]**2+phi_xt_imag[i]**2
    phi_xt_prob.append(res5)
    total+=res5

test=0
for i in range(len(phi_xt_prob)):
    phi_xt_prob[i]=phi_xt_prob[i]/total
    test+=phi_xt_prob[i]

assert test>=0.999 and test<=1.001

fig=plt.figure()
ax1=fig.add_subplot(311)
ax1.plot(phi_xt_real[:100])
ax1.set_title("real part")
ax2=fig.add_subplot(312)
ax2.plot(phi_xt_imag[:100])
ax2.set_title("image part")
ax3=fig.add_subplot(313)
ax3.plot(phi_xt_prob[:100])
ax3.set_title("probability")
plt.show()

np.savetxt("wave_t5000_box.txt", phi_xt)
np.save('wave_t5000_box.npy', phi_xt)


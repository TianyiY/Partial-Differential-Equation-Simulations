import math
import cmath
import matplotlib.pyplot as plt
import numpy as np

x=[i for i in range(0, 1000)]

m=9.1e-31    # particle mass - kg
planck_constant=6.6e-34    # J*s
h_bar=planck_constant/2/math.pi

A=0.15e-9   # m
B=0.19e-9   # m

L=2.*math.pi*0.053e-9    # orbit length - m
n=1.   # energy level

E=(2.*math.pi*h_bar*n/L)**2/(2.*m)  # definite energy
k=math.sqrt(2.*m*E/h_bar**2)
p=k*h_bar   # definite momentum

kx=[]
for i in range(len(x)):
    res1=x[i]*k
    kx.append(res1)

I=complex(0, 1)

Ikx=[]
for i in range(len(kx)):
    res2=I*kx[i]
    Ikx.append(res2)

phi_x=[]
for i in range(len(Ikx)):
    res3=A*cmath.exp(Ikx[i])+B*cmath.exp(-Ikx[i])
    phi_x.append(res3)

time_start=0

E_time_start=cmath.exp(-I*E*time_start/h_bar)

phi_xt_start=[]
for i in range(len(phi_x)):
    res4=E_time_start*phi_x[i]
    phi_xt_start.append(res4)

phi_xt_start_real=[]
for i in range(len(phi_xt_start)):
    res5=phi_xt_start[i].real
    phi_xt_start_real.append(res5)

phi_xt_start_imag=[]
for i in range(len(phi_xt_start)):
    res6=phi_xt_start[i].imag
    phi_xt_start_imag.append(res6)

phi_xt_start_total_prob=0
for i in range(len(phi_xt_start)):
    total1=phi_xt_start[i]*phi_xt_start[i].conjugate()
    total1=total1.real
    phi_xt_start_total_prob+=total1

phi_xt_start_prob=[]
test1=0
for i in range(len(phi_xt_start)):
    res7=(phi_xt_start[i]*phi_xt_start[i].conjugate()).real/phi_xt_start_total_prob
    phi_xt_start_prob.append(res7)
    test1+=res7

assert test1>=0.999 and test1<=1.001

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot(phi_xt_start_real[:100])
ax1.set_title("real part")
ax2 = fig.add_subplot(312)
ax2.plot(phi_xt_start_imag[:100])
ax2.set_title("image part")
ax3 = fig.add_subplot(313)
ax3.plot(phi_xt_start_prob[:100])
ax3.set_title("probabilities")
plt.show()


wave_init=np.array(phi_xt_start)

np.savetxt("wave_t0_free.txt", wave_init)
np.save('wave_t0_free.npy', wave_init)

time_end=1000

E_time_end=cmath.exp(-I*E*time_end/h_bar)

phi_xt_end=[]
for i in range(len(phi_x)):
    res8=E_time_end*phi_x[i]
    phi_xt_end.append(res8)

phi_xt_end_real = []
for i in range(len(phi_xt_end)):
    res9 = phi_xt_end[i].real
    phi_xt_end_real.append(res9)

phi_xt_end_imag = []
for i in range(len(phi_xt_end)):
    res10 = phi_xt_end[i].imag
    phi_xt_end_imag.append(res10)

phi_xt_end_total_prob = 0
for i in range(len(phi_xt_end)):
    total2 = phi_xt_end[i] * phi_xt_end[i].conjugate()
    total2 = total2.real
    phi_xt_end_total_prob += total2

phi_xt_end_prob = []
test2 = 0
for i in range(len(phi_xt_end)):
    res11 = (phi_xt_end[i] * phi_xt_end[i].conjugate()).real / phi_xt_end_total_prob
    phi_xt_end_prob.append(res11)
    test2 += res11

assert test2 >= 0.999 and test2 <= 1.001

fig2 = plt.figure()
ax1 = fig2.add_subplot(311)
ax1.plot(phi_xt_end_real[:100])
ax1.set_title("real part")
ax2 = fig2.add_subplot(312)
ax2.plot(phi_xt_end_imag[:100])
ax2.set_title("image part")
ax3 = fig2.add_subplot(313)
ax3.plot(phi_xt_end_prob[:100])
ax3.set_title("probabilities")
plt.show()

wave_end=np.array(phi_xt_end)

np.savetxt("wave_t1000_free.txt", wave_end)
np.save('wave_t1000_free.npy', wave_end)
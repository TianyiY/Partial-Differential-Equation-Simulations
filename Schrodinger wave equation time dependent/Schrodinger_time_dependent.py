import numpy as np
import math
import matplotlib.pyplot as plt

wave_t0=np.load('wave_init.npy')

m=9.1e-31    # particle mass - kg
planck_constant=6.6e-34    # J*s
h_bar=planck_constant/2/math.pi
V=0   #potential
I=complex(0, 1)
coef1=-h_bar/2./m/I
coef2=V*I/h_bar

time=5000
step_time=1
steps=int(time/step_time)

size=len(wave_t0)
wave_new=np.copy(wave_t0)

for t in range(steps):
    for i in range(1, size-1):
        temp=wave_t0[i]
        wave_new[i]=temp+step_time*(-coef1*(wave_t0[i-1]+wave_t0[i+1]-2.*temp)-coef2*temp)

    # switch for the next iteration
    wave_new, wave_t0=wave_t0, wave_new

wave_real=[]
for i in range(len(wave_t0)):
    res1=wave_t0[i].real
    wave_real.append(res1)

wave_imag=[]
for i in range(len(wave_t0)):
    res2=wave_t0[i].imag
    wave_imag.append(res2)

wave_prob=[]
total=0
for i in range(len(wave_t0)):
    res3=wave_real[i]**2+wave_imag[i]**2
    wave_prob.append(res3)
    total+=res3

test=0
for i in range(len(wave_prob)):
    wave_prob[i]=wave_prob[i]/total
    test+=wave_prob[i]

assert test>=0.999 and test<=1.001

fig=plt.figure()
ax1=fig.add_subplot(311)
ax1.plot(wave_real[:100])
ax1.set_title("real part")
ax2=fig.add_subplot(312)
ax2.plot(wave_imag[:100])
ax2.set_title("image part")
ax3=fig.add_subplot(313)
ax3.plot(wave_prob[:100])
ax3.set_title("probability")
plt.show()

np.savetxt("wave_step5000_box.txt", wave_t0)
np.save('wave_step5000_box.npy', wave_t0)
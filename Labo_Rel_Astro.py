# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 15:04:05 2021

@author: Iván
"""

import matplotlib.pyplot as plt

import numpy as np

import urllib.request

import pandas as pd

from scipy.signal import argrelextrema


#Tarea 1


targetURL = "https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-H.txt"

i = 0

t = []

h = []

for line in urllib.request.urlopen(targetURL):
    
    if i == 0:
        
        column_names = line
        
        i+=1
        
    else:
        
        t.append(float(line.split()[0]))
        
        h.append(float(line.split()[1]))
        
output_df = pd.DataFrame({column_names[:17]: t, column_names[27:-1]: h})

output_df.to_csv('H1f.dat')

        
targetURL = "https://www.gw-openscience.org/GW150914data/P150914/fig1-waveform-H.txt"

i = 0

t_fmod = []

h_fmod = []

for line in urllib.request.urlopen(targetURL):
    
    if i == 0:
        
        column_names = line
        
        i+=1
        
    else:
        
        t_fmod.append(float(line.split()[0]))
        
        h_fmod.append(float(line.split()[1]))
        
output_df = pd.DataFrame({column_names[:17]: t_fmod, column_names[27:-1]: h_fmod})

output_df.to_csv('H1fmod.dat')


targetURL = "https://www.gw-openscience.org/GW150914data/P150914/fig2-unfiltered-waveform-H.txt"

i = 0

t_tmod = []

h_tmod = []

for line in urllib.request.urlopen(targetURL):
    
    if i == 0:
        
        column_names = line
        
        i+=1
        
    else:
        
        t_tmod.append(float(line.split()[0]))
        
        h_tmod.append(float(line.split()[1])*1e-21)
        
output_df = pd.DataFrame({column_names[:17]: t_tmod, column_names[27:-1]: h_tmod})

output_df.to_csv('H1tmod.dat')


#Tarea 2


plt.figure()
        
plt.plot(t, h, color='red', label='Filtered data')
plt.plot(t_fmod, h_fmod, color='blue', label='GR filtered')
plt.xlim(0.36, 0.44)
plt.title('Comparison between signals')
plt.xlabel('Time [s]')
plt.ylabel(r'Strain [$10^{-21}$]')
plt.legend()
plt.grid(True)

plt.savefig('tarea2_IvánVillegas.pdf')


#Tarea 3


f = []

t_3 = []

for i in range(1, len(h_tmod)-1):
    
    delta_t = t_tmod[i+1] - t_tmod[i]
    
    f.append(np.sqrt(-(h_tmod[i+1]-2*h_tmod[i]+h_tmod[i-1])/(h_tmod[i]*delta_t**2))/(2*np.pi))
    
    t_3.append(t_tmod[i])

plt.figure()    

plt.plot(t_3, f, color='black', label='Frequency')
plt.ylim(0, 500)
plt.title('Frequency squared against time')
plt.xlabel('Time [s]')
plt.ylabel(r'$\nu$ [Hz]')
plt.legend()
plt.grid(True)

plt.savefig('tarea3-1_IvánVillegas.pdf')


H_ = np.array(h_tmod)

max_H = argrelextrema(H_, np.greater)
min_H = argrelextrema(H_, np.less)

print('')

print(max_H[0])

print(min_H[0])

print('')

t_3_2 = []

f_2 = []

h_max = []

for i in range(1, len(h_tmod)-1):
        
    if i <= len(f)-1:

        if (h_tmod[i+1] < h_tmod[i] and h_tmod[i-1] < h_tmod[i]) or (h_tmod[i+1] > h_tmod[i] and h_tmod[i-1] > h_tmod[i]):
            
            t_3_2.append(t_tmod[i])
            
            f_2.append(f[i])
            
            h_max.append(h_tmod[i])
            
plt.figure()

plt.plot(t_3_2, f_2, color='green', label='Frequency')
#plt.xlim(0.25, 0.4)
plt.ylim(0, 250)
plt.title('Approximation of frequency')
plt.xlabel('Time [s]')
plt.ylabel(r'$\nu$ [Hz]')
plt.legend()
plt.grid(True)

plt.savefig('tarea3-2_IvánVillegas.pdf')
            


#Tarea 4


f_ = []

t_4 = []

for i in range(1, len(f_2)-1):
    
    delta_t_4 = t_3_2[i+1] - t_3_2[i-1]
    
    f_.append((f_2[i+1]-f_2[i-1])/(delta_t_4))
    
    t_4.append(t_3_2[i])
    

plt.figure()

plt.plot(t_4, f_, color='black')
plt.title(r'$\frac{d\nu}{dt}$ against time')
plt.xlabel('Time [s]')
plt.ylabel(r'$\dot{\nu}$ [Hz$^2$]')
plt.grid(True)


M_c= []

c = 3e8

G = 6.673848e-11

M_S = 2e30

suma_xy = 0

suma_x = 0

suma_y = 0

suma_xx = 0

for i in range(1, len(f_)+1):
    
    M_c.append(((c**3/G)*((5/96)*np.pi**(-8/3)*f_2[i]**(-11/3)*f_[i-1])**(3/5))/M_S)
    
    if t_4[i-1] <= 0.4:
    
        suma_xy = suma_xy + (t_4[i-1]*M_c[i-1])
        
        suma_x = suma_x + t_4[i-1]
        
        suma_y = suma_y + M_c[i-1]
        
        suma_xx = suma_xx + t_4[i-1]**2

a = (len(f_)*suma_xy-suma_x*suma_y)/(len(f_)*suma_xx-suma_x**2)

b = (suma_y-a*suma_x)/len(f_)

ajuste = a*0.2+b

print('')

print(f'M_c(0.2)=a·0.2+b={ajuste}')

print('')

M_C = []

for _ in range(0, len(f_)):
    
    M_C.append(ajuste)

plt.figure()

plt.plot(t_4, M_c, color='red', label='Chirp mass')
plt.plot(t_4, M_C, color='Blue', label='Estimated chirp mass')
plt.title(r'$M_c$ against time')
plt.xlabel('Time [s]')
plt.ylabel(r'$M_c$ [$M_\odot$]')
plt.legend()
plt.grid(True)

plt.savefig('tarea4_IvánVillegas.pdf')         


#Tarea 5


m_2 = np.linspace(0.01, 150, 50)

m_1 = []

m = []

for i in range(0, len(m_2)):
    
    m_1.append(np.sqrt(ajuste**5/m_2[i]**3))
    
    m.append(2**(1/5)*ajuste)

m_tot = 2*(2**(1/5)*ajuste)


plt.figure()

plt.plot(m_2, m_1, color='red', label=r'$m_1>>m_2$')
plt.plot(m_2, m, color='blue', label=r'$m_1\approx m_2$')
plt.xlim(0.01, 150)
plt.ylim(0.01, 150)
plt.title(r'Mass coomparisson')
plt.xlabel(r'$m_2$ [$M_\odot$]')
plt.ylabel(r'$m_1$ [$M_\odot$]')
plt.legend()
plt.grid(True)

plt.savefig('tarea5_IvánVillegas.pdf')        
 

print('')

print(f'Masa total = {m_tot}')

print('')


#Tarea 6


D_L = []

suma_xy = 0

suma_x = 0

suma_y = 0

suma_xx = 0

for i in range(1, len(f_)+1):
    
    D_L.append((10*c*f_[i-1]/(96*np.pi**2*abs(h_max[i])*f_2[i]**3))*3.24e-23)
    
    if t_4[i-1] <= 0.4:
    
        suma_xy = suma_xy + (t_4[i-1]*D_L[i-1])
        
        suma_x = suma_x + t_4[i-1]
        
        suma_y = suma_y + D_L[i-1]
        
        suma_xx = suma_xx + t_4[i-1]**2

a = (len(f_)*suma_xy-suma_x*suma_y)/(len(f_)*suma_xx-suma_x**2)

b = (suma_y-a*suma_x)/len(f_)

ajuste = a*0.2+b

print('')

print(f'D_L(0.2)=a·0.2+b={ajuste}')

print('')

D_l = []

for _ in range(0, len(f_)):
    
    D_l.append(ajuste)

plt.figure()

plt.plot(t_4[:-7], D_L[:-7], color='green', label='Light distance')
plt.plot(t_4, D_l, color='blue', label='Estimated light distance')
plt.title(r'$D_L$ against time')
plt.xlabel('Time [s]')
plt.ylabel(r'$D_L$ [Mpc]')
plt.legend()
plt.grid(True)

plt.savefig('tarea6_IvánVillegas.pdf')
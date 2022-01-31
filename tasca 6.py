# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:54:36 2021

@author: guill
"""

import numpy as np
import matplotlib.pyplot as plt



Reference_wavelength = 0.530  # micròmetres

def Index_TiO2(x):
    return (5.913+0.2441/(x**2-0.0803))**.5

def Index_MgF2(x):
    return (1+0.48755108/(1-(0.04338408/x)**2)+0.39875031/(1-(0.09461442/x)**2)+2.3120353/(1-(23.793604/x)**2))**.5

def Index_BK7(x): 
    return (1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5

#%%

plt.figure(figsize=(8, 4))

plt.subplot(1,3,1, title = 'TiO2')
wl = np.linspace(0.43, 1.53, 1000)
plt.plot(wl, Index_TiO2(wl), color = 'blue')
plt.xlabel("Longitud d'ona (micròmetres)")
plt.ylabel("Índex de Refracció (adimensional)")

plt.subplot(1,3,2, title = 'MgF2')
wl = np.linspace(0.2, 7, 7000)
plt.plot(wl, Index_MgF2(wl), color = 'green')
plt.xlabel("Longitud d'ona (micròmetres)")

plt.subplot(1,3,3, title = 'BK7')
wl = np.linspace(0.3, 2.5, 2000)
plt.plot(wl, Index_BK7(wl), color = 'red')
plt.xlabel("Longitud d'ona (micròmetres)")

plt.tight_layout()
plt.show()


#%%

def Indexs(x):
    return Index_TiO2(x), Index_MgF2(x), Index_BK7(x)  # Funcions 

def Matriu_M(wl, Medi): # Medi de la capa = (TiO2: 0, MgF2: 1, BK7: 2)
    Medi = int(Medi)    # Forçem a enter ja que ha de servir de índex per a la funció Indexs(x)
    gruix = Reference_wavelength / 4 / Indexs(Reference_wavelength)[Medi] # Gruix de la capa
    Fase = 2 * np.pi / wl * Indexs(wl)[Medi] * gruix    # Per comoditat: es correspon al que hi ha en els sinus i cosinus de la matriu M
    
    M = np.array([[ np.cos(Fase)                          ,  1j / Indexs(wl)[Medi] * np.sin(Fase) ],
                  [ 1j * Indexs(wl)[Medi] * np.sin(Fase)  ,  np.cos(Fase)                         ]])
    
    return M  # Matriu de la capa

#%%
# Sistema format per cinc dobles capes, de TiO2 i MgF2 de gruix lambda quarts, i substrat de vidre BK7
all_wl = np.linspace(0.43, 1, 2000)
R = np.array([])
T = np.array([])

for n in all_wl:
    M1 = Matriu_M(n, 0)  # Matriu pel TiO2
    M2 = Matriu_M(n, 1)  # Matriu pel MgF2
    
    MT = np.linalg.matrix_power(M1 @ M2 , 8)
    
    BC = MT @ np.array([[1],[Index_BK7(n)]])
    
    B = BC[0]
    C = BC[1]
    
    Ri = np.abs((B-C)/(B+C)) ** 2
    R = np.append(R, Ri)
    Ti = 1-Ri
    T = np.append(T, Ti)

plt.figure(figsize=(7, 4))
plt.plot(all_wl, R, color = 'orange', label = 'reflectància')
plt.plot(all_wl, T, color = 'cornflowerblue', label = 'transmitància')
plt.xlabel("Longitud d'ona (micròmetres)")
plt.ylabel("Reflectància / Transmitància (adimensional)")
plt.legend()
plt.show()

# Càlcul dels gruixos de les capes (micròmetres)
Gruix_capa_TiO2 = Reference_wavelength / 4 / Indexs(Reference_wavelength)[0]
Gruix_capa_MgF2 = Reference_wavelength / 4 / Indexs(Reference_wavelength)[1]

print('El gruix de la capa de diòxid de Titani és de', round(Gruix_capa_TiO2, 5), 'micròmetres. I el de la capa de Fluorur de Magnesi, de', round(Gruix_capa_MgF2, 5), '''micròmetres. Es pot comprovar que si alterem l'ordre de les capes, s'altera el resultat.''')





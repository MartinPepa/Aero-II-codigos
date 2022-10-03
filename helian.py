# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:09:29 2022

@author: Luis Martín Paredes Ross
"""
#%% Importacion de librerias

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import csv
from scipy import integrate

#%% Funciones

def lectura_helice(archivo_helice):
    with open(archivo_helice, 'rt', encoding = 'utf8') as f:
        rows        = csv.reader(f)
        fila_1      = next(rows)
        fila_2      = next(rows)
        titulo      = str(fila_2[0])
        Ndim1       = int(fila_2[1])
        encabezado  = next(rows)
        xD          = []
        csR         = []
        beta        = []
        for i, row in enumerate(rows):
            xD.append(float(row[0]))
            csR.append(float(row[1]))
            beta.append(float(row[2]))
        return titulo, Ndim1, xD, csR, beta
    
def lectura_curva_polar(archivo_polar):
    with open(archivo_polar, 'rt', encoding = 'utf8') as f:
        rows        = csv.reader(f)
        fila_1      = next(rows)
        fila_2      = next(rows)
        titulo      = str(fila_2[0])
        Ndim2       = int(fila_2[1])
        encabezado  = next(rows)
        alfa        = []
        cl          = []
        cd          = []
        for i, row in enumerate(rows):
            alfa.append(float(row[0]))
            cl.append(float(row[1]))
            cd.append(float(row[2]))
        return titulo, Ndim2, alfa, cl, cd

def lectura_simulacion(parametros):
    with open(parametros, 'rt', encoding = 'utf8') as f:
        rows = csv.reader(f)
        datos = []
        for i, row in enumerate(rows):
            datos.append(row[1])
        titulo  = datos[0]
        AJ0     = float(datos[1])
        AJf     = float(datos[2])
        NJ      = float(datos[3])
        rpm     = float(datos[4])
        h       = float(datos[5])
        B       = float(datos[6])
        D       = float(datos[7])
        X0      = float(datos[8])
        DeltaX  = float(datos[9])
        NX      = int(datos[10])
        Beta75i = float(datos[11])
        Beta75f = float(datos[12])
        NBeta   = int(datos[13])
        AlfaMin = float(datos[14])
        AlfaPer = float(datos[15])
        AlfaMax = float(datos[16])
        Impres  = float(datos[17])
        return titulo, AJ0, AJf, NJ, rpm, h, B, D, X0, DeltaX, NX, Beta75i,\
            Beta75f, NBeta, AlfaMin, AlfaPer, AlfaMax, Impres

#%% Lectura y bajada de datos

titulo, AJ0, AJf, NJ, rpm, h, B, D, X0, DeltaX, NX, Beta75i, Beta75f, NBeta,\
    AlfaMin, AlfaPer, AlfaMax, Impres = lectura_simulacion('parametros.csv')
titulo1, NDim1, xD, csR, beta = lectura_helice('helice.csv')
titulo2, NDim2, alfa, cl, cd = lectura_curva_polar('polar.csv')

#%% Splines y gráficos

cs1 = CubicSpline(xD, csR, bc_type = 'natural')
cs2 = CubicSpline(xD, beta, bc_type = 'natural')

# xs = np.arange(0.2,1.0,0.05)

fig, ax = plt.subplots(dpi=400)

ax.plot(xD, csR    , 'o', label = 'cuerda/radio')
ax.plot(xD, beta   , 'x', label = 'beta [rad]')
ax.plot(xD, cs1(xD), label = 'spline csr')
ax.plot(xD, cs2(xD), label = 'spline beta')

ax.set_xlim(0.0,1.0)
ax.legend(loc='upper right',ncol=1)

plt.show()

cs3 = CubicSpline(alfa, cl, bc_type = 'natural')
cs4 = CubicSpline(alfa, cd, bc_type = 'natural')

fig, bx = plt.subplots(dpi=400)

bx.plot(alfa, cl, 'o', label = 'cl-alfa')
bx.plot(alfa, cd, 'x', label = 'cd-alfa')
bx.plot(alfa, cs3(alfa),label = 'spline cl')
bx.plot(alfa, cs4(alfa),label = 'spline cd')

bx.legend(loc='upper left',ncol=1)
plt.show()

fig, cx = plt.subplots(dpi=400)
cx.plot(cd, cl, 'x', label = 'polar')
plt.show()

#%% Calculos preliminares

if AlfaPer > AlfaMax:
    AlfaPer = AlfaMax

R       = D/2
Omega   = np.pi*rpm/30
v_t     = Omega*R
rho     = 0.1249*(1-2.25577/100000*h)**4.25575
Temp    = 15.0 - 0.065*h
mu      = 1.7894e-5 * ((Temp/15.0)**0.75)
nu      = mu/rho
AlfMax  = AlfaMax*np.pi/180

Beta75    = float(cs2(0.75))

AJtov   = D*rpm/60.0
vei     = AJ0*AJtov
DelVe   = (AJf-AJ0)/NJ*AJtov
nve     = NJ+1
DelBe0  = Beta75i - Beta75/np.pi*180.0
DelBei  = (Beta75f-Beta75i)/NBeta
NDelBe  = NBeta+1

C3 = []
cX3 = 0.0
for i in range(NDim1+1):
    C = csR[i]*R
    cX3 = C*xD[i]**3
    C3.append(cX3)
    C = 0
    cX3 = 0

Integral = integrate.simpson(xD,C3)
Fact_Act = 1e5*Integral/(16*D)

alfa_rad = []
for i in range(NDim2):
    alfa_rad.append(alfa[i]*np.pi/180)

X           = []
cuerda_x    = []
beta_x      = []
sigma_x     = []
for i in range(NX):
    X.append(X0)
    if X0 >= 1.0:
        X[i] = 1.0
        NX = i
        break
    cuerda_x.append(cs1(X0)*R)
    beta_x.append(cs2(X0))
    sigma_x.append(B*cuerda_x[i]/(np.pi*R))
    X0 += DeltaX
#%% Titulo A





#%% 





#%% Titulo B





#%% Rutina de calculo iterativo





# Resultados





#%% 

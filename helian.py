# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:09:29 2022

@author: Luis Martín Paredes Ross
"""
#%% Import todo

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import csv

#%% Lectura y bajada de datos

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
        NX      = float(datos[10])
        Beta75i = float(datos[11])
        Beta75f = float(datos[12])
        NBeta   = float(datos[13])
        AlfaMin = float(datos[14])
        AlfaPer = float(datos[15])
        AlfaMax = float(datos[16])
        Impres  = float(datos[17])
        return titulo, AJ0, AJf, NJ, rpm, h, B, D, X0, DeltaX, NX, Beta75i,\
            Beta75f, NBeta, AlfaMin, AlfaPer, AlfaMax, Impres

titulo, AJ0, AJf, NJ, rpm, h, B, D, X0, DeltaX, NX, Beta75i, Beta75f, NBeta,\
    AlfaMin, AlfaPer, AlfaMax, Impres = lectura_simulacion('parametros.csv')
titulo1, NDim1, xD, csR, beta = lectura_helice('helice.csv')
titulo2, NDim2, alfa, cl, cd = lectura_curva_polar('polar.csv')

#%%

if AlfaPer > AlfaMax:
    AlfaPer = AlfaMax
    R       = D/2
    Omega   = pi*rpm/30
    v_t     = Omega*R
    rho     = 0.1249*(1-2.25577/100000*h)^(4.25575)
    Temp    = 15.0 - 0.065*h
    mu      = 1.7894e-5 * ((Temp/15.0)^0.75)
    nu      = mu/rho
    AlfMax  = AlfaMax*pi/180
    splin(ndim1,xd,beta,z,a)
    be75    = spfun(ndim1,xd,beta,z,0.75)
    
    ajtov   = D*rpm/60.0
    vei     = aji*ajtov
    nve     = nj+1
    delbe0  = bta75i - be75/pi*180.0
    delbei  = (bta75f-bta75i)/nbta
    ndelbe  = nbta+1
#%%







#%%






#%%







# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:35:30 2022

@author: Luis Martín Paredes Ross
"""

#%% Importación de librerias

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
import csv
from scipy import integrate

#%% Leida de datos

def lectura_csv(archivo_csv):
    with open(archivo_csv, 'rt', encoding = 'utf8') as f:
        rows = csv.reader(f)
        encabezado = next(rows)
        enc_1 = encabezado[0]
        enc_2 = encabezado[1]
        enc_3 = encabezado[2]
        xR = []
        csR = []
        beta = []
        Dim = 0
        for i, row in enumerate(rows):
            try: #type(float(row[0])) is float:
                xR.append(float(row[0]))
                csR.append(float(row[1]))
                beta.append(float(row[2]))
                Dim += 1
            except:        
                #row = next(rows)
                titulo_1 = row[0]
                parametros = []
                for i, row in enumerate(rows):
                    parametros.append(row[1])
                Vasc        = float(parametros[0])
                Omega       = float(parametros[1])
                R           = float(parametros[2])
                x0          = float(parametros[3])
                B           = int(  parametros[4])
                NX          = int(  parametros[5])
                rho         = float(parametros[6])
                Va          = float(parametros[7])
                theta_0     = float(parametros[8])
                theta_f     = float(parametros[9])
                delta_theta = float(parametros[10])
                return enc_1, enc_2, enc_3, xR, csR, beta, Dim, titulo_1,\
                    Vasc, Omega, R, x0, B, NX, rho, Va, theta_0, theta_f, delta_theta

#%% Cálculos preliminares

enc_1, enc_2, enc_3, xR, csR, beta, Dim, titulo_1,\
    Vasc, Omega, R, x0, B, NX, rho, Va, theta_0, theta_f,\
        delta_theta = lectura_csv('Hover-input.csv')


#%% Cálculo de tabla


#%% Cálculos finales
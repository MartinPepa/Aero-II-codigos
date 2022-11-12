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
        next(rows)
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

enc_1, enc_2, enc_3, xR, csR, beta, Dim, titulo_1,Vasc, Omega,\
    R, x0, B, NX, rho, Va, theta_0, theta_f, delta_theta = lectura_csv('Hover-input.csv')


#%% Splines y gráficos

XMach   = [x/10 for x in range(3,10)]
AlfMay  = [float(x) for x in range(14,0,-2)]
alfa    = [-0.3490, -0.2792, -0.2443, -0.2094, -0.1745, -0.1396, -0.1047,
        -0.0698, -0.0349, 0.000, 0.0349, 0.0698, 0.1047, 0.1396, 0.1745,
        0.2094, 0.2443, 0.2792, 0.3141, 0.3490]

CLp3 = [-1.08, -1.46, -1.37, -1.25, -1.15, -0.92, -0.67, -0.42, -0.20, 0.05,
        0.30, 0.50, 0.73, 0.96, 1.21, 1.46, 1.50, 1.17, 1.12, 1.10]
CLp4 = [-1.17, -1.21, -1.34, -1.21, -1.10, -0.92, -0.67, -0.42, -0.20, 0.05,
        0.30, 0.50, 0.73, 0.96, 1.21, 1.44, 1.21, 1.12, 1.08, 1.06]
CLp5 = [-0.96, -1.13, -1.25, -1.15, -1.00, -0.85, -0.65, -0.42, -0.20, 0.05,
        0.30, 0.54, 0.77, 1.00, 1.21, 1.25, 1.21, 1.12, 1.08, 1.06]
CLp6 = [-1.17, -1.13, -1.12, -1.10, -1.08, -0.96, -0.75, -0.46, -0.21, 0.06,
        0.37, 0.64, 0.87, 1.00, 1.06, 1.08, 1.12, 1.13, 1.135, 1.14]
CLp7 = [-0.83, -0.86, -0.87, -0.86, -0.85, -0.79, -0.67, -0.54, -0.33, 0.08,
        0.42, 0.71, 0.81, 0.89, 1.00, 1.12, 1.21, 1.33, 1.42, 1.50]
CLp8 = [-0.62, -0.62, -0.62, -0.62, -0.62, -0.584, -0.58, -0.50, -0.37,
        0.09, 0.48, 0.54, 0.58, 0.60, 0.62, 0.67, 0.71, 0.75, 0.79, 0.83]
CLp9 = [-0.50, -0.50, -0.50, -0.50, -0.50, -0.50, -0.50, -0.31, -0.08,
        -0.05, 0.12, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19]

CDp3 = [0.070, 0.060, 0.055, 0.040, 0.038, 0.0353, 0.0176, 0.0106,
        0.0094, 0.008, 0.008, 0.008, 0.009, 0.010, 0.0128, 0.017,
        0.0235, 0.0553, 0.056, 0.057]
CDp4 = [0.071, 0.061, 0.058, 0.049, 0.045, 0.040, 0.033, 0.0117,
        0.008, 0.008, 0.008, 0.008, 0.0086, 0.010, 0.0164, 0.0211,
        0.034, 0.040, 0.045, 0.050]
CDp5 = [0.072, 0.068, 0.064, 0.060, 0.058, 0.0547, 0.0376, 0.0188,
        0.0081, 0.008, 0.008, 0.008, 0.0128, 0.0186, 0.026, 0.0376,
        0.040, 0.042, 0.045, 0.050]
CDp6 = [0.073, 0.070, 0.068, 0.065, 0.062, 0.060, 0.047, 0.025,
        0.0117, 0.008, 0.0082, 0.0128, 0.0294, 0.054, 0.056, 0.060,
        0.065, 0.070, 0.072, 0.075]
CDp7 = [0.074, 0.072, 0.070, 0.068, 0.064, 0.060, 0.054, 0.032,
        0.0182, 0.0086, 0.0141, 0.0443, 0.0506, 0.059, 0.062, 0.065,
        0.068, 0.071, 0.073, 0.076]
CDp8 = [0.075, 0.070, 0.065, 0.060, 0.055, 0.045, 0.0411, 0.0306,
        0.0188, 0.0164, 0.0390, 0.050, 0.055, 0.060, 0.063, 0.066,
        0.069, 0.072, 0.074, 0.077]
CDp9 = [0.076, 0.071, 0.066, 0.062, 0.056, 0.050, 0.045, 0.027,
        0.020, 0.020, 0.022, 0.0306, 0.056, 0.061, 0.064, 0.067,
        0.070, 0.073, 0.075, 0.078]

fig, ax = plt.subplots(dpi=400)

ax.plot(alfa, CLp3,'o', label = 'p3')
ax.plot(alfa, CLp4,'s', label = 'p4')
ax.plot(alfa, CLp5,'^', label = 'p5')
ax.plot(alfa, CLp6,'x', label = 'p6')
ax.plot(alfa, CLp7,'*', label = 'p7')
ax.plot(alfa, CLp8,'+', label = 'p8')
ax.plot(alfa, CLp9,'d', label = 'p9')
ax.legend(loc='upper left',ncol=1)
plt.show()


spline_cuerda = CubicSpline(xR, csR, bc_type = 'natural')
spline_beta = CubicSpline(xR, beta, bc_type = 'natural')
spline_alfamay = CubicSpline(XMach, AlfMay, bc_type = 'natural')

spline_clp3 = CubicSpline(alfa, CLp3, bc_type = 'natural')
spline_clp4 = CubicSpline(alfa, CLp4, bc_type = 'natural')
spline_clp5 = CubicSpline(alfa, CLp5, bc_type = 'natural')
spline_clp6 = CubicSpline(alfa, CLp6, bc_type = 'natural')
spline_clp7 = CubicSpline(alfa, CLp7, bc_type = 'natural')
spline_clp8 = CubicSpline(alfa, CLp8, bc_type = 'natural')
spline_clp9 = CubicSpline(alfa, CLp9, bc_type = 'natural')

spline_cdp3 = CubicSpline(alfa, CDp3, bc_type = 'natural')
spline_cdp4 = CubicSpline(alfa, CDp4, bc_type = 'natural')
spline_cdp5 = CubicSpline(alfa, CDp5, bc_type = 'natural')
spline_cdp6 = CubicSpline(alfa, CDp6, bc_type = 'natural')
spline_cdp7 = CubicSpline(alfa, CDp7, bc_type = 'natural')
spline_cdp8 = CubicSpline(alfa, CDp8, bc_type = 'natural')
spline_cdp9 = CubicSpline(alfa, CDp9, bc_type = 'natural')

#%% Cálculo de tabla

theta = theta_0

while theta <= theta_f:
    print('|  N  | X(N)  |  ZMach  |   AlfMax   |')
    tol = 0.002
    
    #DN          = NX + 0.000001
    Delta_X     = (1.0 - x0)/NX
    theta       = theta/57.3
    VT          = Omega*R
    lambda_g    = Vasc/VT
    alfa_i_0    = 0.017
    NX          = NX + 1
    
    alfa_a  = np.zeros(NX)
    alfa_i  = np.zeros(NX)
    arg     = np.zeros(NX)
    beta_a  = np.zeros(NX)
    c       = np.zeros(NX)
    Cd      = np.zeros(NX)
    Cl      = np.zeros(NX)
    dkt     = np.zeros(NX)
    dkp     = np.zeros(NX)
    F       = np.zeros(NX)
    phi     = np.zeros(NX)
    sigma   = np.zeros(NX)
    VeVT    = np.zeros(NX)
    X       = np.zeros(NX)
    wt_VT   = np.zeros(NX)
    wa      = np.zeros(NX)
    wa_VT   = np.zeros(NX)
    
    for i in range(NX):
        X[i] = x0
        
        if i == NX:
            X[i]        = 1.0
        c[i]        = spline_cuerda(X[i])
        sigma[i]    = B*c[i]/np.pi
        phi[i]      = np.arctan(lambda_g/X[i])
        beta1       = spline_beta(X[i])
        beta_a[i]   = beta1 + theta
        VeVT0       = np.sqrt(X[i]**2+lambda_g**2)
        k           = 0
        while k < 60:
            alfa_i[i]   = alfa_i_0
            VeVT[i]     = VeVT0
            ZMach       = VeVT[i]/(Va/VT)
            
            if ZMach > 0.3:
                AlfaMax = spline_alfamay(ZMach)
            else:
                AlfaMax = 14.0
            
            AlfaMax = AlfaMax/57.3
            arg[i] = phi[i] + alfa_i[i]
            alfa_a[i] = beta_a[i] -arg[i]
            if (alfa_a[i] - AlfaMax) > 0:
                arg[i] = beta_a[i] - AlfaMax
                alfa_i_0 = arg[i] - phi[i] #################### SOLUCIONAR ACA
                continue
            if ZMach <= 0.3:
                Cl[i] = spline_clp3(alfa_a[i])
            elif (ZMach > 0.3) and (ZMach <= 0.4):
                Cl3 = spline_clp3(alfa_a[i])
                Cl4 = spline_clp4(alfa_a[i])
                Cl[i] = (Cl4-Cl3)*10.0*(ZMach-0.3)+Cl3
            elif (ZMach > 0.4) and (ZMach <= 0.5):
                Cl4 = spline_clp4(alfa_a[i])
                Cl5 = spline_clp5(alfa_a[i])
                Cl[i] = (Cl5-Cl4)*10.0*(ZMach-0.4)+Cl4
            elif (ZMach > 0.5) and (ZMach <= 0.6):
                Cl5 = spline_clp5(alfa_a[i])
                Cl6 = spline_clp6(alfa_a[i])
                Cl[i] = (Cl6-Cl5)*10.0*(ZMach-0.5)+Cl5
            elif (ZMach > 0.6) and (ZMach <= 0.7):
                Cl6 = spline_clp6(alfa_a[i])
                Cl7 = spline_clp7(alfa_a[i])
                Cl[i] = (Cl7-Cl6)*10.0*(ZMach-0.6)+Cl6
            elif (ZMach > 0.7) and (ZMach <= 0.8):
                Cl7 = spline_clp7(alfa_a[i])
                Cl8 = spline_clp8(alfa_a[i])
                Cl[i] = (Cl8-Cl7)*10.0*(ZMach-0.7)+Cl7
            elif (ZMach > 0.8) and (ZMach <= 0.9):
                Cl8 = spline_clp8(alfa_a[i])
                Cl9 = spline_clp9(alfa_a[i])
                Cl[i] = (Cl9-Cl8)*10.0*(ZMach-0.8)+Cl8
            else:
                Cl[i] = spline_clp9(alfa_a[i])
            f = (1.0 - X[i])*B/2/(X[i]*np.tan(arg[i]))
            
            if (f - 50.0) >= 0:
                F[i] = 1.0
                if (X[i]-1.0) < 0:
                    wt_VT[i] = sigma[i]*VeVT[i]*Cl[i]/(8.0*X[i]*F[i])
                    if wt_VT[i] >= 0.0:
                        wa[i] = np.sqrt(lambda_g**2 + 4.0*wt_VT[i]*(X[i]-wt_VT[i]))
                    else:
                        wa[i] = np.sqrt(abs(lambda_g**2
                                             + 4.0*wt_VT[i]*(X[i]-wt_VT[i])))
                        print('Warning: Solución compleja.')
                elif (X[i]-1.0) == 0:
                     alfa_a[i] = 0.0
                else:
                     alfa_a[NX] =0.0
            elif (f - 50.0) < 0:
                G = np.exp(-f)
                if (X[i]-1.0) < 0:
                    F[i] = 0.63662*np.arccos(G)
                    wt_VT[i] = sigma[i]*VeVT[i]*Cl[i]/(8.0*X[i]*F[i])
                    if wt_VT[i] >= 0.0:
                        wa[i] = np.sqrt(lambda_g**2 + 4.0*wt_VT[i]*(X[i]-wt_VT[i]))
                    else:
                        wa[i] = np.sqrt(abs(lambda_g**2
                                             + 4.0*wt_VT[i]*(X[i]-wt_VT[i])))
                        print('Warning: Solución compleja.')
                elif (X[i]-1.0) == 0:
                     alfa_a[i] = 0.0
                else:
                     alfa_a[NX-1] =0.0
            wa_VT[i] = (-lambda_g + wa[i])/2
            VEVT_0 = np.sqrt((X[i]-wt_VT[i])**2.0 + (lambda_g+wa_VT[i])**2)
            alf = (lambda_g + wa_VT[i]) / (X[i]-wt_VT[i])
            alfii = np.arctan(alf) - phi[i]
            Error = abs(alfii - alfa_i[i])
            
            if Error <= tol:
                break
            
            k += 1
            if k == 60:
                print(f'ALFAA(I) >= ALFMAX EN I= {i:>3d}. También número de itera\
                      ciones mayor a 60')
            alfa_i_0 = (alfii + alfa_i[i])/2.0
        
        if ZMach <= 0.3:
            Cd[i] = spline_cdp3(alfa_a[i])
        elif (ZMach > 0.3) and (ZMach <= 0.4):
            Cd3 = spline_cdp3(alfa_a[i])
            Cd4 = spline_cdp4(alfa_a[i])
            Cd[i] = (Cd4 - Cd3)*10.0*(ZMach-0.3)+Cd3
        elif (ZMach > 0.4) and (ZMach <= 0.5):
            Cd4 = spline_cdp4(alfa_a[i])
            Cd5 = spline_cdp5(alfa_a[i])
            Cd[i] = (Cd5 - Cd4)*10.0*(ZMach-0.4)+Cd4
        elif (ZMach > 0.5) and (ZMach <= 0.6):
            Cd5 = spline_cdp5(alfa_a[i])
            Cd6 = spline_cdp6(alfa_a[i])
            Cd[i] = (Cd6 - Cd5)*10.0*(ZMach-0.5)+Cd5
        elif (ZMach > 0.6) and (ZMach <= 0.7):
            Cd6 = spline_cdp6(alfa_a[i])
            Cd7 = spline_cdp7(alfa_a[i])
            Cd[i] = (Cd7 - Cd6)*10.0*(ZMach-0.6)+Cd6
        elif (ZMach > 0.7) and (ZMach <= 0.8):
            Cd7 = spline_cdp7(alfa_a[i])
            Cd8 = spline_cdp8(alfa_a[i])
            Cd[i] = (Cd8 - Cd7)*10.0*(ZMach-0.7)+Cd7
        elif (ZMach > 0.8) and (ZMach <= 0.9):
            Cd8 = spline_cdp8(alfa_a[i])
            Cd9 = spline_cdp9(alfa_a[i])
            Cd[i] = (Cd9 - Cd8)*10.0*(ZMach-0.8)+Cd8
        else:
            Cd[i] = spline_cdp9(alfa_a[i])
        
        if X[i] == 1.0:
            alfa_a[NX-1] = 0.0
            Cl[NX-1] = 0.0
            continue
        kt1     = sigma[i]*(VeVT[i]**2)
        kt2     = Cl[i]*np.cos(arg[i])-Cd[i]*np.sin(arg[i])
        dkt[i]  = kt1*kt2
        kp1     = kt1*X[i]
        kp2     = Cl[i]*np.sin(arg[i])+Cd[i]*np.cos(arg[i])
        dkp[i]  = kp1*kp2
        x0      = X[i] + Delta_X
        print(f'{i:>8d} {X[i]:>8.2f} {ZMach:>8.2f} {AlfaMax:>8.2f}')
    print(f'{NX:>8d} {X[NX-1]:>8.2f} {ZMach:>8.2f} {AlfaMax:>8.2f}')
    Beta1 = spline_beta(X[NX-1])
    beta_a[NX-1] = Beta1 + theta
    VR_VT = np.sqrt(X[NX-1]**2 + lambda_g**2)
    phi[NX-1] = np.arctan(lambda_g/X[NX-1])
    alfa_i[NX-1] = beta_a[NX-1] - phi[NX-1]
    wt_VT[NX-1] = VR_VT*np.sin(alfa_i[NX-1])*np.sin(beta_a[NX-1])
    wa_VT[NX-1] = VR_VT*np.sin(alfa_i[NX-1])*np.cos(beta_a[NX-1])
    VeVT[NX-1] = np.sqrt((X[NX-1]-wt_VT[NX-1])**2 + (lambda_g+wa_VT[NX-1])**2)
    dkt[NX-1] = -sigma[NX-1]*(VeVT[NX-1]**2)*Cd[NX-1]*np.sin(beta_a[NX-1])
    dkp[NX-1] = sigma[NX-1]*(VeVT[NX-1]**2)*Cd[NX-1]*np.cos(beta_a[NX-1])
    
    
    #%% Cálculos finales
    
    #A1 = X[0]
    #B1 = X[NX]
    
    tkt = integrate.simpson(dkt, X, even = 'first')
    pkp = integrate.simpson(dkp, X, even = 'first')
    
    T = np.pi/2*rho*(R**2)*(VT**2)*tkt
    P = np.pi/2*rho*(R**2)*(VT**3)*pkp/76.05
    ETA = T*Vasc/(P*76.05)
    n = Omega/(6.2832)
    AJ = Vasc/(n*2*R)
    
    CT = T/(rho*(VT**2)*np.pi*R**2)
    CP = P*76.05/(rho*(VT**3)*np.pi*(R**2))
    
    eta1 = AJ*CT/CP
    Be75 = spline_beta(0.75)
    Beta75 = Be75 + theta
    
    
    #%% Impresión de resultados
    
    print('\n'+'----------'*6)
    print('\nResultados del análisis\n')
    print(f'Vasc = {Vasc:>6.2f} [m/s]   |   Omega = {Omega:>6.2f} [1/s]')
    print(f' R   = {R:>6.2f} [m]   |   B = {B:>6d}')
    print(f' Beta75 = {Beta75:>6.2f} [rad] (Be75 + theta)')
    print(f' Theta = {theta:>6.2f} [rad]\n')
    
    print(f'\nCT = {CT:>5.4f}   |   CP = {CP:>5.4f}')
    print(f'J = {AJ:5.4f}   |   Efic = {eta1:>5.4f}\n')
    print(f'| Tracción [kg] = {T:>8.2f} |\n| Potencia [HP] = {P:>8.2f} |\n| Eficiencia [-] = {ETA:>5.4} |')
          
    if ETA < 0:
        print('Eficiencia negativa.')
        print('##########'*8)
    elif ETA == 0:
        PID = (T**1.5)/(np.sqrt(2*np.pi*rho*R**2))
        Mer = PID/(P*76.05)
        print(f' M [-] = {Mer:>5.4f}')
    else:
        print('No hay Factor de Mérito si asciende.')
    
    print('\n   Alfa   |   CL   |   CD   |   Vind   |')
    for i in range(NX):
        print(f'{alfa_a[i]:>6.4f} {Cl[i]:>5.4f} {Cd[i]:>5.4f} {wa_VT[i]:>5.4f}')
    print('##########'*8)
    
    theta += delta_theta/57.3

#%%
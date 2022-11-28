# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:09:29 2022

@author: Luis Martín Paredes Ross
"""
#%% Importacion de librerias

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
import csv
from scipy import integrate
plt.close(fig = 'all')

#%% Funciones

def lectura_helice(archivo_helice):
    with open(archivo_helice, 'rt', encoding = 'utf8') as f:
        rows        = csv.reader(f)
        fila_1      = next(rows)
        fila_2      = next(rows)
        titulo      = str(fila_2[0])
        encabezado  = next(rows)
        xD          = []
        csR         = []
        beta        = []
        for i, row in enumerate(rows):
            xD.append(float(row[0]))
            csR.append(float(row[1]))
            beta.append(float(row[2]))
            Ndim1   = i+1
        return titulo, Ndim1, xD, csR, beta
    
def lectura_curva_polar(archivo_polar):
    with open(archivo_polar, 'rt', encoding = 'utf8') as f:
        rows        = csv.reader(f)
        fila_1      = next(rows)
        fila_2      = next(rows)
        titulo      = str(fila_2[0])
        encabezado  = next(rows)
        alfa        = []
        cl          = []
        cd          = []
        for i, row in enumerate(rows):
            alfa.append(float(row[0]))
            cl.append(float(row[1]))
            cd.append(float(row[2]))
            NDim2   = i+1
        return titulo, NDim2, alfa, cl, cd

def lectura_simulacion(parametros):
    with open(parametros, 'rt', encoding = 'utf8') as f:
        rows = csv.reader(f)
        datos = []
        for i, row in enumerate(rows):
            datos.append(row[1])
        titulo  = datos[0]
        AJ0     = float(datos[1])
        AJf     = float(datos[2])
        NJ      = int(  datos[3])
        rpm     = float(datos[4])
        h       = float(datos[5])
        B       = int(  datos[6])
        D       = float(datos[7])
        X0      = float(datos[8])
        DeltaX  = float(datos[9])
        NX      = int(  datos[10])
        Beta75i = float(datos[11])
        Beta75f = float(datos[12])
        NBeta   = int(  datos[13])
        AlfaMin = float(datos[14])
        AlfaPer = float(datos[15])
        AlfaMax = float(datos[16])
        Impres  = int(  datos[17])
        return titulo, AJ0, AJf, NJ, rpm, h, B, D, X0, DeltaX, NX, Beta75i,\
            Beta75f, NBeta, AlfaMin, AlfaPer, AlfaMax, Impres

#%% Lectura y bajada de datos

titulo, AJ0, AJf, NJ, rpm, h, B, D, X0, DeltaX,\
    NX, Beta75i, Beta75f, NBeta, AlfaMin, AlfaPer,\
    AlfaMax, Impres = lectura_simulacion('parametros.csv')
    
titulo1, NDim1, xD, csR, beta = lectura_helice('helice.csv')

titulo2, NDim2, alfa, cl, cd = lectura_curva_polar('polar.csv')

#%% Splines y gráficos

cs1 = CubicSpline(xD,  csR, bc_type = 'natural')
cs2 = CubicSpline(xD, beta, bc_type = 'natural')
cs3 = CubicSpline(alfa, cl, bc_type = 'natural')
cs4 = CubicSpline(alfa, cd, bc_type = 'natural')

# xs = np.arange(0.2,1.0,0.05)

fig, ax = plt.subplots(dpi=400)

ax.plot(xD, csR    , 'o', label = 'cuerda/radio')
ax.plot(xD, beta   , 'x', label = 'beta [rad]')
ax.plot(xD, cs1(xD), label = 'spline csr')
ax.plot(xD, cs2(xD), label = 'spline beta')
ax.set_xlim(0.0,1.0)
ax.set_ylim(0.0,0.5)
ax.legend(loc='upper right',ncol=1)
ax.set_xlabel('r/R [adim]')
ax.set_ylabel('beta [rad] | c/R [adim]')
plt.grid( visible = True, which = 'both', axis = 'both')
plt.show()

fig, bx = plt.subplots(dpi=400)

bx.plot(alfa, cl, 'o', label = 'cl-alfa')
bx.plot(alfa, cd, 'x', label = 'cd-alfa')
bx.plot(alfa, cs3(alfa),label = 'spline cl')
bx.plot(alfa, cs4(alfa),label = 'spline cd')
bx.set_xlim(-10.0,30.0)
bx.set_ylim(-0.2,1.4)
bx.set_xlabel('\u03B1 [º]')
bx.set_ylabel('cl | cd')
plt.grid(visible = True, which = 'both', axis = 'both')
bx.legend(loc='upper left',ncol=1)
plt.show()

fig, cx = plt.subplots(dpi=400)

cx.plot(cd, cl, 'x', label = 'polar')
cx.set_xlim(0.0, 0.35)
cx.set_ylim(0.0, 1.4)
cx.set_xlabel('Cd')
cx.set_ylabel('Cl')
plt.grid(visible = True, which = 'both', axis = 'both')
plt.show()

#%% Calculos preliminares

if AlfaPer > AlfaMax:
    AlfaPer = AlfaMax
tol     = 0.0001
R       = D/2
Omega   = np.pi*rpm/30
VT      = Omega*R
rho     = 0.1249*(1-2.25577/100000*h)**4.25575 # Unidad tecnica
Temp    = 15.0 - 0.065*h
mu      = 1.7894e-5 * ((Temp/15.0)**0.75)
nu      = mu/rho
AlfMax  = AlfaMax*np.pi/180

Be75    = float(cs2(0.75))

AJtoV   = D*rpm/60.0 # Factor para convertir AJ a Velocidad
VEi     = AJ0*AJtoV
DelVe   = (AJf-AJ0)/NJ*AJtoV
NVe     = NJ+1
DelBe0  = Beta75i - Be75*180/np.pi
DelBei  = (Beta75f-Beta75i)/NBeta
NDelBe  = NBeta+1

C3 = []
cX3 = 0.0
for i in range(NDim1):
    C   = csR[i]*R
    cX3 = C*xD[i]**3
    C3.append(cX3)
    C   = 0
    cX3 = 0

# xd = [x/100 for x in range(15,105,5)]

# spline_1 = CubicSpline(xD, C3, bc_type = 'natural')
# valor = spline_1(0.15)
# Integral = spline_1.integrate(0.15, 1.0, extrapolate='True')
Integral = integrate.simpson(C3,xD,even='first')
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

print('__________'*6,'\n')
print('Resultados del programa Helian para análisis de hélices\n')
print('----------'*5)
print(f'Título:\n\n{titulo:<}')
print('----------'*5)
print(f'Hélice:\n\n{titulo1:<}\nCant. de puntos: {NDim1:<d}')
print(f'Factor de actividad AF: {Fact_Act:<.2f}')
print('----------'*5)
print(f'Perfil:\n\n{titulo2:<}\nCant. de puntos: {NDim2:<d}')
print('----------'*5)
print('Parámetros de la corrida:\n')
print(f'AJ0: {AJ0:>10.2f} | AJf: {AJf:>10.2f} | NJ: {NJ:>9d} |')
print(f'RPM: {rpm:>10.2f} | h [m]: {h:>8.2f} |')
print(f'Nº de palas: {B:>2d} | D [m]: {D:>8.2f} |')
print(f'X0: {X[0]:>11.2f} | DeltaX: {DeltaX:>7.2f} | NX: {NX:>9d} |')
print(f'Beta75(i): {Beta75i:>4.1f} | Beta75(f): {Beta75f:>4.1f} |\
 NBeta: {NBeta:>6d} |')
print(f'AlfaMin: {AlfaMin:>6.1f} | AlfaPer: {AlfaPer:>6.1f} |\
 AlfaMax: {AlfaMax:>3.1f} |')
print(f'Impresión: {Impres:>4d} |')
print('----------'*5)

enc_01 = '##########'*6 + '\n'
enc_02 = 'Resultados del programa Helian para análisis de hélices\n\
Versión en Python realizada por Martín Paredes\nNoviembre 2022\n'
enc_03 = '##########'*6 + '\n'
enc_04 = '----------'*6 + '\n' + f'Título:\n\n{titulo:<}'
enc_05 = '----------'*6
enc_06 = f'Hélice:\n\n{titulo1:<}\nCant. de puntos: {NDim1:<d}'
enc_07 = f'Factor de actividad AF: {Fact_Act:<.2f}'
enc_08 = '----------'*6
enc_09 = f'Perfil:\n\n{titulo2:<}\nCant. de puntos: {NDim2:<d}'
enc_10 = '----------'*6
enc_11 = 'Parámetros de la corrida:\n'
enc_12 = f'AJ0: {AJ0:>10.2f} | AJf: {AJf:>10.2f} | NJ: {NJ:>9d} |'
enc_13 = f'RPM: {rpm:>10.2f} | h [m]: {h:>8.2f} |'
enc_14 = f'Nº de palas: {B:>2d} | D [m]: {D:>8.4f} |'
enc_15 = f'X0: {X[0]:>11.2f} | DeltaX: {DeltaX:>7.2f} | NX: {NX:>9d} |'
enc_16 = f'Beta75(i): {Beta75i:>4.1f} | Beta75(f): {Beta75f:>4.1f} |\
 NBeta: {NBeta:>6d} |'
enc_17 = f'AlfaMin: {AlfaMin:>6.1f} | AlfaPer: {AlfaPer:>6.1f} |\
 AlfaMax: {AlfaMax:>3.1f} |'
enc_18 = f'Impresión: {Impres:>4d} |'
enc_19 = '----------'*6 + '\n'


lista_encabezado = [enc_01, enc_02, enc_03, enc_04, enc_05, enc_06, enc_07,
                    enc_08, enc_09, enc_10, enc_11, enc_12, enc_13, enc_14,
                    enc_15, enc_16, enc_17, enc_18, enc_19]
for i in lista_encabezado:
    print(i)

with open('Resultados-Helian.txt', 'w', encoding = 'utf-8') as f:
    f.write('\n'.join(lista_encabezado))
with open('Resultados-Helian.csv', 'w', newline = '', encoding = 'utf-8') as g:
    enc_csv = csv.writer(g)
    enc_csv.writerow(['J','Cp','Eta'])

#%% Ciclo grande - Recorre los ángulos 'NBeta' ángulos Beta

fig, grafico = plt.subplots(dpi=500)
# levels = np.linspace(0.0,1.0,10)

for i in range(NDelBe):
    V0      = VEi + DelVe
    if V0 == 0.0:
        continue
    DelBe_G = DelBe0
    DelBe_R = DelBe_G * np.pi/180
    Beta75  = Be75 + DelBe_R
    Be75_G  = Beta75 * 180/np.pi
    
    parametros = [[],[],[]]
    
    #%% Titulo B
    
    seguir = 0
    # seguir = input(f'Continuarmos? Presionar "Enter"')
    print(f'\nDelBe [º]: {DelBe_G:^4.1f} | Be75 [º]: {Be75_G:^4.1f}\n')
    with open('Resultados-helian.txt', 'a', encoding = 'utf-8') as f:
        f.write(f'\nDelBe [º]: {DelBe_G:^4.2f} | Be75 [º]: {Be75_G:^4.2f}\n\n')
    if Impres == 1:
        print('|      J     |     CT     |     CP     | Eficiencia |   Alfa(1)  |  Alfa(N-1) |    Cli     |\n')
        with open('Resultados-Helian.txt', 'a', encoding = 'utf-8') as f:
            f.write('|      J     |     CT     |     CP     | Eficiencia |   Alfa(1)  |  Alfa(N-1) |    Cli     |\n')
    #%% Ciclo interno de iteraciones
    
    for j in range(NVe-1):
        pa = ' '
        pb = ' '
        lambda_g = V0/VT
        alfa_i_0 = 0.0000001
        
        alfa_a      = np.zeros(NX)
        alfa_a_g    = np.zeros(NX)
        alfa_i      = np.zeros(NX)
        alfa_i_g    = np.zeros(NX)
        arg         = np.zeros(NX)
        beta_a      = np.zeros(NX)
        Cd          = np.zeros(NX)
        Cl          = np.zeros(NX)
        Cl_x2       = np.zeros(NX)
        dkp         = np.zeros(NX)
        dkt         = np.zeros(NX)
        F           = np.zeros(NX)
        indice      = np.zeros(NX)
        phi         = np.zeros(NX)
        Re          = np.zeros(NX)
        Ve          = np.zeros(NX)
        Ve_Vt       = np.zeros(NX)
        wa          = np.zeros(NX)
        wa_VT       = np.zeros(NX)
        wt_VT       = np.zeros(NX)


        for k in range(NX):
            check = 0
            VE_VT_0     = np.sqrt( X[k]**2 + lambda_g**2)
            phi[k]      = np.arctan( lambda_g/X[k] ) 
            beta_a[k]   = beta_x[k]+DelBe_R
            
            if X[k] == X[NX-1]:
                fin         = 1
                alfa_a_g[k] = -5.07
                alfa_a[k]   = alfa_a_g[k]*np.pi/180
                Cl[k]       = 0
                Cl_x2[k]    = 0
                Cd[k]       = cs4(alfa_a[k])
                VR_VT       = np.sqrt(X[k]**2+lambda_g**2)
                phi[k]      = np.arctan(lambda_g/X[k])
                alfa_i[k]   = beta_a[k]-phi[k]
                alfa_i_g[k] = alfa_i[k]*180/np.pi
                wt_VT[k]    = VR_VT*np.sin(alfa_i[k])*np.sin(beta_a[k])
                wa_VT[k]    = VR_VT*np.sin(alfa_i[k])*np.cos(beta_a[k])
                Ve_Vt[k]    = np.sqrt((X[k]-wt_VT[k])**2+(lambda_g+wa_VT[k])**2)
                dkt[k]      = -sigma_x[k]*Ve_Vt[k]**2*Cd[k]*np.sin(beta_a[k])
                dkp[k]      = sigma_x[k]*Ve_Vt[k]**2*Cd[k]*np.cos(beta_a[k])
                
                if Impres == 2:
                        Re[k] = cuerda_x[k]*Ve_Vt[k]*VT/nu
                        Ve[k] = Ve_Vt[k]*VT
                if fin == 1:
                    continue
            
            Error       = 1.0
            iteracion   = 0
            while (Error > tol) and (iteracion < 100):
            # for l in range(100):
                alfa_i[k]   = alfa_i_0
                alfa_i_g[k] = alfa_i[k]*180/np.pi
                Ve_Vt[k]    = VE_VT_0
                arg[k]      = phi[k] + alfa_i[k]
                alfa_a[k]   = beta_a[k] - phi[k] - alfa_i[k]
                alfa_a_g[k] = alfa_a[k]*180/np.pi
                
                if ( alfa_a_g[k] > AlfaMax ) or ( alfa_a_g[k] < AlfaMin ):
                    n       = Omega/(2*np.pi)
                    AJ      = V0/(n*D)
                    print(f' {AJ:^12.2f}-----','   Alguna estación se encuentra fuera de la curva cl-alfa   ','-----')
                    with open('Resultados-Helian.txt', 'a', encoding = 'utf-8') as f:
                        f.write(f' {AJ:^12.2f}  -----    Alguna estación se encuentra fuera de la curva cl-alfa     -----')
                        f.write('\n')
                    
                    check = 1
                    break
                
                Cl[k] = cs3(alfa_a_g[k])
                f = B/2*(1-X[k])/(X[k]*np.sin(arg[k]))
                
                if f < 50.0:
                    G = np.exp(-f)
                    F[k] = 2/np.pi*np.arccos(G)
                else:
                    F[k] = 1.0
                    
                wt_VT[k]    = sigma_x[k]*Ve_Vt[k]*Cl[k]/(8.0*X[k]*F[k])
                wa[k]       = np.sqrt(lambda_g**2 + 4.0*wt_VT[k]*(X[k]-wt_VT[k]))
                wa_VT[k]    = 0.5*(-lambda_g+wa[k])
                VE_VT_0     = np.sqrt((X[k]-wt_VT[k])**2+(lambda_g+wa_VT[k])**2)
                alf         = (lambda_g+wa_VT[k])/(X[k]-wt_VT[k])
                alfa_i1     = np.arctan(alf) - phi[k]
                Error       = np.abs(alfa_i1-alfa_i[k])
                indice[k]   = iteracion
                
                if Error <= tol:
                    Cd[k]       = cs4(alfa_a[k])
                    kt1         = sigma_x[k]*Ve_Vt[k]**2
                    kt2         = Cl[k]*np.cos(arg[k])-Cd[k]*np.sin(arg[k])
                    dkt[k]      = kt1*kt2
                    kp1         = kt1*X[k]
                    kp2         = Cl[k]*np.sin(arg[k])+Cd[k]*np.cos(arg[k])
                    dkp[k]      = kp1*kp2
                    Cl_x2[k]    = Cl[k]*X[k]*X[k]
                    if Impres == 2:
                        Re[k] = cuerda_x[k]*Ve_Vt[k]*VT/nu
                        Ve[k] = Ve_Vt[k]*VT
                    break
                else:
                    alfa_i_0 = (alfa_i1+alfa_i[k])/2
                
                if iteracion == 99:
                    print(f'Alfa_a no converge para x ={X[k]:>.2f}, V ={V0:>.2f}, Beta75 = {Be75_G:>.2f}')
                iteracion += 1
            if check == 1:
                break
            if ( alfa_a_g[k] > AlfaPer) or ( alfa_a_g[k] < alfa[0] ):
                pa = '('
                pb = ')'
            #%% Resultados
            
        if check == 1:
            V0 += DelVe
            continue
        
        A1 = X[0]
        B1 = X[NX-1]
        
        cs5 = CubicSpline(X, Cl_x2, bc_type = 'natural')
        Integral_cs5 = integrate.simpson(Cl_x2, X, even = 'first')
        CL_I = 3*Integral_cs5
        
        cs6 = CubicSpline(X, dkt, bc_type = 'natural')
        Integral_dkt = integrate.simpson(dkt, X, even = 'first')
        
        cs7 = CubicSpline(X, dkp, bc_type = 'natural')
        Integral_dkp = integrate.simpson(dkp, X, even = 'first')
        
        T       = 1.5708*rho*(R**2)*(VT**2)*Integral_dkt
        P       = 1.5708*rho*(R**2)*(VT**3)*Integral_dkp/76.05
        ETA     = T*V0/P/76.05
        n       = Omega/(2*np.pi)
        AJ      = V0/(n*D)
        CT      = 3.8759*Integral_dkt
        CP      = 12.1761*Integral_dkp
        ETA_1   = AJ*CT/np.abs(CP)
        
        V0 += DelVe
        
        if ETA_1 == 0.0:
            PID = T**1.5/(R*np.sqrt(6.2832*rho))
            ETA_1 = PID/(P*76.05)
        
        if Impres == 1:
            resultados = [pa,f'{AJ:^11.2f} {CT:^12.4f} {CP:^12.4f} {ETA_1:^12.4f} {alfa_a_g[0]:^12.2f} {alfa_a_g[NX-2]:^12.2f} {CL_I:^11.4f}',pb]
            print(' '.join(resultados))
            with open('Resultados-Helian.txt', 'a', encoding = 'utf-8') as f:
                f.write(' '.join(resultados))
                f.write('\n')
            if pa == ' ':
                parametros[0].append(AJ)
                parametros[1].append(CP)
                parametros[2].append(ETA_1)
                with open('Resultados-Helian.csv', 'a', newline = '', encoding = 'utf-8') as g:
                    fila = csv.writer(g)
                    fila.writerow([AJ, CP, ETA_1])
            if ETA_1 < 0.:
                print('**********'*6)
                # DelBe0 += DelBei
                # break
                continue
            # continue
        else:
            print(pa, f'{AJ:^12.2f}{CT:^12.2f}{CP:^12.2f}{ETA_1:^12.2f}\
                  {alfa_a_g[0]:^12.2f}{alfa_a_g[NX-2]:^12.2f}\
                      {CL_I:^12.2f}',pb)
            print(f'({X[k]:^.2f}, {alfa_a_g[k]:^.2f}, {cl[k]:^.2f}, \
                  {cd[k]:^.2f}, {alfa_i_g[k]:^.2f}, {Ve[k]:^.2f},\
                  {Re[k]:^.2f}, k=1, {NX:^d})')
            if ETA_1 < 0.:
                print('**********'*6)
                # DelBe0 += DelBei
                continue
            # break
    DelBe0 += DelBei
    
    for i, rend in enumerate(parametros[2]):
        if rend <= 0.1:
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[0/255,0/255,255/255]])
        elif (rend > 0.1) and (rend <= 0.2):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[0/255,128/255,255/255]])
        elif (rend > 0.2) and (rend <= 0.3):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[0/255,255/255,255/255]])
        elif (rend > 0.3) and (rend <= 0.4):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[0/255,255/255,128/255]])
        elif (rend > 0.4) and (rend <= 0.5):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[0/255,255/255,0/255]])
        elif (rend > 0.5) and (rend <= 0.6):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[128/255,255/255,0/255]])
        elif (rend > 0.6) and (rend <= 0.7):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[255/255,255/255,0/255]])
        elif (rend > 0.7) and (rend <= 0.8):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[255/255,128/255,0/255]])
        elif (rend > 0.8) and (rend <= 0.9):
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[255/255,0/255,0/255]])
        else:
            grafico.scatter(parametros[0][i], parametros[1][i], c = [[204/255,0/255,0/255]])
    grafico.set_xlim(0.0, AJf)
    grafico.set_ylim(0.0, max(parametros[1]))
    grafico.set_xlabel('\u03B7')
    grafico.set_ylabel('Cp')
    plt.grid(visible = True, which = 'both', axis = 'both')
plt.show()


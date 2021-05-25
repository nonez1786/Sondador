#==============================================================================
#Titulo: Parametrizar_Radio_Enlace
#Autor: Zenon Saavedra 
#==============================================================================     
"""
Parametriza las magnitudes presentes en el proceso Transmision-Recepcion:

    Atenuacion:    
    Clutter: Potencia, fBragg, PDF
    Ruido: Potencia, PDF
    Objetivo: RCS, fD, PDF, retardo, Rango_oblicuo, Rango_terrestre      
        
Created on Fri Mar 19 10:46:53 2021
@author: Zenon Saavedra
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import Funciones
import Generador_Senal
import math 


# Constantes
#==============================================================================
C = 300000000.0 # Velocidad de la luz en el vacio

# Parametros
#==============================================================================

#
Latitud_Tx = 10 # Latitud Geografica [º decimales]
Longitud_Tx = 10 # Longitud Geografica [º decimales]
Hora = 14 # Hora 24 hs
Fecha = "13-05-2021" # ddmmaa
dia,mes,anio = Fecha.split("-")
Rz = 10 # Nivel de Actividad Solar


# Transmisor
Ptx = 500e3 # [W]
Gtx = 30 # [dBi]
Rtx = 50 # [ohm]
Plarizacion = "V" # V o H
Fc = 2e6 # [Hz]
AB = 500e3 # [Hz]
T_pulso = 2e-6 # [s]
N = 1 # Num de integraciones
num_samples =1 #int(T_pulso*fs_banda_base)
f_start = 0 # [Hz]
f_stop = AB # [Hz]
PRF = 50 # [Hz]

# Receptor
Grx = 20 #[dBi]
Rrx = 50 # [ohm]

# Clutter
Estado_Mar = 1
Gain_Clutter = 20 # [dB]
PDF_Clutter = 'Ninguna'# ['Rayleigh','K','Lognormal']
Media_Bragg = 0  # PDF Lognormal
Sigma_Bragg = 0  # PDF Lognormal
Rayleigh_escala_Bragg = 0 # PDF Rayleigh
fB = 0,5 # [Hz]

# Noise
P_Noise = 20 # [W]
PDF_Noise = 'Ninguna' #['Lognormal']
Sigma_Noise = 1
Media_Noise = 1

# Objetivo
Gain_RCS = 20 # [dB]
PDF_RCS = 'Ninguna' #['Swerling1','Swerling2','Swerling3','Swerling4']
fD = 10 # [Hz]
Latitud_Objetivo = 10 # Latitud Geografica [º decimales]
Longitud_Objetivo = 10 # Longitud Geografica [º decimales]
Altitud_Objetivo = 0 # Altitud geografica [m]


# Atenuacion 
Atte = -100 # [dB]

# Parametros de muestro
FsBB = 10e3 # [Hz]
FsRF = 40e6 # [Hz]



class Posicion_Geo:
    def __init__(Pos,Latitud,Longitud,Altitud):
        Pos.Latitud = Latitud
        Pos.Longitud = Longitud
        Pos.Altitud = Altitud
        

class Objetivo:
    def __init__(Obj,RCS,Retardo,fD,Posicion):
        Obj.RCS = RCS
        Obj.Retardo = Retardo
        Obj.fD = fD
        Obj.Posicion = Posicion
        
        



def Trazo_Rayos(fc, Lat_Tx, Lon_Tx, Fecha, Rz, Ang_Elev, Ang_Azi):
    """
    Determina el camino de propagacion de la O.E y una serie de parametros
    relacionados con ello
        
    float Trazo_Rayos(float fc,float Lat_Tx, float Lon_Tx,string Fecha,float Rz)
        Entrada:
            fc: Frec. de Portadora[Hz]
            Lat_Tx: Latitud del Tx [º decimales]
            Lon_Tx: Longitud del Tx [º decimales]
            Fecha: dia mes anio [ddmmaa]
            Rz: Numero de manchas solares
            Ang_Elev: 
            Ang_Azi:
        Salida: 
            Retardo: Ratardo de la O.E (ida) [s]
            Rango_Terrestre: Distancia medida sobre la tierra [m]
            Rango_Oblicuo: Distancia seguida por la O:E [m]
            Atte_Des: Atenuacion con desviación 
            Lat_Obj: Latitud del Objetivo [º decimales]
            Lon_Obj: Longitud del Objetivo [º decimales]
            Alt_Obj: Altitud del Objetivo [m]
    """
    
    Retardo = 1
    Rango_Terrestre = 1 
    Rango_Oblicuo = 1
    Atte_Des = 1
    Lat_Obj = 1
    Lon_Obj =1
    Alt_Obj = 1
    
            
    return (Retardo, Rango_Terrestre, Rango_Oblicuo, 
            Atte_Des, Lat_Obj, Lon_Obj,Alt_Obj)


#------------------------------------------------------------------------------
def RCS():
	print("RCS")
      
#------------------------------------------------------------------------------
def Atenuacion(fc, Rango_Oblicuo, Atte_Des):
    """
    Determina el nivel de atenuacion sufrido por la onda durante la 
    propagacion por el medio (camino de ida).
        
    float Atenuacion(float fc, float Rango_Oblicuo)
        Entrada:
            fc: Frecuencia de portadora [Hz]
            Rango_Oblicuo: Rango Oblicuo seguido por la O. de Radio [Hz]
        Salida: 
            Atte: nivel de atenuacion en el camino de propagacion [dB]
    """ 
    Atte_Geo = 32.5 + 20*math.log10(fc/(1e6)) + 20*math.log10(Rango_Oblicuo/1e3)
    Atte_No_Des = 1
    Atte_Des = Atte_Des
    
    Atte = Atte_Geo + Atte_No_Des + Atte_Des
    
    print("Atenuacion:",Atte, "\n")
    return Atte 
#------------------------------------------------------------------------------

def Potencia_Eco(PTx, Gtx, Atte, RCS, Grx):
    """
    Determina el nivel de potencia de la señal de Eco.
        
    float Potencia_Eco(float PTx,, float Gtx, float RCS, float Atte, float Grx)
        Entrada:
            PTx: [dBm]
            Gtx: [dBi]
            RCS: [dBm2]
            Atte: [dB]
            Grx: [dBi]
        Salida: 
            Peco: nivel de potencia de eco en la entrada del receptor [dBm]       
    """
    Peco = PTx + Gtx - Atte + RCS - Atte - Grx
    
    print("Potencia de Eco: ",Peco, "\n")
    return Peco
	
#------------------------------------------------------------------------------
def Potencia_Clutter(Estado_Mar, fc, Rango_Terrestre, Delta_Azi, Delta_Tau, PTx, Gtx, Atte, Grx):
    """
    Determina el nivel de Potencia de la señal de Clutter y la Frec Bragg.
    
    float float Potencia_Clutter(int Estado_Mar, float fc, float Rango_Terrestre, float Delta_Azi,
                                 float Delta_Tau, float PTx, float Gtx, float Atte, float Grx):
        Entrada:
            Estado_Mar: Los estados son: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
            fc: Frec. de Portadora[Hz]
            Rango_Terrestre: [m]
            Delta_Azi: Ancho de banda del Haz de la Antena Rx [º]
            Delta_Tau: Ancho del sub pulso   1/AB [s]
            PTx: Potencia Transmitida [dBm]
            Gtx: Ganancia de la antena Rx [dB]
            Atte: Atenuacion en el camino de propagacion (Ida) [dB] 
            Grx: Ganancia de la antena Rx [dB]
        Salida: 
            PClutter: nivel de potencia de clutter en la entrada del receptor [dBm]       
            fB: Frecuencia de Bragg  [Hz]
    """
    
    Coef_Dispersion = ([-50,-42,-34,-26,-18,-10,-2,6,10,15]) # Coeficiente de dispersion [dB]
                                                             # Estado del mar 0-9
    Res_Rango_Cruzado = Delta_Azi * Rango_Terrestre
    Res_Rango = (C * Delta_Tau/2)
    Celda_Res = Res_Rango * Res_Rango_Cruzado  # Dimension de Celda de resolucion [m^2]
    Celda_Res_dB = 10 * math.log10(Celda_Res)  # [dBm^2] 
           
    Coef = (Coef_Dispersion[Estado_Mar])
    RCS_Clutter = Coef + Celda_Res_dB          # RCS del area celda_res [dBsm]

    PClutter = PTx + Gtx - Atte + RCS_Clutter - Atte - Grx
    
    fB = 0.102 * (fc/1e6) # Frec de Bragg
    
    return (PClutter,fB)
    
#------------------------------------------------------------------------------
def Potencia_Ruido(fc, AB, Fecha, Lat_Tx, Lon_Tx, Rz, Tipo_Rui_Hum):
    """
    Determina el nivel de potencia de la senal de Ruido.
        
    float Potencia_Ruido(float float, string Fecha, float Lat_Tx, float Lon_Tx, float Rz)
        Entrada:
            fc: Frec. de portadora [Hz]
            AB: Ancho de Banda [Hz]
            Fecha: dia mes anio [ddmmaa]
            Lat_Tx: Latitud geo. del transmisor [º decimal]
            Lon_Tx: Longitud geo. del transmisor[º decimal]
            Rz: Numero de Manchaz solares[dB]
            Tipo_Rui_Hum:
        Salida: 
            Pn: nivel de potencia de ruido en la entrada del receptor [dBm]       
            Sigma: Desviacion estandar pa un PDF Log_Normal

    """
    Fa_Atm = 0 ; Sigma_Atm = 0
    Fa_Gal = 0 ; Sigma_Gal = 0
    Fa_Hum = 0 ; Sigma_Hum = 0
    Fa = 0 ; Sigma = 0
    
    # Ruido Atmosferico
    f = fc/1e6
    if (f > 0.01 & f <= 20):
        frec = ([0.010, 0.020, 0.040, 0.060,0.080, 0.10, 0.20, 0.40, 0.60, 0.80, 1, 2, 4, 6, 8, 10, 20])
        FA = ([144, 128, 105, 90, 82, 72, 48, 20, 10, 5, 0, 6, 20, 28, 29, 28, 0])
        Fa_Atm = np.interp(frec,FA,f)
        Sigma_Atm = 6  # Este valor debe ser dertermina desde tablas en funcion de
                   # la ubicacion, hora del dia y estacion del año
        Fa = np.append(Fa,Fa_Atm)
        Sigma = np.append(Sigma,Sigma_Atm)
    
    # Ruido Galactico               
    if (f >= 5 & f < 100):    
        Fa_Gal = 52 - 23*math.log10(f)
        Sigma_Gal = 1.56    
        Fa = np.append(Fa,Fa_Gal)
        Sigma = np.append(Sigma,Sigma_Gal)
        
    # Ruido Humano          
    if (f >= 0.3 & f < 250):
        if(Tipo_Rui_Hum == 'Comercial'):
            c = 76.8
            d = 27.7
            sigm = 8.4
        elif(Tipo_Rui_Hum == 'Residencial'):           
            c = 72.5
            d = 27.7
            sigm = 5.8
            
        elif(Tipo_Rui_Hum == 'Rural'):           
            c = 67.2
            d = 27.7
            sigm = 6.8        
            
        elif(Tipo_Rui_Hum == 'Rural Tranquila'):           
            c = 52
            d = 23
            sigm = 6.8        
        else:
            print("\n Error: opcion No valida. Ruido Humano\n") 
        Fa_Hum = c - d*math.log10(f)
        Sigma_Hum = sigm   
        
        Fa = np.append(Fa,Fa_Hum)
        Sigma = np.append(Sigma,Sigma_Hum)

    
    n = np.size(Fa)
    alpha = 0 ;  beta = 0  ;  c = 10/(math.log10(10))
    for i in np.arange(n):
        alpha = alpha + math.exp(Fa(i)/c + (Sigma(i))**2/(2*c**2))
        beta = beta + alpha**2*(math.exp(((Sigma(i))**2)/(c**2))-1)
    
    Desv = c*math.sqrt(math.log10(1+(beta/(alpha**2))))
    Fma = c*(math.log10(alpha)-Desv**2/(2*c**2))
    
    
#    [Fna_dB desv_dB] = Fma_total(frec); [dB dB]
    BdB = 10*math.log10(AB) # Ancho de banda en dBHz
    k = 1.38e-23     # constante de Bolzman
    T0 = 290         # Temperatura en Kº
    # Pn[dBW]=Fna[db] + B[dB] - 204[dB]
    # (n0 = k*T0*B Densidad de Ruido Blanco)   10log(kT0)=204 dB_Kel
    Pn = Fma + BdB - 204  #potendia de ruido en dBW
        
    print("Potencia de Ruido: ",Pn, "\n")
    return (Pn , Desv)
	
    

#==============================================================================

# Funcion Principal
def main():
    print('\n=============================================================')
    print('                        INICIO                                 ')
    print('=============================================================\n')
    PDF_RCS = Funciones.RCS_pdf('Swerling4')
    print(PDF_RCS)
    
    Atenuacion(2e6,2000e3,1)
    
    Posicion_Tx = Posicion_Geo(30,40,0)
    
        
    Objetivo1 = Objetivo(10,11,12,Posicion_Tx)
    print("Latitud: ",Objetivo1.Posicion.Latitud)
    
    P_Clutter,fB = Potencia_Clutter(9, 10e6,10e3, 1,1e-3,100, 20, 100, 25)
    print("Frec. Bragg: ", fB)
    
# INICIO    
if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:10:49 2024

@author: Asus
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from sympy.solvers import solve
import math as mt
from sympy import Symbol

"variables"
ar= 1
ad=0.23
T= 198
y= 62
b3=10.2
b2=9.4
k61=0.0333
k56=0.0862
k51= 0.0862
k43=0.00114
k42=0.00189
k34=0.714
k31= 18*(140**(-b3))
k24=0.0164
k23=0.0781
k21= 58*(730**(-b2))
k15= 147
k13=0.0311
k12=0.0931
N_mean= 324.014
a1= -2.4*10**(-7)
b1=7.2*10**(-4)
c1=-2.1*10**(-4)
c0=278



"intervalo de tiempo"
t=np.linspace(1950,2100,151)
t_acum=np.linspace(1850,2100,2100-1850+1)

"2000-2010 y 2010-2019"

"funciones de F"
Fr=0
def Fd(t:int,t0:int)->float:
    R=0.3+0.01*(t-t0)
    return R

def Ff(t:int,t0:int, coef:float)->float:  
    r= 1.4+ (coef)*(t-t0-100)
    return r


"PUNTO 2.1; t0:1850 (tomar como epoca preindustrial)"

"-valores reales"
reales=[315.98,316.91,317.64,318.45,318.99,319.62,320.04,321.37,322.18,323.05,324.62,325.68,326.32,327.46,329.68,330.19,331.12,332.03,333.84,335.41,336.84,338.76,340.12,341.48,343.15,344.85,346.35,347.61,349.31,351.69,353.2,354.45,355.7,356.54,357.21,358.96,360.97,362.74,363.88,366.84,368.54,369.71,371.32,373.45,375.98,377.7,379.98,382.09,384.03,385.83,387.64,390.1,391.85,394.06,396.74,398.87,401.01,404.41,406.76,408.72,411.66,414.24]
tiempito=[1959,1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020]
t_real=np.linspace(1959,2020,61)


defo=[]
fosiles=[]
total=[]

i=0
ndefo=0
nfosiles=0
sumadefo=0
sumafo=0
tut=0

"emisiones totales acumuladas"

while i<251:
    t0=1850
    coe= 4.6/40  
    sumadefo+= (Fd(t_acum[i],t0))
    sumafo+= (Ff(t_acum[i],t0,coe))
    defo.extend([sumadefo])
    fosiles.extend([sumafo])    
    total.extend([sumadefo+sumafo])
    i+=1


"emisiones sin acumulación"
defo2=[]
fosiles2=[]
total2=[]

sumadefo=0
sumafo=0
tut=0
i=0
while i<151:
    t0=1850
    coe= 4.6/40  
    sumadefo= (Fd(t[i],t0))
    sumafo= (Ff(t[i],t0,coe))   
    total2.extend([sumadefo+sumafo])
    i+=1

plt.plot(t, total2)
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel(" PgC", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('E1: Emisiones totales acumuladas', fontdict={'fontweight':'bold'})
plt.legend(["Fd","Ff","emisiones totales"])
plt.grid()
plt.show()

"-valores del modelo"
def fun(A,t):
    M1,M2,M3,M4,M5,M6,M7,G=A
    t0=1850
    coe= 4.6/40
    uno= -(k12+k13)*M1-k15*G*((M1-y)/(M1+T))+k21*(M2**(b2))+ k31*(M3**b3)+k51*M5+k61*M6+Ff(round(t),t0,coe)+Fd(round(t),t0)-Fr
    dos= k12*M1-(k23+k24)*M2-k21*(M2**b2)+k42*M4
    tres= k13*M1+k23*M2-k34*M3-k31*(M3**b3)+k43*M4
    cuatro= k24*M2+k34*M3-(k42+k43)*M4
    cinco= k15*G*((M1-y)/(M1+T))-(k51+k56)*M5-Fd(round(t),t0)+Fr
    seis= k56*M5-k61*M6
    siete=-Ff(round(t),t0,coe)
    ge=-((ad*Fd(round(t),t0)-ar*Fr)/f0[4])    
    return [uno,dos,tres,cuatro,cinco,seis,siete,ge]

"definir valores iniciales"
f0= [612,730,140,37000,580,1500,5300,1]

"solución ecuaciones"
sol= odeint(fun,f0,t)
sol2= odeint(fun,f0,t)
sol_atm=sol[:,0]



for m in range(0,len(sol_atm)):
    sol_atm[m]/=2.16
    sol_atm[m]+=20

cortado=sol_atm[9:62+9]

"Gráfica M1 Busines as usual"
plt.plot(t, sol_atm[:])  
plt.plot(tiempito, cortado[:])  
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel(" [ppm CO2]", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('E1: Concentraciones M1', fontdict={'fontweight':'bold'})
plt.legend(["M1"])
plt.grid()
plt.show()

"FORZAMIENTO"
def forz(c):
    cal=(a1*((c-c0)**2)+b1*abs(c-c0)+c1*N_mean+5.36)*np.log(c/c0)
    return cal

def temp(f):
    cal= (f*0.52)
    return cal


forza=[]
tempe=[]

for i in range(0, len(sol_atm)):
    forza.extend([forz(sol_atm[i])])
    tempe.extend([temp(forza[i])])
    

nuevo=[]
nuevo=sol_atm[70]


print(forza[30],"forzamiento", tempe[30], "temepratura",sol_atm[30], "concentración")

plt.plot(t, tempe, color="purple")  
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel("Temperatura °C", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('Temperatura vs tiempo: E1', fontdict={'fontweight':'bold'})
plt.grid()
plt.show()

plt.plot(sol_atm, tempe,color="green")  
plt.xlabel("atm Co2 (ppm)", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel("Temperatura °C", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('Temperatura vs M1: E1', fontdict={'fontweight':'bold'})
plt.grid()
plt.show()

plt.plot(t, forza, color="purple")  
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel("Forzamiento  Wm^-2 ", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('Forzamiento vs tiempo: E1', fontdict={'fontweight':'bold'})
plt.grid()
plt.show()

plt.plot(sol_atm, forza, color="green")  
plt.xlabel("atm CO2", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel(" Forzamiento Wm^-2 ", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('Forzamiento vs M1: E1', fontdict={'fontweight':'bold'})
plt.grid()
plt.show()

plt.plot(tempe, forza, color="orange")  
plt.xlabel("Temperatura °C", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel(" Forzamiento Wm^-2 ", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('Forzamiento vs Temperatura: E1', fontdict={'fontweight':'bold'})
plt.grid()
plt.show()


"Punto c: estimar el tiempo, limite concentración: 1100 Gt CO2 acum del 2011 al 2050"

"-Tiempo en llegar al límite de emisiones"
c=2011
sumadefo=0
sumafo=0
tot=0
s=True
while s:   
    t0=1850
    coe= 4.6/40     
    sumadefo= (Fd(c,t0))/2.16
    sumafo= (Ff(c,t0,coe))/2.16 
    tot+=sumadefo+sumafo
    if tot>=1100/7.92:
       s=False       
       print("Límite en el año:",c)
       print("Límite:",tot)
       print("tiempo en llegar desde el 2011:",c-2011,"años")
    c+=1

"E4"
F=[]
mu=[]
m_ppm=sol_atm[80]
m_1=sol2[80,0]
#print(m_ppm,m_1)
#print(t[80])
s_m2=sol2[80:,1]
s_m3=sol2[80:,2]
s_m4=sol2[80:,3]
s_m5=sol2[80:,4]
s_m6=sol2[80:,5]
s_m7=sol2[80:,6]
s_g=sol2[80:,7]
s_t=t[80:]

for m in range(0,len(s_m2)):
   mu.extend([m_1/2.16])
   
#print(mu[0],"este")
    

for m in range(0,len(s_m2)):
    cal=(k12+k13)*m_1+k15*s_g[m]*((m_1-y)/(m_1+T))-k21*(s_m2[m]**(b2))-k31*(s_m3[m]**b3)-k51*s_m5[m]-k61*s_m6[m]
    F.extend([cal])
   
"Para la gráfica de emisiones totales sin acumulación"
plt.plot(s_t, F[:])
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel("PgC", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('E4 emisiones totales', fontdict={'fontweight':'bold'})
plt.legend(["Fd+Ff"])
plt.grid()
plt.show()


"Gráfica emisiones totales acumuladas"
acumE4=0
E4=[]
for m in range(0,len(F)):
    acumE4+=F[m]
    E4.extend([acumE4])
    
t_E4=np.linspace(2030,2100,2100-2030+1)
plt.plot(t_E4, E4[:])
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel("PgC", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('E4 Emisiones totales acumuladas', fontdict={'fontweight':'bold'})
plt.legend(["Fd","Ff","Emisiones totales"])
plt.grid()
plt.show()    

"Gráfica concentración M1 E4 tarea 3 "
plt.plot(s_t, mu[:])
plt.xlabel("Año", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.ylabel("ppm CO2", fontdict = {'fontsize':13, 'fontweight':'bold'})
plt.title('E4 concnetraciones M1', fontdict={'fontweight':'bold'})
plt.legend(["M1"])
plt.grid()
plt.show()


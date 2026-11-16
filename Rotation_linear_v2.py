# -*- coding: utf-8 -*-
import re
import math as mt
import scipy.constants as scon
from scipy.integrate import solve_ivp
from sympy.physics.wigner import wigner_3j
import matplotlib.pyplot as plt
#Definition and Data entry------
datafile=open("data.dat",'r')
lines=datafile.readlines()
temp=[]
for line in lines:
    line=line.strip()
    line=re.split(r'[\s, ;-]+',line)
    temp.append(line[0])
datafile.close()

Molecule_name=temp[0]
B=float(temp[1])
alpha_par=float(temp[2])
alpha_per=float(temp[3])
spin_even=float(temp[4])
spin_odd=float(temp[5])
temperature=float(temp[6])
J_MAX_user=int(temp[7])
pulse_width=float(temp[8])
eval_time=float(temp[9])
wave_length=float(temp[10])
I0 = float(temp[11])

del temp
##definition and data entry------
#
##Standardization of data------
J_MAX = J_MAX_user
delta_alpha=(alpha_par-alpha_per)
part_limit = 0.00001            #partition function limit
j=0
temp_part=1.                     #any number bigger than part_limit
partition_function=0.
J_0=0
tpulsed=pulse_width*10**(-15)/(2.41888432*10**(-17))
b=B/(8065.54*27.2114)                            #cm-1 to hartree
tfree=eval_time*10**(-15)/(2.41888432*10**(-17))
tcen=2*tpulsed                                   #pulse peak
tend=4*tpulsed                                   #field is taken as zero after this
c=[]
tpts=3000
time=[((tend-tcen)*(tpts-(i+1))+(tfree*((i+1)-1)))/(tpts-1) for i in range(tpts)]   #Arithematic progression, t=0 is the pulse peak
w=2*mt.pi*3*10**(14)/0.8
E0=mt.sqrt(I0*1e12/3.50945e16)                   #I0 in TW/cm2, E0 in a.u.
V_per=-0.25*alpha_per*E0**2
V_diff=-0.25*delta_alpha*(E0)**2
a=2*(mt.log(2)/(tpulsed)**2)                     #pulse_width is FWHM of intensity
norm=0.
kT = temperature/315774.64                       #K to hartree
top_pop=0.




##standardization of data------
##
#ODE solution Module
def tdse_field(j_0,j_max,m):
    global jmin,vc,top_pop
    if ((j_max-j_0)%2!=0):                       #J changes by 2 only, so same parity as j_0
        j_max -=1
    if ((j_0-m)%2==0):                           #lowest J possible for this m
        jmin=m
    else:
        jmin=m+1
    n=int((j_max-jmin)*0.5+1)
    cjmi=[]
    for cjminit in range (n):
        if (jmin+2*cjminit==j_0):
            cjmi.append(complex(1))
        else:
            cjmi.append(complex(0))
    vc=[]                                        #cos2theta matrix, calculated once per j_0,m
    for i in range (n):
        vc.append([])
        for j in range (n):
            if (i==j) | (abs(i-j)==1):
                vc[i].append(vcos(jmin+2*i,jmin+2*j,m))
            else:
                vc[i].append(0.)
#    print(j_0,j_max,m,cjmi,len(cjmi))
    res = solve_ivp(fcjm, [0,tend], cjmi, rtol=1e-8, atol=1e-10)
    cfin=[res.y[ffj][-1] for ffj in range(n)]
    avg=0.
    for ffj in range(n):
        avg += abs(cfin[ffj])**2
    if (abs(cfin[-1])**2>top_pop):
        top_pop=abs(cfin[-1])**2
    d=[]
    for t in range(tpts):
        cft=[]
        for ffj in range(n):                     #field free evolution
            cft.append(cfin[ffj]*mt.e**(complex(0,-1)*(b*(jmin+2*ffj)*((jmin+2*ffj)+1))*(time[t]-(tend-tcen))))
        dt=0.
        for ffj in range(n):#evaluating final cos2theta
            for ffj1 in range(n):
                if  ((ffj==ffj1) | (abs(ffj-ffj1)==1)) :
                    dt += (complex.conjugate(cft[ffj1])*cft[ffj]*vc[ffj][ffj1])
        d.append(complex(dt).real)
    return(d,avg)

#Function ODE Solver
def fcjm(t,cjm):
    f=[]
    env=(mt.e**(-a*(t-tcen)**2))**2
    for i in range (len(cjm)):
        sum=0.
        for j in range (len(cjm)):
            if (i==j) | (abs(i-j)==1):
                sum += cjm[j] * vc[i][j]
        f.append(complex(0,-1)*((b*(jmin+2*i)*(jmin+2*i+1)+V_per*env)*cjm[i]+V_diff*env*sum))
    return f

#function ode solver

##VCOS2THETA using wigner-3j matrices
def vcos(x,y,m):
    vcos2theta=((-1)**m)*(2/3)*(mt.sqrt((2*x+1)*(2*y+1)))*float(wigner_3j(x,2,y,0,0,0)*wigner_3j(x,2,y,-m,0,m))#wigner function
    if (x==y):
        vcos2theta += (1/3)
    return vcos2theta
##vcos2theta using wigner-3j matrices

##Estimating J_MAX and Partition function
while (temp_part>part_limit):
   temp_part= mt.e**(-j*(j+1)*b/kT)*(2*j+1)
   if (j%2==0):
       g_i=spin_even
   else:
       g_i=spin_odd
   partition_function += g_i*temp_part
   j+=1
   if (j==J_MAX_user):
       print("J out of range, select a bigger J_MAX")
       exit()
max_J=j
#estimating J_MAX and partition function

#Verifying J_MAX
print("J_MAX_user,J_MAX,k")
k=0
sum_w=0
max_w = 0.999
while (sum_w<max_w):
   if (k%2==0):
       g_i=spin_even
   else:
       g_i=spin_odd
   sum_w += g_i* mt.e**(-k*(k+1)*b/kT)*(2*k+1)/partition_function
   k+= 1
print(J_MAX_user,max_J,k)
##verifying J_MAX
#
complete_data=[]
for v in range(tpts):
    complete_data.append(0)

#ODE solution and field free rotation incorporated
for j in range(J_0,max_J+1):
   j_0=j
   if (j%2==0):
       g_i=spin_even
   else:
       g_i=spin_odd
   if (g_i==0): continue
   print(j)
   thermal_prob=g_i*mt.e**(-j*(j+1)*b/kT)/partition_function
   for m in range(j+1):         #m varies from 0 to J...
       temp_data,temp_norm=tdse_field(j,J_MAX,m)
       if m==0:                 #m actually varies from -J to J we'll multiply twice
           mult=1
       else:
           mult=2
       for i in range (len(temp_data)):
           complete_data[i] += mult*thermal_prob*temp_data[i]
       norm += mult*thermal_prob*temp_norm

#ode solution and field free rotation incorporated

if (top_pop>1e-6):
    print("Population reached the last J, select a bigger J_MAX",top_pop)

for i in range (len(complete_data)):
    complete_data[i]=complete_data[i]/norm

#Write to a dat file
for i in range (len(time)):
    time[i]= time[i] *(2.41888432*10**(-17))*10**(12)      #ps

output_file=open(Molecule_name+"_at_"+str(temperature)+"K_,"+str(I0)+"_TWcm2.dat","w")
for i in range (len(complete_data)):
    output_file.writelines([str(time[i]),"   ",str(complete_data[i]),"\n"])

output_file.close()
#write to a dat file

#Plotter
plt.ylabel(r'$  \langle cos^{2}\theta \rangle $')
plt.xlabel('Time (ps)')
plt.title(Molecule_name)
plt.plot(time,complete_data)
plt.show()
#plotter


#ode solution module



"""
Spyder Editor

This is a temporary script file.
"""


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
    line=re.split('[\s, ;-]+',line)
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
b=B/(8065.54*27.2)
tfree=eval_time*10**(-15)/(2.41888432*10**(-17))
c=[]
tpts=100
time=[(tpulsed*(tpts-(i+1))+(tfree*((i+1)-1)))/(tpts-1) for i in range(tpts+1)]   #Arithematic progression
w=2*mt.pi*3*10**(14)/0.8
E0=mt.sqrt(I0*28*(8*scon.pi/137))*1e-3 
V_per=-0.25*alpha_per*E0**2
V_diff=-0.25*delta_alpha*(E0)**2                 #0.002027139589724497#10**(-7)#mt.sqrt(I0)*8.5*10**(-10)*0.8*scon.m_e*scon.c*5.3e-11/(scon.e*27.2)#1.2769*10**(-6)
a=2*(mt.log(2)/(tpulsed)**2)
alper=alpha_per
aldel=delta_alpha
norm=0.
temperature = temperature/315774.64




##standardization of data------
##
#ODE solution Module
def tdse_field(j_0,j_max,m):
#    print('loop1')
    for m in range (j_0+1):
        if (j_0%2==0):
            if (j_max%2!=0):
                j_max -=1
        else:
            if (j_max%2==0):
                j_max -=1
        n=int((j_max-m)*0.5+1)
        if (n==0): break 
        if (j_0+2*(n-1)<=j_max): 
#            print('loop2',"   ",n,m,j_max,j_0)
            cjmi=[]
            for cjminit in range (n):
    #            print('loop3',n,m,j_0,j_max)           
                if (cjminit==0): 
                    cjmi.append(complex(1))
                else:
                    cjmi.append(complex(0))
#            print(j_0,j_max,m,cjmi,len(cjmi))
#            print("ODE Solver")
#            print(j_0,j_max,m)
            res = solve_ivp(fcjm, [0,2*tpulsed], cjmi)
            tt[j_0]=res.y
            avg=0.
            for t in range(tpts):
                for ffj in range(len(res.y)):#evaluating final cos2theta                
                    for ffj1 in range(len(res.y)):                
                        if  ((ffj==ffj1) |abs(ffj-ffj1)==1) :
                            tempdat= ((res.y[ffj][-1])*mt.e**(complex(0,-1)*(b*(j_0+2*ffj)*((j_0+2*ffj)+1))*(time[t])))
                            tempdat1=complex.conjugate((res.y[ffj1][-1])*mt.e**(complex(0,-1)*(b*(j_0+2*ffj1)*((j_0+2*ffj1)+1))*(time[t])))
                            avg +=(tempdat*tempdat1)
#                            print(avg)
                            d[t] += (tempdat1*tempdat*vcos(j_0+2*ffj,j_0+2*ffj1,m))    
    return(d,avg)
    
#Function ODE Solver
def fcjm(t,cjm): 
    f=[]
    for i in range (len(cjm)):
        sum=0.
        for j in range (len(cjm)):
            if (i==j) | (abs(i-j)==1):
                sum += cjm[j] * vcos(j_0+(2*i),j_0+(2*j),m)
#            print(len(cjm),'hhhhhhhhhhhhhhhhhhhhhhh',jinit)
#                print(j_0+2*i,j_0+2*j)
        f.append(complex(0,1)*((b*(j_0+2*i)*(j_0+2*i+1)-alper*((E0*mt.e**(-a*(t-tpulsed)**2))**2))*cjm[i]+aldel*((E0*mt.e**(-a*(t-tpulsed)**2))**2)*sum))
    return f
    
#function ode solver
    
##VCOS2THETA using wigner-3j matrices
def vcos(x,y,m):        
    vcos2theta=((-1)**m)*(2/3)*(mt.sqrt((2*x+1)*(2*y+1)))*wigner_3j(x,2,y,0,0,0)*wigner_3j(x,2,y,-m,0,m)#wigner function
    if (x==y):
        vcos2theta += (1/3)
    return vcos2theta    
##vcos2theta using wigner-3j matrices

##Estimating J_MAX and Partition function
while (temp_part>part_limit):
   temp_part= mt.e**(-j*(j+1)*B/(scon.k*temperature))*(2*j+1)
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
   sum_w += g_i* mt.e**(-k*(k+1)*B/(scon.k*temperature))*(2*k+1)
   k+= 1
print(J_MAX_user,max_J,k)
##verifying J_MAX
#
d=[]
complete_data=[]
for v in range(tpts):
    d.append(0),complete_data.append(0)
    


tt=[]
for v in range(J_MAX+1):
    tt.append(0)
#ODE solution and field free rotation incorporated   
for j in range(J_0,J_MAX+1):
   j_0=j
   print(j)
   if (j%2==0):
       g_i=spin_even
   else:
       g_i=spin_odd
   for m in range(j+1):         #m varies from 0 to J...
       if m==0:                 #m actually varies from -J to J we'll multiply twice                        
           temp_data,temp_norm=tdse_field(j,J_MAX,m)
       else:
           temp_data,temp_norm=tdse_field(j,J_MAX,m)
       for i in range (len(temp_data)):
           complete_data[i] += temp_data[i]
           norm += temp_norm    
       
#ode solution and field free rotation incorporated
          

#Plotter
plt.ylabel(r'$  \langle cos^{2}\theta \rangle $')
plt.xlabel('Time (ps)')
plt.title(Molecule_name)
plt.plot([complex(i/norm).real for i in (complete_data)])
#plotter
            
#Write to a dat file
for i in range (len(time)):
    time[i]= time[i] *(2.41888432*10**(-17))

    
output_file=open(Molecule_name+"_at_"+str(temperature)+"K_,"+str(I0)+"_Wcm2.dat","w")
for i in range (len(time)-1):
    output_file.writelines([str(time[i]),str(complete_data[i])])

output_file.close()
#write to a dat file
output_file=open(Molecule_name+"_at_"+str(temperature)+"K_,"+str(I0)+"_Wcm2.dat","w")
for i in range (len(complete_data)):
    output_file.writelines([str(time[i]),str(complex(complete_data[i]).real)])
    output_file.writelines('/n')
    

output_file.close()
#write to a dat file
            
       
#ode solution module       
       
    

"""
Spyder Editor

This is a temporary script file.
"""


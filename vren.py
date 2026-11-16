from scipy.integrate import solve_ivp
from sympy.physics.wigner import wigner_3j
import math as mt
import matplotlib.pyplot as plt 
import scipy.constants as scon

def fcjm(t,cjm): 
    f=[]
    for i in range (len(cjm)):
        sum=0.
        for j in range (len(cjm)):
            if (i==j) | (abs(i-j)==1):
                sum += cjm[j] * vcos(j0+(2*i),j0+(2*j),m)
#            print(len(cjm),'hhhhhhhhhhhhhhhhhhhhhhh',jinit)
#                print(j0+2*i,j0+2*j)
        f.append(complex(0,1)*((b*(j0+2*i)*(j0+2*i+1)-alper*((E0*mt.e**(-a*(t-tpulsed)**2))**2))*cjm[i]+aldel*((E0*mt.e**(-a*(t-tpulsed)**2))**2)*sum))
    return f

def vcos(x,y,m):        
    vcos2theta=((-1)**m)*(2/3)*(mt.sqrt((2*x+1)*(2*y+1)))*wigner_3j(x,2,y,0,0,0)*wigner_3j(x,2,y,-m,0,m)#wigner function
    if (x==y):
        vcos2theta += (1/3)
    return vcos2theta

j0=0
max_J=10
tpulsed=10*10**(-15)/(2.41888432*10**(-17))
b=1.99824/(8065.54*27.2)
tfree=30000*10**(-15)/(2.41888432*10**(-17))
c=[]
tpts=100
time=[(tpulsed*(tpts-(i+1))+(tfree*((i+1)-1)))/(tpts-1) for i in range(tpts+1)]   #Arithematic progression
w=2*mt.pi*3*10**(14)/0.8
I0=0.8e15
E0=10**(-7)#mt.sqrt(I0)*8.5*10**(-10)*0.8*scon.m_e*scon.c*5.3e-11/(scon.e*27.2)#1.2769*10**(-6)
a=2*(mt.log(2)/(tpulsed)**2)
alper=9.8
aldel=5.2
avg=0
#d=[[[]]]
#for u in range(jmax):
#    d.append([[]])
#    for v in range(jmax):
#        d[u].append([])
#        for w in range(tpts):
#            d[u][v].append(0)
d=[]
for v in range(tpts):
    d.append(0)

tt=[]
for v in range(max_J+1):
    tt.append(0)


        
for jinit in range(j0,max_J+1):
    j0 = jinit
    jmax=max_J
#    print('loop1')
    for m in range (j0+1):
        if (j0%2==0):
            if (jmax%2!=0):
                jmax -=1
        else:
            if (jmax%2==0):
                jmax -=1
        n=int((jmax-m)*0.5+1)
        if (n==0): break 
        if (j0+2*(n-1)<=jmax): 
#            print('loop2',"   ",n,m,jmax,j0)
            cjmi=[]
            for cjminit in range (n):
    #            print('loop3',n,m,j0,jmax)           
                if (cjminit==0): 
                    cjmi.append(complex(1))
                else:
                    cjmi.append(complex(0))
#            print(j0,jmax,m,cjmi,len(cjmi))
            print("ODE Solver")
            print(j0,jmax,m)
            res = solve_ivp(fcjm, [0,2*tpulsed], cjmi)
            tt[j0]=res.y
                        
            for t in range(tpts):
                for ffj in range(len(res.y)):#evaluating final cos2theta                
                    for ffj1 in range(len(res.y)):                
                        if  ((ffj==ffj1) |abs(ffj-ffj1)==1) :
                            tempdat= ((res.y[ffj][-1])*mt.e**(complex(0,-1)*(b*(jinit+2*ffj)*((jinit+2*ffj)+1))*(time[t])))
                            tempdat1=complex.conjugate((res.y[ffj1][-1])*mt.e**(complex(0,-1)*(b*(jinit+2*ffj1)*((jinit+2*ffj1)+1))*(time[t])))
                            avg +=(tempdat*tempdat1)
#                            print(avg)
                            d[t] += (tempdat1*tempdat*vcos(jinit+2*ffj,jinit+2*ffj1,m))
#                            print(jinit+ffj,jinit+ffj1,m,complex(tempdat*tempdat1*vcos(jinit+2*ffj,jinit+2*ffj1,m)).real)
                
#    plt.plot([complex(i).real for i in (d)])
#            print(jinit)
###        
####for i in range (jmax):
####    (c[i]**2)*vcos(i,i)
###    
plt.plot([complex(i/avg).real for i in (d)])

###print(res)

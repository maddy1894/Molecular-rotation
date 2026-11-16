from scipy.integrate import solve_ivp
from sympy.physics.wigner import wigner_3j
import math as mt
import matplotlib.pyplot as plt 

def fcjm(t,cjm): 
    f=[]
#    m=0
#    j0=0
#    jmax=10
#    sum=0.
    for i in range (len(cjm)):
        sum=0.
        for j in range (len(cjm)):
            sum += cjm[j] * vcos(jinit+i,jinit+j,m)
#            print(len(cjm),'hhhhhhhhhhhhhhhhhhhhhhh',jinit)
#            print(i,j,cjm[j])
        f.append(complex(0,1)*(b*(jinit+i)*(jinit+i+1)-alper*((E0*mt.e**(-a*(t-tpulsed)**2))**2)*cjm[i]+aldel*((E0*mt.e**(-a*(t-tpulsed)**2))**2)*sum))
    return f

def vcos(x,y,m):        
    vcos2theta=((-1)**m)*(2/3)*(mt.sqrt((2*x+1)*(2*y+1)))*wigner_3j(x,2,y,0,0,0)*wigner_3j(x,2,y,-m,0,m)#wigner function
    if (x==y):
        vcos2theta += (1/3)
#        print(m,'This is m')
    return vcos2theta

j0=0
jmax=5
tpulsed=100000*10**(-15)
b=9.1085e-06
tfree=20000*10**(-15)
c=[]
tpts=200
time=[(tpulsed*(tpts-(i+1))+(tfree*((i+1)-1)))/(tpts-1) for i in range(tpts+1)]   #Arithematic progression
E0=1.2769*10**(-4)
a=8.1112*10**(-14)
alper=9.8
aldel=5.2
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
for v in range(jmax):
    tt.append(0)

        
for jinit in range(j0,jmax+1):

    for t in range(tpts):
        for ffj in range(len(res.y)):#evaluating final cos2theta
            tempdat= ((res.y[ffj][-1])*mt.e**(complex(0,-1)*(b*(jinit+ffj)*((jinit+ffj)+1))*(time[t])))
            for ffj1 in range(len(res.y)):                
                tempdat1=complex.conjugate((res.y[ffj1][-1])*mt.e**(complex(0,-1)*(b*(jinit+ffj1)*((jinit+ffj1)+1))*(time[t])))
                d[t] += tempdat1*vcos(jinit+ffj,jinit+ffj1,m)*tempdat
#                print(jinit+ffj,jinit+ffj1,m,tempdat,tempdat1,vcos(jinit+ffj,jinit+ffj1,m))
    print(jinit)
        
#for i in range (jmax):
#    (c[i]**2)*vcos(i,i)
    
plt.plot([complex(i).real for i in (d)])
#print(res)
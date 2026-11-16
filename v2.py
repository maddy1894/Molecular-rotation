from scipy.integrate import solve_ivp
from sympy.physics.wigner import wigner_3j
import math as mt
import matplotlib.pyplot as plt 

def fcjm(t,cjm): 
    f=[]
#    m=0
#    j0=0
#    jmax=10
    sum=0.
    for i in range (len(cjm)):
        for j in range (len(cjm)):
            sum += cjm[j] * vcos(jinit+i,jinit+j,m)
#            print(len(cjm),'hhhhhhhhhhhhhhhhhhhhhhh',jinit)
#            print(i,j,cjm[j])
        if (abs(i-j)%2==0):
            f.append(complex(0,1)*(cjm[i]+12*sum))
    return f

def vcos(x,y,m):        
    vcos2theta=((-1)**m)*(2/3)*(mt.sqrt((2*x+1)*(2*y+1)))*wigner_3j(x,2,y,0,0,0)*wigner_3j(x,2,y,-m,0,m)#wigner function
    if (x==y):
        vcos2theta += (1/3)
#        print(m,'This is m')
    return vcos2theta

j0=0
jmax=5
tpulsed=0.1*10**(-15)
b=10.
tfree=2000*10**(-15)
c=[]
tpts=200
time=[(tpulsed*(tpts-(i+1))+(tfree*((i+1)-1)))/(tpts-1) for i in range(tpts+1)]   #Arithematic progression
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

        
for jinit in range(j0,jmax+1):
    j0 = jinit
    print('loop1')
#    for m in range (j0+1)
    m=j0
    print('loop2')
    cjmi=[]
    for cjminit in range (j0,jmax+1):
        print('loop3',m,j0,jmax)
        if (cjminit==j0): 
            cjmi.append(complex(1))
        else:
            cjmi.append(complex(0))
#                print(j0,m,cjmi, len(cjmi),i,cjmi[i])
    res = solve_ivp(fcjm, (0, tpulsed), cjmi)
    for t in range(tpts):
        for ffj in range(len(res.y)):#evaluating final cos2theta
            for ffj1 in range(len(res.y)):
                if (abs(ffj-ffj1)%2==0):
                    tempdat= ((res.y[ffj][-1])*mt.e**(complex(0,-1)*(b*(jinit+ffj)*((jinit+ffj)+1))*(time[t])))
                    tempdat1=complex.conjugate((res.y[ffj1][-1])*mt.e**(complex(0,-1)*(b*(jinit+ffj1)*((jinit+ffj1)+1))*(time[t])))
                    d[t] += tempdat1*vcos(jinit+ffj,jinit+ffj1,m)*tempdat
                    print(jinit+ffj,jinit+ffj1,m)
                    
#for i in range (jmax):
#    (c[i]**2)*vcos(i,i)
    
plt.plot([complex(i).real for i in (d)])
#print(res)
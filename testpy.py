from scipy.integrate import odeint,solve_ivp
import matplotlib.pyplot as plt 
def f(y,tspan):
    return(2*y[0]+3*y[1] , 3*y[0]+2*y[1])
    
#tspan=(0,20)
numpoints = 370
stoptime=20
# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
tspan = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

#for i in range(20001):
#    tspan.append(.001*i)    
#tspan = [ (0.001*i)  for i in range(20001)]                             # set interval of integration
y0=(0,1)
sol = odeint(f,y0,tspan);
print(len(tspan),len(sol),sol)
#plt.plot(sol)
x=[]
y=[]
print(sol[100],sol[100][1])
for i in range(len(sol)):
      x.append(sol[i][0])
      y.append(sol[i][1])
#plt.plot(x,tspan,'ro')
#plt.plot(x)
def f1(tspan,y):
    return(2*y[0]+3*y[1] , 3*y[0]+2*y[1])
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

res = solve_ivp(f1, [0,0.00000000001], y0)
print(res)
print(len(res))
print(res.y)
plt.plot(res.y[1])
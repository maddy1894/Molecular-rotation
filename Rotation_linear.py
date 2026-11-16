# -*- coding: utf-8 -*-
import re
import math as mt
from scipy.integrate import solve_ivp
from sympy.physics.wigner import wigner_3j
import matplotlib.pyplot as plt

class Rotation_linear:

    #Definition and Data entry------
    def __init__(self,filename):
        datafile=open(filename,'r')
        lines=datafile.readlines()
        temp=[]
        for line in lines:
            line=line.strip()
            line=re.split(r'[\s, ;-]+',line)
            temp.append(line[0])
        datafile.close()

        self.Molecule_name=temp[0]
        self.B=float(temp[1])
        self.alpha_par=float(temp[2])
        self.alpha_per=float(temp[3])
        self.spin_even=float(temp[4])
        self.spin_odd=float(temp[5])
        self.temperature=float(temp[6])
        self.J_MAX_user=int(temp[7])
        self.pulse_width=float(temp[8])
        self.eval_time=float(temp[9])
        self.wave_length=float(temp[10])
        self.I0 = float(temp[11])
        del temp
        self.standardize()
    ##definition and data entry------

    ##Standardization of data------
    def standardize(self):
        self.J_MAX = self.J_MAX_user
        self.delta_alpha=(self.alpha_par-self.alpha_per)
        self.part_limit = 0.00001            #partition function limit
        self.partition_function=0.
        self.J_0=0
        self.tpulsed=self.pulse_width*10**(-15)/(2.41888432*10**(-17))
        self.b=self.B/(8065.54*27.2114)                            #cm-1 to hartree
        self.tfree=self.eval_time*10**(-15)/(2.41888432*10**(-17))
        self.tcen=2*self.tpulsed                                   #pulse peak
        self.tend=4*self.tpulsed                                   #field is taken as zero after this
        self.tpts=3000
        tpts=self.tpts
        self.time=[((self.tend-self.tcen)*(tpts-(i+1))+(self.tfree*((i+1)-1)))/(tpts-1) for i in range(tpts)]   #Arithematic progression, t=0 is the pulse peak
        self.E0=mt.sqrt(self.I0*1e12/3.50945e16)                   #I0 in TW/cm2, E0 in a.u.
        self.V_per=-0.25*self.alpha_per*self.E0**2
        self.V_diff=-0.25*self.delta_alpha*(self.E0)**2
        self.a=2*(mt.log(2)/(self.tpulsed)**2)                     #pulse_width is FWHM of intensity
        self.norm=0.
        self.kT = self.temperature/315774.64                       #K to hartree
        self.top_pop=0.
        self.complete_data=[]
        for v in range(tpts):
            self.complete_data.append(0)
    ##standardization of data------

    def spin(self,j):
        if (j%2==0):
            g_i=self.spin_even
        else:
            g_i=self.spin_odd
        return g_i

    ##VCOS2THETA using wigner-3j matrices
    def vcos(self,x,y,m):
        vcos2theta=((-1)**m)*(2/3)*(mt.sqrt((2*x+1)*(2*y+1)))*float(wigner_3j(x,2,y,0,0,0)*wigner_3j(x,2,y,-m,0,m))#wigner function
        if (x==y):
            vcos2theta += (1/3)
        return vcos2theta
    ##vcos2theta using wigner-3j matrices

    #Function ODE Solver
    def fcjm(self,t,cjm):
        f=[]
        b=self.b
        jmin=self.jmin
        vc=self.vc
        env=(mt.e**(-self.a*(t-self.tcen)**2))**2
        for i in range (len(cjm)):
            sum=0.
            for j in range (len(cjm)):
                if (i==j) | (abs(i-j)==1):
                    sum += cjm[j] * vc[i][j]
            f.append(complex(0,-1)*((b*(jmin+2*i)*(jmin+2*i+1)+self.V_per*env)*cjm[i]+self.V_diff*env*sum))
        return f
    #function ode solver

    #ODE solution Module
    def tdse_field(self,j_0,j_max,m):
        b=self.b
        time=self.time
        if ((j_max-j_0)%2!=0):                       #J changes by 2 only, so same parity as j_0
            j_max -=1
        if ((j_0-m)%2==0):                           #lowest J possible for this m
            jmin=m
        else:
            jmin=m+1
        self.jmin=jmin
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
                    vc[i].append(self.vcos(jmin+2*i,jmin+2*j,m))
                else:
                    vc[i].append(0.)
        self.vc=vc
#        print(j_0,j_max,m,cjmi,len(cjmi))
        res = solve_ivp(self.fcjm, [0,self.tend], cjmi, rtol=1e-8, atol=1e-10)
        cfin=[res.y[ffj][-1] for ffj in range(n)]
        avg=0.
        for ffj in range(n):
            avg += abs(cfin[ffj])**2
        if (abs(cfin[-1])**2>self.top_pop):
            self.top_pop=abs(cfin[-1])**2
        d=[]
        for t in range(self.tpts):
            cft=[]
            for ffj in range(n):                     #field free evolution
                cft.append(cfin[ffj]*mt.e**(complex(0,-1)*(b*(jmin+2*ffj)*((jmin+2*ffj)+1))*(time[t]-(self.tend-self.tcen))))
            dt=0.
            for ffj in range(n):#evaluating final cos2theta
                for ffj1 in range(n):
                    if  ((ffj==ffj1) | (abs(ffj-ffj1)==1)) :
                        dt += (complex.conjugate(cft[ffj1])*cft[ffj]*vc[ffj][ffj1])
            d.append(complex(dt).real)
        return(d,avg)
    #ode solution module

    ##Estimating J_MAX and Partition function
    def partition(self):
        j=0
        temp_part=1.                     #any number bigger than part_limit
        while (temp_part>self.part_limit):
            temp_part= mt.e**(-j*(j+1)*self.b/self.kT)*(2*j+1)
            self.partition_function += self.spin(j)*temp_part
            j+=1
            if (j==self.J_MAX_user):
                print("J out of range, select a bigger J_MAX")
                exit()
        self.max_J=j
        #estimating J_MAX and partition function

        #Verifying J_MAX
        print("J_MAX_user,J_MAX,k")
        k=0
        sum_w=0
        max_w = 0.999
        while (sum_w<max_w):
            sum_w += self.spin(k)* mt.e**(-k*(k+1)*self.b/self.kT)*(2*k+1)/self.partition_function
            k+= 1
        print(self.J_MAX_user,self.max_J,k)
        ##verifying J_MAX

    #ODE solution and field free rotation incorporated
    def run(self):
        self.partition()
        for j in range(self.J_0,self.max_J+1):
            g_i=self.spin(j)
            if (g_i==0): continue
            print(j)
            thermal_prob=g_i*mt.e**(-j*(j+1)*self.b/self.kT)/self.partition_function
            for m in range(j+1):         #m varies from 0 to J...
                temp_data,temp_norm=self.tdse_field(j,self.J_MAX,m)
                if m==0:                 #m actually varies from -J to J we'll multiply twice
                    mult=1
                else:
                    mult=2
                for i in range (len(temp_data)):
                    self.complete_data[i] += mult*thermal_prob*temp_data[i]
                self.norm += mult*thermal_prob*temp_norm

        if (self.top_pop>1e-6):
            print("Population reached the last J, select a bigger J_MAX",self.top_pop)

        for i in range (len(self.complete_data)):
            self.complete_data[i]=self.complete_data[i]/self.norm

        self.time_ps=[]
        for i in range (len(self.time)):
            self.time_ps.append(self.time[i] *(2.41888432*10**(-17))*10**(12))      #ps
    #ode solution and field free rotation incorporated

    #Write to a dat file
    def write(self):
        output_file=open(self.Molecule_name+"_at_"+str(self.temperature)+"K_,"+str(self.I0)+"_TWcm2.dat","w")
        for i in range (len(self.complete_data)):
            output_file.writelines([str(self.time_ps[i]),"   ",str(self.complete_data[i]),"\n"])
        output_file.close()
    #write to a dat file

    #Plotter
    def plot(self):
        plt.ylabel(r'$  \langle cos^{2}\theta \rangle $')
        plt.xlabel('Time (ps)')
        plt.title(self.Molecule_name)
        plt.plot(self.time_ps,self.complete_data)
        plt.show()
    #plotter


molecule=Rotation_linear("data.dat")
molecule.run()
molecule.write()
molecule.plot()

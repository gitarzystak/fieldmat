# -*- coding: utf-8 -*-
def layerm(n1,n2,thi,pol):
    tht=arcsin(sin(thi)*n1/n2);
    if pol==1: #s-pol
        r_21 = ((n2*cos(tht))-(n1*cos(thi)))/((n2*cos(tht))+(n1*cos(thi)))
        r_12 = ((n1*cos(thi))-(n2*cos(tht)))/((n2*cos(tht))+(n1*cos(thi)))
        t_21 = (2*(n2*cos(tht)))/((n2*cos(tht))+(n1*cos(thi)))
        t_12 = (2*(n1*cos(thi)))/((n2*cos(tht))+(n1*cos(thi)))        
        r=(n2*cos(thi)-n1*cos(tht))/(n2*cos(thi)+n1*cos(tht));
    else: # p-pol
        r_21 = ((n1*cos(tht))-(n2*cos(thi)))/((n1*cos(tht))+(n2*cos(thi)))
        r_12 = ((n2*cos(thi))-(n1*cos(tht)))/((n2*cos(thi))+(n1*cos(tht)))
        t_21 = (2*(n2*cos(tht)))/((n2*cos(thi))+(n1*cos(tht)))
        t_12 = (2*(n1*cos(thi)))/((n1*cos(tht))+(n2*cos(thi)))
        r=(n1*cos(thi)-n2*cos(tht))/(n2*cos(tht)+n1*cos(thi));
    #m=array([[1, r],[r,1]]);
    #return 1/(1-r)*m
    m=array([[1,-r_12],[r_21,(t_21*t_12-r_21*r_12)]])
    return 1/t_21*m

def prop(l,th,d,n):
    k=2*pi/l;
    fi=k*cos(th)*d*n;
    m=[[exp(1j*fi),0],[0,exp(-1j*fi)]]
    return array(m)

class multil:
    def __init__(self, name):
        self.name = name
        self.layers = []
        self.source = source(400,1100,1,1,0,600)
        self.refl = []
        self.trans=[]
        self.globm = []
        self.fpoints = []
        self.phpoints = []

    def add_layer(self, name,thick,n):
        self.layers.append(layer(name,thick,n))

    def set_source(self,lmin,lmax,step,pol,th,cwl):
        self.source=source(lmin,lmax,step,pol,th,cwl)

    def get_layerinfo(self):
        for i in range(len(self.layers)):
            print('Type: '+self.layers[i].name)
            print('Refractive index: '+str(self.layers[i].n))
            print('Thickness: '+str(self.layers[i].thick))
            

    def calc_tr(self):
        source=self.source
        layers=self.layers
        layers.reverse()
        print('not ready yet')
        N=len(source.waves)
        R=[0]*N
        T=[0]*N
        n=0
        for l in source.waves: #pętla po falach
            m=eye(2)
            smax=len(layers); #liczba warstw
            for s in arange(1,smax-1,1):
		#obliczanie n1 i n2
                if not len(self.layers[s-1].ndata):
                    n1=layers[s-1].n
                else:
                    n1=interp(l,layers[s-1].wl,layers[s-1].ndata)+1j*interp(l,layers[s-1].wl,layers[s-1].kdata)
                if not len(layers[s].ndata):
                    n2=layers[s].n
                else:
                    n2=interp(l,layers[s].wl,layers[s].ndata)+1j*interp(l,layers[s].wl,layers[s].kdata)
		#obliczanie macierzy
                d2=layers[s].thick
                gr1=layerm(n1,n2,0,1);
                pr1=prop(l,0,d2,n2);
                m=pr1.dot(gr1).dot(m)
	    #ostatni etap
            n1=n2
            n2=layers[smax-1].n #ostatnia nie może być nk
            grn=layerm(n1,n2,0,1);
            #pr1=prop(l,0,d2,n2);
            m=grn.dot( m)
            t=1/m[0,0];
            r=m[1,0];
            R[n]=abs(r*t)**2; #może nie działać dla padania innego niż prostopadłe
            T[n]=abs(t)**2
            n=n+1
        self.refl=R
        self.trans=T
        #print(l)
        self.globm=m

    def calc_field(self):
        if not len(self.globm):
            self.calc_tr()
        t=1/self.globm[0,0]
        r=self.globm[1,0]*t
        #v=[1/(1+r),r/(1+r)]
        v=[t,0]
	#print(dot(self.globm,v))

        layers=self.layers
#	layers.reverse() calc_tr odwraca, nie ma po co drugi raz
        smax=len(layers); #liczba warstw

        field=zeros([smax,2],dtype="complex")
        field[0,:]=v
        fi=[0]*(smax)
        l=self.source.cwl
        th=self.source.thi

        m=eye(2)
	
        for s in arange(1,smax-1,1):
            if not len(layers[s-1].ndata):
                n1=layers[s-1].n
            else:
                n1=interp(l,layers[s-1].wl,layers[s-1].ndata)+1j*interp(l,layers[s-1].wl,layers[s-1].kdata)
            if not len(layers[s].ndata):
                n2=layers[s].n
            else:
                n2=interp(l,layers[s].wl,layers[s].ndata)+1j*interp(l,layers[s].wl,layers[s].kdata)
            d2=layers[s].thick
            #print(layers[s-1].name)
            gr1=layerm(n1,n2,0,1);
            pr1=prop(l,0,d2,n2);
            k=2*pi/l;
            th=arcsin(sin(th)*n1/n2);
            fi[s]=k*cos(th)*n2;
            m=dot(gr1,m) #odtąd trzeba kombinować
            v1=dot(m,v)
            m=dot(pr1,m)
            field[s,:]=v1
            #field[s]=abs(v1[0])**2+abs(v1[1])**2
	    #print(v1)
        n1=n2
        n2=layers[smax-1].n #ostatnia nie może być nk
        th=arcsin(sin(th)*n1/n2);
        grn=layerm(n1,n2,0,1);
	#pr1=prop(l,0,d2,n1);
        m=grn.dot(m)
        v1=dot(m,v)
	#print(v1)
        field[smax-1,:]=v1
        fi[0]=2*pi/l*cos(self.source.thi)*layers[0].n #Starting point for phase
        fi[smax-1]=2*pi/l*cos(th)*n2
	#field.reverse() #get 
        #field[smax-1]=abs(v1[0])**2+abs(v1[1])**2
        self.fpoints=field
	#print(array(m)-array(self.globm))
	#print(abs(v1[1]/v1[0])**2)
        self.phpoints=fi
        print(self.fpoints)
        print(self.phpoints)

    def plot_r(self):
        plt.plot(10**9*self.source.waves,self.refl)
        plt.show()

    def plot_field(self):
        thv=[0]*(len(self.layers))
        d=[0]*(len(self.layers))
        thv[0]=self.layers[0].thick

        x0=0
        for i in range(len(self.layers)):
            x=linspace(0,self.layers[i].thick,100)
            if i==0:
                fi=exp(1j*self.phpoints[i]*(x-self.layers[i].thick))
            else:
                fi=exp(1j*self.phpoints[i]*x)
            field=self.fpoints[i,0]*fi+self.fpoints[i,1]/fi
            plt.plot((x+x0)*10**9,abs(field))
            x0=x0+self.layers[i].thick
	    
        plt.show()

    
class layer:
    def __init__(self, name,thick,n):
        self.name = name
        self.thick=thick
        self.n=n
        self.filename='sodank.txt'
        self.wl=[]
        self.ndata=[]
        self.kdata=[]

    def get_nk(self):
        if self.name!='nk':
            print('no nk data available for this layer')
        else:
            n=[]
            k=[]
            wl=[]
            # Open file
            f = open(self.filename, 'r')

            # Read and ignore header lines
            header1 = f.readline()
            header2 = f.readline()
            header3 = f.readline()
            # Read file
            for line in f:
                line = line.strip()
                columns = line.split()
                wl.append(float(columns[0]))
                n.append(float(columns[1]))
                k.append(float(columns[2]))
            f.close()
        self.wl=array(wl)*10**-9
        self.ndata=array(n)
        self.kdata=array(k)

    def set_filename(self,filename):
        self.filename=filename
        self.get_nk()

class source:
    def __init__(self, lmin, lmax, step,pol,thi,cwl):
        self.lmin = 10**-9*lmin
        self.lmax = 10**-9*lmax
        self.step = 10**-9*step
        self.pol=pol
        self.thi=thi
        self.waves=arange(self.lmin, self.lmax+self.step, self.step)
        self.cwl=cwl*10**-9

from numpy import cos, arcsin, sin, linspace, exp, dot, eye, zeros, pi, array, arange, interp, real
import matplotlib.pyplot as plt

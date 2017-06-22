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
        self.source = source(400,1100,1,1,0)
        self.refl = []
        self.trans=[]
        self.globm = []
        self.cwl=6*10**-7

    def add_layer(self, name,thick,n):
        self.layers.append(layer(name,thick,n))

    def set_source(self,lmin,lmax,step,pol,th):
        self.source=source(lmin,lmax,step,pol,th)

    def get_layerinfo(self):
        for i in range(len(self.layers)):
            print('Type: '+self.layers[i].name)
            print('Refractive index: '+str(self.layers[i].n))
            print('Thickness: '+str(self.layers[i].thick))
            

    def calc_tr(self):
        source=self.source
        layers=self.layers
        print('not ready yet')
        N=len(source.waves)
        R=[0]*N
        T=[0]*N
        n=0
        for l in source.waves:
            '''
            np=self.layers[0].n
            n1=self.layers[1].n
            n2=self.layers[2].n
            d1=self.layers[1].thick
            d2=self.layers[2].thick
            gr0=layerm(np,n1,0,1);
            grn=layerm(n1,np,0,1);
            gr1=layerm(n1,n2,0,1);
            pr2=prop(l,0,d2,n2);
            pr1=prop(l,0,d1,n1);
            gr2=layerm(n2,n1,0,1);
            '''
            m=eye(2)
            smax=len(layers); #liczba warstw
            for s in arange(1,smax-1,1):
                if not len(self.layers[s-1].ndata):
                    n1=self.layers[s-1].n
                else:
                    n1=interp(l,self.layers[s-1].wl,self.layers[s-1].ndata)+1j*interp(l,self.layers[s-1].wl,self.layers[s-1].kdata)
                if not len(self.layers[s].ndata):
                    n2=self.layers[s].n
                else:
                    n2=interp(l,self.layers[s].wl,self.layers[s].ndata)+1j*interp(l,self.layers[s].wl,self.layers[s].kdata)
                d2=self.layers[s].thick
                gr1=layerm(n1,n2,0,1);
                pr1=prop(l,0,d2,n2);
                m=pr1.dot(gr1).dot(m)
            n1=n2
            n2=self.layers[smax-1].n #ostatnia nie może być nk
            grn=layerm(n1,n2,0,1);
            m=grn.dot(pr1).dot( m)
            t=1/m[0,0];
            r=m[1,0];
            R[n]=abs(r*t)**2; #może nie działać dla padania innego niż prostopadłe
            T[n]=abs(t)**2
            n=n+1
        self.refl=R
        self.trans=T
        self.globm=m

    def calc_field(self):
        if len(self.globm):
            self.calc_tr()
        t=1/self.globm[0,0]
        r=self.globm[1,0]*t
        v=[1 r]
        m=eye(2)
        smax=len(layers); #liczba warstw
        field=[0]*(smax+1)
        for s in arange(1,smax-1,1):
            if not len(self.layers[s-1].ndata):
                n1=self.layers[s-1].n
            else:
                n1=interp(l,self.layers[s-1].wl,self.layers[s-1].ndata)+1j*interp(l,self.layers[s-1].wl,self.layers[s-1].kdata)
            if not len(self.layers[s].ndata):
                n2=self.layers[s].n
            else:
                n2=interp(l,self.layers[s].wl,self.layers[s].ndata)+1j*interp(l,self.layers[s].wl,self.layers[s].kdata)
            d2=self.layers[s].thick
            gr1=layerm(n1,n2,0,1);
            pr1=prop(l,0,d2,n2);
            m=pr1.dot(gr1).dot(m) #odtąd trzeba kombinować
            v1=dot(m,v)
            field[s]=v1[0]+v1[1]
        n1=n2
        n2=self.layers[smax-1].n #ostatnia nie może być nk
        grn=layerm(n1,n2,0,1);
        m=grn.dot(pr1).dot( m)
        v1=dot(m,v)
        field[smax]=v1[0]+v1[1]
        
    def plot_r(self):
        plt.plot(10**9*self.source.waves,self.refl)
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
    def __init__(self, lmin, lmax, step,pol,thi):
        self.lmin = 10**-9*lmin
        self.lmax = 10**-9*lmax
        self.step = 10**-9*step
        self.pol=pol
        self.thi=thi
        self.waves=arange(self.lmin, self.lmax+self.step, self.step)

from numpy import cos, arcsin, sin, linspace, exp, dot, eye, zeros, pi, array, arange, interp
import matplotlib.pyplot as plt

#initialize
d=multil('Figo')

#configuration
l=6e-7;
d1=0.25*l/1.5;
d2=0.25*l/1;
d.add_layer('const',d2,1)
d.add_layer('nk',d1,1)
d.layers[1].set_filename('sodank.txt')
d.add_layer('const',d2,1)
d.add_layer('const',d1,1.5)
d.add_layer('const',d2,1)
d.add_layer('const',d1,1.5)
d.add_layer('const',d2,1)
d.add_layer('const',d1,1.5)
d.add_layer('const',d2,1)
d.add_layer('const',d1,1.5)
d.add_layer('const',d2,1)
d.add_layer('const',d1,1.5)
d.add_layer('const',d2,1)
d.add_layer('const',d1,1.5)
d.add_layer('const',d2,1)


#setup & run simulation
d.get_layerinfo()
d.set_source(400,800,0.5,1,0)
d.calc_tr()
d.plot_r()

#initialize
from fieldmat import  *
d=multil('Figo') # give your multilayer a sweet name

#configuration
l=6e-7;
d1=0.25*l/1.5;
d2=0.25*l/1;
d.add_layer('const',d1,1)
d.add_layer('const',d1,1)
d.add_layer('nk',d2,1)
d.add_layer('nk2',d1,1.5)
d.add_layer('nk3',d1,1)
d.add_layer('nk3',d1,1)

d.get_layerinfo()
d.set_source(400,800,0.5,1,0,800.5)
d.calc_tr()
d.calc_field()
d.plot_r()
d.plot_field()

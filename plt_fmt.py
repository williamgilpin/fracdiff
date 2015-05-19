# Plot formatting

from numpy import *
from scipy import *
from matplotlib.pyplot import *

rcParams['font.family'] = 'Myriad Pro'
rcParams['font.weight'] = 'bold'
# # These are the "Tableau 20" colors as RGB.  
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
#tableau20 = [(0.59375, 0.257813, 0.886719),(0.800781, 0.0625, 0.460938), (0.278431, 0.788235, 0.478431), (0.917647, 0.682353, 0.105882), (0.372549, 0.596078, 1.), (0.8, 0.8, .8), (1.0, .3882, .2784)];

tableaugray = [(123,102,210),(207,207,207),(165,172,175),(143,135,130),(96,99,106),(65,68,81)]

tab_purp = [(255,192,218),(216,152,186),(180,177,155),(153,86,136),(139,124,110),\
         (220,95,189),(95,90,65),(171,106,213),(219,212,197)]

tab_og = [(255,127,15), (172,217,141),(60,183,204),(50,162,81),(255,217,74),(184,90,13),\
(152,217,228),(134,180,169),(57,115,124),(204,201,77),(130,133,59)]

tab_br = [(255,182,176),(240,39,32),(181,200,226),(172,97,60),(44,105,176),(233,195,155),\
(221,201,180),(181,223,253),(172,135,99),(107,163,214),(244,115,122),(189,10,54)]

tab_extra = [(208,152,238), (166,153,232), (219,212,197)]

all_colors = tableau20 + tableaugray + tab_purp + tab_og + tab_br + tab_extra

for i in range(len(all_colors)):  
    r, g, b = all_colors[i]  
    all_colors[i] = (r, g, b)
    all_colors[i] = (r / 255., g / 255., b / 255.)  

bcolor1 = (53 / 255., 56 / 255., 57 / 255.)
bcolor2 = (59 / 255., 60 / 255., 54 / 255.)

rcParams['axes.color_cycle'] = all_colors
rcParams['axes.labelcolor'] = bcolor1
rcParams['axes.edgecolor'] = bcolor2
rcParams['xtick.color'] = bcolor2
rcParams['ytick.color'] = bcolor2

rcParams['legend.frameon'] = False

rcParams['axes.facecolor'] = 'none'

rcParams['lines.linewidth'] = 2

rcParams['figure.facecolor'] = 'none'
rcParams['savefig.facecolor'] = 'none'

rcParams['image.cmap'] = 'gray'
rcParams['image.interpolation'] = 'nearest'

rcParams['figure.figsize'] = (6, 6)

# This slows everything down because of equation rendering. It also breaks fonts
# rcParams['text.usetex'] = True

#rcParams['figure.savefig.bbox']= 'none'

# for ii in xrange(50):
#     plot([rand(), rand()],[rand(), rand()])
#savefig("out.pdf")
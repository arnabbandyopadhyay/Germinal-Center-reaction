#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 13:27:53 2019

@author: abp19
"""
import numpy as np
import random

from matplotlib import pyplot as plt
from matplotlib import animation

mu, sigma = 7.1, 3 # mean and standard deviation
box_width_dz=500
box_width_lz=500
class Bcell:
    """ B cell objects are the main unit of this simulation. They consist of
    the following information:
        - their key sequence (nkey codons given as number code)
        - their germline ancestors key sequence
        - their current affinity to the antigen sequence
        - their ancestors affinity to the antigen sequence
        - their origin: recruited into this process from naive or pre-existing
            memory pool?
        - their overall mutation count
        - their family ID (for grouping cells into clones)
        - their time of birth (for turnover of naive pool)
        - the most recent time they entered a GC (relevant for discarding
            cells from GC that have competed for survival signals without
            success for too long)
        - the time this cell or, if it was born inside the GC, its ancestor,
            entered the GC - relevant for knowing whether AID is aready acting
            and mutations can happen, as the enzyme first needs to be
            upregulated upon GC entry
        - whether or not affinity maturation (affinity improvement as compared
            to the germline ancestor) is blocked by detrimental mutations in
            the non-key region
    """
    def __init__(self, birthtime, **kwargs):
        self.birthtime=birthtime
         
        self.nextdiv=self.birthtime+np.random.randint(7, 14)#np.random.normal(mu, sigma)
        
        if 'generation' in kwargs: self.generation = kwargs['generation']
        else:                  self.generation = 0
        
        if 'selected' in kwargs: self.selected = kwargs['selected']
        else:                  self.selected = 0
        
        self.tselected=0
        
        if 'LZentrytime' in kwargs: self.LZentrytime = kwargs['LZentrytime']
        else:                  self.LZentrytime = 0
        
        self.xpos=round(np.random.random(),2)*box_width_dz
        self.ypos=round(np.random.random(),2)*box_width_dz
        
        self.xvel=round(2*(np.random.random()-0.5),2)*box_width_dz
        self.yvel=round(2*(np.random.random()-0.5),2)*box_width_dz
        
    @classmethod    
    def clone(cls, b):
        return cls(birthtime=tnow, generation = b.generation, selected=b.selected, LZentrytime=b.LZentrytime, nextdiv=tnow+np.random.randint(7, 14))
    def gen(self):
        self.generation += 1

    def proliferate(self):


        nl=[Bcell.clone(self),Bcell.clone(self)]

        return nl
#        self.nextdiv=self.birthtime+7
#                 birthtime, mutations=0, family=None, GCentrytime=None,
#                 AIDstart=None, block=False):
#    uID = 0  # unique ID counter
#    ufamID = 0  # unique family ID counter
#
#    def __init__(self, sequence, sequence0, affinity, affinity0, origin,
#                 birthtime, mutations=0, family=None, GCentrytime=None,
#                 AIDstart=None, block=False):
#        self.ID = self.uID
#        Bcell.uID += 1
#
#        if family is None:  # new family ID only if new family created
#            self.family = self.ufamID
#            Bcell.ufamID += 1
#        else:
#            self.family = family
#
#        self.sequence = sequence  # current sequence
#        self.sequence0 = sequence0  # sequence of the original ancestor
#
#        self.affinity = affinity  # highest affinity
#        self.affinity0 = affinity0  # affinity of original ancestor
#        self.origin = origin  # was first ancestor naive or unspecific memory?
#        self.mutations = mutations  # counter of bp mutations
#
#        self.birthtime = birthtime  # time the family was produced from the BM
#        self.GCentrytime = GCentrytime  # time the cell entered the waitlist
#        self.AIDstart = AIDstart  # time since ancestor entered GC site
#        self.block = block


def proliferate(lists):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    nl=[]
    for i in range(len(lists)):
        nc1=Bcell(birthtime=tnow)
        nl.append(nc1)
        nc2=Bcell(birthtime=tnow)
        nl.append(nc2)
    
    
    return nl

def proliferate2(b):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    b.gen()
    nl=[b.clone(b),b.clone(b)]

    return nl

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(1, 1000), ylim=(1, 500))
plt.axvline(x=500,linewidth=4,color='r')

tfh_x_coord=np.random.randint(500,1000,200)
tfh_y_coord=np.random.randint(0,500,200)
plt.plot(tfh_x_coord,tfh_y_coord,'ro',ms=1)

fdc_x_coord=np.random.randint(500,1000,200)
fdc_y_coord=np.random.randint(0,500,200)
antg_x_coord=[round(np.random.normal(k, 5),0) for k in fdc_x_coord for b in range(6)]
antg_y_coord=[round(np.random.normal(k, 5),0) for k in fdc_y_coord for b in range(6)]
plt.plot(antg_x_coord,antg_y_coord,'bo',ms=1)


particles, = ax.plot([], [], 'bo', ms=1)
#line, = ax.plot([], [], lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.7, 0.05, '', transform=ax.transAxes)


Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)

# initialization function: plot the background of each frame
def init():
    particles.set_data([], [])
    time_text.set_text('')
    return particles, time_text


#lists=[Bcell(birthtime=0),Bcell(birthtime=0)]
#for tnow in range(10):
#    birthtime=tnow
#    print(tnow,len(lists))
#    nlists=lists
#    for i in range(len(lists)):
#        print(len(nlists))
#        b=lists[i]
#        if b.nextdiv<tnow:
#            nlists.append(proliferate(lists))
#        else:
#            nlists.append(Bcell(birthtime))
#            nlists.append(Bcell(birthtime))
#        
#        
#        
#bt=[b.nextdiv for b in lists]
dt=1
generation=0
lists=[Bcell(birthtime=0),Bcell(birthtime=0)]
#lists=[Bcell(birthtime=0) for b in range(144)]
pltlist=[lists[0::]]
LZ=[]
r=[]
#trange=np.arange(1, 250, 0.1).tolist()
for tnow in range(1,250):
    
    birthtime=tnow
#    print(tnow,len(lists),[b.generation for b in lists])
#    print(tnow,len(lists),[[b.birthtime,b.nextdiv] for b in lists])
#    dlist=[b.nextdiv for b in lists]
#    dlist2=[i for i,x in enumerate(dlist) if x <tnow]
    nl2=[]
    rc=[]

#    nlists=lists
    for i in range(len(lists)):
        nl=[]
#        print(i)
#        print(len(lists))
        b=lists[i]
        x1=b.xvel*dt
        y1=b.yvel*dt
        b.xpos += b.xvel*dt
        b.ypos += b.yvel*dt
        
        if abs(b.xpos) > box_width_dz:
            b.xpos = -b.xvel
            
        if abs(b.ypos) > box_width_dz:
            b.ypos = -b.yvel
        
        if b.selected==0:
            if (b.nextdiv<=tnow) and (b.generation<6):
                rc.append(b)
#            print('y')
                generation=b.generation+1
                nl=proliferate2(b)
        elif b.selected==1:
            if (b.nextdiv<=tnow) and (b.generation<13):
                rc.append(b)
#            print('y')
                generation=b.generation+1
                nl=proliferate2(b)
        nl2+=nl
        
        
    lists+=nl2

    
    [lists.remove(x) for x in rc]
            
#    if len(nl)>0:
#        lists.append(nl)
    
 
### 2 new cells comes in per hour, but this version starts with 144 cells so no cells comes in 

#    lists.append(Bcell(birthtime))
#    lists.append(Bcell(birthtime))
    
#    pltlist2.append(lists)
    pltlist.append(lists[0::])
#    print(lists)
#    print(pltlist)
    
    gen4=[b for b in lists if b.generation>=6 and b.selected==0]
    
    LZ_cells=random.choices(gen4, k=np.random.binomial(len(gen4),0.5))
    if len(LZ_cells)>0:
        LZ_cells=set(LZ_cells) # making unique set
#        print(tnow,len(lists),len(LZ_cells))
        [lists.remove(x) for x in LZ_cells]
        
        for x in LZ_cells:
            x.xpos += 500+x.xvel*dt
            x.ypos += 500+x.yvel*dt
            
        if abs(x.xpos) > box_width_lz+box_width_dz:
            x.xpos = -x.xvel
            
        if abs(x.ypos) > box_width_lz+box_width_dz:
            x.ypos = -x.yvel
            
        for b in LZ_cells:
            b.LZentrytime=tnow
#        
    lzcells=[b for b in LZ_cells]
    
    gen5=[b for b in lists if b.generation>6 and b.selected ==1]
#    [lists.remove(x) for x in gen5]
    
    LZ_selected_cells=random.choices(gen5, k=np.random.binomial(len(gen5),0.5))
    if len(LZ_selected_cells)>0:
        LZ_selected_cells=set(LZ_selected_cells) # making unique set
#        print(tnow,len(lists),len(LZ_cells))
        [lists.remove(x) for x in LZ_selected_cells]
        
        for x in LZ_selected_cells:
            x.xpos += 500+x.xvel*dt
            x.ypos += 500+x.yvel*dt
            
        if abs(x.xpos) > box_width_lz+box_width_dz:
            x.xpos = -x.xvel
            
        if abs(x.ypos) > box_width_lz+box_width_dz:
            x.ypos = -x.yvel
            
        for b in LZ_selected_cells:
            b.LZentrytime=tnow
##        
#    lz_selected_cells=[b for b in LZ_selected_cells] ## didn't update any list assuming leaving from DZ
    
    
#    print([b.LZentrytime for b in lzcells])
    appcells=[b for b in LZ if tnow > b.LZentrytime+6]
    ap_cells=random.choices(appcells, k=np.random.binomial(len(appcells),0.5))
    ap_cells=set(ap_cells)
    [LZ.remove(x) for x in ap_cells]
    LZ+=lzcells
    
    
### those cells need to be selected which are less than 6 hr as I implemented appcells    
    selected_cells=random.choices(LZ, k=np.random.binomial(len(LZ),0.5))
    selected_cells=set(selected_cells)
    for b in selected_cells:
        b.selected=1
    selected_cells=list(selected_cells)
    
    # 10 percent cell leaves GC as early plasma cell
    
    ps_cells=random.choices(selected_cells, k=np.random.binomial(len(selected_cells),0.2))
    ps_cells=set(ps_cells)
    
    [LZ.remove(x) for x in selected_cells]
    
    [selected_cells.remove(x) for x in ps_cells]
    
    lists+=selected_cells
        
#    print(tnow,len(gen4))
    print([tnow,len(lists)+len(LZ_cells)])
#    print(len(lists),len(LZ))
    
    
    r.append([tnow,len(lists),len(LZ)])
#






    
xcoord=[[b.xpos for b in pltlist[item]] for item in range(len(pltlist))]  
ycoord=[[b.ypos for b in pltlist[item]] for item in range(len(pltlist))]   




      
def animate(i):
    particles.set_data(xcoord[i], ycoord[i])
    time_text.set_text(time_template % (i))
    return particles, time_text

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=1, blit=True, repeat=False)

anim.save('gc.mp4',writer=writer)

plt.show()

r=np.array(r)


fig=plt.figure()
plt.subplot(2,1,1)
np.seterr(divide='ignore', invalid='ignore')
plt.plot(r[:,0], r[:,1]/r[:,2])
plt.plot([0,300],[2,2],[0,300],[3,3],'r')
plt.axis([0,300,0,10])
plt.ylabel('DZ/LZ ratio')

plt.subplot(2,1,2)
#labels=(,'LZ Bcells','All Bcells')
p1=plt.plot(r[:,0], r[:,1], label='DZ Bcells')
p2=plt.plot(r[:,0],r[:,2], label='LZ Bcells')
p3=plt.plot(r[:,0],r[:,1]+r[:,2], label='All Bcells')
plt.legend()
plt.show()

fig.savefig('gc.pdf')
#plt.plot(r[:,0],r[:,1],r[:,2])
#plt.plot(r[:,0],r[:,1]/r[:,2])


    
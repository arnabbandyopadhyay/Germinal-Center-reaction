#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 13:27:53 2019

@author: abp19

to optimize parameters 2 cells/hr entry
"""
import numpy as np
import random
import math
import time
from matplotlib import pyplot as plt
from matplotlib import animation
#from pandas import DataFrame
import csv


start = time.time()

mu, sigma = 7.1, 3 # mean and standard deviation
mu2, sigma2 = 7, 12
box_width_dz=500
box_width_lz=500
opseq=[3,3,3,3]
mum, sigm = 5, 1 # hr from LZ to DZ and vice versa
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
         
        self.nextdiv=self.birthtime+np.random.randint(mu2, sigma2)#np.random.normal(mu, sigma)
        
        if 'generation' in kwargs: self.generation = kwargs['generation']
        else:                  self.generation = 0
        
        if 'selected' in kwargs: self.selected = kwargs['selected']
        else:                  self.selected = 0
        
        if 'tfhcontact' in kwargs: self.tfhcontact = kwargs['tfhcontact']
        else:                  self.tfhcontact = 0
        
        
        if 'LZentrytime' in kwargs: self.LZentrytime = kwargs['LZentrytime']
        else:                  self.LZentrytime = 0
        
        if 'sequence' in kwargs: self.sequence = kwargs['sequence']
        else: self.sequence = np.random.randint(0,10,4).tolist()
        
        if 'affinity' in kwargs: self.affinity = kwargs['affinity']
        else: self.affinity = math.exp(-hamdist(opseq,self.sequence)**2/2.8**2)
        
        if 'antigen_collect' in kwargs: self.antigen_collect = kwargs['antigen_collect']
        else: self.antigen_collect = 0
        
        if 'fdc_selected' in kwargs: self.fdc_selected = kwargs['fdc_selected']
        else: self.fdc_selected = 0
        
        if 'tmigration' in kwargs: self.tmigration = kwargs['tmigration']
        else: self.tmigration = tnow + np.random.normal(mum,sigm)
        
        if 'antg_clt' in kwargs: self.antg_clt = kwargs['antg_clt']
        else: self.antg_clt = 0
        
        if 'tfhsignal' in kwargs: self.tfhsignal = kwargs['tfhsignal']
        else: self.tfhsignal = 0
        
        if 'tfh_sig_t' in kwargs: self.tfh_sig_t = kwargs['tfh_sig_t']
        else: self.tfh_sig_t = 0
        
        if 'tfh_contact_start' in kwargs: self.tfh_contact_start = kwargs['tfh_contact_start']
        else: self.tfh_contact_start = 0
        
        if 'tfh_contact_end' in kwargs: self.tfh_contact_end = kwargs['tfh_contact_end']
        else: self.tfh_contact_end = 0
        
        if 'deleting_time' in kwargs: self.deleting_time = kwargs['deleting_time']
        else: self.deleting_time = 600
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = round(np.random.random(),2)*box_width_dz
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos = round(np.random.random(),2)*box_width_dz
        
        
        
        
#        self.xpos=round(np.random.random(),2)*box_width_dz
#        self.ypos=round(np.random.random(),2)*box_width_dz
        
        self.xvel=round(2*(np.random.random()-0.5),2)*box_width_dz
        self.yvel=round(2*(np.random.random()-0.5),2)*box_width_dz
        
    @classmethod    
    def clone(cls, b):
        return cls(birthtime=tnow, generation = b.generation, selected=b.selected, tfhcontact = b.tfhcontact, 
                   LZentrytime=b.LZentrytime, sequence=b.sequence, affinity=b.affinity, antigen_collect=b.antigen_collect, 
                   fdc_selected=b.fdc_selected, tmigration = tnow + np.random.normal(mum,sigm), 
                   antg_clt = b.antg_clt, tfhsignal = b.tfhsignal, tfh_sig_t = b.tfh_sig_t,
                   tfh_contact_start = b.tfh_contact_start, tfh_contact_end = b.tfh_contact_end, 
                   deleting_time = b.deleting_time,
                   nextdiv=tnow+np.random.randint(mu2, sigma2), xpos = round(np.random.random(),2)*box_width_dz,
                   ypos = round(np.random.random(),2)*box_width_dz)
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

def hamdist(str1, str2):
    
   """Count the # of differences between equal length strings str1 and str2"""
        
   diffs = 0
   for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
   return diffs

def mutate(seq):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    nl=[]
    for i in range(len(seq)):
        seq2=seq[:]
        seq2[i]=seq[i]+1
        nl.append(seq2)
        seq2=seq[:]
        seq2[i]=seq[i]-1
        nl.append(seq2)
    nxtpos1=random.choice(nl)
    nl.pop(nl.index(nxtpos1))
    nxtpos2=random.choice(nl)
    affinity1=math.exp(-hamdist(opseq,nxtpos1)**2/2.8**2)
    affinity2=math.exp(-hamdist(opseq,nxtpos2)**2/2.8**2)
    return nxtpos1,nxtpos2, affinity1, affinity2



def antigen_collect(lists):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    bc=[b for b in lists if tnow < b.antg_clt]
    for i in range(len(bc)):
        b=bc[i]
        dist=[math.sqrt((b.xpos-i[0])**2 + (b.ypos-i[1])**2) for i in zip(fdc_x_coord,fdc_y_coord)]
        min_dist=min(dist)
        if min_dist < 10:
            antigen_collect=np.random.binomial(1,b.affinity)
            if antigen_collect ==1:
                b.antigen_collect +=1
                b.fdc_selected = 1
        
    return lists


def antigen_process(lists):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    bc=[b for b in lists if tnow > b.antg_clt]
    for i in range(len(bc)):
        b=bc[i]
        dist=[math.sqrt((b.xpos-i[0])**2 + (b.ypos-i[1])**2) for i in zip(fdc_x_coord,fdc_y_coord)]
        min_dist=min(dist)
        if min_dist < 10:
            antigen_collect=np.random.binomial(1,b.affinity)
            if antigen_collect ==1:
                b.antigen_collect +=1
                b.fdc_selected = 1
        
    return lists




def contain_points(check_x, check_y, center_x, center_y, radius):
    if (check_x - center_x)**2 + (check_y - center_y)**2 < radius^2:
        return('True')
        
    else:
        return('False')

def t_help(lists):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    for i in range(n_tfh):
        nl=[]
        for j in range(len(lists)):
            b=lists[j]
            if b.tfhcontact==1 and b.tfh_sig_t < tnow:
                b.tfhcontact=0
                b.tfh_sig_t=0
            elif b.tfhcontact==1 and b.tfh_sig_t > tnow:
                b.tfhcontact=1
                b.tfh_sig_t=b.tfh_sig_t
                nl.append(b)
                    
            elif b.tfhcontact==0 and contain_points(b.xpos,b.ypos,tfh_x_coord[i], tfh_y_coord[i],radius)=='True':
                b.tfhcontact=1
                b.tfh_sig_t=tnow + 0.5
                
                nl.append(b)
                
        k=[b.antigen_collect for b in nl]
        if len(k)>0:
            
            pos=k.index(max(k))
            nl[pos].tfhsignal += timestep  # find bcell with max antigen, find position in list index and change tfhsignal to +1
        
        
        
    return lists



def proliferate2(b):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    b.gen()
    if tnow < 72:
        
        nl=[b.clone(b),b.clone(b)]
    else:
        k=b.clone(b)
        mu_prob=abs(0.5-(0.5-0)**k.affinity)
        chance_mutation=np.random.binomial(1,mu_prob)
        if chance_mutation == 1:
            nxtpos1, nxtpos2, aff1, aff2 = mutate(k.sequence) 
            k.sequence=nxtpos1
            k.affinity=aff1
            nl1=k.clone(k)
            k.sequence=nxtpos2
            k.affinity=aff2
            nl2=k.clone(k)
            nl=[nl1, nl2]
        else:
            nl=[b.clone(b),b.clone(b)]

    return nl

def LZposition(lists):
    
    for x in lists:
        if x.tfhcontact==1 and x.tfh_sig_t > tnow:
            x.xpos=x.xpos
            x.ypos=x.ypos
        else:
            x.xpos = 500 + round(np.random.random(),2)*box_width_lz
            x.ypos = round(np.random.random(),2)*box_width_lz
            
            if abs(x.xpos) > box_width_lz+box_width_dz or abs(x.xpos) < box_width_lz:
                x.xpos = round(np.random.random(),2)+box_width_lz
            
            if abs(x.ypos) > box_width_lz:
                x.ypos = round(np.random.random(),2)*box_width_lz
        
    return lists
            
    

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(1, 1000), ylim=(1, 500))
plt.axvline(x=500,linewidth=4,color='r')

n_tfh=200
tfh_x_coord=np.random.randint(500,1000,n_tfh)
tfh_y_coord=np.random.randint(0,500,n_tfh)
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
radius=10
dt=1
generation=0
tnow=0
timestep=1
lists=[Bcell(birthtime=0),Bcell(birthtime=0)]
#lists=[Bcell(birthtime=0) for b in range(144)]
pltlist=[lists[0::]]
LZ=[]
pot_bc=[]
outputcells=[]
r=[]
bb=[]
datatoexport=[]

#trange=np.arange(1, 250, 0.1).tolist()
for tnow in np.arange(1,200,timestep):
    print(tnow, len(lists),len(LZ))

    bbc=[]
    
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
#        x1=b.xvel*dt
#        y1=b.yvel*dt
#        b.xpos += b.xvel*dt
#        b.ypos += b.yvel*dt
#        
#        if abs(b.xpos) > box_width_dz:
#            b.xpos = -b.xvel
#            
#        if abs(b.ypos) > box_width_dz:
#            b.ypos = -b.yvel
        
        if b.selected==0:
            if (b.nextdiv<=tnow) and (b.generation<6):
                rc.append(b)
#            print('y')
                generation=b.generation+1
                nl=proliferate2(b)
        elif b.selected==1:
            if (b.nextdiv<=tnow) and (b.generation<20):
                rc.append(b)
#            print('y')
                generation=b.generation+1
                nl=proliferate2(b)
        nl2+=nl
        
        
    lists+=nl2

    
    [lists.remove(x) for x in rc]
            
#    if len(nl)>0:
#        lists.append(nl)
    if tnow > 72:
        if (round(tnow,2) % 3) == 0:
            lists.append(Bcell(birthtime))  ### 1 cell per 3 hour entry in GC
            print('enter 1 cell per 3 hr')
#    lists.append(Bcell(birthtime))
    else:
        if (round(tnow,2) % 1) == 0:
            lists.append(Bcell(birthtime))
            lists.append(Bcell(birthtime))
            print('enter 2 cells')
#        lists.append(Bcell(birthtime))

#    pltlist2.append(lists)
    pltlist.append(lists[0::])
#    print(lists)
#    print(pltlist)
    
    gen4=[b for b in lists if b.generation>=6 and b.selected == 0]
    
#    LZ_cells=random.choices(gen4, k=np.random.binomial(len(gen4),0.05))
    LZ_cells = [b for b in gen4 if tnow > b.tmigration]
    if len(LZ_cells)>0:
#        LZ_cells=set(LZ_cells) # making unique set
#        print(tnow,len(lists),len(LZ_cells))
        [lists.remove(x) for x in LZ_cells]
        
#        for x in LZ_cells:
#            x.xpos = 500 + round(np.random.random(),2)*box_width_lz
#            x.ypos = round(np.random.random(),2)*box_width_lz
#            
#        if abs(x.xpos) > box_width_lz+box_width_dz or abs(x.xpos) < box_width_lz:
#            x.xpos = round(np.random.random(),2)+box_width_lz
#            
#        if abs(x.ypos) > box_width_lz:
#            x.ypos = round(np.random.random(),2)*box_width_lz
        
#        LZ_cells=LZposition(LZ_cells)
            
        for b in LZ_cells:
            b.LZentrytime=tnow
            b.antg_clt=tnow+0.7
            b.tfh_contact_start = b.antg_clt + np.random.normal(3,0.1)
            b.tfh_contact_end = b.tfh_contact_start + np.random.normal(3,0.1)
            
    
#### antigen collection
            
            
#        
    lzcells=[b for b in LZ_cells]
    
    gen5=[b for b in lists if b.generation>8 and b.selected ==1]
#    gen6=set(random.choices(gen5, k=np.random.binomial(len(gen5),0.5)))
#    [lists.remove(x) for x in gen5]
#    
##    LZ_selected_cells=random.choices(gen5, k=np.random.binomial(len(gen5),0.003))
    LZ_selected_cells = [b for b in gen5 if tnow > b.tmigration]
#    
    if len(LZ_selected_cells)>0:
        LZ_selected_cells=set(LZ_selected_cells) # making unique set
#        print(tnow,len(lists),len(LZ_cells))
        [lists.remove(x) for x in LZ_selected_cells]
        
    for b in LZ_selected_cells:
        b.LZentrytime=tnow
        b.antg_clt=tnow+0.7
        b.tfh_contact_start = b.antg_clt + np.random.normal(3,0.1)
        b.tfh_contact_end = b.tfh_contact_start + np.random.normal(3,0.1)
        
#    print([[b,b.generation,b.tmigration,b.LZentrytime] for b in LZ_selected_cells])
#        
##        for x in LZ_selected_cells:
##            x.xpos += 500+x.xvel*dt
##            x.ypos += 500+x.yvel*dt
##            
##        if abs(x.xpos) > box_width_lz+box_width_dz:
##            x.xpos = -x.xvel
##            
##        if abs(x.ypos) > box_width_lz+box_width_dz:
##            x.ypos = -x.yvel
#            
#        for b in LZ_selected_cells:
#            b.LZentrytime=tnow
#            b.antg_clt=tnow+0.7
##        
#    lz_selected_cells=[b for b in LZ_selected_cells] ## didn't update any list assuming leaving from DZ
    
    
#    print([b.LZentrytime for b in lzcells])
#    appcells=[b for b in LZ if tnow > b.LZentrytime+6]
#    appcells=[b for b in LZ if tnow > b.LZentrytime+6]
#    ap_cells=random.choices(appcells, k=np.random.binomial(len(appcells),0.01))
#    ap_cells=set(ap_cells)
    
#    [LZ.remove(x) for x in ap_cells]
    LZ+=lzcells
    LZ+=LZ_selected_cells

    
    LZ=LZposition(LZ)
#    bb.append([[b,b.xpos,b.ypos,b.tfhcontact,b.tfh_sig_t] for b in LZ])
#    bb.append(LZ[:])
    

    LZ=antigen_collect(LZ)
    
#    appcells=[b for b in LZ if b.antigen_collect==0 and b.antg_clt < tnow]
    for b in LZ:
        if b.antigen_collect == 0 and b.antg_clt < tnow:
            b.deleting_time = b.antg_clt + np.random.normal(6,0.01)
            
    
#    print([['id=',b,'generation=',b.generation,'selected=',b.selected,
#            'birthtime=',b.birthtime,'migrationtime=',b.tmigration,'LZentrytime=',b.LZentrytime,
#            'fdcselected=',b.fdc_selected,'antigen collection time=',b.antg_clt,
#            'antigen collected=',b.antigen_collect,'tfh signal=',b.tfhsignal,
#            'tfh sig time=',b.tfh_sig_t,'tfh contact=',b.tfhcontact,
#            'tfh contact start=',b.tfh_contact_start,'tfh_contact_end=',
#            b.tfh_contact_end,'xpos=',b.xpos,'ypos=',b.ypos] for b in LZ])
    
    datatoexport.append([[b,b.generation,b.selected,b.birthtime,b.tmigration,b.LZentrytime,b.fdc_selected,b.antg_clt
                          ,b.antigen_collect,b.tfhsignal,b.tfh_sig_t,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end,b.xpos,b.ypos] for b in LZ])
    
    
            
    appcells = [b for b in LZ if b.antigen_collect == 0 and b.deleting_time < tnow]
    
    [LZ.remove(x) for x in appcells]
    
    pbc=[b for b in LZ if b.antigen_collect > 0 and b.tfh_contact_start < tnow and b.tfh_contact_end > tnow] ### how about implement this in t_help function, like LZposition and antigen collect function
    
    pbc=t_help(pbc)
    pot_bc.append(pbc)
    
    
    selected_cells=[b for b in LZ if b.tfhsignal > 50 and b.tfh_contact_end < tnow]

    
    unselected_cells=[b for b in LZ if b.tfhsignal < 50 and b.tfh_contact_end < tnow]
    
    for b in selected_cells:
        b.selected=1
#    selected_cells=list(selected_cells)
    
    # 10 percent cell leaves GC as early plasma cell
    
    ps_cells=random.choices(selected_cells, k=np.random.binomial(len(selected_cells),0.01))
    ps_cells=set(ps_cells)
    outputcells+=ps_cells
    
    print([len(selected_cells),len(ps_cells)])
    
    [LZ.remove(x) for x in selected_cells]
    
    [LZ.remove(x) for x in unselected_cells]
    
    [selected_cells.remove(x) for x in ps_cells]
    
    lists+=selected_cells
    
    
#### those cells need to be selected which are less than 6 hr as I implemented appcells    
#    selected_cells=random.choices(LZ, k=np.random.binomial(len(LZ),0.002))
#    selected_cells=set(selected_cells)
#    for b in selected_cells:
#        b.selected=1
#    selected_cells=list(selected_cells)
#    
#    # 10 percent cell leaves GC as early plasma cell
#    
#    ps_cells=random.choices(selected_cells, k=np.random.binomial(len(selected_cells),0.0001))
#    ps_cells=set(ps_cells)
#    outputcells+=ps_cells
#    
#    [LZ.remove(x) for x in selected_cells]
#    
#    [selected_cells.remove(x) for x in ps_cells]
#    
#    lists+=selected_cells
#        
##    print(tnow,len(gen4))
#    print([tnow,len(lists),len(LZ)])
##    print(len(lists),len(LZ))
#    
#    
    r.append([tnow,len(lists),len(LZ),len(outputcells)])
#
    
xcoord=[[b.xpos for b in pltlist[item]] for item in range(len(pltlist))]  
ycoord=[[b.ypos for b in pltlist[item]] for item in range(len(pltlist))]   

#df = DataFrame(datatoexport, columns= ['id', 'generation', 'selected','birthtime','migrationtime'
#                                           ,'LZentrytime','fdcselected','antigen collection time',
#                                           'antigen collected','tfh signal','tfh sig time','tfh contact',
#                                           'tfh contact start','tfh contact end','xpos','ypos'])
    
#export_excel = df.to_excel (r'export_dataframe.xlsx', index = None, header=True)



with open("test.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(datatoexport)


def animate(i):
    particles.set_data(xcoord[i], ycoord[i])
    time_text.set_text(time_template % (i))
    return particles, time_text

#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=50000, interval=10000, blit=True, repeat=False)
#
#anim.save('gc.mp4',writer=writer)

plt.show()

r=np.array(r)


fig=plt.figure()
plt.subplot(2,1,1)
np.seterr(divide='ignore', invalid='ignore')
plt.plot(r[:,0]/24, r[:,1]/r[:,2])
plt.plot([0,20],[2,2],[0,20],[3,3],'r')
plt.axis([0,20,0,10])
plt.ylabel('DZ/LZ ratio')

plt.subplot(2,1,2)
#labels=(,'LZ Bcells','All Bcells')
p1=plt.plot(r[:,0]/24, r[:,1], label='DZ Bcells')
p2=plt.plot(r[:,0]/24,r[:,2], label='LZ Bcells')
p3=plt.plot(r[:,0]/24,r[:,1]+r[:,2], label='All Bcells')
p4=plt.plot(r[:,0]/24,r[:,3], label='output Bcells')
plt.legend()
plt.show()

fig.savefig('gc.pdf')
#plt.plot(r[:,0],r[:,1],r[:,2])
#plt.plot(r[:,0],r[:,1]/r[:,2])

end = time.time()

print(end-start)

    
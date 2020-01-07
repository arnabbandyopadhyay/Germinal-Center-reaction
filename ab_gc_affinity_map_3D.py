#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 12:20:26 2019

@author: abp19
"""

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
import itertools
from mpl_toolkits.mplot3d import Axes3D 
#import antigen
#import multiprocessing
#pool = multiprocessing.Pool(processes=2)
from multiprocessing import Pool
pool=Pool(2)

start = time.time()

mu, sigma = 7.1, 3 # mean and standard deviation
mu2, sigma2 = 7, 12
box_width_dz=500
box_width_lz=500
opseq=[3,3,3,3]
mum, sigm = 5, 1 # hr from LZ to DZ and vice versa
class Bcell:
    """ Class of B cell that contains the followings:
        - current generation
        - whether it is selected or not
        - whether it is in contact with Tfh
        - LZ entry time
        - current sequence
        - current affinity
        - amount of antigen collected
        - fdc selected
        - antigen collection time
        - tfh signal aquired
        - tfh signaling time
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
        else:                  self.LZentrytime = 600
        
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
        else: self.tfh_contact_end = 600
        
        if 'tfh_attached' in kwargs: self.tfh_attached = kwargs['tfh_attached']
        else: self.tfh_attached = 'Null'
        
        if 'ndiv' in kwargs: self.ndiv = kwargs['ndiv']
        else: self.ndiv = np.random.randint(4,6)
        
        if 'deleting_time' in kwargs: self.deleting_time = kwargs['deleting_time']
        else: self.deleting_time = 600
        
        if 'xpos' in kwargs: self.xpos = kwargs['xpos']
        else: self.xpos = round(np.random.random(),2)*box_width_dz
        
        if 'ypos' in kwargs: self.ypos = kwargs['ypos']
        else: self.ypos = round(np.random.random(),2)*box_width_dz
        
        if 'zpos' in kwargs: self.zpos = kwargs['zpos']
        else: self.zpos = round(np.random.random(),2)*box_width_dz
        
        
        
        
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
                   ndiv=b.ndiv, deleting_time = b.deleting_time,
                   nextdiv=tnow+np.random.randint(mu2, sigma2), xpos = round(np.random.random(),2)*box_width_dz,
                   ypos = round(np.random.random(),2)*box_width_dz, zpos = round(np.random.random(),2)*box_width_dz)
    def gen(self):
        self.generation += 1
#        self.ndiv = self.generation + pmin + round((pmax-pmin)*((b.antigen_collect**2)/((b.antigen_collect**2)+(kp**2))))

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


class Tcell:
    """ B cell objects are the main unit of this simulation. They consist of
            the non-key region
    """
    def __init__(self):
        self.xpos=np.random.randint(0,490)
         
        self.ypos=np.random.randint(-490,490)
        self.zpos=np.random.randint(-490,490)
        
        self.bcontact=[]




n_tfh=200
tlist = [Tcell() for t in range(n_tfh)]


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

def distance(ll):
    dist=math.sqrt((b.xpos-ll[0])**2 + (b.ypos-ll[1])**2)
    return dist


def antigen_collect_map(b):
    dist=[math.sqrt((b.xpos-i[0])**2 + (b.ypos-i[1])**2 + (b.zpos-i[2])**2) for i in zip(antg_x_coord,antg_y_coord,antg_z_coord)]

    min_dist=min(dist)
    if min_dist < radius_fdc:
        pos=dist.index(min_dist)
        antigen_collect=np.random.binomial(1,b.affinity)
        if antigen_collect ==1 and antg_load[pos]>0:
            antg_load[pos] -=1
            b.antigen_collect +=1
            b.fdc_selected = 1
    return b


def antigen_collect_2(lists,antg_load):
    dist=[[math.sqrt((b.xpos-i[0])**2 + (b.ypos-i[1])**2 + (b.zpos-i[2])**2) for i in zip(antg_x_coord,antg_y_coord, antg_z_coord)] for b in lists]
    min_dist=[min(dist[i]) for i in range(len(lists))]
    dist2=[i for i in min_dist if i<radius_fdc]
    posl=[[i for i,val in enumerate(min_dist) if val==j] for j in dist2]
    posl2=[list(t) for t in set(tuple(element) for element in posl)] ## deleting duplicates
    pos2=[k for i in posl2 for k in i]

    for i in pos2:
        b=lists[i]
        antigen_collect=np.random.binomial(1,b.affinity)
        if antigen_collect ==1 and antg_load[dist[i].index(min(dist[i]))]>0:
            antg_load[dist[i].index(min(dist[i]))] -=1
            b.antigen_collect +=1
            b.fdc_selected = 1
    return lists, antg_load 



def antigen_collect(lists,antg_load):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    bc=[b for b in lists if b.LZentrytime <= tnow and tnow <= b.antg_clt]
#    for i in range(len(bc)):
#        b=bc[i]
#        dist=[math.sqrt((b.xpos-i[0])**2 + (b.ypos-i[1])**2) for i in zip(antg_x_coord,antg_y_coord)]
#        min_dist=min(dist)
#        if min_dist < radius_fdc:
#            pos=dist.index(min_dist)
#            antigen_collect=np.random.binomial(1,b.affinity)
#            if antigen_collect ==1 and antg_load[pos]>0:
#                antg_load[pos] -=1
#                b.antigen_collect +=1
#                b.fdc_selected = 1
    bc=list(map(antigen_collect_map2,bc))
#   pool.close()
        
    return lists, antg_load


#xpos=[[b.xpos-5+i,b.ypos-5+j] for i in range(11) for j in range(11)]
#antg_xy=[[i,j] for i,j in zip(antg_x_coord,antg_y_coord)]
#res=[]
#for _a in xpos:
#  if _a in antg_xy:
#    res.append(_a)
#    antg_xy.remove(_a)
#pp=[]



def intersection(lst1, lst2): 
    tup1 = map(tuple, lst1) 
    tup2 = map(tuple, lst2)  
    return list(map(list, set(tup1).intersection(tup2))) 


def antigen_collect_map2(b):
    xpos=[[b.xpos-5+i,b.ypos-5+j, b.zpos-5+k] for i in range(11) for j in range(11) for k in range(11)]
        
#    res=[]
#    for _a in xpos:
#        if _a in antg_xy:
#            res.append(_a)
#            break
    res=intersection(xpos,antg_xy)
    if len(res)>0:
        pos=antg_xy.index(random.choices(res,k=1)[0])
        antigen_collect=np.random.binomial(1,b.affinity)
        if antigen_collect ==1 and antg_load[pos]>0:
            antg_load[pos] -=1
            b.antigen_collect +=1
            b.fdc_selected = 1
    return b
#def antigen_collect(lists,antg_load):
#    """ Returns lists for keeping track of free (outside of GCs) memory and
#    naive B cells as well as a list of lists of B cells waiting for surivival
#    signals in each GC. """
#    bc=[b for b in lists if b.LZentrytime <= tnow and tnow <= b.antg_clt]
##    for i in range(len(bc)):
##        b=bc[i]
##        xpos=[[b.xpos-5+i,b.ypos-5+j] for i in range(11) for j in range(11)]
##        
##        res=[]
##        for _a in xpos:
##            if _a in antg_xy:
##                res.append(_a)
##                break
##        if len(res)>0:
##            pos=antg_xy.index(res[0])
##            antigen_collect=np.random.binomial(1,b.affinity)
##            if antigen_collect ==1 and antg_load[pos]>0:
##                antg_load[pos] -=1
##                b.antigen_collect +=1
##                b.fdc_selected = 1
#    bc=list(pool.map(antigen_collect_map2,bc))
#    pool.close()
#        
#    return lists, antg_load






def contain_points(check_x, check_y, check_z, center_x, center_y, center_z, radius):
    if (check_x - center_x)**2 + (check_y - center_y)**2 + (check_z - center_z)**2 < radius**2:
        return('True')
        
    else:
        return('False')

def t_help(lists):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    pbc=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end ] ### how about implement this in t_help function, like LZposition and antigen collect function
    
    for i in range(n_tfh):
        nl=[]
        for j in range(len(pbc)):
            b=pbc[j]
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


def t_help2(lists):
    """ Returns lists for keeping track of free (outside of GCs) memory and
    naive B cells as well as a list of lists of B cells waiting for surivival
    signals in each GC. """
    pbc1=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==1]
    for k in range(len(pbc1)):
        b=pbc1[k]
        if b.tfhcontact==1 and b.tfh_sig_t < tnow:
            b.tfhcontact=0
            b.tfh_sig_t=0
        elif b.tfhcontact==1 and b.tfh_sig_t > tnow:
            b.tfhcontact=1
            b.tfh_sig_t=b.tfh_sig_t
            
    pbc=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==0 ] ### how about implement this in t_help function, like LZposition and antigen collect function

    
    for i in range(n_tfh):
        poslist=[]
        dist=[[b,math.sqrt((b.xpos-tfh_x_coord[i])**2 + (b.ypos-tfh_y_coord[i])**2)] for b in pbc if math.sqrt((b.xpos-tfh_x_coord[i])**2 + (b.ypos-tfh_y_coord[i])**2) < 50]
#        mdist=[l for l in dist if l<20]
        if len(dist)>0:
            dist2=[x[0] for x in dist]
            for j in range(len(dist2)):
                b=dist2[j]
                b.tfhcontact=1
                b.tfh_sig_t=tnow+0.1
                poslist.append(b)
            
            [pbc.remove(x) for x in poslist]
            
            
            k=[b.antigen_collect for b in poslist]
            pos=k.index(max(k))
            poslist[pos].tfhsignal += timestep
                   
    pbc2=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and tnow >= b.tfh_contact_end]
    
    for k in range(len(pbc2)):
        b=pbc2[k]
        b.tfhcontact=0
        b.tfh_sig_t=0
        
        
    return lists

def thelp3(lists):
    
    """ 
    bcells that already collected antigen and are within proper time and contact less
    """
    pbc=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==0 ]
    """ 
    tcells that can accomodate another b cell
    """
    at=[t for t in tlist if len(t.bcontact) < 6]
    
    
    for i in range(len(at)):
        t=at[i]
        canaccomodate=6-len(t.bcontact)
        dist=[[b,math.sqrt((b.xpos-t.xpos)**2 + (b.ypos-t.ypos)**2)] for b in pbc if math.sqrt((b.xpos-t.xpos)**2 + (b.ypos-t.ypos)**2) < 20]
#        mdist=[l for l in dist if l<20]
        
        if len(dist)>0:
            dist2 = sorted(dist,key=lambda random: random[1])
            for j in range(min(canaccomodate,len(dist))):
                t.bcontact.append(dist2[j][0])
                dist2[j][0].tfhcontact=1
                dist2[j][0].tfh_sig_t=tnow+0.1
                dist2[j][0].tfh_attached=t
                #
                
    """ 
    bcells that already collected antigen and are within proper time and are under contact
    """
    pbc1=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==1]
    
    
    for k in range(len(pbc1)):
        b=pbc1[k]
        if b.tfhcontact==1 and b.tfh_sig_t < tnow:
            b.tfhcontact=0
            b.tfh_sig_t=0
            tc=b.tfh_attached
            tc.bcontact.remove(b)
            b.tfh_attached='null'
                        
            
    """
    Tcell signaling
    """
    for tcell in tlist:
        antgbcell=sorted([[b,b.antigen_collect] for b in tcell.bcontact],key=lambda random: random[1],reverse=True)
        if len(antgbcell) > 0:
            antgbcell[0][0].tfhsignal += timestep
                
    
    """ 
    bcells that already collected antigen and are out of time
    """
    pbc2=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and tnow >= b.tfh_contact_end]
    
    for k in range(len(pbc2)):
        b=pbc2[k]
        if b.tfhcontact==1:
            b.tfhcontact=0
            b.tfh_sig_t=0
            tc=b.tfh_attached
            tc.bcontact.remove(b)
            b.tfh_attached='null'
            
    return lists

def module1(t):
    global pbc
#    pbc=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==0 ]
    canaccomodate=6-len(t.bcontact)
    dist=[[b,math.sqrt((b.xpos-t.xpos)**2 + (b.ypos-t.ypos)**2 + (b.zpos-t.zpos)**2)] for b in pbc if math.sqrt((b.xpos-t.xpos)**2 + (b.ypos-t.ypos)**2 + (b.zpos-t.zpos)**2) < 20]
#        mdist=[l for l in dist if l<20]
        
    if len(dist)>0:
        dist2 = sorted(dist,key=lambda random: random[1])
        for j in range(min(canaccomodate,len(dist))):
            t.bcontact.append(dist2[j][0])
            dist2[j][0].tfhcontact=1
            dist2[j][0].tfh_sig_t=tnow+0.1
            dist2[j][0].tfh_attached=t
    return t


def module2(b):
    if b.tfhcontact==1 and b.tfh_sig_t < tnow:
        b.tfhcontact=0
        b.tfh_sig_t=0
        tc=b.tfh_attached
        tc.bcontact.remove(b)
        b.tfh_attached='null'
    return b

def module3(tcell):
    antgbcell=sorted([[b,b.antigen_collect] for b in tcell.bcontact],key=lambda random: random[1],reverse=True)
    if len(antgbcell) > 0:
        antgbcell[0][0].tfhsignal += timestep
    return tcell

def module4(b):
    if b.tfhcontact==1:
        b.tfhcontact=0
        b.tfh_sig_t=0
        tc=b.tfh_attached
        tc.bcontact.remove(b)
        b.tfh_attached='null'
    return b
    

def thelp33(lists,tlist):
    
    """ 
    bcells that already collected antigen and are within proper time and contact less
    """
    pbc=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==0 ]
    """ 
    tcells that can accomodate another b cell
    """
    at=[t for t in tlist if len(t.bcontact) < 6]
    
    
    for i in range(len(at)):
        t=at[i]
        canaccomodate=6-len(t.bcontact)
        dist=[[b,math.sqrt((b.xpos-t.xpos)**2 + (b.ypos-t.ypos)**2 + (b.zpos-t.zpos)**2)] for b in pbc if math.sqrt((b.xpos-t.xpos)**2 + (b.ypos-t.ypos)**2 + (b.zpos-t.zpos)**2) < 20]
#        mdist=[l for l in dist if l<20]
        
        if len(dist)>0:
            dist2 = sorted(dist,key=lambda random: random[1])
            for j in range(min(canaccomodate,len(dist))):
                t.bcontact.append(dist2[j][0])
                dist2[j][0].tfhcontact=1
                dist2[j][0].tfh_sig_t=tnow+0.1
                dist2[j][0].tfh_attached=t
                #
                
#    at=list(map(module1,at))
    
    """ 
    bcells that already collected antigen and are within proper time and are under contact
    """
    pbc1=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start <= tnow and tnow < b.tfh_contact_end and b.tfhcontact==1]
    
    
#    for k in range(len(pbc1)):
#        b=pbc1[k]
#        if b.tfhcontact==1 and b.tfh_sig_t < tnow:
#            b.tfhcontact=0
#            b.tfh_sig_t=0
#            tc=b.tfh_attached
#            tc.bcontact.remove(b)
#            b.tfh_attached='null'
#    pool = multiprocessing.Pool(processes=2)
    pbc1=list(map(module2,pbc1))
#    pool.close()
                        
            
    """
    Tcell signaling
    """
#    for tcell in tlist:
#        antgbcell=sorted([[b,b.antigen_collect] for b in tcell.bcontact],key=lambda random: random[1],reverse=True)
#        if len(antgbcell) > 0:
#            antgbcell[0][0].tfhsignal += timestep
#    pool = multiprocessing.Pool(processes=2)
    tlist=list(map(module3,tlist))
#    pool.close()
                
    
    """ 
    bcells that already collected antigen and are out of time
    """
    pbc2=[b for b in lists if b.antigen_collect > 0 and b.antg_clt < tnow and tnow >= b.tfh_contact_end]
    
#    for k in range(len(pbc2)):
#        b=pbc2[k]
#        if b.tfhcontact==1:
#            b.tfhcontact=0
#            b.tfh_sig_t=0
#            tc=b.tfh_attached
#            tc.bcontact.remove(b)
#            b.tfh_attached='null'
#    pool = multiprocessing.Pool(processes=2)
    pbc2=list(map(module4,pbc2))
#    pool.close()
            
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


def dzdynamics(b):
    global nl2
    if b.selected==0:
        if (b.nextdiv<=tnow) and (b.generation<b.ndiv):
            rc.append(b)
#            print('y')
#                generation=b.generation+1
            nl=proliferate2(b)
            nl2+=nl

            
        
    elif b.selected==1:
        if (b.nextdiv<=tnow) and (b.generation<b.ndiv):
            rc.append(b)
#            print('y')
#                generation=b.generation+1
            nl=proliferate2(b)
            nl2+=nl

    

def DZdynamics(lists):
    nl2=[]
    rc=[]
    nl=[]
#    for i in range(len(lists)):
#        
#        nl=[]
##        print(i)
##        print(len(lists))
#        b=lists[i]
##        x1=b.xvel*dt
##        y1=b.yvel*dt
##        b.xpos += b.xvel*dt
##        b.ypos += b.yvel*dt
##        
##        if abs(b.xpos) > box_width_dz:
##            b.xpos = -b.xvel
##            
##        if abs(b.ypos) > box_width_dz:
##            b.ypos = -b.yvel
#        
#        if b.selected==0:
#            if (b.nextdiv<=tnow) and (b.generation<b.ndiv):
#                rc.append(b)
##            print('y')
##                generation=b.generation+1
#                nl=proliferate2(b)
#        elif b.selected==1:
#            if (b.nextdiv<=tnow) and (b.generation<b.ndiv):
#                rc.append(b)
##            print('y')
##                generation=b.generation+1
#                nl=proliferate2(b)
#        nl2+=nl
    list(map(dzdynamics,lists))
#    pool.close()
    nl2+=nl
    return nl2,rc

def LZposition(lists):
    
    for x in lists:
        if x.tfhcontact==1:
            x.xpos=x.xpos
            x.ypos=x.ypos
            x.zpos=x.zpos
        else:
            xpos = np.random.randint(0,490)
            ypos = np.random.randint(-490,490)
            zpos = np.random.randint(-490,490)
            while math.sqrt((xpos)**2 + (ypos)**2 + (zpos)**2)<=490:
                x.xpos=xpos
                x.ypos=ypos
                x.zpos=zpos
                
                xpos = np.random.randint(0,490)
                ypos = np.random.randint(-490,490)
                zpos = np.random.randint(-490,490)
            
#            if abs(x.xpos) > box_width_lz+box_width_dz or abs(x.xpos) < box_width_lz:
#                x.xpos = round(np.random.random(),2)+box_width_lz
#            
#            if abs(x.ypos) > box_width_lz:
#                x.ypos = round(np.random.random(),2)*box_width_lz
        
    return lists

def LZposition2(x):
    if x.tfhcontact==1:
        x.xpos=x.xpos
        x.ypos=x.ypos
        x.zpos=x.zpos
    else:
        xpos = np.random.randint(0,490)
        ypos = np.random.randint(-490,490)
        zpos = np.random.randint(-490,490)
        while math.sqrt((xpos)**2 + (ypos)**2 + (zpos)**2)<=490:
            x.xpos=xpos
            x.ypos=ypos
            x.zpos=zpos
                
            xpos = np.random.randint(0,490)
            ypos = np.random.randint(-490,490)
            zpos = np.random.randint(-490,490)
        
    return(x)

def DZposition(x):
    
    xpos = np.random.randint(-490,0)
    ypos = np.random.randint(-490,490)
    zpos = np.random.randint(-490,490)
    while math.sqrt((xpos)**2 + (ypos)**2 + (zpos)**2)<=490:
        x.xpos=xpos
        x.ypos=ypos
        x.zpos=zpos
                
        xpos = np.random.randint(-490,0)
        ypos = np.random.randint(-490,490)
        zpos = np.random.randint(-490,490)
        
    return(x)
    
def fun(b):
    b.antigen_collect=0
    b.LZentrytime=tnow
    b.antg_clt=tnow+0.7
    b.tfh_contact_start = b.antg_clt + np.random.normal(0.5,0.1)
    b.tfh_contact_end = b.tfh_contact_start + np.random.normal(5,0.1)
    b.tfhsignal=0
    return b

def fun2(b):
    b.selected=1
    b.nextdiv=tnow+np.random.randint(1, 3)
    b.tmigration=600
    b.LZentrytime=600
    b.ndiv=b.generation + pmin+round((pmax-pmin)*((b.antigen_collect**2)/((b.antigen_collect**2)+(kp**2))))
    b.antigen_collect=0
    b.tfhsignal=0
    return b

def fun3(b):
    b.LZentrytime=tnow
    b.antg_clt=tnow+0.7
    b.tfh_contact_start = b.antg_clt + np.random.normal(0.5,0.1)
    b.tfh_contact_end = b.tfh_contact_start + np.random.normal(5,0.1)
    return b

def fun4(b):
    if b.antigen_collect == 0 and b.antg_clt < tnow:
        b.deleting_time = b.antg_clt + np.random.normal(6,0.1)
    return b

def fun5(b):
    if b.tfh_contact_end < tnow:
        b.tfhcontact=0
    return b
            
def antigen_distribution(ntotal,ndiv):
    a=round(ndiv/2)
    a1=round(ntotal/a)
    ar=[]
    for i in range(a):
        a2=random.randint(1,a1)
        ar+=[a2,a1-a2]
    return ar

    

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(1, 1000), ylim=(1, 500))
plt.axvline(x=500,linewidth=4,color='r')

n_tfh=200
#tfh_x_coord=np.random.randint(500,1000,n_tfh)
#tfh_y_coord=np.random.randint(0,500,n_tfh)

#antigen_num=1000000
tlist = [Tcell() for t in range(n_tfh)]

#tt=np.array([[t.xpos,t.ypos] for t in tlist]);plt.plot(tt[:,0],tt[:,1],'ro',ms=1)

#plt.plot(tfh_x_coord,tfh_y_coord,'ro',ms=1)
n_fdc=200
antg_p_fdc=3000
arm_fdc=6
fdc_x_coord=np.random.randint(0,490,n_fdc)
fdc_y_coord=np.random.randint(-490,490,n_fdc)
fdc_z_coord=np.random.randint(-490,490,n_fdc)
antg_x_coord=[round(np.random.normal(k, 5),0) for k in fdc_x_coord for b in range(6)]
antg_y_coord=[round(np.random.normal(k, 5),0) for k in fdc_y_coord for b in range(6)]
antg_z_coord=[round(np.random.normal(k, 5),0) for k in fdc_z_coord for b in range(6)]
#antg_x_coord=[k+5*math.cos(b) for k in fdc_x_coord for b in range(60,361,60)]
#antg_y_coord=[k+5*math.sin(b) for k in fdc_y_coord for b in range(60,361,60)]
antg_xy=[[i,j,k] for i,j,k in zip(antg_x_coord,antg_y_coord,antg_z_coord)]
n_antg=500

#antg_load=[n_antg for n in antg_x_coord]

antg_load=list(itertools.chain(*[antigen_distribution(antg_p_fdc,arm_fdc) for i in range(n_fdc)]))

#plt.plot(antg_x_coord,antg_y_coord,'bo',ms=1)


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
radius=5
radius_fdc=5
dt=1
generation=0
tnow=0
timestep=0.01
pre=2 # precession for timestep
pmin=1
pmax=6
kp=12
lists=[Bcell(birthtime=0),Bcell(birthtime=0)]
#lists=[Bcell(birthtime=0) for b in range(144)]
pltlist=[lists[0::]]
LZ=[]
pot_bc=[]
outputcells=[]
r=[]
bb=[]
appopcells=[]
affy=[]
antgcollect=[]
affy2=[]
ct=[]
#trange=np.arange(1, 250, 0.1).tolist()
for tnow in np.arange(1,500,timestep):
#    pool = multiprocessing.Pool(processes=2)
    l890 = time.time()
    if (round(tnow,pre) % 10) == 0:
        print(tnow, len(lists),len(LZ),np.mean([b.antigen_collect for b in LZ]),len(selected_cells),len(ps_cells),len(outputcells))
#    print(tnow, len(lists))
    
    lists=list(map(DZposition,lists))
    
    bbc=[]
    
    birthtime=tnow
#    print(tnow,len(lists),[b.generation for b in lists])
#    print(tnow,len(lists),[[b.birthtime,b.nextdiv] for b in lists])
#    dlist=[b.nextdiv for b in lists]
#    dlist2=[i for i,x in enumerate(dlist) if x <tnow]
#    nl2=[]
#    rc=[]

#    nlists=lists
#    for i in range(len(lists)):
#        nl=[]
##        print(i)
##        print(len(lists))
#        b=lists[i]
##        x1=b.xvel*dt
##        y1=b.yvel*dt
##        b.xpos += b.xvel*dt
##        b.ypos += b.yvel*dt
##        
##        if abs(b.xpos) > box_width_dz:
##            b.xpos = -b.xvel
##            
##        if abs(b.ypos) > box_width_dz:
##            b.ypos = -b.yvel
#        
#        if b.selected==0:
#            if (b.nextdiv<=tnow) and (b.generation<6):
#                rc.append(b)
##            print('y')
#                generation=b.generation+1
#                nl=proliferate2(b)
#        elif b.selected==1:
#            if (b.nextdiv<=tnow) and (b.generation<18):
#                rc.append(b)
##            print('y')
#                generation=b.generation+1
#                nl=proliferate2(b)
#        nl2+=nl
        
#    nl2,rc=DZdynamics(lists)
    rc=[]
    nl2=[]
    list(map(dzdynamics,lists))
    

#    pool.close()

#    print('dz dynamics =',l865-l863)
    
#    nl2+=nl
    lists+=nl2

    
    [lists.remove(x) for x in rc]
            
#    if len(nl)>0:
#        lists.append(nl)
    if tnow > 72:
#        print('y')
        if (round(tnow,pre) % 3) == 0:
            print('y')
#            lists.append(Bcell(birthtime))  ### 1 cell per 3 hour entry in GC
#    lists.append(Bcell(birthtime))
    else:
        if (round(tnow,pre) % 1) == 0:
            lists.append(Bcell(birthtime))
            lists.append(Bcell(birthtime))
            lists.append(Bcell(birthtime))

#    pltlist2.append(lists)
            
#    pltlist.append(lists[0::])
    
    
#    print(lists)
#    print(pltlist)
    
    LZ_cells=[b for b in lists if b.generation==b.ndiv and b.selected == 0 and tnow > b.tmigration]
    
#    LZ_cells=random.choices(gen4, k=np.random.binomial(len(gen4),0.05))
#    LZ_cells = [b for b in gen4 if tnow > b.tmigration]
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
            
#        for b in LZ_cells:
#            b.LZentrytime=tnow
#            b.antg_clt=tnow+0.7
#            b.tfh_contact_start = b.antg_clt + np.random.normal(0.5,0.1)
#            b.tfh_contact_end = b.tfh_contact_start + np.random.normal(5,0.1)
        
        LZ_cells=list(map(fun3,LZ_cells))
#        pool.close()
        
#        print('l_919 =',l921-l919)
        LZ+=LZ_cells
    l1003 = time.time()          
#### antigen collection
            
            
#        
#    lzcells=[b for b in LZ_cells]
    l1009 = time.time()   
    
    gen5=[b for b in lists if b.generation==b.ndiv and b.selected ==1 and tnow > b.tmigration]
#    gen6=[b for b in gen5 if tnow > b.tmigration]
#    gen6=set(random.choices(gen5, k=np.random.binomial(len(gen5),0.1)))
#    print([[b,b.birthtime, b.tmigration,b.LZentrytime,b.generation] for b in gen6])
#    [lists.remove(x) for x in gen5]
#    
##    LZ_selected_cells=random.choices(gen5, k=np.random.binomial(len(gen5),0.003))
#    LZ_selected_cells = set(random.choices(gen6, k=np.random.binomial(len(gen6),0.5)))
    
    LZ_selected_cells = gen5
    
#    print([[b,b.birthtime, b.tmigration,b.LZentrytime,b.generation] for b in LZ_selected_cells])
#    
    if len(LZ_selected_cells)>0:
#        LZ_selected_cells=set(LZ_selected_cells) # making unique set
#        print(tnow,len(lists),len(LZ_cells))
        [lists.remove(x) for x in LZ_selected_cells]
        
#        for b in LZ_selected_cells:
#            
#            b.antigen_collect=0
#            b.LZentrytime=tnow
#            b.antg_clt=tnow+0.7
#            b.tfh_contact_start = b.antg_clt + np.random.normal(0.5,0.1)
#            b.tfh_contact_end = b.tfh_contact_start + np.random.normal(5,0.1)
#            b.tfhsignal=0
        LZ_selected_cells=list(map(fun,LZ_selected_cells))
#        pool.close()
            
        LZ+=LZ_selected_cells
#    print([[b,b.birthtime, b.tmigration,b.LZentrytime,b.generation] for b in LZ_selected_cells])
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
#    LZ+=LZ_cells
#    LZ+=LZ_selected_cells

#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in lists if b.generation ==6 and b.tmigration < tnow])
#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in lists if b.generation ==9 and b.tmigration < tnow])
#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in lists if b.generation ==10 and b.tmigration < tnow])
#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in lists if b.generation ==11 and b.tmigration < tnow])
#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in lists if b.generation ==12 and b.tmigration < tnow])
#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in lists if b.generation ==13 and b.tmigration < tnow])

#    LZ=LZposition(LZ)

#    pool = multiprocessing.Pool(processes=2)
    LZ=list(map(LZposition2,LZ))
#    pool.close()
    
#    print('lz position =',l999-l997)
#    bb.append([[b,b.xpos,b.ypos,b.tfhcontact,b.tfh_sig_t] for b in LZ])
#    bb.append(LZ[:])
    
    
    bc=[b for b in LZ if b.LZentrytime <= tnow and tnow <= b.antg_clt]
#    pool = multiprocessing.Pool(processes=2)
#    bc=list(map(antigen_collect_map,bc))
#    pool.close()


    LZ, antg_load =antigen_collect(LZ,antg_load)
    
#    print('antigencollect =',l1006-l1004)
#    print(sum(antg_load))
    
#    print([[b,b.birthtime, b.tmigration, b.LZentrytime,b.antg_clt, b.affinity,b.antigen_collect] for b in LZ])
    
#    appcells=[b for b in LZ if b.antigen_collect==0 and b.antg_clt < tnow]
#    for b in LZ:
#        if b.antigen_collect == 0 and b.antg_clt < tnow:
#            b.deleting_time = b.antg_clt + np.random.normal(6,0.1)
    
#    pool = multiprocessing.Pool(processes=2)
    LZ=list(map(fun4,LZ))
#    pool.close()
    
#    print('l1017 =',l1018-l1016)
    
    appcells = [b for b in LZ if b.antigen_collect == 0 and b.deleting_time < tnow]
    
    appopcells+=appcells
    [LZ.remove(x) for x in appcells]
    
    l1117 = time.time()
#    pbc=[b for b in LZ if b.antigen_collect > 0 and b.antg_clt < tnow and b.tfh_contact_start < tnow and b.tfh_contact_end > tnow] ### how about implement this in t_help function, like LZposition and antigen collect function
#    print(len([b.antigen_collect for b in LZ if b.antigen_collect > 0]))
    l1120 = time.time()
#    pool = multiprocessing.Pool(processes=2)
    LZ=thelp33(LZ,tlist)
    
#    print('thelp =',l1030-l1028)
#    LZ=thelp3(LZ)
#    for b in LZ:
#        if b.tfh_contact_end < tnow:
#            b.tfhcontact=0
    
#    ct.append([l865-l863,l921-l919,l999-l997,l1006-l1004,l1018-l1016,l1030-l1028])
    
    LZ=list(map(fun5,LZ))
#    pool.close()
    
#    pot_bc.append(pbc)
    
    affy.append(LZ)
    affy2.append([tnow,np.mean([b.affinity for b in LZ])])
    antgcollect.append([[b.generation,b.affinity,b.selected, b.antigen_collect,b.nextdiv,b.tmigration,b.ndiv] for b in LZ])
            
        
    selected_cells=[b for b in LZ if b.tfhsignal > 0.5 and b.tfh_contact_end < tnow]
#    print(len(selected_cells))
    
    unselected_cells=[b for b in LZ if b.tfhsignal <= 0.5 and b.tfh_contact_end < tnow]
    appopcells+=unselected_cells
    
#    selected_cells=[b for b in LZ if b.antigen_collect >6 and b.tfh_contact_end < tnow]
#    print(len(selected_cells))
#    
#    unselected_cells=[b for b in LZ if b.antigen_collect <= 6 and b.tfh_contact_end < tnow]
#    
    
    
#    print([[b.selected,b.generation,b.birthtime,b.nextdiv,b.tmigration,b.LZentrytime,b.antigen_collect,b.antg_clt,b.tfhsignal,b.tfhcontact,b.tfh_contact_start,b.tfh_contact_end] for b in unselected_cells if b.tfhsignal>0])
    
#    print([len(lists),len(LZ_cells),len(LZ_selected_cells),len(appcells),len(selected_cells),len(unselected_cells)])
    
    [LZ.remove(x) for x in unselected_cells]
    
#    for b in selected_cells:
#        b.selected=1
#        b.nextdiv=tnow+np.random.randint(1, 3)
#        b.tmigration=600
#        b.LZentrytime=600
#        b.ndiv=b.generation + pmin+round((pmax-pmin)*((b.antigen_collect**2)/((b.antigen_collect**2)+(kp**2))))
#        b.antigen_collect=0
#        b.tfhsignal=0
#    pool = multiprocessing.Pool(processes=2)
    selected_cells=list(map(fun2,selected_cells))
#    pool.close()
    
#    selected_cells=list(selected_cells)
    
    # 10 percent cell leaves GC as early plasma cell
    ps_cells=[]
    for s in selected_cells:
        randn=np.random.random()
        if randn<=0.3:
            ps_cells.append(s)
        
#    ps_cells=random.choices(selected_cells, k=np.random.binomial(len(selected_cells),0.1))
#    ps_cells=set(ps_cells)
    outputcells+=ps_cells
#    print([len(selected_cells),len(ps_cells)])
    
    [LZ.remove(x) for x in selected_cells]
    
    [selected_cells.remove(x) for x in ps_cells]
    
    lists+=selected_cells
    
#    pool.close()
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
    r.append([tnow,len(lists),len(LZ),len(outputcells),len(appopcells)])
    
    l1220 = time.time()
    ct.append([l1003-l890,l1117-l1009,l1220-l1120])
#
    
xcoord=[[b.xpos for b in pltlist[item]] for item in range(len(pltlist))]  
ycoord=[[b.ypos for b in pltlist[item]] for item in range(len(pltlist))]  
zcoord=[[b.zpos for b in pltlist[item]] for item in range(len(pltlist))]   


def animate(i):
    particles.set_data(xcoord[i], ycoord[i], zcoord[i])
    time_text.set_text(time_template % (i))
    return particles, time_text

#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=50000, interval=10000, blit=True, repeat=False)
#
#anim.save('gc.mp4',writer=writer)

#plt.show()



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x=[b.xpos for b in lists]
y=[b.ypos for b in lists]
z=[b.zpos for b in lists]
ax.scatter(x, y, z, c='r')
#    ax.scatter(xs, ys, zs, marker=m)
x1=[b.xpos for b in LZ]
y1=[b.ypos for b in LZ]
z1=[b.zpos for b in LZ]
ax.scatter(x1, y1, z1, c='g')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Now set color to white (or whatever is "invisible")
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')

# Bonus: To get rid of the grid as well:
ax.grid(False)

plt.show()
fig.savefig('gc2.pdf')









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
#p5=plt.plot(r[:,0]/24,r[:,4], label='appop Bcells')
plt.legend()


#plt.show()

fig.savefig('gc.pdf')
#plt.plot(r[:,0],r[:,1],r[:,2])
#plt.plot(r[:,0],r[:,1]/r[:,2])

affy2=np.array(affy2)

fig=plt.figure()
#labels=(,'LZ Bcells','All Bcells')
plt.plot(affy2[:,0], affy2[:,1], label='Affinity')
plt.legend()

#plt.show()

fig.savefig('affinity.pdf')

end = time.time()

print(end-start)

#max([antgcollect[b][i][0] for b in range(1000,1500) for i in range(len(antgcollect[b]))])
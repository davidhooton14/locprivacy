# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 14:44:45 2019

@author: David Hooton
"""

import math, random, scipy, numpy as np
from scipy import special, stats
from paths import *
from levenshtein import find_all_matches, lookup_path

def GeoLBS(x,y,eps):
    # adds Laplacian random noise to a point in R^2
    theta = random.uniform(0,2*math.pi)
    p = random.random()
    r = np.real(-(1/eps)*(1+scipy.special.lambertw((p-1)/eps,-1)))
    x += r*math.cos(theta)
    y += r*math.sin(theta)
    return((x,y))

def locMap(Q):
    # creates a mapping from integer states to coordinate positions
    # eg for l = locMap(Q6), l[1] == (1,0); l[44] == (4,4)
    n = np.shape(Q)[0]
    if n != np.shape(Q)[1]: return False
    if abs(math.sqrt(n) % 1) > 1e-8: return False 
    m = round(math.sqrt(n))
    locList = [(0,0)]*n
    for i in range(0,n):
        x = i % m
        y = (i // (m))
        locList[i] = (x,y)
    return locList

def mapToInt(x,y,n):
    # pseudo-inverse of locMap, converts coordinates to their integer states
    # coordinate grid must be square of width n
    # eg maptoInt(4,4,10) = 44
    return int(y*n + x)

def pathLBS(path,l,eps,Q,k,ll,S): # Q,k,ll,S included for compatibility with pathLBS2
    # adds Polar Laplacian random noise to each point in the path independently
    try:
        path = [int(x) for x in path.split()]
    except (ValueError,TypeError):
        path = [int(x) for x in UCToPath(path).split()]
    pathCoords = [l[xy] for xy in path]
    
    hiddenCoords = [GeoLBS(*xy,eps) for xy in pathCoords]
    return hiddenCoords

def pathLBSFixed(path,l,eps,Q,k,ll,S): # Q,k,ll,S included for compatibility with pathLBS2
    # generates one Polar Laplacian variate and adds it to all points in the path
    try:
        path = [int(x) for x in path.split()]
    except (ValueError,TypeError):
        path = [int(x) for x in UCToPath(path).split()]
    pathCoords = [l[xy] for xy in path]
    
    theta = random.uniform(0,2*math.pi)
    p = random.random()
    r = np.real(-(1/eps)*(1+scipy.special.lambertw((p-1)/eps,-1)))
    hiddenCoords = [(x+r*math.cos(theta),y+r*math.sin(theta)) for (x,y) in pathCoords]
    return coordSnap(hiddenCoords,max([xy[0] for xy in l])+1)

def pathLBS2(path,l,eps,Q,k=-1,ll=-1,S=-1,verbose=False):
    
    # set all required variables
    # can be set externally to avoid needlessly recreating them
    try:
        if ll==-1: ll = toLinkedList(Q)   
    except TypeError: pass # it is late in the project and there are workarounds galore
    n = max([xy[0] for xy in l])+1
    try: 
        if S==-1: S = getJump(Q)
    except ValueError: pass
    
    # generate one Polar Laplacian variate and add it to all points in the path
    hiddenCoords = pathLBSFixed(path,l,eps,Q,k,ll,S)
    hpUC = pathToUC(" ".join([str(mapToInt(*xy,n)) for xy in hiddenCoords]))
    
    # find the most probable path within levenshtein distance k of hpUC
    bestProb = 0
    if k==-1: k = max(1,round(len(hpUC)/3))
    if verbose: count=0
    for j in range(0,k+1):
        for x in find_all_matches(hpUC,j,lookup_path,ll):
            if verbose: count +=1
            if verbose and count%1000==0: print(str(count))
            x = [int(s) for s in UCToPath(x).split()]
            p = int(all([x[i] in range(0,n**2) for i in range(0,len(x))]))
            if not p: continue
            for i in range(0,len(x)-1):
                if x[i+1] not in ll[x[i]]: 
                    p = 0
                    break
                else: p *= S[x[i],x[i+1]]
            if not p: continue
            if p > bestProb:
                bestProb = p
                hiddenCoords = [l[s] for s in x]
        
    return hiddenCoords
    
def coordSnap(coordPath,n):
    # snap each point in its path to the nearest point on the coordinate grid
    # coordinate grid must be square with width n
    coordPath = [[round(xy[0]),round(xy[1])] for xy in coordPath]
    for i in range(0,len(coordPath)): 
        if coordPath[i][0] < 0: coordPath[i][0] = 0
        if coordPath[i][0] > n-1: coordPath[i][0] = n-1
        if coordPath[i][1] < 0: coordPath[i][1] = 0
        if coordPath[i][1] > n-1: coordPath[i][1] = n-1
    return coordPath

def geoAttack(hiddenCoords,Q,eps,l=-1,ll=-1,nreturn=-1,k=-1,sort=True,returnlength=-1):
    # finds the most probable path at edit distance k from the hidden path
    
    # set the location map and linked list if not already set
    if l==-1: l = locMap(Q)
    if ll==-1: ll = toLinkedList(Q)
    # n: width of location map
    n = max([xy[0] for xy in l])+1
    # path processing
    path = [str(mapToInt(x,y,n)) for [x,y] in coordSnap(hiddenCoords,n)]
    path = pathToUC(" ".join(path))
    # k: max edit distance to consider
    # if k not set, choose a sensible value
    if k==-1: k = max(round(len(path)/2),1)
    # other return lengths currently don't work as
    # probability calculations for these still yet to be implemented
    if returnlength==-1: returnlength = len(hiddenCoords)
    
    
    # find all possible paths within edit distance k of the released path
    cp = [x for x in find_all_matches(path,k,lookup_path,ll) if len(x)==returnlength]
    # sort by probability
    S = getJump(Q)
    cpp = {x:pathProb(x,hiddenCoords,S,l,eps) for x in cp}
    cpp = sorted(cpp.items(),key = lambda d: d[1],reverse=True)
    # release the nreturn most probable close paths, or all if nreturn not set
    if nreturn==-1: return cpp
    return cpp[0:nreturn]


    





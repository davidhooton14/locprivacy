# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:28:49 2019

@author: David Hooton
"""

import math, numpy as np

def getJump(Q):
    S = np.zeros(np.shape(Q))
    for i in range(0,np.shape(Q)[1]):
        if(Q[i,i]>=0):  S[i,i] = 1
        else:  
            S[i,:] = -Q[i,:]/Q[i,i]
            S[i,i] = 0
    return S
    
def toLinkedList(Q):
    ll = [[]]*np.shape(Q)[0]
    for i in range(0,np.shape(Q)[0]):
        ll[i] = np.arange(0,np.shape(Q)[1])[Q[i,] > 0]
    return ll

def pathsByProbs(pps,N):
    pps2 = {}
    if N == "all":
        pps2 = pps
    else:
        for k,v in pps.items():
            if getattr(v,"d","") == N: pps2[k] = v
    pps2 = sorted(pps2.items(),key = lambda d: d[1].p,reverse=True)
    return pps2

def pathsByProbs2(pps,N):
    return sorted(pps.items(),key = lambda d: d[1],reverse=True)

def pathToUC(st):
    if st == "root": return "root"
    tmpstr = st.split(" ")
    tmpstr = [chr(int(x)+65) for x in tmpstr]
    return "".join(tmpstr)

def UCToPath(st):
    if st == "root": return "root"
    tmpstr = [str(ord(x)-65) for x in st]
    return " ".join(tmpstr)

def isValidPath(path,Q):
    # needs integer path
    n = np.shape(Q)[0]
    path = [int(x) for x in path.split()]
    if len(path) == 1: return path[0] in range(0,n)
    for i in range(0,len(path)-1):
        if path[i] in range(0,n) and path[i+1] in range(0,n):
            if Q[path[i],path[i+1]] <= 0:
                return False
        else: return False
    return True

def pathProb(path,hp,S,l,eps):
    if not isValidPath(UCToPath(path),S): return 0
    path = [int(x) for x in UCToPath(path).split()]
    p = 1
    for i in range(0,len(path)):
        if i < len(path)-1: p *= S[path[i],path[i+1]]
        (xh,yh) = hp[i]
        (xp,yp) = l[path[i]]
        d = math.sqrt((xh-xp)**2 + (yh-yp)**2)
        p *= math.exp(-d*eps)*eps**2/(2*math.pi)
    return p

def pathProb2(path,Q):
    S = getJump(Q)
    path = [int(x) for x in UCToPath(path).split()]
    p = 1
    for i in range(0,len(path)-1):
        p *= S[path[i],path[i+1]]
    return p

def avgPathDist(path1,path2):
    # requires coordinate paths
    if not len(path1) == len(path2): return False
    d = 0
    for i in range(0,len(path1)):
        (x1,y1),(x2,y2) = path1[i],path2[i]
        d += math.sqrt((x1-x2)**2 + (y1-y2)**2)
    return d/len(path1)

def simulateChain(Q,N,outstyle="list",a=None):
    n = np.shape(Q)[0]
    S = getJump(Q)
    # if no initial probabilities given, start state is uniformly random
    if a is None: a = [1/n]*n
    state = np.random.choice(range(0,len(a)),1,p=a)[0]
    outpath = [state]
    for i in range(0,N-1):
        state = np.random.choice(range(0,n),1,p=S[outpath[i],])[0]
        outpath.append(state)
    if outstyle=="list": return outpath
    outpath = " ".join([str(x) for x in outpath])
    if outstyle=="numeric": return outpath
    if outstyle=="UC": return pathToUC(outpath)
    return "invalid outstyle"




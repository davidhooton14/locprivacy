# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:24:51 2019

@author: david
"""

import math, random, itertools,numpy as np, pandas, sys, statistics as st
from datetime import datetime
from paths import *
from geo import *
from levenshtein import *

random.seed(444)

# example matrices
Q6 = np.zeros([100,100])
for i in range(0,100):
    Q6[i,i] = -4
    if i%10 == 9:
        Q6[i,i] = -3
    else: Q6[i,i+1] = 1
    if i%10 == 0:
        Q6[i,i] = -3
    else: Q6[i,i-1] = 1
    if i >= 90:
        Q6[i,i] = -3
    else: Q6[i,i+10] = 1
    if i <= 9:
        Q6[i,i] = -3
    else: Q6[i,i-10] = 1
Q6[0,0] = Q6[9,9] = Q6[90,90] = Q6[99,99] = -2

Q7 = np.zeros([100,100])
for i in range(0,100):
    Q7[i,i] = -6
    if i%10 == 9:
        Q7[i,i] = -4
    else: Q7[i,i+1] = 2
    if i%10 == 0:
        Q7[i,i] = -4
    else: Q7[i,i-1] = 2
    if i >= 90:
        Q7[i,i] = -5
    else: Q7[i,i+10] = 1
    if i <= 9:
        Q7[i,i] = -5
    else: Q7[i,i-10] = 1
Q7[0,0] = Q7[9,9] = Q7[90,90] = Q7[99,99] = -3

QRome = pandas.read_csv("data/QRome.csv").to_numpy()[:,1:]

def experiment(Q,_eps,_pathlength,_niters,LBS,filewrite=True,filePrefix="",truePath=None,verbose=True):
    start = datetime.now()
    # create location map, linked list, jump matrix from Q - the rates matrix
    l = locMap(Q)
    ll = toLinkedList(Q)
    S = getJump(Q)
    # Q must be square
    n = int(math.sqrt(np.shape(Q)[0]))
    
    # uniform distribution on the non-absorbing states
    a = (QRome.diagonal()!=1)/sum(QRome.diagonal()!=1) 
    # simulate path as a random walk throughout the network
    if truePath==None: truePath = simulateChain(Q,_pathlength,"numeric",a)
    if verbose: print("\n",truePath)
    
    # generate _niters paths hidden by geo-indistinguishability
    hiddenPaths = []        
    for i in range(0,_niters): hiddenPaths.append(LBS(truePath,l,_eps,Q,-1,ll,S))
    
    approxPaths = []
    closePathsAll = []
    for hp in hiddenPaths[:]:
        hpstart = datetime.now()
        hpUC = pathToUC(" ".join([str(mapToInt(x,y,n)) for [x,y] in coordSnap(hp,n)]))
        
        # get all paths within the required distance
        k = 0
        if verbose: print(UCToPath(hpUC))
        closePaths = geoAttack(hp,Q,_eps,l,ll,-1,k)
        # increase edit distance until we get a "reasonable" selection of paths
        while (len(closePaths)<500 or k < len(hpUC)/2) and k < len(hpUC)-2:
            k += 1
            if verbose: print(k)
            closePaths = closePaths + geoAttack(hp,Q,_eps,l,ll,-1,k)
        # if no paths were found even with a large edit distance, continue loop
        if len(closePaths)==0:
            if verbose: print("problematic hidden path")
            hiddenPaths.remove(hp)
            continue
        
        # sort paths by likelihood, select the best ones
        if verbose: print("sorting")
        closePaths.sort(key=lambda x:x[1],reverse=True)
        closePathsAll.append(closePaths)
        approxPaths.append(closePaths[0])
        if verbose: print(datetime.now()-hpstart)
    
    # compute the distance from the true path to the hidden paths and our approximations
    trueCoords = [l[int(x)] for x in truePath.split()]
    approxCoords = [[l[int(x)] for x in UCToPath(ap[0]).split()] for ap in approxPaths]
    approxDists = [avgPathDist(trueCoords,ac) for ac in approxCoords]
    origDists = [avgPathDist(trueCoords,hp) for hp in hiddenPaths]
    
    # write output
    if filewrite:
        out1 = open("out/rawdists/"+filePrefix+"_eps"+str(_eps)+"_len"+str(_pathlength)+".csv","w+")
        out1.write("truePath, origDist, approxDist, hiddenPath, approxPath\n")
        for i in range(0,len(origDists)):
            out1.write(truePath+", "+str(origDists[i])+", "+str(approxDists[i])+", "+str(hiddenPaths[i])+", "+UCToPath(approxPaths[i][0])+"\n")
        out1.close()
     
    if filewrite and len(approxDists)>0:
        out2 = open("out/summary2.csv","a+")
        out2.write(filePrefix+", "+str(_eps)+", "+str(_pathlength)+", "+str(_niters)+", "+str(st.mean(approxDists))+", "+str(st.mean(origDists))+", "+truePath+", "+str(datetime.now()-start)+"\n")
        out2.close()
    
if __name__ == "__main__":
    niters = int(sys.argv[1])
    LBSs = [pathLBS,pathLBSFixed,pathLBS2]
    
    out2 = open("out/summary2.csv","w+")
    out2.write("Q_LBS, eps, pathLength, niters, approxDist, origDist, truePath, itertime\n")
    out2.close()
    
#    for eps,pathLen,i in itertools.product([0.3,0.4,0.5,0.6,0.8,1],range(5,12,2),range(0,3)):
#        experiment(Q6,eps,pathLen,niters,LBSs[i],True,"Q6_LBS"+str(i))
#    for eps,pathLen,i in itertools.product([0.3,0.4,0.5,0.6,0.8,1],range(5,12,2),range(0,3)):
#        experiment(Q7,eps,pathLen,niters,LBSs[i],True,"Q7_LBS"+str(i))
    a = (QRome.diagonal()!=1)/sum(QRome.diagonal()!=1) 
    for pathLen, _ in itertools.product(range(6,7),range(0,3)):
        truePathR = simulateChain(QRome,pathLen,"numeric",a)
        for eps,i in itertools.product([0.3,0.4,0.5,0.6,0.8,1],range(0,3)):
            experiment(QRome,eps,pathLen,niters,LBSs[i],True,"QRome_LBS"+str(i),truePathR)
        

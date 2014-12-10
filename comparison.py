
import igraph
import time
import random
from igraph import *
import copy
import math
#from scipy.sparse import *
import time
import numpy as np
#from reachabilityMatrix import reachability
from igraph._igraph import OUT
from numpy.random import sample
from shortestpath import reachabilityMulti, transformLinkWeight
from shortestpath import dijkstra, floydWarshall
import getopt,sys    
from datetime import datetime
from sys import getsizeof
    
T=0.001 

def reachabilitySample(seed,g,N,rSample, ind, rounds):
    
    
    for i in xrange(rounds):
        gm=g.copy()
        edgeList=gm.es()
        probList=np.random.random(len(edgeList))
        removeList=[]
        for i in xrange(len(edgeList)):
            if probList[i] > edgeList[i]["weight"]:
                removeList.append(i)
        
        gm.delete_edges(removeList)
        ret = gm.subcomponent(seed, OUT)
        for i in xrange(N):
            if i in ret and i!=seed:
                rSample[ind,i]+=1
    
    rSample[ind,]=1.0*rSample[ind,]/rounds
    
    
def MIP(seed, g, N, rMIP, ind):
    #AM = g.get_adjacency()
    (distance, previous)=dijkstra(g, seed, None, 1)
    
    P={}
    for i in xrange(N):
        if seed==i:
            P[i]=1
        else:
            P[i]=0
        
    step=0
    notFinished=True
    while notFinished==True:
        step+=1
        notFinished=False
        for i in xrange(N):
            if i!=seed and step<=distance[i]+1:
                P[i]=0
                inNeighbors=g.neighbors(i,"in")
                for node in inNeighbors:
                    P[i] = P[i] + P[node]*g.es[g.get_eid(node,i)]["weight"]
                notFinished=True
                
    for i in xrange(N):
        if seed!=i:
            rMIP[ind,i]=P[i]
        
def ISP(seed, g, N, rISP, ind):
    #AM = g.get_adjacency()
    (distance, previous)=dijkstra(g, seed, None, 1)
    
    P={}
    for i in xrange(N):
        if seed==i:
            P[i]=1
        else:
            P[i]=0
        
    step=0
    notFinished=True
    while notFinished==True:
        step+=1
        notFinished=False
        for i in xrange(N):
            if i!=seed and step<=distance[i]+1:
                P[i]=0
                inNeighbors=g.neighbors(i,"in")
                for node in inNeighbors:
                    P[i] = 1- (1-P[i])*(1-P[node]*g.es[g.get_eid(node,i)]["weight"])
                notFinished=True
                
    for i in xrange(N):
        if seed!=i:
            rISP[ind,i]=P[i]    
    
    
    
def DLI(seed, g, N, rDLI, vertexName, inNeighbors):
    d=0.85
    P={}
    for i in xrange(N):
        P[i]=0
    
    found=False
    while found==False:
        prev_P = dict(P)
        for i in xrange(N):
            if vertexName[i]==str(seed):
                P[i]=1
            else:
                #neighbors = g.neighbors(vertexName[i], "in")
                neighbors = inNeighbors[i]
                P[i]=0
                for n in neighbors:
                    P[i] = P[i] + P[n]*g.es[g.get_eid(n,i)]["weight"]
        
        diff = math.sqrt(sum((P[k] - prev_P[k])**2 for k in P.keys()))
        if diff<T:
            found=True
    #print P
    for i in xrange(N):
        if seed!=i:
            rDLI[seed,i]=P[i] 
    

def SP1M(seed, g, N, rSP1M, ind):
    #AM = g.get_adjacency()
    
    (distance, previous)=dijkstra(g, seed, None, 1)
    
    P={}
    for i in xrange(N):
        P[i]={}
        P[i][distance[i]]=0
        P[i][distance[i]+1]=0
        
    P[seed][0]=1
    parents=set()
    parents.add(seed)
    step=0
    while len(parents)>0:
        step+=1
        tmp=set()
        plus1=set()
        for p in parents:
            children = g.neighbors(p, "out")
            for c in children:
                if step in P[c]:
                    tmp.add(c)
                    if step-1 in P[c]:
                        plus1.add(c)
                        
                    P[c][step]=1-(1-P[c][step])*(1-g.es[g.get_eid(p,c)]["weight"]*P[p][step-1])
        
        #for c in plus1:
        #    P[c][step]=1-(1-P[c][step-1])*(P[c][step])
         
        parents=tmp
        
    for i in xrange(N):
        if seed!=i:
            for k in P[i]:
                rSP1M[ind,i]=1-(1-rSP1M[ind,i])*(1-P[i][k])
    
        

def MIA(g, N, rMIA, theta):
    MIIA={}
    tmp=g.copy()
    transformLinkWeight(tmp)
    #WM=np.asmatrix(tmp.get_adjacency(attribute="weight").data)
    
    for i in xrange(N):
        MIIA[i]=set()
    
    for i in xrange(N):

        (distances, previous)=dijkstra(g, i, None, 2)    
        for j in xrange(N):
            if j!=i:
                if distances[j]<-math.log(theta):
                    ind=j
                    while ind!=i:
                        MIIA[j].add(tmp.get_eid(previous[ind], ind))
                        ind=previous[ind]
    
    return MIIA


def findAP(P, t, g, edges, vertexName):
    ret=0
    inNeighbors=g.neighbors(vertexName[t], "in")
    for node in inNeighbors:
        if P[node]==1:
            ap=1
        elif g.get_eid(node,t) in edges:
            ap=findAP(P,node,g,edges,vertexName)
        else:
            ap=0
        ret = 1-(1-ret)*(1-ap*g.es[g.get_eid(node,t)]["weight"])
    return ret
    

def reachabilityMIA(seed, g, N, MIIA, rMIA, ind, vertexName):
    P={}
    for i in xrange(N):
        if i==seed:
            P[i]=1
        else:
            P[i]=0
            
    for i in xrange(N):
        if seed != i:
            tmp=P.copy()
            P[i]=findAP(tmp,i,g,MIIA[i], vertexName)
        
    for i in xrange(N):
        if i!=seed:
            rMIA[ind, i]=P[i]
    

def SSS(seed, g, N, vertexName, inNeighbors, rSSS, ind):
    P={}
    for i in xrange(N):
        P[i]=0
    
    found=False
    while found==False:
        prev_P = dict(P)
        for i in xrange(N):
            if vertexName[i]==str(seed):
                P[i]=1
            else:
                #neighbors = g.neighbors(vertexName[i], "in")
                neighbors = inNeighbors[i]
                pt=1
                for n in neighbors:
                    pt = pt*(1-P[n]*g.es[g.get_eid(n,i)]["weight"])
                P[i] = 1-pt
        
        #print P
        diff = math.sqrt(sum((P[k] - prev_P[k])**2 for k in P.keys()))
        if diff<T:
            found=True
    #print P
    for i in xrange(N):
        if seed!=i:
            rSSS[ind,i]=P[i]
    
    
    
    
def main(argv):    
    optlist,args=getopt.getopt(argv,'')
    #g=Graph(directed=True)
#     for i in xrange(3):
#         g.add_vertices(str(i))
#     g.add_edges([(0,1), (1,2), (2,1)])
#     g.es[0]["weight"]=0.5
#     g.es[1]["weight"]=0.5
#     g.es[2]["weight"]=0.5
    
#     for i in xrange(4):
#         g.add_vertices(str(i))
#     g.add_edges([(0,1), (1,2), (2,3), (3,1), (2,0), (3,0)])
#     g.es[0]["weight"]=0.5
#     g.es[1]["weight"]=0.4
#     g.es[2]["weight"]=0.6
#     g.es[3]["weight"]=0.2
#     g.es[4]["weight"]=0.3
#     g.es[5]["weight"]=0.1

#     for i in xrange(5):
#         g.add_vertices(str(i))
#     g.add_edges([(0,1), (0,2), (1,2), (0,3), (1,4), (3,4), (4,0), (4,2)])
#     g.es[0]["weight"]=0.5
#     g.es[1]["weight"]=0.4
#     g.es[2]["weight"]=0.6
#     g.es[3]["weight"]=0.2
#     g.es[4]["weight"]=0.3
#     g.es[5]["weight"]=0.1
#     g.es[6]["weight"]=0.3
#     g.es[7]["weight"]=0.2

#     for i in xrange(7):
#         g.add_vertices(str(i))
#     g.add_edges([(0,1), (1,2), (2,3), (3,1), (1,6), (0,4), (4,5), (5,6), (6,0)])
#     g.es[0]["weight"]=0.5
#     g.es[1]["weight"]=0.4
#     g.es[2]["weight"]=0.6
#     g.es[3]["weight"]=0.2
#     g.es[4]["weight"]=0.3
#     g.es[5]["weight"]=0.1
#     g.es[6]["weight"]=0.3
#     g.es[7]["weight"]=0.2
#     g.es[8]["weight"]=0.2


#     for i in xrange(10):
#         g.add_vertices(str(i))
#     g.add_edges([(0,5), (0,8), (1,5), (1,8), (1,9), (3,6), (4,0), (4,6), (4,8), (5,7), (5,8), (5,9), (5,3), \
#                  (6,0), (6,9), (7,4), (7,8), (8,5), (8,2), (9,5)])

    #g = Graph.Read_Lgl("/Users/biru/workspace/Reachability/slashdot.lgl")    
    #N=g.vcount()
    #E=g.ecount()
    N = int(args[0])
    if N==1: # digg network
        g = Graph.Read_Lgl("/Users/biru/workspace/Reachability/digg.lgl")
        N=g.vcount()
        E=g.ecount()
    elif N==2:
        g = Graph.Read_Lgl("/Users/biru/workspace/Reachability/slashdot.lgl")
        N=g.vcount()
        E=g.ecount()
    else:
        E=N*8
        g = Graph.Static_Power_Law(N,E,2.45,2.45)    
    
    print N, E
        
    
    link_weight = float(args[1])
    
    
    for i in xrange(N):
        g.vs[i]["name"]=str(i)
    for i in xrange(E):
        g.es[i]["weight"]=np.random.uniform(0,link_weight)
    
#     WM=np.asmatrix(g.get_adjacency(attribute="weight").data)
#     col_sums = WM.transpose().sum(axis=1)
#     for i in xrange(N):
#         WM[:,i] = WM[:,i]/col_sums[i]
    
#     for i in xrange(N):
#         for j in xrange(N):
#             if WM[i,j]!=0:
#                 g.es[g.get_eid(i,j)]["weight"]=WM[i,j]
    
    
    #print g
    #print "graph done"    
    
        
    vertexName={}
    vertexIndex={}
    for i in xrange(N):
        vertexName[i]=g.vs[i]["name"]
        vertexIndex[g.vs[i]["name"]]=i
        
    inNeighbors={}
    for i in xrange(N):
        inNeighbors[i] = g.neighbors(vertexName[i], "in")

    seeds=set()
    L=min(10,N)
    while len(seeds)<L:
        t = random.randint(0,N-1)
        seeds.add(t)
    
    rMIP=np.zeros([L, N])
    tMIP=0
    rMIA=np.zeros([L,N])
    tMIA=0
    rSP1M=np.zeros([L, N])
    tSP1M=0
    rSSS=np.zeros([L, N])
    tSSS=0
    rISP=np.zeros([L, N])
    tISP=0
    rSample=np.zeros([L, N])
    tSample=0
    rDLI=np.zeros([L, N])
    tDLI=0
    


    
    #print "start time", (time.strftime("%H:%M:%S"))
    tSample = datetime.now()
    ind=0
    for seed in seeds:
        reachabilitySample(seed,g,N,rSample, ind, 10000)
        ind+=1
        #print seed, (time.strftime("%H:%M:%S"))
    tSample = datetime.now() - tSample    
    #print "sample done", (time.strftime("%H:%M:%S"))
 
    
    tMIA = datetime.now()
    ind=0
    MIIA=MIA(g, N, rMIA, T)
    print "get MIIA"
    for seed in seeds:
        reachabilityMIA(seed,g,N,MIIA,rMIA, ind, vertexName)
        ind+=1
    tMIA = datetime.now() - tMIA
        #print seed, (time.strftime("%H:%M:%S"))
    #print "MIA done", (time.strftime("%H:%M:%S"))
     
    tSP1M = datetime.now()
    ind=0
    for seed in seeds:
        SP1M(seed,g,N,rSP1M, ind)
        ind+=1
        #print seed, (time.strftime("%H:%M:%S"))
    tSP1M = datetime.now() -tSP1M
    #print "SP1M done", (time.strftime("%H:%M:%S"))
    
    tMIP = datetime.now()
    ind=0
    for seed in seeds:
        MIP(seed,g,N,rMIP, ind)
        ind+=1
        #print seed, (time.strftime("%H:%M:%S"))
    tMIP = datetime.now() - tMIP
    #print "SSS done", (time.strftime("%H:%M:%S"))
    
#     for seed in xrange(N):
#         DLI(seed,g,N,rDLI, vertexName, inNeighbors)
#         print seed, (time.strftime("%H:%M:%S"))
#     print "DLI done", (time.strftime("%H:%M:%S")) 
     
    tSSS = datetime.now()
    ind=0
    for seed in seeds:
        SSS(seed, g, N, vertexName, inNeighbors, rSSS, ind)
        ind+=1
        #print seed, (time.strftime("%H:%M:%S"))
    tSSS = datetime.now() - tSSS
    #print "iterate done", (time.strftime("%H:%M:%S"))
    
    tISP = datetime.now()
    ind=0
    for seed in seeds:
        ISP(seed, g, N, rISP, ind)
        ind+=1
        #print seed, (time.strftime("%H:%M:%S"))
    tISP = datetime.now() - tISP
    #print "DINDP1 done", (time.strftime("%H:%M:%S"))
    
    
    
#     for seed in xrange(N):
#         reachabilityIterate(seed, g, N, vertexName, inNeighbors, rGS, 1)
#         print seed, (time.strftime("%H:%M:%S"))
#     print "iterate done", (time.strftime("%H:%M:%S"))
            
#     for seed in xrange(N):    
#         reachabilityLinear(seed, g, N, rLinear)
#         print seed, (time.strftime("%H:%M:%S"))
#     print "linear done", (time.strftime("%H:%M:%S"))
    

    
    #print "Sample:", rSample    
    #print "MIA:", rMIA
    #print "SP1M:", rSP1M
    #print "SSS:", rSSS
    #print "DLI:", rDLI
    #print "INDP:", rINDP
    #print "INDP2:", rINDP2
    #print "Matrix:", rMatrix
    
    diffs=[np.linalg.norm(rSample-rMIA), np.linalg.norm(rSample-rSP1M), np.linalg.norm(rSample-rMIP), np.linalg.norm(rSample-rSSS), np.linalg.norm(rSample-rISP)]
#     for i in xrange(N):
#         print i, np.linalg.norm(rSample[i,]-rMatrix[i,])
    print diffs
    print tSample, tMIA, tSP1M, tMIP, tSSS, tISP
    #print reachability(0, M, N), reachability(1, M, N), reachability(2, M, N), reachability(3, M, N)#, reachability(4, M, N) 

if __name__=="__main__":
    main(sys.argv[1:])    

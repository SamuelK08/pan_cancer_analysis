# -*- coding: utf-8 -*-
"""
Created on Thu May 22 15:03:58 2025

@author: Eli
"""

import rpy2
import os
os.environ['R_HOME'] = "C:/Program Files/R/R-4.3.3"
import rpy2.robjects as ro
import networkx as nx
import pandas as pd
import numpy as np
import itertools as itr
import condor
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import random
random.seed(5151)
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors2 = list(mcolors.XKCD_COLORS.values())
random.shuffle(colors2)
colors += colors2




rSource = ro.r('''
               source('analyzeModuleFunctions.R')
               ''')
r_getImportantPaths = ro.globalenv["getImportantPaths"]

def getModules(cancerType):
    G = nx.read_edgelist(cancerType+"_EdgeList2.txt")
    LCC = sorted(list(nx.connected_components(G)),key=len,reverse=True)[0]
    G = nx.subgraph(G, LCC)
    mRNAs = []
    miRs = []
    for n in G:
        if("hsa" in n):
            miRs.append(n)
        else:
            mRNAs.append(n)
    edgeList = []
    for miR in miRs:
        for mRNA in G.neighbors(miR):
            edgeList.append([miR,mRNA])
    edgeList = pd.DataFrame(edgeList)
    import random
    random.seed(0)
    co = condor.condor_object(dataframe=edgeList,silent=True)
    co.initial_community()
    co.brim()
    
    modDict = {}
    for i in range(len(co.tar_memb)):
        if(co.tar_memb.iloc[i]["community"] not in modDict):
            modDict[co.tar_memb.iloc[i]["community"]] = [co.tar_memb.iloc[i]["tar"][4:]]
        else:
            modDict[co.tar_memb.iloc[i]["community"]].append(co.tar_memb.iloc[i]["tar"][4:])
    for i in range(len(co.reg_memb)):
        if(co.reg_memb.iloc[i]["community"] not in modDict):
            modDict[co.reg_memb.iloc[i]["community"]] = [co.reg_memb.iloc[i]["reg"][4:]]
        else:
            modDict[co.reg_memb.iloc[i]["community"]].append(co.reg_memb.iloc[i]["reg"][4:])
    modDict = {k:modDict[k] for k in sorted(modDict)}

    f = open(cancerType+"_Modules_Bipartite.txt",'w')
    for m in modDict.values():
        temp = sorted(m)
        string = ", ".join(temp)
        string += "\n"
        f.write(string)
    f.close()

def getModuleFunctions(cancerType):
    f = open(cancerType+"_Modules_Bipartite.txt","r")
    modules = f.readlines()
    f.close()
    modules = [m.rstrip().split(", " ) for m in modules]
    writer = pd.ExcelWriter(cancerType+"_moduleFunctions.xlsx",engine="xlsxwriter")
    for moduleNum in range(1,len(modules)+1):
        
        pathways =  pd.read_csv(cancerType+"_Module"+str(moduleNum)+"_pathwayAll_FisherResults.csv",index_col=0)
        if(len(pathways) == 0):
            resDF = pd.DataFrame(index=["Functions","Average FE"]).transpose()
            resDF.to_excel(writer,sheet_name="Module #"+str(moduleNum))
            continue
        relevantPaths = list(pathways[(pathways["FEs"] > 1) & (pathways["padj"] < 0.05)]["allPathNames.inDF."])
        
        uniquePaths = set(relevantPaths)
        
        for i in range(len(modules)):
            if(i != moduleNum-1):
                pathwaysOR2 =  pd.read_csv(cancerType+"_Module"+str(i+1)+"_pathwayAll_FisherResults.csv",index_col=0)
                if(len(pathwaysOR2) > 0):
                    relevantPaths2 = list(pathwaysOR2[(pathwaysOR2["FEs"] > 1) & (pathwaysOR2["padj"] < 0.05)]["allPathNames.inDF."])
                    uniquePaths = uniquePaths - set(relevantPaths2)
                
        uniquePaths = list(uniquePaths)
        uniquePathFEs = list(pathways.loc[[t in uniquePaths for t in pathways.loc[:,"allPathNames.inDF."]]]["FEs"])
        
        uniqueG = nx.Graph()
        for i in range(len(uniquePaths)):
            uniqueG.add_node(uniquePaths[i],weight=uniquePathFEs[i])
        
        
        for i in range(len(uniquePaths)):
            g1 = pathways[pathways["allPathNames.inDF."] == uniquePaths[i]]["pathwayGeneLists"].values[0]
            g1 = g1.split(', ')
            for j in range(i):
                g2 = pathways[pathways["allPathNames.inDF."] == uniquePaths[j]]["pathwayGeneLists"].values[0]
                g2 = g2.split(', ')

                numCommon = len(set(g1).intersection(set(g2)))
                if(numCommon > 0):
                    uniqueG.add_edge(uniquePaths[i],uniquePaths[j],weight=numCommon)
        
        comms = nx.community.greedy_modularity_communities(uniqueG,weight="weight")
        comms = [list(set(c)) for c in sorted(comms,key=len)]
        avgFEs = []
        for i in range(len(comms)):
            FEs = 0
            for n in comms[i]:
                FEs += uniqueG.nodes[n]['weight']
            avgFEs.append(FEs/len(comms[i]))
        
        resDF = pd.DataFrame([comms,avgFEs],index=["Functions","Average FE"]).transpose()
        resDF.to_excel(writer,sheet_name="Module #"+str(moduleNum))
    writer.close()
    
#cancerType = "BRCACommunity2"
#getModules(cancerType)
#r_getImportantPaths(cancerType)
#getModuleFunctions(cancerType)


for moduleNum in [14,15,16,5,8]:
    print(moduleNum)
    cancerType="BRCACommunity"+str(moduleNum)
    getModules(cancerType)
    r_getImportantPaths(cancerType)
    getModuleFunctions(cancerType)
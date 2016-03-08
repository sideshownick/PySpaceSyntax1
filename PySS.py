"""
PySS.py module 

requires: numpy matplotlib and networkx
* recommended to install all-in-one "anaconda python" distribution

Created on Fri Mar  4 16:06:49 2016

@author: Nick McCullen
"""
from __future__ import division

def getjunctions(ways):
    osmid_list=ways.keys()
    Nroads=len(ways)
    junctions=[]
    #need to have list of connections for each junction what directs to where...
    for i in range(Nroads-1):
        osmid = osmid_list[i]
        tag1, ref1 = ways[osmid]
        for osm2 in osmid_list[i+1:]:
            tag2, ref2 = ways[osm2]
            common_nodes = list(set(ref1).intersection(ref2))
            for junc in common_nodes:
                if junc not in junctions:
                    junctions.append(junc)
    return junctions


def make_segments(coords, ways, intersections, arcangle, onewayroads=False):
    from numpy import arctan, sign, rad2deg
    #define dictionary for nodes and their accessible segments 
    tempdict={} #use a temporary dict so as to trim it at the end
    axdict={} #dictionary to store axial lines {AID: [SEGIDs]}
    segdict={} #dictionary to store segment data {SID: [SEGENDS]}
    for j,osm_id in enumerate(ways):
        #refs are full list of waypoint references
        tags, refs = ways[osm_id]
        if onewayroads and 'oneway' in tags and tags['oneway']=='yes':
            oneway=True
        else: oneway=False
        #onlyconsider out-sements at this stage, start 0 for initial waypoint
        #segID refers to the one following on from current waypoint
        n=0
        axID=str(osm_id)+".A"+str(n)
        segID=str(osm_id)+".S"+str(n)
        axdict[axID]=[segID]
        
        if refs[0] in tempdict:
            tempdict[refs[0]].append(segID) #record where each node leads
        else: tempdict[refs[0]]=[segID]
        segID0=segID #remember segment ID before moving on
        
        totalangle=0.0
        
        segdict[segID]=[refs[0]] #start waypoint    
        
        x1, y1 = coords[refs[0]]
        x2, y2 = coords[refs[1]]

        if x2-x1!=0: oldangle=arctan((y2-y1)/(x2-x1))
        else: oldangle=90*sign(y2-y1)
        
        #start from waypoint 1 and check whether to segment or not 
        #iterate through the waypoint refs from zero (in case 0 is intersection)
        for i in range(0,len(refs)-1):
            #always looking forward in case initial point is a junction
            x2,y2 = coords[refs[i+1]]  
            x1,y1 = coords[refs[i]] 
            #set angle to initial orientation
            if x2-x1!=0: newangle=arctan((y2-y1)/(x2-x1))
            else: newangle=90*sign(y2-y1)
            totalangle+=newangle-oldangle
            oldangle=newangle

            if abs(rad2deg(totalangle)) >= arcangle:# and i<len(refs):   
                #put latest waypoint into the segment's list of points
                segdict[segID].append(refs[i]) 
                #in all cases can exit forward--->
                #remember previous id and increment the segment id
                segID0=segID
                n+=1
                segID=str(osm_id)+".S"+str(n) #new segment ID for following segment
                segdict[segID]=[refs[i]] #start next segment list with same point
                axID=str(osm_id)+".A"+str(n) #new axline ID for following axline
                axdict[axID]=[segID] #put the new segment in the new axline entry
                
                #put following seg into dictionary {junction: [segments]}
                if refs[i] in tempdict:
                    tempdict[refs[i]].append(segID)
                else: tempdict[refs[i]]=[segID]
                #two way streets can return backwards<---
                if not oneway:    
                    tempdict[refs[i]].append(segID0) #previous (current) segment
                    
                totalangle=0.0    #reset total angle to zero    

            elif refs[i] in intersections:# and i<len(refs):
                segdict[segID].append(refs[i]) #last waypoint for previous segment                
                #in all cases can exit forward--->
                #remember previous and increment the segment id
                segID0=segID
                n+=1
                segID=str(osm_id)+".S"+str(n) #following segment ID
                axdict[axID].append(segID) #put the next segment in the current axline entry
                segdict[segID]=[refs[i]] #start next segment list with same point
                
                if refs[i] in tempdict:
                    tempdict[refs[i]].append(segID)
                else: tempdict[refs[i]]=[segID]
                #two way streets can return backwards<---
                if not oneway:    
                    tempdict[refs[i]].append(segID0)
                #DO NOT INCREMENT AXIAL LINES HERE
        #check if last waypoint is an intersection
        i=-1
        segdict[segID].append(refs[i]) #last waypoint for previous segment
        if refs[i] in intersections:# and i<len(refs):            
            if refs[i] in tempdict:
                if not oneway:    
                    tempdict[refs[i]].append(segID)
                
        #make sure can return down segment if not one-way
        if refs[i] in tempdict:
            if not oneway: 
                tempdict[refs[i]].append(segID) #previous (current) segment
        else: 
            if not oneway:
                tempdict[refs[i]] = [segID] 
        #record dead-ends

    linkdict={}
    #remove single segment junctions by writing a new list from the temp one
    for link in tempdict:
        if len(tempdict[link])>=0:
            linkdict[link]=tempdict[link]
    return linkdict, axdict, segdict
    
def do_integration(mapname, axlines, segments, ways, links, onewayroads=False):
    import numpy as np
    import networkx as nx

    #output file, rename as appropriate
    outfile=open(mapname+"-integration.csv", "w")
    outfile.write("AxLineID, n, TotalDepth, MeanDepth, Integration\n")

    #initialise a networkx "directed graph" (network)
    G=nx.DiGraph()
    # this will be a bipartite graph with axial lines as start/end nodes 
    # and internal links between segments. 
    # Correct pathlengths between axlines comes from only giving a weighting to
    # links that connect different axial lines.

    #1. for each axline link axID (weight 0) to each segID within it, and
    #2. look at other accessible segments from junctions and weight them 1
    for AID in sorted(axlines):
        segIDs = axlines[AID] #lists the segments in each axline
        #2. for each segend in each axline 
        for SID in sorted(segIDs):
            G.add_edge(AID, SID, weight=0.0) #from subgraph of axline nodes to segments
            G.add_edge(SID, AID, weight=0.0) #back from segment subgraph to axlines
            
            ##NEED TO CHECK IF THE CURRENT SEGMENT IS ONE-WAY      
            wayID = int(SID.split(".")[0])
            tags, refs = ways[wayID]
            oneway=False
            if onewayroads and 'oneway' in tags and tags['oneway']=='yes':
                oneway=True
            
            #(i) look where the segment leads
            end1, end2 = segments[SID] #get the segment ends
            startlinks = links[end1] #segments reachable from the current seg start point
            endlinks = links[end2] #segments reachable from the current seg end point
            # * Note: if onewayroads=True in make_segments() then 
            #  only accessible segments are considered already
            
            #in any case can exit from the end waypoint:
            for toSID in endlinks:
                wt=0.0 
                # SEE IF toSID IS NOT IN THE SAME AXLINE
                if toSID not in segIDs:
                    wt=1.0
                G.add_edge(SID, toSID, weight=wt)
                #return edges are considered in their own axial line routine 
            
            #if it's oneway can not exit via the start waypoint links
            if not oneway:
                for toSID in startlinks:
                    wt=0.0 
                    # SEE IF toSID IS NOT IN THE SAME AXLINE
                    if toSID not in segIDs:
                        wt=1.0
                    G.add_edge(SID, toSID, weight=wt)
                    #return edges are considered in their own axial line routine 
                    

    #MAIN PART OF THE SPACE SYNTAX CALCULATION FOLLOWS
    
    #calculate all shortest (weighted) routes between all segments
    SPL=nx.shortest_path_length(G, weight="weight")
    # * returns a dictionary with the origin as the key
    #   and a destion dictionary as the given entry as a key
    # each dest dict entry is the shortest path length

    #thin out keys to only include axial lines as calculation nodes
    segIDs = SPL.keys()
    axialIDs = []
    for ID in segIDs:
        tag=ID.split(".")
        if tag[1][0]=="A":
            axialIDs.append(ID)

    totalpaths=len(axialIDs)
    #calculate Ingegration
    for orig in axialIDs:
        TD = 0.0 #start total depth at zero
        n = 1 #start n at 1 for self only
        paths = SPL[orig]
        for dest in axialIDs:
            if dest in paths and dest!=orig:
                Length = SPL[orig][dest]
                TD += Length #add to total depth
                n += 1 #increment n
        if n>1: MD=TD/(n-1.0) 
        else: MD=np.inf
        if n>1 and TD>n-1:
            MD=TD/(n-1.0) 
            RA=(MD-1.0)#*2.0/(n-2.0)
            DN=(n*(np.log2((n+2.0)/3.0)-1.0)+1.0)/(n-1.0)#*2.0/(n-2.0)
            RRA=RA/DN
            I=1.0/RRA /(totalpaths/(n-1)) #new divisor to account for fragmentation
        else: I=0

        outfile.write(",    ".join(map(str, [orig, n, TD, MD, I]))+"\n")
         
        axlines[orig] = {"segments":axlines[orig], "Integration": I}
    return axlines


def plotmap(mapname, axlines, segments, coords, figtype="png"):
    """plotmap(mapname, axlines, segments, coords)
    plot axial map with integration data as colour values
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig,axes=plt.subplots(figsize=(12,9))
    
    valary=[]
    for ax in axlines:
        I = axlines[ax]["Integration"]
        valary.append(I)

    cmap = cm.gnuplot
    norm=plt.Normalize(vmin=max(0, min(valary)), vmax=max(valary))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(valary)

    for i,a in enumerate(axlines):
        cval=float(axlines[a]["Integration"])
        for seg in axlines[a]["segments"]:
            n1, n2 = segments[seg]
            x1, y1 = coords[n1]
            x2, y2 = coords[n2]
            axes.plot([x1,x2], [y1,y2], color=sm.to_rgba(cval), linewidth=3)
    axes.set_xlabel("Longitude")
    axes.set_ylabel("Latitude")
    cb=plt.colorbar(sm, ax=axes)
    cb.set_label("Integration")

    plt.savefig(mapname+"-integration."+figtype)

        
def savedata(data, fname, delimiter="  ;  "):
    with open(fname, "w") as of:
        for key in data:
            values = [key]+map(str, data[key])
            of.write(delimiter.join(values) + "\n")
    
###for testing purposes###
#this test is executed only if the current file is run as a standalone  
if __name__=="__main__":
    import osmparser as osm
    mapname="whs"
    dl="  ;  "
    arcangle=45
    oneway=False
    restrict_to = "highway"
    exclude=['footway', 'steps', 'pedestrian', 'cycleway', 'path', 'service']

    osm.parsefile(mapname, delimiter=dl, restrict=restrict_to, excludeway=exclude)

    coords=osm.getcoords(mapname)
    ways=osm.getways(mapname)
    
    ###find all junctions and make list of IDs
    junctions=getjunctions(ways)
    links, axlines, segments = make_segments(coords, ways, junctions, arcangle, onewayroads=oneway)
    axlines=do_integration(mapname, axlines, segments, ways, links, onewayroads=oneway)
    plotmap(mapname, axlines, segments, coords)

    
    
    
        

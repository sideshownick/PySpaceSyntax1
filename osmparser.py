#parser.py osm openstreetmaps parser module
#osmid: openstreetmaps-ID
"""
Created on Mon Feb 29 17:22:32 2016

@author: Nick McCullen 
"""

mapname="whs"
restrict_to="highway"
exclude=['footway', 'steps', 'pedestrian', 'cycleway', 'path', 'service']

def parsefile(mapname, delimiter="  ;  ", restrict=False, excludeway=[]):
    """parser.parsefile(<MAPNAME>, delimiter="  ;  ")
    parses <MAPNAME>.osm xml file to extract data
    saves three csv files with delimiters "  ;  ": 
        <MAPNAME>-coords.csv
            - contains: osmid  ;  lon  ;  lat  ;  {attributes}
            - the last entry is a dictionary of attibutes if any exist
        <MAPNAME>-ways.csv
            - contains: osmid  ;  {attributes}  ;  [osmid1, osmid2,...]
            - third entry are the node IDs making up the way
        <MAPNAME>-relations.csv
            -contains: osmid  ;  {attributes}  ;  [ways]  ;  [nodes]
            -ways and nodes are lists of the way and node IDs in the relation
    optional arguments restrict="waytype" and excludeway=[] to narrow search
    """
             
    import xml.etree.ElementTree

    e = xml.etree.ElementTree.parse(mapname+'.osm').getroot()
    
    coordfilename = mapname+"-coords.csv"        
    with open(coordfilename, "w") as cfile:
        for coords in e.findall('node'):
            osmid = coords.get("id")
            lat = coords.get("lat")
            lon = coords.get("lon")
            attribs={}
            for obj in coords:
                att = obj.attrib
                if len(att) > 0:
                    key = att["k"]
                    val = att["v"]
                    attribs[key]=val
            cfile.write(delimiter.join([osmid,lon,lat,str(attribs)]) + "\n")
        
    wayfilename = mapname+"-ways.csv"
    with open(wayfilename, "w") as wfile:
        for way in e.findall('way'):
            osmid = way.get("id")
            attribs={}
            refs=[]
            for obj in way:
                att = obj.attrib
                if "k" in att: 
                    key = att["k"]
                    val = att["v"]
                    attribs[key]=val
                if "ref" in att:
                    refs.append(att["ref"])
            if restrict==False:
                wfile.write(delimiter.join([osmid,str(attribs),str(refs)]) + "\n")
            elif restrict in attribs and attribs[restrict] not in excludeway:
                wfile.write(delimiter.join([osmid,str(attribs),str(refs)]) + "\n")

    relfilename = mapname+"-relations.csv"
    with open(relfilename, "w") as rfile:
        for rel in e.findall('relation'):
            osmid = rel.get("id")
            attribs={}
            ways=[]
            nodes=[]
            for obj in rel:
                att = obj.attrib
                if "k" in att: 
                    key = att["k"]
                    val = att["v"]
                    attribs[key]=val
                if "type" in att and att["type"]=="way":
                    ways.append(att["ref"])
                if "type" in att and att["type"]=="node":
                    nodes.append(att["ref"])
            rfile.write(delimiter.join([osmid,str(attribs),str(ways),str(nodes)]) + "\n\n")




def getcoords(mapname, delimiter="  ;  "):
    """getcoords(<MAPNAME>, delimiter="  ;  ")
    get coordinate data from file <MAPNAME>-coords.csv
    returns: dictionary of coordinates in format: {'OSMID': (lon, lat)}
    """
    coords={}
    for line in file(mapname+"-coords.csv"):
        line=line.split(delimiter) 
        osmid, lon, lat, attribs = map(eval, line)
        coords[str(osmid)] = (lon, lat)
    return coords



def getways(mapname, delimiter="  ;  "):
    """getways(<MAPNAME>, delimiter="  ;  ")
    get data on ways from file
    returns: dictionary of ways in format: {"OSMID": [{attributes}, [refs]]}
    """
    ways={}
    for line in file(mapname+"-ways.csv"):
        line=line.split(delimiter) #split the line based on the delimiter dl
        osmid, attribs, refs = map(eval, line) #literal evaluation of each entry
        ways[osmid] = [attribs, refs]
    return ways





if __name__=="__main__":
    parsefile(mapname, restrict=restrict_to, excludeway=exclude)


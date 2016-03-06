### Process OpenStreetMaps map and generate SpaceSyntax Integration Maps and Data
## the modules osmparser.py and axmap.py need to be in the same working directory
import osmparser as osm
import PySS as ss


#######################################
### LIST OF PARAMETERS AND SETTINGS ###
##
## Base name for <NAME>.osm map file 
## and generated <NAME>[-ways,-coords,-relations].csv files
mapname="tinymap" 
##
## Delimiter for the .csv files 
## (default tries to be distinct from anything in any text fields)
dl="  ;  "
## * Ensure to specify "semicolon" (;) when loading these into spreadsheets
##
## Maximum angle a road-way can turn before splitting to a new Axial Line
arcangle=45
## 30 (degrees) seems about right
##
## Consider one-way roads or not (True, False)
oneway=True
#######################################

######################################################################
### SETTINGS FOR PARSING THE <NAME>.osm FILE
##
## FIRST RUN ON NEW .osm MAP SET THIS TO True TO GENERATE .csv FILES
## TO USE GENERATED <NAME>-ways.csv and <NAME>-coords.csv SET TO False
FirstRun=True
##
## restrict to "highway" 
restrict_to = "highway"
## * others are "waterway" etc, set to False to see all in .csv file (very expensive to compute below!)
## Don't inlcude this list of minor ways in case of traffic analysis
exclude=['footway', 'steps', 'pedestrian', 'cycleway', 'path', 'service']
#######################################################################


#####################################
## main functions are called below ##
#####################################


## parse osm file and output to csv text files with delimiter dl above
if FirstRun: 
    osm.parsefile(mapname, delimiter=dl, restrict=restrict_to, excludeway=exclude)

#get coordinate lookup dictionary and dictionary of ways from .csv files
coords=osm.getcoords(mapname)
ways=osm.getways(mapname)
    
###find all junctions and make list of IDs
junctions=ss.getjunctions(ways)


###create axial lines and append to end of entries in ways 
links, axlines, segments = ss.make_segments(coords, ways, junctions, arcangle)
ss.savedata(links, mapname+"-links.csv", delimiter=dl)
ss.savedata(axlines, mapname+"-axlines.csv", delimiter=dl)
ss.savedata(segments, mapname+"-segments.csv", delimiter=", ")

###calculate integration measure and save to file mapname-integration.csv
axlines=ss.do_integration(mapname, axlines, segments, ways, links, onewayroads=oneway)

###plot graph and save to file mapname-integration.png
ss.plotmap(mapname, axlines, segments, coords)

    
    

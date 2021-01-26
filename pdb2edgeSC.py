import sys
import os
import numpy
import math
from scipy import spatial
import re
import copy
import string
import time 

t0 = time.time()
#Pi-Cation Strength
#Arg-Tyr = -2.8 kCal/mol = -11.7 kJ/mol
#Arg-Phe = -3.1 kCal/mol = -13 kJ/mol
#Arg-Trp = -4.1 kCal/mol = -17.2 kJ/mol
#Lys-Tyr = -3.8 kCal/mol = -15.9 kJ/mol
#Lys-Phe = -2.7 kCal/mol = -11.3 kJ/mol
#Lys-Trp = -3.7 kCal/mol = -15.48 kJ/mol
#Disulfide bond energy -60 kcal/mol = -251 kJ/mol
def getResidue(line):
    # Input line of pdb file (mostly it relates to ATOM lines) and return info of this line, in format of HIS-26A
    residuename = line[17:20];residuename=residuename.lstrip().strip()
    chainname = line[21];chainname=chainname.lstrip().strip()
    if chainname == "":
        chainname = "X"
    residueposition = line[22:26]; residueposition=residueposition.lstrip().strip()
    return(residuename+"-"+residueposition+chainname)

def getXYZ(line):
    # Input pdb line and return coordinates for each atom
    xcoord = line[30:38]; xcoord=xcoord.lstrip().strip()
    ycoord = line[38:46]; ycoord=ycoord.lstrip().strip()
    zcoord = line[46:54]; zcoord=zcoord.lstrip().strip()
    return([float(xcoord),float(ycoord),float(zcoord)])

def getAtom(line):
    # Input pdb line and return name of atom
    atomname = line[12:16]; atomname=atomname.lstrip().strip()
    return(atomname)

def getBfactor(line):
    # Input pdb line and return bfactor 
    return(line[60:66]).lstrip()

def COM(residue):
    # Input: name of residues in format HIS-26A-CA, same format in other functions, and calculate the center of mass for acids
    if residue[0:3] == "LYS":
        if not residue+"-NZ" in coords:
            return("NA")
        else:
            return(numpy.array(coords[residue+"-NZ"]))
    elif residue[0:3] == "ASP":
        coord1 = coords[residue+"-CG"]
        coord2 = coords[residue+"-OD1"]
        coord3 = coords[residue+"-OD2"]
        return(numpy.array([(coord1[0]+coord2[0]+coord3[0])/3,(coord1[1]+coord2[1]+coord3[1])/3,(coord1[2]+coord2[2]+coord3[2])/3]))
    elif residue[0:3] == "GLU":
        coord1 = coords[residue+"-CD"]
        coord2 = coords[residue+"-OE1"]
        coord3 = coords[residue+"-OE2"]  
        return(numpy.array([(coord1[0]+coord2[0]+coord3[0])/3,(coord1[1]+coord2[1]+coord3[1])/3,(coord1[2]+coord2[2]+coord3[2])/3]))
    elif residue[0:3] == "ARG":
        if not residue+"-CZ" in coords:
            return("NA")
        else:
            return(coords[residue+"-CZ"])
    elif residue[0:3] == "HIS":
        coord1 = coords[residue+"-ND1"]
        coord2 = coords[residue+"-CE1"]
        coord3 = coords[residue+"-NE2"]  
        return(numpy.array([(coord1[0]+coord2[0]+coord3[0])/3,(coord1[1]+coord2[1]+coord3[1])/3,(coord1[2]+coord2[2]+coord3[2])/3]))
    elif residue[0:3] == "TYR" or residue[0:3] == "PHE":
        coord1 = coords[residue+"-CG"]
        coord2 = coords[residue+"-CD1"]
        coord3 = coords[residue+"-CD2"]
        coord4 = coords[residue+"-CE1"]
        coord5 = coords[residue+"-CE2"]
        coord6 = coords[residue+"-CZ"]
        return(numpy.array([(coord1[0]+coord2[0]+coord3[0]+coord4[0]+coord5[0]+coord6[0])/6,(coord1[1]+coord2[1]+coord3[1]+coord4[1]+coord5[1]+coord6[1])/6,(coord1[2]+coord2[2]+coord3[2]+coord4[2]+coord5[2]+coord6[2])/6]))
    elif residue[0:3] == "TRP":
        coord1 = coords[residue+"-CG"]
        coord2 = coords[residue+"-CD1"]
        coord3 = coords[residue+"-CD2"]
        coord4 = coords[residue+"-NE1"]
        coord5 = coords[residue+"-CE2"]
        coord6 = coords[residue+"-CE3"]
        coord7 = coords[residue+"-CZ2"]
        coord8 = coords[residue+"-CZ3"]
        coord9 = coords[residue+"-CH2"]
        return(numpy.array([(12*coord1[0]+12*coord2[0]+12*coord3[0]+14*coord4[0]+12*coord5[0]+12*coord6[0]+12*coord7[0]+12*coord8[0]+12*coord9[0])/(12*8+14),(12*coord1[1]+12*coord2[1]+12*coord3[1]+14*coord4[1]+12*coord5[1]+12*coord6[1]+12*coord7[1]+12*coord8[1]+12*coord9[1])/(12*8+14),(12*coord1[2]+12*coord2[2]+12*coord3[2]+14*coord4[2]+12*coord5[2]+12*coord6[2]+12*coord7[2]+12*coord8[2]+12*coord9[2])/(12*8+14)]))
    else:
        return("NONE")

def dist3D(coord1,coord2):
    # Return euclidean distance b/w coords of 2 atoms
    return(numpy.linalg.norm(coord1-coord2))

def isPi(residue):
    # Input: residue; check if residue is PHE, TYR or TRP for pipi interactions
    if residue[0:3]=="PHE" or residue[0:3]=="TYR":
        if residue.split("-")[2] in ["CG","CD1","CD2","CE1","CE2","CZ","OH"]:
            return(True)
        else:
            return(False)
    elif residue[0:3]=="TRP":
        if residue.split("-")[2] in ["NE1","CG","CD1","CD2","CE2","CE3","CZ2","CZ3","CH2"]:
            return(True)
        else:
            return(False)
    else:
        return(False)

def normalVec(residue):
    # Calculate and return coords of normal vector to center of mass for TRP, TYR, PHE acids
    if residue[0:3]=="TRP":
        origin = COM(residue)
        a = numpy.array([coords[residue+"-CG"][0]-origin[0],coords[residue+"-CG"][1]-origin[1],coords[residue+"-CG"][2]-origin[2]])
        b = numpy.array([coords[residue+"-CD1"][0]-origin[0],coords[residue+"-CD1"][1]-origin[1],coords[residue+"-CD1"][2]-origin[2]])
        normalVec = numpy.cross(a,b)
    elif residue[0:3]=="TYR" or residue[0:3]=="PHE":
        origin = COM(residue)
        a = numpy.array([coords[residue+"-CG"][0]-origin[0],coords[residue+"-CG"][1]-origin[1],coords[residue+"-CG"][2]-origin[2]])
        b = numpy.array([coords[residue+"-CD1"][0]-origin[0],coords[residue+"-CD1"][1]-origin[1],coords[residue+"-CD1"][2]-origin[2]])        
        normalVec = numpy.cross(a,b)
    return(normalVec)

def normalize(v):
    # normalization function
    norm=numpy.linalg.norm(v)
    if norm==0:
        return v
    return v/norm

def normalVecPIPI(residue):
    # Calculate coords of origin as mean value of coords in the ring 
    # Calculate coords of cross product of different combinations of ring atoms and take mean value, 
    # this will be vector which perpendicular to ring plane
    x=0;y=0;z=0
    if residue[0:3]=="TRP":
        #If resolution is bad, skip
        if not residue+"-CD2" in coords and not residue+"-CE2" in coords and not residue+"-CE3" in coords:
            return(False)
        #Calculate normal for hexagon
        origin1 = [mean([coords[residue+"-CD2"][0],coords[residue+"-CE2"][0],coords[residue+"-CE3"][0],coords[residue+"-CZ2"][0],coords[residue+"-CZ3"][0],coords[residue+"-CH2"][0]]),mean([coords[residue+"-CD2"][1],coords[residue+"-CE2"][1],coords[residue+"-CE3"][1],coords[residue+"-CZ2"][1],coords[residue+"-CZ3"][1],coords[residue+"-CH2"][1]]),mean([coords[residue+"-CD2"][2],coords[residue+"-CE2"][2],coords[residue+"-CE3"][2],coords[residue+"-CZ2"][2],coords[residue+"-CZ3"][2],coords[residue+"-CH2"][2]])]
        a = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-CE3"][0],coords[residue+"-CD2"][1]-coords[residue+"-CE3"][1],coords[residue+"-CD2"][2]-coords[residue+"-CE3"][2]]))
        b = normalize(numpy.array([coords[residue+"-CE2"][0]-coords[residue+"-CE3"][0],coords[residue+"-CE2"][1]-coords[residue+"-CE3"][1],coords[residue+"-CE2"][2]-coords[residue+"-CE3"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CE2"][0]-coords[residue+"-CZ2"][0],coords[residue+"-CE2"][1]-coords[residue+"-CZ2"][1],coords[residue+"-CE2"][2]-coords[residue+"-CZ2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CE3"][0]-coords[residue+"-CZ2"][0],coords[residue+"-CE3"][1]-coords[residue+"-CZ2"][1],coords[residue+"-CE3"][2]-coords[residue+"-CZ2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CE3"][0]-coords[residue+"-CZ3"][0],coords[residue+"-CE3"][1]-coords[residue+"-CZ3"][1],coords[residue+"-CE3"][2]-coords[residue+"-CZ3"][2]]))
        b = normalize(numpy.array([coords[residue+"-CZ2"][0]-coords[residue+"-CZ3"][0],coords[residue+"-CZ2"][1]-coords[residue+"-CZ3"][1],coords[residue+"-CZ2"][2]-coords[residue+"-CZ3"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CZ2"][0]-coords[residue+"-CH2"][0],coords[residue+"-CZ2"][1]-coords[residue+"-CH2"][1],coords[residue+"-CZ2"][2]-coords[residue+"-CH2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CZ3"][0]-coords[residue+"-CH2"][0],coords[residue+"-CZ3"][1]-coords[residue+"-CH2"][1],coords[residue+"-CZ3"][2]-coords[residue+"-CH2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CZ3"][0]-coords[residue+"-CD2"][0],coords[residue+"-CZ3"][1]-coords[residue+"-CD2"][1],coords[residue+"-CZ3"][2]-coords[residue+"-CD2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CH2"][0]-coords[residue+"-CD2"][0],coords[residue+"-CH2"][1]-coords[residue+"-CD2"][1],coords[residue+"-CH2"][2]-coords[residue+"-CD2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CH2"][0]-coords[residue+"-CE2"][0],coords[residue+"-CH2"][1]-coords[residue+"-CE2"][1],coords[residue+"-CH2"][2]-coords[residue+"-CE2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-CE2"][0],coords[residue+"-CD2"][1]-coords[residue+"-CE2"][1],coords[residue+"-CD2"][2]-coords[residue+"-CE2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        normalVec1 = numpy.array([x/6,y/6,z/6])
        #Calculate normal for pentagon
        x=0;y=0;z=0
        origin2 = [mean([coords[residue+"-CD2"][0],coords[residue+"-CE2"][0],coords[residue+"-NE1"][0],coords[residue+"-CD1"][0],coords[residue+"-CG"][0]]),mean([coords[residue+"-CD2"][1],coords[residue+"-CE2"][1],coords[residue+"-NE1"][1],coords[residue+"-CD1"][1],coords[residue+"-CG"][1]]),mean([coords[residue+"-CD2"][2],coords[residue+"-CE2"][2],coords[residue+"-NE1"][2],coords[residue+"-CD1"][2],coords[residue+"-CG"][2]])]
        a = normalize(numpy.array([coords[residue+"-CG"][0]-coords[residue+"-CD2"][0],coords[residue+"-CG"][1]-coords[residue+"-CD2"][1],coords[residue+"-CG"][2]-coords[residue+"-CD2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CD1"][0]-coords[residue+"-CD2"][0],coords[residue+"-CD1"][1]-coords[residue+"-CD2"][1],coords[residue+"-CD1"][2]-coords[residue+"-CD2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CD1"][0]-coords[residue+"-NE1"][0],coords[residue+"-CD1"][1]-coords[residue+"-NE1"][1],coords[residue+"-CD1"][2]-coords[residue+"-NE1"][2]]))
        b = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-NE1"][0],coords[residue+"-CD2"][1]-coords[residue+"-NE1"][1],coords[residue+"-CD2"][2]-coords[residue+"-NE1"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-CE2"][0],coords[residue+"-CD2"][1]-coords[residue+"-CE2"][1],coords[residue+"-CD2"][2]-coords[residue+"-CE2"][2]]))
        b = normalize(numpy.array([coords[residue+"-NE1"][0]-coords[residue+"-CE2"][0],coords[residue+"-NE1"][1]-coords[residue+"-CE2"][1],coords[residue+"-NE1"][2]-coords[residue+"-CE2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-NE1"][0]-coords[residue+"-CG"][0],coords[residue+"-NE1"][1]-coords[residue+"-CG"][1],coords[residue+"-NE1"][2]-coords[residue+"-CG"][2]]))
        b = normalize(numpy.array([coords[residue+"-CE2"][0]-coords[residue+"-CG"][0],coords[residue+"-CE2"][1]-coords[residue+"-CG"][1],coords[residue+"-CE2"][2]-coords[residue+"-CG"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CE2"][0]-coords[residue+"-CD1"][0],coords[residue+"-CE2"][1]-coords[residue+"-CD1"][1],coords[residue+"-CE2"][2]-coords[residue+"-CD1"][2]]))
        b = normalize(numpy.array([coords[residue+"-CG"][0]-coords[residue+"-CD1"][0],coords[residue+"-CG"][1]-coords[residue+"-CD1"][1],coords[residue+"-CG"][2]-coords[residue+"-CD1"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        normalVec2 = numpy.array([x/6,y/6,z/6])
        normalVec = [normalVec1,normalVec2]
        origin = [origin1,origin2]
    elif residue[0:3]=="TYR" or residue[0:3]=="PHE":
        if not residue+"-CG" in coords and not residue+"-CD1" in coords and not residue+"-CD2" in coords:
            return(False)
        origin = [mean([coords[residue+"-CG"][0],coords[residue+"-CD1"][0],coords[residue+"-CD2"][0],coords[residue+"-CE1"][0],coords[residue+"-CE2"][0],coords[residue+"-CZ"][0]]),mean([coords[residue+"-CG"][1],coords[residue+"-CD1"][1],coords[residue+"-CD2"][1],coords[residue+"-CE1"][1],coords[residue+"-CE2"][1],coords[residue+"-CZ"][1]]),mean([coords[residue+"-CG"][2],coords[residue+"-CD1"][2],coords[residue+"-CD2"][2],coords[residue+"-CE1"][2],coords[residue+"-CE2"][2],coords[residue+"-CZ"][2]])]
        a = normalize(numpy.array([coords[residue+"-CG"][0]-coords[residue+"-CD2"][0],coords[residue+"-CG"][1]-coords[residue+"-CD2"][1],coords[residue+"-CG"][2]-coords[residue+"-CD2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CD1"][0]-coords[residue+"-CD2"][0],coords[residue+"-CD1"][1]-coords[residue+"-CD2"][1],coords[residue+"-CD1"][2]-coords[residue+"-CD2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CD1"][0]-coords[residue+"-CE1"][0],coords[residue+"-CD1"][1]-coords[residue+"-CE1"][1],coords[residue+"-CD1"][2]-coords[residue+"-CE1"][2]]))
        b = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-CE1"][0],coords[residue+"-CD2"][1]-coords[residue+"-CE1"][1],coords[residue+"-CD2"][2]-coords[residue+"-CE1"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-CE2"][0],coords[residue+"-CD2"][1]-coords[residue+"-CE2"][1],coords[residue+"-CD2"][2]-coords[residue+"-CE2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CE1"][0]-coords[residue+"-CE2"][0],coords[residue+"-CE1"][1]-coords[residue+"-CE2"][1],coords[residue+"-CE1"][2]-coords[residue+"-CE2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CE1"][0]-coords[residue+"-CZ"][0],coords[residue+"-CE1"][1]-coords[residue+"-CZ"][1],coords[residue+"-CE1"][2]-coords[residue+"-CZ"][2]]))
        b = normalize(numpy.array([coords[residue+"-CE2"][0]-coords[residue+"-CZ"][0],coords[residue+"-CE2"][1]-coords[residue+"-CZ"][1],coords[residue+"-CE2"][2]-coords[residue+"-CZ"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CE2"][0]-coords[residue+"-CG"][0],coords[residue+"-CE2"][1]-coords[residue+"-CG"][1],coords[residue+"-CE2"][2]-coords[residue+"-CG"][2]]))
        b = normalize(numpy.array([coords[residue+"-CZ"][0]-coords[residue+"-CG"][0],coords[residue+"-CZ"][1]-coords[residue+"-CG"][1],coords[residue+"-CZ"][2]-coords[residue+"-CG"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CZ"][0]-coords[residue+"-CD1"][0],coords[residue+"-CZ"][1]-coords[residue+"-CD1"][1],coords[residue+"-CZ"][2]-coords[residue+"-CD1"][2]]))
        b = normalize(numpy.array([coords[residue+"-CG"][0]-coords[residue+"-CD1"][0],coords[residue+"-CG"][1]-coords[residue+"-CD1"][1],coords[residue+"-CG"][2]-coords[residue+"-CD1"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        normalVec = numpy.array([x/6,y/6,z/6])
    elif residue[0:3]=="HIS":
        if not residue+"-CG" in coords and not residue+"-CND1" in coords and not residue+"-CD2" in coords:
            return(False)
        origin = [mean([coords[residue+"-CG"][0],coords[residue+"-ND1"][0],coords[residue+"-CD2"][0],coords[residue+"-CE1"][0],coords[residue+"-NE2"][0]]),mean([coords[residue+"-CG"][1],coords[residue+"-ND1"][1],coords[residue+"-CD2"][1],coords[residue+"-CE1"][1],coords[residue+"-NE2"][1]]),mean([coords[residue+"-CG"][2],coords[residue+"-ND1"][2],coords[residue+"-CD2"][2],coords[residue+"-CE1"][2],coords[residue+"-NE2"][2]])]
        a = normalize(numpy.array([coords[residue+"-CG"][0]-coords[residue+"-CD2"][0],coords[residue+"-CG"][1]-coords[residue+"-CD2"][1],coords[residue+"-CG"][2]-coords[residue+"-CD2"][2]]))
        b = normalize(numpy.array([coords[residue+"-ND1"][0]-coords[residue+"-CD2"][0],coords[residue+"-ND1"][1]-coords[residue+"-CD2"][1],coords[residue+"-ND1"][2]-coords[residue+"-CD2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-ND1"][0]-coords[residue+"-CE1"][0],coords[residue+"-ND1"][1]-coords[residue+"-CE1"][1],coords[residue+"-ND1"][2]-coords[residue+"-CE1"][2]]))
        b = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-CE1"][0],coords[residue+"-CD2"][1]-coords[residue+"-CE1"][1],coords[residue+"-CD2"][2]-coords[residue+"-CE1"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CD2"][0]-coords[residue+"-NE2"][0],coords[residue+"-CD2"][1]-coords[residue+"-NE2"][1],coords[residue+"-CD2"][2]-coords[residue+"-NE2"][2]]))
        b = normalize(numpy.array([coords[residue+"-CE1"][0]-coords[residue+"-NE2"][0],coords[residue+"-CE1"][1]-coords[residue+"-NE2"][1],coords[residue+"-CE1"][2]-coords[residue+"-NE2"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-CE1"][0]-coords[residue+"-CG"][0],coords[residue+"-CE1"][1]-coords[residue+"-CG"][1],coords[residue+"-CE1"][2]-coords[residue+"-CG"][2]]))
        b = normalize(numpy.array([coords[residue+"-NE2"][0]-coords[residue+"-CG"][0],coords[residue+"-NE2"][1]-coords[residue+"-CG"][1],coords[residue+"-NE2"][2]-coords[residue+"-CG"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        a = normalize(numpy.array([coords[residue+"-NE2"][0]-coords[residue+"-ND1"][0],coords[residue+"-NE2"][1]-coords[residue+"-ND1"][1],coords[residue+"-NE2"][2]-coords[residue+"-ND1"][2]]))
        b = normalize(numpy.array([coords[residue+"-CG"][0]-coords[residue+"-ND1"][0],coords[residue+"-CG"][1]-coords[residue+"-ND1"][1],coords[residue+"-CG"][2]-coords[residue+"-ND1"][2]]))
        normalVec = numpy.cross(a,b)
        if normalVec[0]<0:
            normalVec = normalVec*-1
        x+=normalVec[0];y+=normalVec[1];z+=normalVec[2]
        normalVec = numpy.array([x/5,y/5,z/5])
    return([origin,normalVec])

def mean(vec):
    return(float(sum(vec))/len(vec))

def isH(atom):
    # Check if it is hydrogne atom
    if "H" in atom.split("-")[2] and not "OH" in atom.split("-")[2] and not "NH" in atom.split("-")[2]:
        return(True)
    else:
        return(False)

def isNeighbor(atom1,atom2):
    # Check if 2 acids are neighbors, ex. HIS-26A and GLY-27A - neighbors, HIS-26A and GLY-28A - not neighbors
    if abs(int(atom1.split("-")[1][0:-1])-int(atom2.split("-")[1][0:-1])) == 1:
        return True
    else:
        return False

def isAcceptor(atom):
    # Check whether atoms is acceptor atoms
    SCacceptors = ["ASN-OD1","ASP-OD1","ASP-OD2","GLN-OE1","GLU-OE1","GLU-OE2","HIS-ND1","HIS-NE2","SER-OG","THR-OG1","TYR-OH","CYS-SG","HOH-O"]
    MCacceptors = ["ARG-O","HIS-O","LYS-O","ASP-O","GLU-O","SER-O","THR-O","ASN-O","GLN-O","CYS-O","SEC-O","GLY-O","PRO-O","ALA-O","VAL-O","ILE-O","LEU-O","MET-O","PHE-O","TYR-O","TRP-O"]
    if atom.split("-")[0]+"-"+atom.split("-")[2] in SCacceptors or atom.split("-")[0]+"-"+atom.split("-")[2] in MCacceptors:
        return(True)
    else:
        return(False)

def isDonor(atom):
    # Check whether atoms is donor atoms
    SCdonors = ["ARG-NE","ARG-NH1","ARG-NH2","ASN-ND2","GLN-NE2","HIS-ND1","HIS-NE2","LYS-NZ","SER-OG","THR-OG1","TRP-NE1","TYR-OH","CYS-SG","HOH-O"]
    MCdonors = ["ARG-N","HIS-N","LYS-N","ASP-N","GLU-N","SER-N","THR-N","ASN-N","GLN-N","CYS-N","SEC-N","GLY-N","PRO-N","ALA-N","VAL-N","ILE-N","LEU-N","MET-N","PHE-N","TYR-N","TRP-N"]
    if atom.split("-")[0]+"-"+atom.split("-")[2] in SCdonors or atom.split("-")[0]+"-"+atom.split("-")[2] in MCdonors:
        return(True)
    else:
        return(False)

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def isAminoAcid(acid):
    if acid in ["ARG","HIS","LYS","ASP","GLU","SER","THR","ASN","GLN","CYS","PRO","GLY","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP","HOH"]+metals:
        return(True)
    else:
        return(False)

def isDNA(acid):
    if acid in ["DA","DG","DC","DT"]:
        return(True)
    else:
        return(False)

def SPHyb(atom):
    # Return the type of hybidization
    atom = atom.split("-")
    if (atom[0]=="SER" and atom[2]=="OG") or (atom[0]=="THR" and atom[2]=="OG1") or (atom[0]=="CYS" and atom[2]=="SG") or (atom[0]=="LYS" and atom[2]=="NZ") or (atom[0]=="HOH" and atom[2]=="O"):
        return("SP3")
    else:
        return("SP2")
# Used in search of hydrogen bonds
hydrogen = {
"ARG-H": "N", "ARG-HE": "NE", "ARG-HH11": "NH1", "ARG-HH12": "NH1", "ARG-HH21": "NH2", "ARG-HH22": "NH2",
"HIS-H": "N", "HIS-HE1": "NE2",
"LYS-H": "N", "LYS-HZ1": "NZ", "LYS-HZ2": "NZ", "LYS-HZ3": "NZ",
"ASP-H": "N",
"GLU-H": "N",
"SER-H": "N", "SER-HG": "OG",
"THR-H": "N", "THR-HG1": "OG1",
"ASN-H": "N", "ASN-HD21": "ND2", "ASN-HD22": "ND2",
"GLN-H": "N", "GLN-HE21": "NE2", "GLN-HE22": "NE2",
"CYS-H": "N", "CYS-HG": "SG",
"SEC-H": "N",
"GLY-H": "N",
"ALA-H": "N",
"VAL-H": "N",
"ILE-H": "N",
"LEU-H": "N",
"MET-H": "N",
"PHE-H": "N",
"TYR-H": "N", "TYR-HH": "OH",
"TRP-H": "N", "TRP-HE1": "NE1",
"HOH-H1": "O", "HOH-H2": "O"}

nextAcid = {
"ARG-N": ["CA"],"ARG-NE": ["CD","CZ"], "ARG-NH1": ["CZ"], "ARG-NH2": ["CZ"], "ARG-O": ["C"],
"HIS-N": ["CA"], "HIS-NE2": ["CG","CE1"], "HIS-ND1": ["CD2","CE1"], "HIS-O": ["C"],
"LYS-N": ["CA"], "LYS-NZ": ["CE"], "LYS-O": ["C"],
"ASP-N": ["CA"], "ASP-OD1": ["CG"], "ASP-OD2": ["CG"], "ASP-O": ["C"],
"GLU-N": ["CA"], "GLU-OE1": ["CD"], "GLU-OE2": ["CD"], "GLU-O": ["C"],
"SER-N": ["CA"], "SER-OG": ["CB"], "SER-O": ["C"],
"THR-N": ["CA"], "THR-OG1": ["CB"], "THR-O": ["C"],
"ASN-N": ["CA"], "ASN-OD1": ["CG"], "ASN-ND2": ["CG"], "ASN-O": ["C"],
"GLN-N": ["CA"], "GLN-OE1": ["CD"], "GLN-NE2": ["CD"], "GLN-O": ["C"],
"CYS-N": ["CA"], "CYS-O": ["C"], "CYS-SG": ["CB"],
"SEC-N": ["CA"], "SEC-O": ["C"],
"GLY-N": ["CA"], "GLY-O": ["C"],
"PRO-N": ["CA","CB"],"PRO-O": ["C"],
"ALA-N": ["CA"], "ALA-O": ["C"],
"VAL-N": ["CA"], "VAL-O": ["C"],
"ILE-N": ["CA"], "ILE-O": ["C"],
"LEU-N": ["CA"], "LEU-O": ["C"],
"MET-N": ["CA"], "MET-O": ["C"], "MET-SD": ["CG","CE"],
"PHE-N": ["CA"], "PHE-O": ["C"],
"TYR-N": ["CA"], "TYR-OH": ["CZ"], "TYR-O": ["C"],
"TRP-N": ["CA"], "TRP-NE1": ["CD1","CE2"], "TRP-O": ["C"],
"HOH-O": ["H1","H2"]
}

print('Define all functions, dicts and import modules', (time.time()-t0)*1000, 'ms')

##open pdb file
pdbfile = open(sys.argv[1],'r').readlines()

##open phi angle file
phi = {}
phifile = open(re.sub(".pdb","",sys.argv[1])+".phi",'r').readlines()
secStructure = open(sys.argv[1].replace(".pdb","")+"_secondaryStructure",'w')
counter = 0
currentSecStructure = ""
for line in phifile:
    if line[0:3]=="ASG":
        prevStruct = line[33:42].strip()
        break
if prevStruct == "Coil" or prevStruct == "Bridge":
    prevStruct = "CoilA"
if prevStruct == "Turn":
    prevStruct = "TurnA"

currentAngle = -1

for i in range(len(phifile)):
    line = phifile[i]
    if line[0:3]=="ASG":
        if line[12:15].strip()[-1] in string.ascii_uppercase:
            continue
        phi[line[5:8]+line[11:15].strip()+line[9]] = line[42:50].strip()
        secStruct = line[33:42].strip()
        if line[5:8] == "PRO" and line[33:42].strip() == "aHelix":
            if i+1 != len(phifile) and phifile[i+1][33:42].strip() != "aHelix":
                secStruct = phifile[i+1][33:42].strip()
        if secStruct=="Bridge" or secStruct=="Coil":
            if (float(line[42:51]) > 0 and currentAngle < 0 and float(line[42:51])!=360) or (float(line[42:51]) < 0 and currentAngle < 0) or (float(line[42:51])==360 and currentAngle < 0):
                if prevStruct == "CoilA":
                    secStruct="CoilA"; prevStruct = "CoilA"
                else:
                    secStruct="CoilB"; prevStruct = "CoilB"
            elif (float(line[42:51]) < 0 and currentAngle > 0) or (float(line[42:51]) > 0 and currentAngle > 0 and float(line[42:51])!=360) or (float(line[42:51])==360 and currentAngle > 0):
                if prevStruct == "CoilA":
                    secStruct="CoilB"; prevStruct = "CoilB"
                else:
                    secStruct="CoilA"; prevStruct = "CoilA"
        if secStruct=="Turn":
            if (float(line[42:51]) > 0 and currentAngle < 0 and float(line[42:51])!=360) or (float(line[42:51]) < 0 and currentAngle < 0) or (float(line[42:51])==360 and currentAngle < 0):
                if prevStruct == "TurnA":
                    secStruct="TurnA"; prevStruct = "TurnA"
                else:
                    secStruct="TurnB"; prevStruct = "TurnB"
            elif (float(line[42:51]) < 0 and currentAngle > 0) or (float(line[42:51]) > 0 and currentAngle > 0 and float(line[42:51])!=360) or (float(line[42:51])==360 and currentAngle > 0):
                if prevStruct == "TurnA":
                    secStruct="TurnB"; prevStruct = "TurnB"
                else:
                    secStruct="TurnA"; prevStruct = "TurnA"
        if ("Coil" in secStruct or "Turn" in secStruct) and (line[5:8]=="GLY" and float(line[42:51]) > 0 and float(line[42:51])!=360):
            phiangle = line[42:51].strip()
            if float(phiangle)==360:
                phiangle=0
            secStructure.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+str(phiangle)+"\n")
            currentSecStructure = secStruct                
            counter+=1
        elif secStruct != currentSecStructure:
            counter+=1
            phiangle = line[42:51].strip()
            if float(phiangle)==360:
                phiangle=0
            secStructure.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+str(phiangle)+"\n")
            currentSecStructure = secStruct
        else:
            secStructure.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+line[42:51].strip()+"\n")
        currentAngle = float(line[42:51])
        if currentAngle==360:
            curentAngle=-1

secStructure.close()

#Version 2 of secondary structure
phifile = open(re.sub(".pdb","",sys.argv[1])+".phi",'r').readlines()
secStructure = open(sys.argv[1].replace(".pdb","")+"_secondaryStructure2",'w')
counter = 0
prevStruct = ""

for i in range(len(phifile)):
    line = phifile[i]
    if line[0:3]=="ASG":
        if line[12:15].strip()[-1] in string.ascii_uppercase:
            continue
        secStruct = line[33:42].strip()
        if secStruct == "aHelix" and secStruct == prevStruct:
            secStructure.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+line[42:51].strip()+"\n")
        else:
            counter+=1
            secStructure.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+line[42:51].strip()+"\n")
        prevStruct = secStruct

secStructure.close()

##Calculate RSA
phifile = open(re.sub(".pdb","",sys.argv[1])+".phi",'r').readlines()
rsafile = open(re.sub(".pdb","",sys.argv[1])+".rsa",'w')
areafile = open("AAarea",'r')
area = {}
for line in areafile:
    line = line.strip().split("\t")
    area[line[0]] = float(line[1])

for line in phifile:
    if line[0:3] == "ASG":
        if line[12:15].strip()[-1] in string.ascii_uppercase:
            continue
        if line[5:8] in area.keys():
            rsafile.write(line[5:8]+line[12:15].strip()+line[9]+"\t"+str(float(line[64:69].strip())/area[line[5:8]])+"\n")

print('secStructure and RSA files created', (time.time()-t0)*1000, 'ms')

##Define metals
metals = ["NA","MG","K","CA","MN","FE","CO","CU","ZN"]

##Load atomic weights
atomicWeightFile = open("atomicWeights",'r')
atomicWeights = {}
for line in atomicWeightFile:
    line = line.strip("\n").split()
    atomicWeights[line[0]] = float(line[1])

##Dictionaries
coords = {}
residues = []
residue = ""
chain = {}
subunits = []
coordsvec=[]
allatoms=[]

##Bfactor
bfactor = open(re.sub(".pdb","",sys.argv[1])+"_Bfactor",'w')

##Unique atoms
uniqueAtomFile = open("uniqueAtoms",'r')
uniqueAtoms = {}
for line in uniqueAtomFile:
    line = line.strip("\n").split()
    uniqueAtoms[line[0]] = line[1].split(",")
uniqueAtoms["HOH"] = ["O","H1","H2"]

##Initialize HETATMs
ligands = []
ligandCoords = {}
centroid = {}
tmpResidue = ""
tmpX = []
tmpY = []
tmpZ = []
totalWeight = 0

##read through pdb file
##lines starting with ATOM are assumed to be amino acids 
##lines starting with HETATM are assumed to be ligands, but if they are amino acids or waters they will also be included as nodes (waters treated specially)

for line in pdbfile:
    #if not line[0:4]=="ATOM" and not (line[0:6]== "HETATM" and line[17:20]=="HOH") and not (line[0:6]=="HETATM" and line[17:20].lstrip().strip() in metals):
    #    continue
    if not line[0:4]=="ATOM" and not line[0:6]=="HETATM":
        continue
    #if line[0:4]=="TER":
        #break
    atomname = line[12:16]; atomname=atomname.lstrip().strip()
    if (atomname == "OXT"):
        continue
    residuename = line[17:20];residuename=residuename.lstrip().strip()
    chainname = line[21];chainname=chainname.lstrip().strip()
    if chainname == "":
        chainname = "X"
    residueposition = line[22:26]; residueposition=residueposition.lstrip().strip()
    if line[26] in string.ascii_uppercase:
        continue
    if (not isAminoAcid(residuename) or not residueposition.isdigit()) and (not isDNA(residuename) and (not line[0:6]=="HETATM")):
        continue
    #Make dictionary of centroids
    if not atomname in atomicWeights.keys():
        atomicWeights[atomname] = 1
    if  residuename in uniqueAtoms:
        if tmpResidue == getResidue(line) and (not atomname in ["C","O","N","CA"] or (atomname == "CA" and residuename == "GLY")):
        #if tmpResidue == getResidue(line) and atomname in uniqueAtoms[residuename] and residuename in uniqueAtoms:
            tmpX.append(getXYZ(line)[0]*atomicWeights[atomname])
            tmpY.append(getXYZ(line)[1]*atomicWeights[atomname])
            tmpZ.append(getXYZ(line)[2]*atomicWeights[atomname])
            totalWeight = totalWeight+atomicWeights[atomname]
        else:
            if len(tmpX) == 0 and (not atomname in ["C","O","N","CA"] or (atomname == "CA" and residuename == "GLY")):
            #if len(tmpX) == 0 and atomname in uniqueAtoms[residuename] and residuename in uniqueAtoms:
                tmpResidue = getResidue(line)
                tmpX = [getXYZ(line)[0]*atomicWeights[atomname]]
                tmpY = [getXYZ(line)[1]*atomicWeights[atomname]]
                tmpZ = [getXYZ(line)[2]*atomicWeights[atomname]]
                totalWeight = atomicWeights[atomname]
            else:
                if (not atomname in ["C","O","N","CA"] or (atomname == "CA" and residuename == "GLY")):
                #if atomname in uniqueAtoms[residuename]:
                    centroid[tmpResidue] = [sum(tmpX)/float(totalWeight),sum(tmpY)/float(totalWeight),sum(tmpZ)/float(totalWeight)]
                    tmpResidue = getResidue(line)
                    tmpX = [getXYZ(line)[0]*atomicWeights[atomname]]
                    tmpY = [getXYZ(line)[1]*atomicWeights[atomname]]
                    tmpZ = [getXYZ(line)[2]*atomicWeights[atomname]]
                    totalWeight = atomicWeights[atomname]
        #Add last residue to centroid dictionary
        if (not atomname in ["C","O","N","CA"] or (atomname == "CA" and residuename == "GLY")):
        #if atomname in uniqueAtoms[residuename]:
            centroid[tmpResidue] = [sum(tmpX)/float(totalWeight),sum(tmpY)/float(totalWeight),sum(tmpZ)/float(totalWeight)]
    if isAminoAcid(residuename) or line[17:20]=="HOH" or (line[17:20].lstrip().strip() in metals) or isDNA(residuename):
        coords[getResidue(line)+"-"+getAtom(line)]=numpy.array(getXYZ(line))
        coordsvec.append(getXYZ(line))
        subunits.append(chainname)
        allatoms.append(getResidue(line)+"-"+getAtom(line))
        if getResidue(line)!=residue:
            residues.append(getResidue(line))
            residue=getResidue(line)
        if getAtom(line) in ["N","O","C","CA","HA2","HA3"]:
            if "GLY" in getResidue(line) and getAtom(line) in ["O","CA","HA2","HA3","N","NH"]:
                chain[getResidue(line)+"-"+getAtom(line)] = "SC"
            else:
                chain[getResidue(line)+"-"+getAtom(line)] = "MC"
        else:
            chain[getResidue(line)+"-"+getAtom(line)] = "SC"
    if line[0:6]=="HETATM" and (not (line[17:20]=="HOH") and not (line[17:20].lstrip().strip() in metals)):
        ligands.append(residuename+"-"+residueposition+chainname)    
        ligandCoords[getResidue(line)+"-"+getAtom(line)]=numpy.array(getXYZ(line))
    if isDNA(residuename):
        ligands.append(residuename+"-"+residueposition+chainname)
        ligandCoords[getResidue(line)+"-"+getAtom(line)]=numpy.array(getXYZ(line))
    if getAtom(line) == "CA":
        bfactor.write(residuename+residueposition+chainname+"\t"+getBfactor(line)+"\n")
    #deal with HETATMs (Nonstandard residues include inhibitors, cofactors, ions, and solvent)



ligands = list(set(ligands))

bfactor.close()
subunits = list(set(subunits))

neighboringAtoms = {}
toremove = []
for atom in allatoms:
    #remove HOH that don't have -H1 and -H2 (the reduction failed)
    if "HOH" in atom:
        specAtom = atom.split("-")[2]
        water = atom.split("-")[1]
        if not ("HOH-"+water+"-H1" in allatoms and "HOH-"+water+"-H2" in allatoms and "HOH-"+water+"-O" in allatoms):
            toremove.append(atom)
            continue
    #Initialize beighboringAtoms vec
    neighboringAtoms[atom]=[]

for atom in toremove:
    allatoms.remove(atom)

if len(allatoms) > 52000:
    # print "There are "+str(len(allatoms))+" atoms. This is too many. Exiting."
    sys.exit()

if len(allatoms) < 17000:
    ##Fast approach
    N=len(allatoms)
    distances = spatial.distance.pdist(numpy.array(coordsvec))
    for i in range(len(allatoms)-1):
        atom1 = allatoms[i]
        distanceVec = distances[i*N-sum(range(0,i+1)):i*N-sum(range(0,i+1))+N-(i+1)]
        otherAtoms = allatoms[i+1:]
        neighbors = [otherAtoms[k] for k in range(len(distanceVec)) if distanceVec[k] < 7]
        neighboringAtoms[atom1] = neighboringAtoms[atom1]+neighbors
        [neighboringAtoms[j].append(atom1) for j in neighbors]
else:
    ##Slow approach
    # print str(len(allatoms))+" atoms in total. Using slower low memory approach."
    for i in range(len(allatoms)-1):
        atom1 = allatoms[i]
        #print atom1+": "+str(numpy.array([[coords[atom1][0],coords[atom1][1],coords[atom1][2]]]))+"  "+str(numpy.array(coordsvec[i+1:]))
        distanceVec = spatial.distance.cdist(numpy.array([[coords[atom1][0],coords[atom1][1],coords[atom1][2]]]),numpy.array(coordsvec[i+1:]))
        otherAtoms = allatoms[i+1:]
        neighbors = [otherAtoms[k] for k in range(len(distanceVec[0])) if distanceVec[0][k] < 7]
        neighboringAtoms[atom1] = neighboringAtoms[atom1]+neighbors
        [neighboringAtoms[j].append(atom1) for j in neighbors]

print('Found neighboringAtoms', (time.time()-t0)*1000, 'ms')

  

##Barlow and Thorton salt bridge (saltMethod = 1)
##Kumar and Nussinov salt bridges (saltMethod = 2)
##ASP and GLU (OEs) can connect to ARG, HIS and LYS
##Only 1 connection allowed between the 2 O's and only 1 connection btwn 2 N's on ARG
##Hydrogen bonds should not be counted if salt bridge is presentas
print("Starting salt bridges.", (time.time()-t0)*1000, 'ms')
saltMethod = 1
out = open(sys.argv[1].replace(".pdb","_saltBridges_Barlow"),'w')
saltBridges = []
saltBridgesTmp = []
# allSalt = open(sys.argv[1].replace(".pdb","_allSalt"),'w')
for atom1 in allatoms:
    if not atom1[0:3] in ["ASP","GLU","ARG","LYS","HIS"] or not chain[atom1] == "SC":
        continue
    for atom2 in neighboringAtoms[atom1]:
        if not atom2[0:3] in ["ASP","GLU","ARG","LYS","HIS"] or not chain[atom2] == "SC":
            continue
        #if we've already evaluated these two atoms, skip
        if allatoms.index(atom2) < allatoms.index(atom1):
            continue
        if not "".join(atom1.split("-")[0:2])=="".join(atom2.split("-")[0:2]) and not isNeighbor(atom1,atom2):
            if (atom1[0:3] in ["ASP","GLU"] and atom2[0:3] in ["ARG","LYS","HIS"]) or (atom2[0:3] in ["ASP","GLU"] and atom1[0:3] in ["ARG","LYS","HIS"]):
                if (("NE" in atom1.split("-")[2] or "NH" in atom1.split("-")[2] or "ND" in atom1.split("-")[2] or "NZ" in atom1.split("-")[2])  and ("OE" in atom2.split("-")[2] or "OD" in atom2.split("-")[2])) or (("NE" in atom2.split("-")[2] or "NH" in atom2.split("-")[2] or "ND" in atom2.split("-")[2] or "NZ" in atom2.split("-")[2])  and ("OE" in atom1.split("-")[2] or "OD" in atom1.split("-")[2])):
                    # allSalt.write(atom1+"\t"+atom2+"\t"+str(dist3D(coords[atom1],coords[atom2]))+"\n")
                    if (saltMethod==1 and dist3D(coords[atom1],coords[atom2]) < 4) or (saltMethod==2 and dist3D(COM("-".join(atom1.split("-")[0:2])),COM("-".join(atom2.split("-")[0:2]))) < 4):
                        saltBridges.append([atom1,atom2])
                        saltBridges.append([atom2,atom1])
                        if not ["".join(atom1.split("-")[0:2]),"".join(atom2.split("-")[0:2])] in saltBridgesTmp and not ["".join(atom2.split("-")[0:2]),"".join(atom1.split("-")[0:2])] in saltBridgesTmp:
                            atom1tmp = atom1
                            atom2tmp = atom2
                            if atom1.split("-")[2]=="NH1" or atom1.split("-")[2]=="NH2":
                                atom1tmp = atom1.split("-")[0]+"-"+atom1.split("-")[1]+"-NH1/2"
                            if atom2.split("-")[2]=="NH1" or atom2.split("-")[2]=="NH2":
                                atom2tmp = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-NH1/2"
                            if atom1.split("-")[2]=="OD1" or atom1.split("-")[2]=="OD2":
                                atom1tmp = atom1.split("-")[0]+"-"+atom1.split("-")[1]+"-OD1/2"
                            if atom2.split("-")[2]=="OD1" or atom2.split("-")[2]=="OD2":
                                atom2tmp = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-OD1/2"
                            if atom1.split("-")[2]=="OE1" or atom1.split("-")[2]=="OE2":
                                atom1tmp = atom1.split("-")[0]+"-"+atom1.split("-")[1]+"-OE1/2"
                            if atom2.split("-")[2]=="OE1" or atom2.split("-")[2]=="OE2":
                                atom2tmp = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-OE1/2"
                            out.write(atom1tmp+"\t"+atom2tmp+"\t20\t"+atom1tmp.split("-")[2]+"\t"+atom2tmp.split("-")[2]+"\n")
                            saltBridgesTmp.append(["".join(atom1.split("-")[0:2]),"".join(atom2.split("-")[0:2])])
                            saltBridgesTmp.append(["".join(atom2.split("-")[0:2]),"".join(atom1.split("-")[0:2])])
out.close()
# allSalt.close()

print('Found salt bridges', (time.time()-t0)*1000, 'ms')

# ------------------------------------------------------------------------------------------------------------------------------
##Hydrogen bonds
#For now, donor can only be N,O,S
#Acceptor is N,O,S
#Water hybridization: https://socratic.org/questions/what-is-the-hybridization-of-h2o
print("Starting hydrogen bonds search", (time.time()-t0)*1000, 'ms')

hydrogenBonds=[]
acceptorBonds = {}
for donor in allatoms:
    dH = donor.split("-")[0]+"-"+donor.split("-")[2]
    if not dH in hydrogen:
        continue
    if 'HOH' in donor: # add to no include water
        continue
    if isDNA(donor.split("-")[0]):
        continue
    dA = donor.split("-")[0]+"-"+donor.split("-")[1]+"-"+hydrogen[dH]
    if not isDonor(dA):
        continue
    for acceptor in neighboringAtoms[donor]:
        if [dA,acceptor] in saltBridges or [acceptor,dA] in saltBridges:
            continue
        if 'HOH' in acceptor: # add to no include water
            continue
        if not "-".join(donor.split("-")[0:2])=="-".join(acceptor.split("-")[0:2]) and not isNeighbor(donor,acceptor):
            if (acceptor.split("-")[0]+"-"+acceptor.split("-")[2] in nextAcid) and isAcceptor(acceptor):
                d1 = dist3D(coords[donor],coords[acceptor])
                if d1 > 3.5:
                    continue
                d2 = dist3D(coords[dA],coords[acceptor])
                if d2 > 4.5:
                    continue
                a = [coords[acceptor][0]-coords[donor][0],coords[acceptor][1]-coords[donor][1],coords[acceptor][2]-coords[donor][2]]
                b = [coords[dA][0]-coords[donor][0],coords[dA][1]-coords[donor][1],coords[dA][2]-coords[donor][2]]
                alpha = numpy.degrees(angle(a,b))
                aNext = nextAcid[acceptor.split("-")[0]+"-"+acceptor.split("-")[2]]
                if len(aNext)==2:
                    neighbor1 = "-".join(acceptor.split("-")[0:2])+"-"+aNext[0]
                    neighbor2 = "-".join(acceptor.split("-")[0:2])+"-"+aNext[1]
                    A1 = [(coords[neighbor1][0]+coords[neighbor2][0])/2,(coords[neighbor1][1]+coords[neighbor2][1])/2,(coords[neighbor1][2]+coords[neighbor2][2])/2]
                else:
                    A1 = coords["-".join(acceptor.split("-")[0:2])+"-"+aNext[0]]                                                                    
                a = [A1[0]-coords[acceptor][0],A1[1]-coords[acceptor][1],A1[2]-coords[acceptor][2]]
                b = [coords[donor][0]-coords[acceptor][0],coords[donor][1]-coords[acceptor][1],coords[donor][2]-coords[acceptor][2]]
                beta = numpy.degrees(angle(a,b))
                a = [A1[0]-coords[acceptor][0],A1[1]-coords[acceptor][1],A1[2]-coords[acceptor][2]]
                b = [coords[dA][0]-coords[acceptor][0],coords[dA][1]-coords[acceptor][1],coords[dA][2]-coords[acceptor][2]]
                gamma = numpy.degrees(angle(a,b))
                if (not "S" in dA.split("-")[2] and not "S" in acceptor.split("-")[2] and d1 < 2.5 and d2 < 3.9 and alpha > 90 and beta > 90 and gamma > 90) or ("S" in dA.split("-")[2] and not "S" in acceptor.split("-")[2] and d1 < 3.0 and d2 < 4.3 and alpha > 90 and beta > 90 and gamma > 90) or (not "S" in dA.split("-")[2] and "S" in acceptor.split("-")[2] and d1 < 3.2 and d2 < 4.1 and alpha > 90 and beta > 90 and gamma > 90) or ("S" in dA.split("-")[2] and "S" in acceptor.split("-")[2] and d1 < 3.3 and d2 < 4.5 and alpha > 90 and beta > 90 and gamma > 90):
                    if not acceptor in acceptorBonds.keys():
                        acceptorBonds[acceptor] = []
                    if acceptor.split("-")[2] in ["OD1","OD2","OE1","OE2","OG","OG1","SG"]:
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 2:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:2]
                    if acceptor.split("-")[2] in ["O","ND1","NE2","OH"] and acceptor.split("-")[0]!="HOH":
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 1:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:1]
                    #if acceptor is a water, allow unlimited bonds
                    if acceptor.split("-")[2] == "O" and acceptor.split("-")[0]=="HOH":
                        acceptorBonds[acceptor].append(d1)
                    #Ef = V*[5*(do/d)^12 - 6*(do/d)^10]*F where do is 2.8, d is the distance between donor and acceptor atom (ie d2), V is 8 (Jacobs et al)
                    #sp3 donor - sp3 acceptor F=cos^2(alpha)cos^2(beta-109.5)
                    #sp3 donor - sp2 acceptor F=cos^2(alpha)cos^2(beta)
                    #sp2 donor - sp3 acceptor F=cos^4(alpha)
                    #sp2 donor - sp2 acceptor F=cos^2(alpha)cos^2(max[alpha,beta])
                    if beta < 90:
                        beta = 180-beta
                    if SPHyb(dA)=="SP3" and SPHyb(acceptor)=="SP3":                        
                        E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(beta-109.5)))**2
                    elif SPHyb(dA)=="SP3" and SPHyb(acceptor)=="SP2":
                        E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(beta)))**2
                    elif SPHyb(dA)=="SP2" and SPHyb(acceptor)=="SP3":
                        E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**4)
                    elif SPHyb(dA)=="SP2" and SPHyb(acceptor)=="SP2":
                        normalVecDonor=""
                        normalVecAcceptor=""
                        #plane donor
                        if dA.split("-")[0] == "HIS" and dA.split("-")[2]=="ND1":
                            origin = coords["-".join(dA.split("-")[0:2])+"-ND1"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CG"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CE1"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CE1"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CE1"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "HIS" and dA.split("-")[2]=="NE2":
                            origin = coords["-".join(dA.split("-")[0:2])+"-NE2"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CD2"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CD2"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CD2"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CE1"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CE1"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CE1"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "ARG" and dA.split("-")[2]=="NE":
                            origin = coords["-".join(dA.split("-")[0:2])+"-NE"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CD"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CD"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CD"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CZ"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CZ"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CZ"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "ARG" and dA.split("-")[2]=="NH1":
                            origin = coords["-".join(dA.split("-")[0:2])+"-NH1"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-HH11"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-HH11"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-HH11"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CZ"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CZ"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CZ"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "ARG" and dA.split("-")[2]=="NH2":
                            origin = coords["-".join(dA.split("-")[0:2])+"-NH2"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-HH21"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-HH21"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-HH21"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CZ"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CZ"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CZ"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "TRP" and dA.split("-")[2]=="NE1":
                            origin = coords["-".join(dA.split("-")[0:2])+"-NE1"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CD1"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CD1"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CD1"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CE2"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CE2"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CE2"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "TYR" and dA.split("-")[2]=="OH":
                            origin = coords["-".join(dA.split("-")[0:2])+"-OH"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CZ"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CZ"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CZ"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CE2"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CE2"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CE2"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "ASN" and dA.split("-")[2]=="ND2":
                            origin = coords["-".join(dA.split("-")[0:2])+"-ND2"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-HD21"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-HD21"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-HD21"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CG"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "GLN" and dA.split("-")[2]=="NE2":
                            origin = coords["-".join(dA.split("-")[0:2])+"-NE2"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-HE21"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-HE21"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-HE21"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CD"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CD"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CD"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[2] == "N":
                            origin = coords["-".join(dA.split("-")[0:2])+"-N"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-CA"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-CA"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-CA"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-H"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-H"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-H"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        elif dA.split("-")[0] == "HOH" and dA.split("-")[2]=="O":
                            origin = coords["-".join(dA.split("-")[0:2])+"-O"]
                            a = numpy.array([coords["-".join(dA.split("-")[0:2])+"-H1"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-H1"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-H1"][2]-origin[2]])
                            b = numpy.array([coords["-".join(dA.split("-")[0:2])+"-H2"][0]-origin[0],coords["-".join(dA.split("-")[0:2])+"-H2"][1]-origin[1],coords["-".join(dA.split("-")[0:2])+"-H2"][2]-origin[2]])
                            normalVecDonor = numpy.cross(a,b)
                        #plane acceptor
                        if acceptor.split("-")[0] == "ASN" and acceptor.split("-")[2]=="OD1":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OD1"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CG"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CB"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CB"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CB"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "ASP" and acceptor.split("-")[2]=="OD1":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OD1"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CG"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-OD2"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-OD2"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-OD2"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "ASP" and acceptor.split("-")[2]=="OD2":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OD2"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CG"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-OD1"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-OD1"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-OD1"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "GLN" and acceptor.split("-")[2]=="OE1":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OE1"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CD"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CD"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CD"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CG"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "GLU" and acceptor.split("-")[2]=="OE1":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OE1"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CD"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CD"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CD"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-OE2"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-OE2"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-OE2"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "GLU" and acceptor.split("-")[2]=="OE2":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OE2"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CD"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CD"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CD"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-OE1"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-OE1"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-OE1"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "HIS" and acceptor.split("-")[2]=="ND1":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-ND1"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CG"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CG"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CG"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CE1"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CE1"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CE1"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "HIS" and acceptor.split("-")[2]=="ND2":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-ND2"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CD2"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CD2"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CD2"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CE1"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CE1"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CE1"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[0] == "TYR" and acceptor.split("-")[2]=="OH":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-OH"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CZ"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CZ"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CZ"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CE1"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CE1"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CE1"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[2] == "O" and acceptor.split("-")[0] != "HOH":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-O"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-C"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-C"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-C"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-CA"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-CA"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-CA"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        elif acceptor.split("-")[2] == "O":
                            origin = coords["-".join(acceptor.split("-")[0:2])+"-O"]
                            a = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-H1"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-H1"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-H1"][2]-origin[2]])
                            b = numpy.array([coords["-".join(acceptor.split("-")[0:2])+"-H2"][0]-origin[0],coords["-".join(acceptor.split("-")[0:2])+"-H2"][1]-origin[1],coords["-".join(acceptor.split("-")[0:2])+"-H2"][2]-origin[2]])
                            normalVecAcceptor = numpy.cross(a,b)
                        if normalVecDonor == "" or normalVecAcceptor == "":
                            continue
                        psi = numpy.degrees(angle(normalVecDonor,normalVecAcceptor))
                        E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(max([beta,psi]))))**2
                    hydrogenBonds.append(donor+"\t"+acceptor+"\t"+chain[dA]+chain[acceptor]+"\t"+str(d1)+"\t"+str(d2)+"\t"+str(alpha)+"\t"+str(beta)+"\t"+str(gamma)+"\t"+dA+"\t"+str(E)+"\t"+donor.split("-")[2]+"\t"+acceptor.split("-")[2]+"\n")
                

with open(sys.argv[1].replace(".pdb","_hb"),'w') as out:
    finalHydrogenBonds = []
    donorBonds = {}
    for interaction in hydrogenBonds:
        donor = interaction.split("\t")[0]
        dH = donor.split("-")[0]+"-"+donor.split("-")[2]
        dA = donor.split("-")[0]+"-"+donor.split("-")[1]+"-"+hydrogen[dH]
        acceptor = interaction.split("\t")[1]
        d1 = interaction.split("\t")[3]
        if not dA in donorBonds:
            donorBonds[dA] = 0
        #if water, allow all interactions
        if dA.split("-")[0]=="HOH" or acceptor.split("-")[0]=="HOH":
            out.write(interaction)
        #if not water, there are restrictions on the number of hydrogen bonds it can have
        if len(acceptorBonds[acceptor])==2 and (d1==str(acceptorBonds[acceptor][0]) or d1==str(acceptorBonds[acceptor][1])):
            if dA.split("-")[2] not in ["NH1","NH2","ND2","NE2","NZ"] and donorBonds[dA]==0:
                out.write(interaction)
                donorBonds[dA]+=1
            elif dA.split("-")[2] in ["NH1","NH2","ND2","NE2"] and donorBonds[dA] < 2:
                out.write(interaction)
                donorBonds[dA]+=1
            elif dA.split("-")[2] in ["NZ"] and donorBonds[dA] < 3:
                out.write(interaction)
                donorBonds[dA]+=1
        elif d1 == str(acceptorBonds[acceptor][0]):
            if dA.split("-")[2] not in ["NH1","NH2","ND2","NE2","NZ"] and donorBonds[dA]==0:
                out.write(interaction)
                donorBonds[dA]+=1
            elif dA.split("-")[2] in ["NH1","NH2","ND2","NE2"] and donorBonds[dA] < 2:
                out.write(interaction)
                donorBonds[dA]+=1
            elif dA.split("-")[2] in ["NZ"] and donorBonds[dA] < 3:
                out.write(interaction)
                donorBonds[dA]+=1
            elif dA.split("-")[2] in ["0"]:
                out.write(interaction)

print("hydrogen bonds created", (time.time()-t0)*1000, 'ms')
print(30*'-')

# ------------------------------------------------------------------------------------------------------------------------------
##pipi version 2
##Planes from https://www.schrodinger.com/kb/1556
print("Starting pi-pi bonds search.", (time.time()-t0)*1000, 'ms')
pipi2=[]
pipiDetailed2=[]
prevRes = ""
for atom1 in allatoms:
    if not atom1[0:3] in ["HIS","TYR","TRP","PHE"]:
        continue
    for atom2 in neighboringAtoms[atom1]:
        ring1 = "NORMAL"; ring2 = "NORMAL"
        #skip if we've seen residue 2 before
        if "-".join(atom2.split("-")[0:2]) == prevRes:
            continue
        if not atom2[0:3] in ["HIS","TYR","TRP","PHE"]:
            continue
        #if residue1 equals residue2, continue
        if "-".join(atom1.split("-")[0:2]) == "-".join(atom2.split("-")[0:2]):
            continue
        x=normalVecPIPI("-".join(atom1.split("-")[0:2]))
        if not x:
            continue
        #sort out TRP
        if atom1[0:3]=="TRP" and atom1.split("-")[2] in ["CG","NE1","CD1"]:
            origin1 = numpy.array(x[0][1]); norm1 = x[1][1]
            ring1 = "PENT"
        elif atom1[0:3]=="TRP" and atom1.split("-")[2] not in ["CG","NE1","CD1"]:
            origin1 = numpy.array(x[0][0]); norm1 = x[1][0]
            ring1 = "HEX"
        else:
            origin1 = numpy.array(x[0]); norm1 = x[1]
        y=normalVecPIPI("-".join(atom2.split("-")[0:2]))
        if not y:
            continue
        if atom2[0:3]=="TRP" and atom2.split("-")[2] in ["CG","NE1","CD1"]:
            origin2 = numpy.array(y[0][1]); norm2 = y[1][1]
            ring2 = "PENT"
        elif atom2[0:3]=="TRP" and atom2.split("-")[2] not in ["CG","NE1","CD1"]:
            origin2 = numpy.array(y[0][0]); norm2 = y[1][0]
            ring2 = "HEX"
        else:
            origin2 = numpy.array(y[0]); norm2 = y[1]
        #if we've seen these residues before
        # was mistake with order of ring1 and ring2
        if atom1.split("-")[0]+atom1.split("-")[1]+"\t"+atom2.split("-")[0]+atom2.split("-")[1]+"\t10\t"+ring1+"\t"+ring2+"\n" in pipi2 or atom2.split("-")[0]+atom2.split("-")[1]+"\t"+atom1.split("-")[0]+atom1.split("-")[1]+"\t10\t"+ring2+"\t"+ring1+"\n" in pipi2:
            continue
        pipiangle = numpy.degrees(angle(norm1,norm2))
        if pipiangle > 90:
            pipiangle = 180-pipiangle
        if (dist3D(origin1,origin2) <= 4.4 and pipiangle<30) or (dist3D(origin1,origin2) <=5.5 and 60 < pipiangle < 120):
            pipi2.append(atom1.split("-")[0]+atom1.split("-")[1]+"\t"+atom2.split("-")[0]+atom2.split("-")[1]+"\t10\t"+ring1+"\t"+ring2+"\n")
            pipiDetailed2.append(atom1.split("-")[0]+atom1.split("-")[1]+"\t"+atom2.split("-")[0]+atom2.split("-")[1]+"\t10\t"+atom1.split("-")[2]+"\t"+atom2.split("-")[2]+"\n")

pipi2=list(set(pipi2))
with open(sys.argv[1].replace(".pdb","_pipi2"),'w') as out:
    for interaction in pipi2:
        out.write(interaction)

print('Found pipi interactions', (time.time()-t0)*1000, 'ms')
print(30*'-')


# ------------------------------------------------------------------------------------------------------------------------------
##Pication version 2
#http://www.pnas.org.ezp-prod1.hul.harvard.edu/content/96/17/9459.full
#https://www.schrodinger.com/kb/1556
print("Running pi-cation bonds search", (time.time()-t0)*1000, 'ms')

pication=[]
for atom1 in allatoms:
    if not atom1[0:3] in ["HIS","PHE","TYR","TRP","ARG","LYS"]:
        continue
    for atom2 in neighboringAtoms[atom1]:
        if not atom2[0:3] in ["HIS","PHE","TYR","TRP","ARG","LYS"]:
            continue
        if not "".join(atom1.split("-")[0:2])=="".join(atom2.split("-")[0:2]) and not isNeighbor(atom1,atom2):
            if atom1[0:3] in ["HIS","PHE","TYR","TRP"] and (atom2[0:3] in ["ARG","LYS"] and atom2.split("-")[2] in ["CZ","NZ"]):
                x = normalVecPIPI("-".join(atom1.split("-")[0:2]))
                if x == False:
                    continue
                #Sort out TRP
                if atom1[0:3]=="TRP" and atom1.split("-")[2] in ["CG","NE1","CD1"]:
                    origin = numpy.array(x[0][1]); norm = x[1][1]
                    ring = "PENT"
                elif atom1[0:3]=="TRP" and atom1.split("-")[2] not in ["CG","NE1","CD1"]:
                    origin = numpy.array(x[0][0]); norm = x[1][0]
                    ring = "HEX"
                else:
                    origin = numpy.array(x[0]); norm = x[1]
                    ring = "NORMAL"
                catvec = numpy.array([coords[atom2][0]-origin[0],coords[atom2][1]-origin[1],coords[atom2][2]-origin[2]])
                catangle = numpy.degrees(angle(norm,catvec))
                if catangle > 90:
                    catangle = 180-catangle
                if dist3D(origin,coords[atom2]) < 6.6 and catangle < 30:
                    pication.append(atom1.split("-")[0]+atom1.split("-")[1]+"\t"+atom2.split("-")[0]+atom2.split("-")[1]+"\t9.4\t"+ring+"\t"+atom2.split("-")[2]+"\n")
            if atom2[0:3] in ["HIS","PHE","TYR","TRP"] and (atom1[0:3] in ["ARG","LYS"] and atom1.split("-")[2] in ["CZ","NZ"]):
                x = normalVecPIPI("-".join(atom2.split("-")[0:2]))
                if x==False:
                    continue
                #Sort out TRP
                if atom2[0:3]=="TRP" and atom2.split("-")[2] in ["CG","NE1","CD1"]:
                    origin = numpy.array(x[0][1]); norm = x[1][1]
                    ring = "PENT"
                elif atom2[0:3]=="TRP" and atom2.split("-")[2] not in ["CG","NE1","CD1"]:
                    origin = numpy.array(x[0][0]); norm = x[1][0]
                    ring = "HEX"
                else:
                    origin = numpy.array(x[0]); norm = x[1]
                    ring = "NORMAL"
                catvec = numpy.array([coords[atom1][0]-origin[0],coords[atom1][1]-origin[1],coords[atom1][2]-origin[2]])
                catangle = numpy.degrees(angle(norm,catvec))
                if catangle > 90:
                    catangle = 180-catangle
                if dist3D(origin,coords[atom1]) < 6 and catangle < 30:
                    pication.append(atom2.split("-")[0]+atom2.split("-")[1]+"\t"+atom1.split("-")[0]+atom1.split("-")[1]+"\t9.4\t"+ring+"\t"+atom1.split("-")[2]+"\n")
                           

pication = list(set(pication))
with open(sys.argv[1].replace(".pdb","_pication2"),'w') as out:
    for interaction in pication:
        out.write(interaction)

print('Found pication interactions', (time.time()-t0)*1000, 'ms')
print(30*'-')

# ------------------------------------------------------------------------------------------------------------------------------
##disulfide
print("Starting disulfide bonds among CYS acids", (time.time()-t0)*1000, 'ms')

disulf=[]
for atom1 in allatoms:
    if not atom1[0:3] == "CYS":
        continue
    for atom2 in neighboringAtoms[atom1]:
        if not atom2[0:3] == "CYS":
            continue
        #if we've already evaluated these two atoms, skip
        if allatoms.index(atom2) < allatoms.index(atom1):
            continue
        if atom1.split("-")[2]=="SG" and atom2.split("-")[2]=="SG" and not isNeighbor(atom1,atom2):
            if dist3D(coords[atom1],coords[atom2]) <= 3:
                disulf.append(atom1.split("-")[0]+atom1.split("-")[1]+"\t"+atom2.split("-")[0]+atom2.split("-")[1]+"\t167\t"+atom1.split("-")[2]+"\t"+atom2.split("-")[2]+"\n")

disulf=list(set(disulf))
with open(sys.argv[1].replace(".pdb","_disulf"),'w') as out:
    for interaction in disulf:
        out.write(interaction)

print('Found disulfide bonds', (time.time()-t0)*1000, 'ms')
print(30*'-')

# ------------------------------------------------------------------------------------------------------------------------------
##vanderwaals
##Energy.4-4
print("Define vdw radii and Starting van der waals.", (time.time()-t0)*1000, 'ms')

radii={}
with open("allatoms_radii",'r') as surfacefile:
    for line in surfacefile:
        line = line.strip("\n").split("\t")
        radii[line[0]] = float(line[1])


vdw={}
scores=[]
vdw2={}
scores2=[]
for atom1 in allatoms:
    if "HOH" in atom1 or atom1.split("-")[0] in metals or isDNA(atom1.split("-")[0]) or not isAminoAcid(atom1.split("-")[0]):
        continue
    for atom2 in neighboringAtoms[atom1]:
        if "HOH" in atom2 or atom2.split("-")[0] in metals or not isAminoAcid(atom2.split("-")[0]):
            continue
        #if we've already evaluated these two atoms, skip
        if "".join(atom2.split("-")[0:2])+"-"+"".join(atom1.split("-")[0:2])+"-"+chain[atom2]+"-"+chain[atom1]+"-"+atom2.split("-")[2]+"-"+atom1.split("-")[2] in vdw:
            continue
        #if not "".join(atom1.split("-")[0:2])=="".join(atom2.split("-")[0:2]) and atom1 in radii and atom2 in radii:
        if not "".join(atom1.split("-")[0:2])=="".join(atom2.split("-")[0:2]) and atom1.split("-")[2] in radii and atom2.split("-")[2] in radii:
            rm=radii[atom1.split("-")[2]]+radii[atom2.split("-")[2]]
            #rm=radii[atom1]+radii[atom2]
            r=dist3D(coords[atom1],coords[atom2])
            if r - rm < .5 and not (isNeighbor(atom1,atom2) and chain[atom1]=="MC" and chain[atom2]=="MC"):
                if not "".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2] in vdw:
                    vdw["".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2]] = []
                E = (-.997*((rm/r)**12-2*(rm/r)**6))*4
                vdw["".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2]].append(E)
                scores.append(sum(vdw["".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2]]))
                #scores.append(E)
                #vdw2
                if ("C" in atom1.split("-")[2] and "C" in atom2.split("-")[2]) or ("C" in atom1.split("-")[2] and atom2.split("-")[2] in ["NE2","OE1","ND2","OD1"] and atom2.split("-")[0] in ["GLN","ASN"]) or (atom1.split("-")[2] in ["NE2","OE1","ND2","OD1"] and atom1.split("-")[0] in ["GLN","ASN"] and "C" in atom2.split("-")[2]):
                    if not "".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2] in vdw2:
                        vdw2["".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2]] = []
                    vdw2["".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2]].append(E)
                    scores2.append(sum(vdw["".join(atom1.split("-")[0:2])+"-"+"".join(atom2.split("-")[0:2])+"-"+chain[atom1]+"-"+chain[atom2]+"-"+atom1.split("-")[2]+"-"+atom2.split("-")[2]]))

with open(sys.argv[1].replace(".pdb","_vdw"),'w') as out:
    if len(scores) > 0:
        maxScore=max(scores)
        for contact in vdw:
            if not (sum(vdw[contact]) < 0 and abs(int(re.sub("[^0-9]","",contact.split("-")[0]))-int(re.sub("[^0-9]","",contact.split("-")[1])))==1):
                out.write(contact+"\t"+str(sum(vdw[contact]))+"\t"+contact.split("-")[4]+"\t"+contact.split("-")[5]+"\n")

with open(sys.argv[1].replace(".pdb","_vdw2"),'w') as out:
    if len(scores2) > 0:
        maxScore=max(scores2)
        for contact in vdw2:
            if not (sum(vdw2[contact]) < 0 and abs(int(re.sub("[^0-9]","",contact.split("-")[0]))-int(re.sub("[^0-9]","",contact.split("-")[1])))==1):
                out.write(contact+"\t"+str(sum(vdw2[contact]))+"\t"+contact.split("-")[4]+"\t"+contact.split("-")[5]+"\n")

print('Found vanderwaals bonds', (time.time()-t0)*1000, 'ms')

# ------------------------------------------------------------------------------------------------------------------------------
##Metals
print("Starting metals", (time.time()-t0)*1000, 'ms')
metalBonds=[]
metalBondsFiltered = []
#MCOdist = {"NA": 2.58, "MG": 2.5, "K": 2.84, "CA": 2.56, "MN": 2.29, "FE": 2.3, "CO": 2.3, "CU": 2.3, "ZN": 2.3}
#GLUASPmonoDist = {"NA": 2.61, "MG": 2.17, "K": 3.21, "CA": 2.59, "MN": 2.29, "FE": 2.29, "CO": 2.29, "CU": 2.53, "ZN": 2.49}
#GLUASPbiDist = {"NA": 2.8, "MG": 2.5, "K": 3.2, "CA": 2.8, "MN": 2.5, "FE": 2.5, "CO": 2.5, "CU": 2.5, "ZN": 2.5}
#ASNdist = {"NA": 2.61, "MG": 2.17, "K": 3.21, "CA": 2.59, "MN": 2.29, "FE": 2.29, "CO": 2.29, "CU": 2.53, "ZN": 2.49}
#GLNdist = {"NA": 2.61, "MG": 2.17, "K": 3.21, "CA": 2.59, "MN": 2.29, "FE": 2.29, "CO": 2.29, "CU": 2.53, "ZN": 2.49}
#SERdist = {"NA": 2.61, "MG": 2.17, "K": 3.21, "CA": 2.59, "MN": 2.29, "FE": 2.29, "CO": 2.29, "CU": 2.53, "ZN": 2.49}
#THRdist = {"NA": 2.61, "MG": 2.17, "K": 3.21, "CA": 2.59, "MN": 2.29, "FE": 2.29, "CO": 2.29, "CU": 2.53, "ZN": 2.49}
#TYRdist = {"NA": 2.41, "MG": 2.02, "K": 2.91, "CA": 2.39, "MN": 2.14, "FE": 2.09, "CO": 2.09, "CU": 2.4, "ZN": 2.04}
#HISdist = {"MN": 2.41, "FE": 2.56, "CO": 2.34, "CU": 2.22, "ZN": 2.13}
#CYSdist = {"MN": 2.5, "FE": 2.4, "CO": 2.55, "CU": 2.5, "ZN": 2.51}
#METdist = {"FE": 2.5, "CU": 2.8}
#LYSdist = {"ZN": 2.3}

MCOdist = {"NA": 2.88, "MG": 2.76, "K": 3.24, "CA": 2.86, "MN": 2.69, "FE": 2.54, "CO": 2.58, "CU": 2.54, "ZN": 2.57}
GLUASPmonoDist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 3.29, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
GLUASPbiDist = {"NA": 3, "MG": 2.7, "K": 3.4, "CA": 3.4, "MN": 2.7, "FE": 2.7, "CO": 2.7, "CU": 2.7, "ZN": 2.7}
ASNdist = {"NA": 2.88, "MG": 2.76, "K": 3.24, "CA": 2.86, "MN": 2.69, "FE": 2.54, "CO": 2.58, "CU": 2.54, "ZN": 2.57}
GLNdist = {"NA": 2.88, "MG": 2.76, "K": 3.24, "CA": 2.86, "MN": 2.69, "FE": 2.54, "CO": 2.58, "CU": 2.54, "ZN": 2.57}
SERdist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 2.89, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
THRdist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 2.89, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
TYRdist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 2.89, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
HISdist = {"MN": 2.71, "FE": 2.66, "CO": 2.64, "CU": 2.52, "ZN": 2.73}
CYSdist = {"MN": 2.85, "FE": 2.8, "CO": 2.75, "CU": 2.65, "ZN": 2.81}
METdist = {"FE": 3, "CU": 3.3}
LYSdist = {"ZN": 2.8}


for atom1 in allatoms:
    if atom1.split("-")[0] in metals:
        metal2atom=[]
        for atom2 in neighboringAtoms[atom1]:
            if "HOH" in atom2 or isDNA(atom2.split("-")[0]):
                continue
            if not atom2.split("-")[0] in metals:
                if atom2.split("-")[2] == "O":
                    #Main chain
                    if dist3D(coords[atom1],coords[atom2]) <= MCOdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "GLU") and atom2.split("-")[2] in ["OE1","OE2"]:
                    #GLU
                    #mono
                    if dist3D(coords[atom1],coords[atom2]) <= GLUASPmonoDist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                    #bi
                    tmpatom21 = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-"+"OE1"
                    tmpatom22 = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-"+"OE2"
                    if (dist3D(coords[atom1],coords[tmpatom21]) <= GLUASPbiDist[atom1.split("-")[0]]) and (dist3D(coords[atom1],coords[tmpatom22]) <= GLUASPbiDist[atom1.split("-")[0]]):
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "ASP") and atom2.split("-")[2] in ["OD1","OD2"]:
                    #ASP
                    #mono
                    if dist3D(coords[atom1],coords[atom2]) <= GLUASPmonoDist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                    #bi
                    tmpatom21 = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-"+"OD1"
                    tmpatom22 = atom2.split("-")[0]+"-"+atom2.split("-")[1]+"-"+"OD2"
                    if (dist3D(coords[atom1],coords[tmpatom21]) <= GLUASPbiDist[atom1.split("-")[0]]) and (dist3D(coords[atom1],coords[tmpatom22]) <= GLUASPbiDist[atom1.split("-")[0]]):
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "ASN") and atom2.split("-")[2] == "OD1":
                    #ASN
                    if dist3D(coords[atom1],coords[atom2]) <= ASNdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "GLN") and atom2.split("-")[2] == "OE1":
                    #GLN
                    if dist3D(coords[atom1],coords[atom2]) <= GLNdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "SER") and atom2.split("-")[2] == "OG":
                    #SER
                    if dist3D(coords[atom1],coords[atom2]) <= SERdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "THR") and atom2.split("-")[2] == "OG1":
                    #THR
                    if dist3D(coords[atom1],coords[atom2]) <= THRdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "TYR") and atom2.split("-")[2] == "OH":
                    #TYR
                    if dist3D(coords[atom1],coords[atom2]) <= TYRdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "HIS") and atom2.split("-")[2] in ["NE2","ND1"] and atom1.split("-")[0] in HISdist:
                    #HIS
                    if dist3D(coords[atom1],coords[atom2]) <= HISdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "CYS") and atom2.split("-")[2] == "SG" and atom1.split("-")[0] in CYSdist:
                    #CYS
                    if dist3D(coords[atom1],coords[atom2]) <= CYSdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "MET") and atom2.split("-")[2] == "SD" and atom1.split("-")[0] in METdist:
                    #MET
                    if dist3D(coords[atom1],coords[atom2]) <= METdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
                if (atom2.split("-")[0] == "LYS") and atom2.split("-")[2] == "NZ" and atom1.split("-")[0] in LYSdist:
                    #LYS
                    if dist3D(coords[atom1],coords[atom2]) <= LYSdist[atom1.split("-")[0]]:
                        metalBonds.append(atom1+"\t"+atom2+"\n")
                        metal2atom.append(atom2)
        origin = coords[atom1]
        for atom21 in metal2atom:
            atom21coords = coords[atom21]
            atom21coords = [atom21coords[0]-origin[0],atom21coords[1]-origin[1],atom21coords[2]-origin[2]]
            for atom22 in metal2atom:
                atom22coords = coords[atom22]
                atom22coords = [atom22coords[0]-origin[0],atom22coords[1]-origin[1],atom22coords[2]-origin[2]]
                if atom21 != atom22:
                    ang = numpy.degrees(angle(atom21coords,atom22coords))
                    if atom21+"\t"+atom22+"\t"+chain[atom21]+chain[atom22]+"\t"+atom1+"\t"+str(ang)+"\n" in metalBondsFiltered:
                        continue
                    if atom22+"\t"+atom21+"\t"+chain[atom22]+chain[atom21]+"\t"+atom1+"\t"+str(ang)+"\n" in metalBondsFiltered:
                        continue
                    if atom1.split("-")[0] in ["NA","K","FE","CU","ZN"]:
                        if 90 < ang < 180:
                            metalBondsFiltered.append(atom21+"\t"+atom22+"\t"+chain[atom21]+chain[atom22]+"\t"+atom1+"\t"+str(ang)+"\n")
                    else:
                        if ang > 98:
                            metalBondsFiltered.append(atom21+"\t"+atom22+"\t"+chain[atom21]+chain[atom22]+"\t"+atom1+"\t"+str(ang)+"\n")


#Write out metal bonds
with open(sys.argv[1].replace(".pdb","_metal"),'w') as out:
    for bond in metalBondsFiltered:
        out.write(bond)

print("metals ended", (time.time()-t0)*1000, 'ms')
print(30*'-')


##DNA
#Omega is set to 90
print("Starting DNA bonds...", (time.time()-t0)*1000, 'ms')
DNAbindingPairs = []

for atom1 in allatoms:
    if isDNA(atom1.split("-")[0]):
        nucleotide = atom1
        for atom2 in neighboringAtoms[nucleotide]:
            if not isDNA(atom2.split("-")[0]) and not "HOH" in atom2:
                #Figure out if atom2 and nucleotide interact
                bound = False
                if ("N" in nucleotide.split("-")[2] and (isDonor(atom2) or isAcceptor(atom2))) or ("O"in nucleotide.split("-")[2] and (isDonor(atom2) or isAcceptor(atom2))) or (("N" in nucleotide.split("-")[2] or "O" in nucleotide.split("-")[2]) and "S" in atom2.split("-")[2]) or ("C" in nucleotide.split("-")[2] and "O" in atom2.split("-")[2] and isAcceptor(atom2)) or ("O" in nucleotide.split("-")[2] and atom2.split("-")[2] in ["NZ","NH1","NH2","OD1","OD2","OE1","OE2"]) or ("C" in nucleotide.split("-")[2] and "C" in atom2.split("-")[2]):
                    for atom3 in neighboringAtoms[nucleotide]:
                        if atom3==atom2 or isH(atom3):
                            continue
                        vector1 = [coords[nucleotide][0]-coords[atom3][0],coords[nucleotide][1]-coords[atom3][1],coords[nucleotide][2]-coords[atom3][2]]
                        vector2 = [coords[atom2][0]-coords[atom3][0],coords[atom2][1]-coords[atom3][1],coords[atom2][2]-coords[atom3][2]]
                        omega = numpy.degrees(angle(vector1,vector2))
                        if omega <= 90:
                            bound = True
                        else:
                            bound = False
                            break
                if bound:
                    DNAbindingPairs.append([nucleotide,atom2])
                        

with open(sys.argv[1].replace(".pdb","_DNA"),'w') as out:
    for pair in DNAbindingPairs:
        out.write(pair[0]+"\t"+pair[1]+"\t"+chain[pair[1]]+"\n")
print('Found DNA bonds', (time.time()-t0)*1000, 'ms')


##Ligands
print("Identify centroid of each ligand and find ligand bonds", time.time()-t0)
ligandCentroids = {}
for ligand in ligands:
    x = []
    y = []
    z = []
    for atom in ligandCoords.keys():
        if ligand in atom:
            x.append(ligandCoords[atom][0])
            y.append(ligandCoords[atom][1])
            z.append(ligandCoords[atom][2])
    meanx = numpy.mean(x)
    meany = numpy.mean(y)
    meanz = numpy.mean(z)
    ligandCentroids[ligand] = [meanx,meany,meanz]

# add with open, and sorted(ligands)
with open(sys.argv[1].replace(".pdb","_ligand"),'w') as out:
    for atom in allatoms:
        if "HOH" in atom:
            continue
        if chain[atom]=="SC":
            tmpdist0 = 150
            for ligand in sorted(ligands):
                tmpdist1 = dist3D(coords[atom],ligandCentroids[ligand])
                if tmpdist1 > tmpdist0:
                    continue
                out.write(atom.split("-")[0]+atom.split("-")[1]+"-"+atom.split("-")[2]+"\t"+ligand+"\t"+str(tmpdist1)+"\n")
                tmpdist0 = float(copy.copy(tmpdist1))
print("Found ligand bonds and write to file", (time.time()-t0)*1000, 'ms')


# add with open, sorted(centroid.keys)
with open(re.sub(".pdb","",sys.argv[1])+"_centroidNetSC",'w') as out:
    for residue1 in sorted(centroid.keys()):
        for residue2 in sorted(centroid.keys()):
            if "HOH" in residue1 and "HOH" in residue2:
                continue
            if residue1==residue2:
                continue
            x=dist3D(numpy.array(centroid[residue1]),numpy.array(centroid[residue2])) 
            if x < 8.5:
                out.write(re.sub("-","",residue1)+"\t"+re.sub("-","",residue2)+"\t"+"SCSC"+"\t"+"1"+"\t"+"CENTROID\t"+"CENTROID\t"+"CENTROID\t"+str(x)+"\n")

# add with open file, sorted(centroid.keys) and sorted(ligands)
with open(re.sub(".pdb","",sys.argv[1])+"_centroidNetLigand",'w') as out:
    for residue1 in sorted(centroid.keys()):
        if "HOH" in residue1:
            continue
        tmpdist0 = 150
        for ligand in sorted(ligands):
            tmpdist1 = dist3D(numpy.array(centroid[residue1]),ligandCentroids[ligand])
            if tmpdist1 > tmpdist0:
                continue
            out.write(re.sub("-","",residue1)+"\t"+ligand+"\t"+str(tmpdist1)+"\n")
            tmpdist0 = float(copy.copy(tmpdist1))

print("Found centroidNets and write to file", time.time()-t0, 's')
print(30*'-')

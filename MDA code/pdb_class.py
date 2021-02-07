import MDAnalysis as mda
import numpy as np
import string
import sys, os, re, copy, math
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import subprocess
from MDAnalysis.analysis.hbonds import HydrogenBondAnalysis as hbonds
from MDAnalysis.lib.distances import calc_bonds as MDdist3D
import time
import parameters as prs 

'''
    Before need to preprocess pdb file
    1. Add hydrogens (phenix.reduce)
    2. Protonate water (phenix.ready_set)
    3. And create phi file (stride)
'''

def SecondaryStructure(file1):
    '''
    @description
     	Read pdb and phi file. Extract phi angles and secondary Structures of each acid from phi file
    @input
    	file1 - pdb file, ex. 4eiy.pdb
    @output 
        secondary structure files, ex. 4eiy_secondaryStructure and 4iey_secondaryStructure2
    '''
    pdbfile = open(file1,'r').readlines() 
    phifile = open(re.sub(".pdb","",file1)+".phi",'r').readlines()
    secStructure = open(file1.replace(".pdb","")+"_secondaryStructure",'w')
    
    phi = {} 
    counter = 0
    currentSecStructure = ""
    prevStruct = ''
    # define prevStruct and change it 
    for line in phifile:
        if line[0:3]=="ASG":
            prevStruct = line[33:42].strip()
            break
    if prevStruct == "Coil" or prevStruct == "Bridge":
        prevStruct = "CoilA"
    if prevStruct == "Turn":
        prevStruct = "TurnA"

    currentAngle = -1
    ''' 
    1. go through phi file and define prevStruct as first ASG line
    2. if it is Coil or Bridge change to CoilA, Turn -> TurnA
    3. iterate in phi file, define residue + id + segid, phi angle and name of secstructure
    4. if residue is PRO and structure is not aHelix and (next line not last and next structure is not aHelix), secStruct=new struct
    5. Rename structure
    '''
    for i in range(len(phifile)):
        line = phifile[i]   
        if line[0:3]=="ASG":
            res = line[5:8]
            resid = line[11:15].strip()
            segid = line[9]
            phi_angle = float(line[42:50].strip())
            structure = line[33:42].strip()
            phi[res+resid+segid] = phi_angle
            secStruct = structure

            if res == 'PRO' and structure != 'aHelix':
                if i+1 != len(phifile) and phifile[i+1][33:42].strip() != "aHelix":
                    secStruct = phifile[i+1][33:42].strip()
            if secStruct=="Bridge" or secStruct=="Coil":
                if ((phi_angle > prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or
                    (phi_angle < prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE) or 
                    (phi_angle== prs.FULL_ANGLE and currentAngle < prs.ZERO_ANGLE)):

                    if prevStruct == "CoilA":
                        secStruct="CoilA"; prevStruct = "CoilA"
                    else:
                        secStruct="CoilB"; prevStruct = "CoilB"

                elif ((phi_angle > prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or
                      (phi_angle < prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE) or  
                      (phi_angle== prs.FULL_ANGLE and currentAngle > prs.ZERO_ANGLE)):
                    
                    if prevStruct == "CoilA":
                        secStruct="CoilB"; prevStruct = "CoilB"
                    else:
                        secStruct="CoilA"; prevStruct = "CoilA"

            if secStruct=="Turn":
                if ((phi_angle > prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or 
                    (phi_angle < prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE) or 
                    (phi_angle==prs.FULL_ANGLE and currentAngle < prs.ZERO_ANGLE)):

                    if prevStruct == "TurnA":
                        secStruct="TurnA"; prevStruct = "TurnA"
                    else:
                        secStruct="TurnB"; prevStruct = "TurnB"

                elif ((phi_angle > prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or 
                      (phi_angle < prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE) or
                      (phi_angle==prs.FULL_ANGLE and currentAngle > prs.ERO_ANGLE)):
                    
                    if prevStruct == "TurnA":
                        secStruct="TurnB"; prevStruct = "TurnB"
                    else:
                        secStruct="TurnA"; prevStruct = "TurnA"

            if ("Coil" in secStruct or "Turn" in secStruct) and (res=="GLY" and phi_angle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE):
                phiangle = phi_angle
                if phiangle==prs.FULL_ANGLE:
                    phiangle=0
                secStructure.write(res+resid+segid+"\t"+secStruct+str(counter)+"\t"+str(phiangle)+"\n")
                currentSecStructure = secStruct                
                counter+=1
            elif secStruct != currentSecStructure:
                counter+=1
                phiangle = line[42:51].strip()
                if float(phiangle)==prs.FULL_ANGLE:
                    phiangle=0
                secStructure.write(res+resid+segid+"\t"+secStruct+str(counter)+"\t"+str(phiangle)+"\n")
                currentSecStructure = secStruct
            else:
                secStructure.write(res+resid+segid+"\t"+secStruct+str(counter)+"\t"+str(phi_angle)+"\n")
            currentAngle = phi_angle
            if currentAngle==prs.FULL_ANGLE:
                curentAngle=-1
    secStructure.close()

    #Version 2 of secondary structure
    secStructure2 = open(file1.replace(".pdb","")+"_secondaryStructure2",'w')
    counter = 0
    prevStruct = ""

    for i in range(len(phifile)):
        line = phifile[i]
        if line[0:3]=="ASG":
            if line[12:15].strip()[-1] in string.ascii_uppercase:
                continue
            secStruct = line[33:42].strip()
            if secStruct == "aHelix" and secStruct == prevStruct:
                secStructure2.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+line[42:51].strip()+"\n")
            else:
                counter+=1
                secStructure2.write(line[5:8]+line[11:15].lstrip()+line[9]+"\t"+secStruct+str(counter)+"\t"+line[42:51].strip()+"\n")
            prevStruct = secStruct

    secStructure2.close()

def RSA(file1):
    '''
    @description
        Read phi file. Take solvent area values from phi file and divide by max solvent area of aminoacid
        Max solvent area taken from Rose 1985 and defined in parametres.py
    @input
        file1 - pdb file, ex. 4eiy.pdb
    @output 
        rsa file, ex. 4eiy_rsa 
    '''
    phifile = open(file1.replace('pdb', 'phi'),'r').readlines()
    with open(file1.replace('pdb', 'rsa'),'w') as rsafile:
        for line in phifile:
            if line[0:3] == "ASG" and line[5:8] in prs.area.keys():
                rsafile.write(line[5:8]+line[12:15].strip()+line[9]+"\t"+str(float(line[64:69].strip())/prs.area[line[5:8]])+"\n")


class AminoAcid:
    '''
    @description
        Define parameters of interactions for each aminoacid 
        - parameters, like coords of atoms involved in certain interaction, normalVecs, etc.
        - pipi, pication, disulf, salt bridge bonds options
    @input
        residue object from mda.Universe.residues 
    '''
    def __init__(self, aminoacid):        
        self.acid = aminoacid  # mda.residue object  
        self.res = aminoacid.resname+'-'+str(aminoacid.resid)+aminoacid.segid
        origin, vec = self.normalVecPIPI(self.acid)
        self.pipi = self.PiPi(self.acid, origin, vec)
        self.pication = self.PiCation(self.acid, origin, vec)
        self.disulf = self.Disulf(self.acid)
        self.sb = self.Salt_Bridge(self.acid)

    def normalVecPIPI(self, acid):
        '''
        @description
            create dictionaries of coords of center of geometry and normals to ring planes of residues
            these dicts used in pipi and pication bonds finding
        @input
            self.acid, mda.residue object
        '''
        origin_dict = {} # center of geometry
        vec_dict = {} # coords of normals
        # select residues and their atoms for calculations
        if acid.resname in ['HIS', 'HISD', 'HISE', 'HISH', 'HIS1']:
            select = [f'resname {acid.resname} and resid {acid.resid} and name CG ND1 CD2 CE1 NE2']
        elif acid.resname in ['TYR', 'PHE']:
            select = [f'resname {acid.resname} and resid {acid.resid} and name CG CD1 CD2 CE1 CE2 CZ']
        elif acid.resname in ['TRP']:
            select = [f'resname {acid.resname} and resid {acid.resid} and name CD2 CE2 CE3 CZ2 CZ3 CH2', 
                    f'resname {acid.resname} and resid {acid.resid} and name CG CD1 CD2 NE1 CE2']
        else:
            return [None, None]
        # calculate coordinates using mda.center_of_geometry and parameters.calc_norm_vecs() function 
        res_coms = [u.select_atoms(sel) for sel in select]
        origins = [res_com.center_of_geometry(compound='residues') for res_com in res_coms]
        vecs = [prs.calc_norm_vecs(acid, res_com, origin) for res_com, origin in zip(res_coms, origins)] 

        # TRP has hex and pent ring planes, so name them; for others - NORMAL
        if acid.resname in ['TRP']:
            for i, name in enumerate(['-HEX', '-PENT']):
                origin_dict[prs.getTotalResidue(res_coms[0], flag=False)[0]+name] = origins[i]
                vec_dict[prs.getTotalResidue(res_coms[0], flag=False)[0]+name] = vecs[i]
        else:
            origin_dict[prs.getTotalResidue(res_coms[0], flag=False)[0]+'-NORMAL'] = origins[0]
            vec_dict[prs.getTotalResidue(res_coms[0], flag=False)[0]+'-NORMAL'] = vecs[0]
        return origin_dict, vec_dict


    def PiPi(self, res, origin_dict, vec_dict):
        '''
        @description
            residues involved in pipi bonds are HIS, TYR, PHE, TRP
            return coords of their center of geometry and normalVecs
        @input
            - self.acid, mda.residue object
            - origin_dict, vec_dict - def normalVecPIPI() output
        '''
        if res.resname in prs.pication_aa:
            residue = res.resname + '-' + str(res.resid) + res.segid
            res_keys = [key for key in list(origin_dict.keys())]
            origin = [origin_dict[key] for key in res_keys]
            norm = [vec_dict[key] for key in res_keys]
            return res_keys, origin, norm
        else:
            return None

    def PiCation(self, res, origin_dict, vec_dict):
        '''
        @description
            residues involved in pi bonding are HIS, TYR, PHE, TRP + cation atoms NZ, CZ of ARG, LYS
            return coords of their center of geometry and normalVecs
        @input
            - self.res, mda.residue object
            - origin_dict, vec_dict - def normalVecPIPI() output
        '''
        if res.resname in prs.pication_aa:
            residue = res.resname+'-'+str(res.resid)+res.segid
            res_keys = [key for key in list(origin_dict.keys())]
            origin = [origin_dict[key] for key in res_keys]
            norm = [vec_dict[key] for key in res_keys]
            return res_keys, origin, norm

        if res.resname in prs.cation_aa:
            cations = res.atoms.select_atoms('name NZ CZ')
            return cations

    def Disulf(self, res):
        '''
        @description
            for CYS return coords of S* atoms
        @input
            - self.res, mda.residue object
        '''
        if res.resname in prs.disulf_aa:
            coord = res.atoms.select_atoms('name SG').positions
            return coord
        else:
            return None

    def Salt_Bridge(self, res):
        '''
        @description
            selection of basic (ARG, LYS, HIS) and acidic (ASP, GLU) residues involved in salt bridges
            return mda.atoms and coords of atoms N* and O*
        @input
            - self.res, mda.residue object
        '''
        sel_basic = f"(resname {' '.join(i for i in prs.basic_salt_aa)}) and (name NH* NZ* NE* ND*)"
        sel_acidic = f"(resname {' '.join(i for i in prs.acidic_salt_aa)}) and (name OE* OD*)"
        if res.resname in prs.basic_salt_aa:
            atoms = res.atoms.select_atoms(sel_basic)
            coords = atoms.positions
        elif res.resname in prs.acidic_salt_aa:
            atoms = res.atoms.select_atoms(sel_acidic)
            coords = atoms.positions
        else:
            return None
        return atoms, coords


def PiPi(res1, res2):
    '''
    @description
    	Find pipi interactions b/w atoms based on distance and angle cutoff
        - check if res1, res2 is not None in pipi bonding
        - for loop in res1, res2, find angle b/w normals and distance b/w centers
        - return list of pipi bonds under conditions (defined in home page)
    @input
        - res1, AminoAcid(mda.residue) class 
        - res2, AminoAcid(mda.residue) class
    @output 
        - list of acids 
    '''

    out = []
    if res1.pipi is None or res2.pipi is None:
        return None
    r1_key, origin1, norm1 = res1.pipi # name of res, origin and normal coords
    r2_key, origin2, norm2 = res2.pipi
    r1 = res1.acid.resname+' '+str(res1.acid.resid)+res1.acid.segid
    r2 = res2.acid.resname+' '+str(res2.acid.resid)+res2.acid.segid
    for i in range(len(r1_key)):
        for j in range(len(r2_key)):
            pipiangle = np.degrees(mda.lib.mdamath.angle(norm1[i][0], norm2[j][0]))
            dist = MDdist3D(origin1[i], origin2[j])
            # condition for angle and distance b/w ring planes
            if pipiangle > prs.RIGHT_ANGLE:
                pipiangle = prs.STRAIGHT_ANGLE-pipiangle
            if (dist <= prs.PIPI_D1 and pipiangle<prs.PIPI_ANG1) or (dist <=prs.PIPI_D2 and prs.PIPI_ANG2 < pipiangle < prs.PIPI_ANG3):
                out.append(''.join(r1.split(' '))+'\t'+''.join(r2.split(' '))+prs.PIPI_CONST+str(r1_key[i].split('-')[2])+'\t'+str(r2_key[j].split('-')[2])+'\n')
    return out

def PiPi_writing(file1):   
    '''
    Write pipi interactions into file (output of def PiPi)
    '''
    with open(file1.replace(".pdb","_pipi2"),'w') as out:
        for i in range(len(acids_class)):
            for j in range(i+1, len(acids_class)):
                result = PiPi(acids_class[i], acids_class[j])
                if result is not None:
                    # for TRP we have HEX and PENT output
                    for res in result:
                        out.write(res)

def PiCation(res1, res2):
    '''
    @description
    	Find pi-cation interactions b/w atoms based on distance and angle cutoff
        - check if res1, res2 is PI residues or contain cation atoms else return None
        - for loop in res1, res2, find angle and distance b/w
        - return list of pication bonds under conditions (defined in home page)
    @input
        - res1, AminoAcid(mda.residue) class 
        - res2, AminoAcid(mda.residue) class
    @return 
        - list of acids
    '''
    out = []
    # Delete this restriction bc return None, end) -> not workfirst return None for None values
#     if res1.pication is None or res2.pication is None:
#         return None
    if len(res1.pication) == 3 and len(res2.pication) == 1:
        r1_keys, origin, norm = res1.pication
        cations = res2.pication
    elif len(res1.pication) == 1 and len(res2.pication) == 3:
        r1_keys, origin, norm = res2.pication
        cations = res1.pication
    else:
        return None

    for i in range(len(r1_keys)):
        for cat in cations:
            catvec = cat.position - origin[i]
            catangle = np.degrees(mda.lib.mdamath.angle(norm[i][0], catvec[0]))
            if catangle > prs.RIGHT_ANGLE: #change 90 to prs.RIGHT_ANGLE
                catangle = prs.STRAIGHT_ANGLE-catangle
            if MDdist3D(origin[i][0], cat.position) < prs.PICAT_DIST and catangle < prs.PICAT_ANG:
                out.append(''.join(r1_keys[i].split('-')[:-1])+'\t'+(cat.resname+str(cat.resid)+cat.segid)+prs.PICATION_CONST+r1_keys[i].split('-')[2]+'\t'+cations.names[0]+'\n')
    return out

def PiCation_writing(file1):
    '''
    Write pipi interactions into file (output of def PiCation)
    '''
    with open(file1.replace(".pdb","_pication2"),'w') as out:
        for i in range(len(acids_class)):
            for j in range(i+1, len(acids_class)):
                result = PiCation(acids_class[i], acids_class[j])
                if result is not None:
                    # for TRP we have HEX and PENT output
                    for res in result:
                        out.write(res)

def SaltBridge(res1, res2):
    '''
    @description
    	Find salt bridges b/w atoms based on distance and angle cutoff
        - check if res1, res2 have salt bridges option
        - pairs, distances via mda.capped_distance (dist cutoff - SALT_distance)
        - for loop in pairs, and add to list
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @output 
        - list of detailed salt bridges
        - list of acids for writing to file  
    '''
    saltBridges = []
    out = []
    saltBridgesTmp = []
    # check not None residues and call coords, atoms involved in salt bridges
    if res1.sb is None or res2.sb is None:
        return None
    else:
        atoms1, coords1 = res1.sb
        atoms2, coords2 = res2.sb
    # check basic and acidic residues
    name1, name2 = res1.acid.resname, res2.acid.resname
    if all(name in prs.basic_salt_aa for name in [name1, name2]) or all(name in prs.acidic_salt_aa for name in [name1, name2]):
        return None
    if np.abs(res1.acid.resid - res2.acid.resid) == 1:
        return None
    # combinations of coords with distance cutoff
    pairs, distances = mda.lib.distances.capped_distance(coords1, coords2, max_cutoff=prs.SALT_DISTANCE, min_cutoff=0)

    if distances.size != 0:
        for k, [i, j] in enumerate(pairs):
            # use prs.rplsAtom, replace NH1 and NH2 to NH1/2, 
            atom1 = atoms1[i].resname+'-'+str(atoms1[i].resid)+atoms1[i].segid+'-'+prs.rplcAtom(atoms1[i].name)
            atom2 = atoms2[j].resname+'-'+str(atoms2[j].resid)+atoms2[j].segid+'-'+prs.rplcAtom(atoms2[j].name)
            # write atom1 & atom2 and atom2 & atom1
            saltBridges.append(atoms1[i].resname+'-'+str(atoms1[i].resid)+atoms1[i].segid+'-'+atoms1[i].name+'\t'+atoms2[j].resname+'-'+str(atoms2[j].resid)+atoms2[j].segid+'-'+atoms2[j].name)
            saltBridges.append(atoms2[j].resname+'-'+str(atoms2[j].resid)+atoms2[j].segid+'-'+atoms2[j].name+'\t'+atoms1[i].resname+'-'+str(atoms1[i].resid)+atoms1[i].segid+'-'+atoms1[i].name)
            # add to out list non repeating residues with same ids
            # ARG-1034A ASP-1073A must appear once
            # out.append(atom1+'\t'+atom2+prs.SALT_CONST+prs.rplcAtom(atoms1[i].name)+'\t'+prs.rplcAtom(atoms2[j].name)+'\t'+'\n')
            if not [''.join(atom1.split('-')[0:2]), ''.join(atom2.split('-')[0:2])] in saltBridgesTmp and not [''.join(atom2.split('-')[0:2]), ''.join(atom1.split('-')[0:2])] in saltBridgesTmp:
                out.append(atom1+'\t'+atom2+prs.SALT_CONST+prs.rplcAtom(atoms1[i].name)+'\t'+prs.rplcAtom(atoms2[j].name)+'\t'+'\n')
                saltBridgesTmp.append([''.join(atom1.split('-')[0:2]), ''.join(atom2.split('-')[0:2])])
                saltBridgesTmp.append([''.join(atom2.split('-')[0:2]), ''.join(atom1.split('-')[0:2])])
        return list(set(out)), list(set(saltBridges))

        
def SB_writing(file1):
    '''
    Write salt bridges into file (output of def SaltBridge)
    '''
    saltBridges = []
    with open(file1.replace(".pdb","_SaltBridges_Barlow"),'w') as out:
        for i in range(len(acids_class)):
            for j in range(i+1, len(acids_class)):
                result = SaltBridge(acids_class[i], acids_class[j])
                if result is not None:
                    # hydrogen bonds that are in salt bridges will be deleted
                    for res in result[1]:
                        saltBridges.append(res)
                    for res in result[0]:
                        out.write(res)
    return saltBridges

def Disulfide(res1, res2):
    '''
    @description
        for atoms of CYS find disulfide bonds based on distance cutoff
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @return 
        string of disulfide bonds
    '''
    # check if residues have disulfide params and return coords of SG atoms
    if res1.disulf is None or res2.disulf is None:  
        return None
    else:   
        coord1, coord2 = res1.disulf, res2.disulf
    
    dist = MDdist3D(coord1, coord2)
    if dist < prs.DISULF_D:
        r1 = res1.res
        r2 = res2.res
        out = r1+'\t'+r2+prs.DISULF_CONST+'SG'+'\t'+'SG'+'\n'
        return out

def Disulfide_writing(file1):
    '''
    Write disulfide bonds into file (output of def Disulfide)
    '''
    with open(file1.replace(".pdb","_disulf1"),'w') as out:
        for i in range(len(acids_class)):
            for j in range(i+1, len(acids_class)):
                result = Disulfide(acids_class[i], acids_class[j])
                if result is not None:
                    out.write(result)

def Hydrogen_calculation():
    '''
    @description
    	Find hydrogen bonds in created hydrogen table by MDAnalysis
        - define saltBridges, if hydrogen atoms in salt bridges, continue
        - iterate through hydrogen table 
        - define donor (X-H) and acceptor (---A) parameters
    @input
    @output 
        - list of acids 
        - write file, ex. 4eiy_hb file 
    '''
    saltBridges = SB_writing(file1)
    hydrogenBonds = []
    acceptorBonds = {}
    for hbond in h.table:
        # here donor is hydrogen atom
        donor = hbond[3]+'-'+str(hbond[4])+segids[hbond[1]]+'-'+hbond[5]
        acceptor = hbond[6]+'-'+str(hbond[7])+segids[hbond[2]]+'-'+hbond[8]
        dH = hbond[3]+'-'+hbond[5]
        if not dH in prs.hydrogen:
            continue
        # replace in donor hydrogen atom
        dA = hbond[3]+'-'+str(hbond[4])+segids[hbond[1]]+'-'+prs.hydrogen[dH]
        if not prs.isDonor(dA):
            continue
        # check whether bond not in saltbridges
        if (dA+'\t'+acceptor) in saltBridges or (acceptor+'\t'+dA) in saltBridges:
            continue
        # don't include neighboring acids
        if donor.split('-')[:2] != acceptor.split('-')[:2] and (np.abs(hbond[4]-hbond[7]) != 1):
            # '-'.join(acceptor.split('-')[::2]) in nextAcids and check isAcceptor
            if (acceptor.split("-")[0]+"-"+acceptor.split("-")[2] in prs.nextAcid) and prs.isAcceptor(acceptor):
                d1 = hbond[9] # H-A; hydrogen-acceptor distance 
                if d1 > 3.5:
                    continue
                d2 = MDdist3D(coords[dA], coords[acceptor]) # donor-acceptor distance
                if d2 > 4.5:
                    continue

                alpha = hbond[10] # DHA angle > 120
                aNext = prs.nextAcid[acceptor.split("-")[0]+"-"+acceptor.split("-")[2]] # next acid close to acceptor
                # acceptor can have some number of bonds
                if len(aNext)==2:
                    neighbor1 = "-".join(acceptor.split("-")[0:2])+"-"+aNext[0]
                    neighbor2 = "-".join(acceptor.split("-")[0:2])+"-"+aNext[1]
                    A1 = (coords[neighbor1]+coords[neighbor2])/2
                else:
                    A1 = coords["-".join(acceptor.split("-")[0:2])+"-"+aNext[0]]

                a = A1 - coords[acceptor]
                b = coords[donor] - coords[acceptor]
                beta = np.degrees(mda.lib.mdamath.angle(a, b)) # angle b/w A1-acceptor-donor (A1-acceptor neighbor)

                a = A1 - coords[acceptor]
                b = coords[dA] - coords[acceptor]
                gamma = np.degrees(mda.lib.mdamath.angle(a, b)) # angle b/w A1-acceptor-dA

                if [(prs.hydrogen[dH]!='S' and hbond[8]!='S' and d1 < 2.5 and d2 < 3.9 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]=='S' and hbond[8]!='S' and d1 < 2.5 and d2 < 3.9 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]!='S' and hbond[8]=='S' and d1 < 3.2 and d2 < 4.1 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]=='S' and hbond[8]=='S' and d1 < 3.3 and d2 < 4.5 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE)]:
                    #hbond[8] - acceptor atom
                    #hydrogen[dH] - hydrogen
                    if not acceptor in acceptorBonds.keys():
                        acceptorBonds[acceptor] = []
                    if hbond[8] in ["OD1","OD2","OE1","OE2","OG","OG1","SG"]: 
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 2:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:2]
                    if hbond[8] in ["O","ND1","NE2","OH"] and hbond[6]!="HOH": 
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 1:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:1] 
                    if hbond[8] == "O" and hbond[6]=="HOH":
                        acceptorBonds[acceptor].append(d1)
                    """
                    Energy for various hybridizations
                    Ef = V*[5*(do/d)^12 - 6*(do/d)^10]*F where do is 2.8, d is the distance between donor and acceptor atom (ie d2), V is 8 
                    sp3 donor - sp3 acceptor F=cos^2(alpha)cos^2(beta-109.5)
                    sp3 donor - sp2 acceptor F=cos^2(alpha)cos^2(beta)
                    sp2 donor - sp3 acceptor F=cos^4(alpha)
                    sp2 donor - sp2 acceptor F=cos^2(alpha)cos^2(max[alpha,beta])
                    """
                    if beta < 90: # prs.STRAIGHT_ANGLE
                        beta = prs.STRAIGHT_ANGLE-beta if beta>prs.RIGHT_ANGLE else beta

		    # define Hydrogen Energy function in parameters.py
                    if prs.SPHyb(dA)=="SP3" and prs.SPHyb(acceptor)=="SP3":
                        E = prs.energy_hb(d2, 'SP3', 'SP3', alpha, beta, 0)
                        # E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(beta-109.5)))**2
                    elif prs.SPHyb(dA)=="SP3" and prs.SPHyb(acceptor)=="SP2":
                        E = prs.energy_hb(d2, 'SP3', 'SP2', alpha, beta, 0)
                        # E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(beta)))**2
                    elif prs.SPHyb(dA)=="SP2" and prs.SPHyb(acceptor)=="SP3":
                        E = prs.energy_hb(d2, 'SP2', 'SP3', alpha, beta, 0)
                        # E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**4) 
                    elif prs.SPHyb(dA)=="SP2" and prs.SPHyb(acceptor)=="SP2":
                        # for sp2-sp2 need to find angle psi b/w normals 
                        normalVecDonor = prs.normalDonorVecToPlane1(dA, coords)
                        if normalVecDonor is None: # change to avoid warning
                            continue
                        normalVecAcceptor = prs.normalAcceptorVecToPlane1(acceptor, coords)
                        if normalVecAcceptor is None: # change to avoid warning
                            continue
                        psi = np.degrees(mda.lib.mdamath.angle(normalVecDonor,normalVecAcceptor))
                        E = prs.energy_hb(d2, 'SP2', 'SP2', alpha, beta, psi)
                        # E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)*((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(max([beta,psi]))))**2
                    interaction = donor+"\t"+acceptor+"\t"+dA+chain[dA]+chain[acceptor]+"\t"+str(d1)+"\t"+str(d2)+"\t"+str(alpha)+"\t"+str(beta)+"\t"+str(gamma)+"\t"+dA+"\t"+str(E)+"\t"+donor.split("-")[2]+"\t"+acceptor.split("-")[2]+"\n"
                    hydrogenBonds.append(donor+"\t"+acceptor+"\t"+chain[dA]+chain[acceptor]+"\t"+str(d1)+"\t"+str(d2)+"\t"+str(alpha)+"\t"+str(beta)+"\t"+str(gamma)+"\t"+dA+"\t"+str(E)+"\t"+donor.split("-")[2]+"\t"+acceptor.split("-")[2]+"\n")
                    
    out = open(file1.replace(".pdb","_hb"),'w')
    finalHydrogenBonds = []
    donorBonds = {}
    for interaction in hydrogenBonds:
        donor = interaction.split("\t")[0]
        dH = donor.split("-")[0]+"-"+donor.split("-")[2]
        dA = donor.split("-")[0]+"-"+donor.split("-")[1]+"-"+prs.hydrogen[dH]
        acceptor = interaction.split("\t")[1]
        d1 = interaction.split("\t")[3]
        if not dA in donorBonds:
            donorBonds[dA] = 0
        #if water, allow all interactions
        if dA.split("-")[0]=="HOH" or acceptor.split("-")[0]=="HOH":
            out.write(interaction)
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
    out.close()

def VdW_calculation():
    '''
    @description
        - select atoms and their coords
        - capped distance with restricted distance (4.25)
        - iterate through atoms and find distance between atoms and radii of these atoms 
        - if it's difference less 0.5 and not from neighboring acids and not from main chain
        - find vdw energy 
    @input
    @return 
        - list of acids 
        - files, ex. 4eiy_vdw and 4eiy_vdw2
    '''

    allatoms = list(prs.getTotalResidue(u) + names) 
    van = u1.atoms.positions
    # capped distnace with distance cutoff based on max vdw radii
    pairs, distances = mda.lib.distances.capped_distance(van, van, max_cutoff=4.25, min_cutoff=0)
    van_atom = np.split(pairs[:, 1], np.cumsum(np.unique(pairs[:, 0], return_counts=True)[1])[:-1])
    scores=[]
    scores2 = []
    vdw1 = {}
    vdw2 = {}
    output1 = []
    # iterate through combinations of atoms and find bonds under condtions
    for i, k in enumerate(van_atom):
        elem1 = allatoms[i]
        res1 = resnames[i]+' '+str(resids[i])+' '+segids[i]
        atom1 = names[i]
        for j in k:
            elem2 = allatoms[j]
            res2 = resnames[j]+' '+str(resids[j])+' ' +segids[j]
            atom2 =names[j]
            if res2+':'+res1+':'+chain[elem2]+":"+chain[elem1]+":"+atom2+":"+atom1 in vdw1:
                continue
            if (not res1==res2 and atom1 in prs.radii and atom2 in prs.radii):
                rm = prs.radii[atom1] + prs.radii[atom2] # sum vdw radii of 2 atoms
                r = MDdist3D(coords[elem1], coords[elem2]) # distance of 2 atoms

                if r - rm < .5 and not (prs.isNeighbor(elem1,elem2) and chain[elem1]=="MC" and chain[elem2]=="MC"):
                    if not  res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2 in vdw1:
                        vdw1[res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2] = []
                    E = (-.997*((rm/r)**12-2*(rm/r)**6))*4
                    # E = prs.energy_vdw(rm , r)
                    vdw1[res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2].append(E)
                    # scores.append(E)
                    scores.append(sum(vdw1[res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2]))

                    if ("C" in atom1 and "C" in atom2) or ("C" in atom1 and atom2 in ["NE2","OE1","ND2","OD1"] and res2.split(" ")[0] in ["GLN","ASN"]) or (atom1 in ["NE2","OE1","ND2","OD1"] and res1.split(" ")[0] in ["GLN","ASN"] and "C" in atom2):
                        if not res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2 in vdw2:
                            vdw2[res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2] = []
                        vdw2[res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2].append(E)
                        scores2.append(sum(vdw1[res1+":"+res2+":"+chain[elem1]+":"+chain[elem2]+":"+atom1+":"+atom2]))

    with open(file1.replace(".pdb","_vdw"),'w') as out:
        for contact in vdw1:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split(':')[0].split(' ')[1]) - int(contact.split(':')[1].split(' ')[1]))==1):
                out.write(''.join(contact.split(' '))+"\t"+str(sum(vdw1[contact]))+"\t"+contact.split(":")[4]+"\t"+contact.split(":")[5]+"\n")
    
    with open(file1.replace(".pdb","_vdw2"),'w') as out:
        for contact in vdw2:
            if not (sum(vdw2[contact]) < 0 and abs(int(contact.split(':')[0].split(' ')[1]) - int(contact.split(':')[1].split(' ')[1]))==1):
                out.write(''.join(contact.split(' '))+"\t"+str(sum(vdw1[contact]))+"\t"+contact.split(":")[4]+"\t"+contact.split(":")[5]+"\n")

def Calc_MetalBonds(met, AA1):
    '''
    @description
    	find metal bonds, distance b/w metal and acid atoms
	for each acid take own distance from dictionary
    @input 
    	met - mda.select of metalls
    	AA1 - aminoacids class
    @output
	list of residues involved in metal bonds
    '''
    metal2atom = []
    acid = AA1.acid
    metal = met.resname+'-'+str(met.resid)+met.segid+'-'+met.name
    r, resid, segid = acid.resname, acid.resid, acid.segid
    res = r + '-' + str(resid) + segid + '-'

    for atom in acid.atoms:
        a = atom.name
        distance = MDdist3D(coords[metal], coords[res+a])
        if distance < 3.5 and r in prs.METALLS_DICT:
            dist, atoms2 = prs.METALLS_DICT[r]
            if met.resname in dist:
                d = dist[met.resname]
                if r == 'ASP' and a in atoms2:
                    if distance < prs.GLUASPmonoDist[met.resname]:
                        metal2atom.append(res+a)
                    tmpatom21 = res+'OD1'
                    tmpatom12 = res+'OD2'
                    if MDdist3D(coords[res+a], coords[tmpatom21]) <= prs.GLUASPbiDist[met.resname] and MDdist3D(coords[res+a], coords[tmpatom12]) <= prs.GLUASPbiDist[met.resname]:
                        metal2atom.append(res+a)
                
                elif r == 'GLU' and a in atoms2:
                    if dist3D(coords[metal], coords[res+a]) < prs.GLUASPmonoDist[met.resname]:
                        metal2atom.append(res+a)
                    tmpatom21 = res+'OE1'
                    tmpatom12 = res+'OE2'
                    if MDdist3D(coords[res+a], coords[tmpatom21]) <= prs.GLUASPbiDist[met.resname] and MDdist3D(coords[res+a], coords[tmpatom12]) <= prs.GLUASPbiDist[met.resname]:
                        metal2atom.append(res+a)   

                elif distance <= d and a in atoms2:
                    metal2atom.append(res+a)
    return metal2atom

def write_metalbonds(met, metal2atom):
	'''
	@description
	    Iterate through the list of residues and find distances b/w metal atom
	@input
	    met - mda.selection of metalls 
	    metal2atom - list of residues involved in metal bonds
	'''
    metal = met.resname+'-'+str(met.resid)+met.segid+'-'+met.name
    metalBondsFiltered = []
    origin = met.position
    for atom21 in metal2atom:
        atom21coords = coords[atom21]
        atom21coords = atom21coords - origin
        for atom22 in metal2atom:
            atom22coords = coords[atom22]
            atom22coords = atom22coords - origin
            if atom21 != atom22:
                ang = np.degrees(mda.lib.mdamath.angle(atom21coords, atom22coords))
                if atom22+"\t"+atom21+"\t"+chain[atom22]+chain[atom21]+"\t"+metal+"\t"+str(ang)+'\n' in metalBondsFiltered:
                    continue
                if atom21+"\t"+atom22+"\t"+chain[atom21]+chain[atom22]+"\t"+metal+"\t"+str(ang)+'\n' in metalBondsFiltered:
                    continue
                if met.name in ["NA","K","FE","CU","ZN"] and prs.RIGHT_ANGLE<ang<prs.STRAIGHT_ANGLE:
                    metalBondsFiltered.append(atom21+"\t"+atom22+"\t"+chain[atom21]+chain[atom22]+"\t"+metal+"\t"+str(ang)+'\n')
                elif ang > prs.RIGHT_ANGLE:
                    metalBondsFiltered.append(atom21+"\t"+atom22+"\t"+chain[atom21]+chain[atom22]+"\t"+metal+"\t"+str(ang)+'\n')
    return metalBondsFiltered

def MetalBonds_calculation1():
    '''
    @description
    	iterate through acid classes and metalls and find metal bonds
    @parameteres
    	input pdb file defined in __main__, ex. 4eiy.pdb
    @return
    	file, ex. 4eiy_metal
    '''
    metals = prs.metals # from parameters
    met_select = 'resname {}'.format(' '.join(list(metals)))
    metalls = u.select_atoms(met_select)
    metal2atom = []
    output = []
    for met in metalls:
        for i in range(len(acids_class)):
            result = Calc_MetalBonds(met, acids_class[i])
            if len(result) != 0:
                for res in result:
                    metal2atom.append(res)
        output.append(write_metalbonds(met, metal2atom))

    out = open(file1.replace(".pdb","_metal"),'w')
    for bond in output:
        for i in bond:
            out.write(i)
    out.close()

def DNA_bonds(file1):
    '''
    @description
    	choose DA, DG, DC, DT nucleotides and among neighbors find bonds
    @input
    	pdb file, return ex. 4eiy_DNA file
    '''
    dna = u.select_atoms('resname DA DG DC DT') # select dna acids
    dna_atoms = prs.getTotalResidue(dna) + dna.atoms.names # dna atoms
    DNAbindingPairs = []
    for nucleotide in dna_atoms:
        for atom2 in neighboringAtoms:
            if not prs.isDNA(atom2.split("-")[0]) and not "HOH" in atom2: # atom2 not be DNA and HOH
                bound = False
                # check if atom2 and nucleotide interact
                if (("N" in nucleotide.split("-")[2] and (prs.isDonor(atom2) or prs.isAcceptor(atom2))) or 
                    ("O"in nucleotide.split("-")[2] and (prs.isDonor(atom2) or prs.isAcceptor(atom2))) or 
                    (("N" in nucleotide.split("-")[2] or "O" in nucleotide.split("-")[2]) and "S" in atom2.split("-")[2]) 
                    or ("C" in nucleotide.split("-")[2] and "O" in atom2.split("-")[2] and prs.isAcceptor(atom2)) or 
                    ("O" in nucleotide.split("-")[2] and atom2.split("-")[2] in ["NZ","NH1","NH2","OD1","OD2","OE1","OE2"]) or 
                    ("C" in nucleotide.split("-")[2] and "C" in atom2.split("-")[2]) ):

                    for atom3 in neighboringAtoms[nucleotide]: # atom3 - neighbor of nucleotide
                        if atom3==atom2 or prs.isH(atom3): # check that atom2!=atom3 and atom3 is hydrogen atom
                            continue
                        vector1 = coords[nucleotide] - coords[atom3]
                        vector2 = coords[atom2] - coords[atom3]
                        omega = np.degrees(mda.lib.mdamath.angle(vector1, vector2))
                        if omega <= prs.RIGHT_ANGLE:
                            bound = True
                        else:
                            bound = False
                            break
                if bound:
                        DNAbindingPairs.append([nucleotide,atom2])

    with open(file1.replace('.pdb', '_DNA'), 'w') as out:
        for nucleotide, atom in DNAbindingPairs:
            out.write(nucleotide+'\t'+atom+'\t'+chain[atom]+'\n')

def Ligand_bonds(file1):
    '''
	@description
	    iterate through protein atoms and ligands, finding distance b/w them and write to file
    	replacing distance to less value try to find the less distance b/w atom and ligand
    '''
    atoms = list(prs.getTotalResidue(protein) + (protein).atoms.names)
    with open(file1.replace(".pdb","_ligand"),'w') as out:
        for atom in atoms:
            if chain[atom] == 'SC':
                tmpdist0 = 150
                for ligand in sorted(ligandCentroids.keys()):
                    tmpdist1 = MDdist3D(coords[atom], ligandCentroids[ligand])
                    if tmpdist1 > tmpdist0:
                        continue
                    out.write(atom.split("-")[0]+atom.split("-")[1]+"-"+atom.split("-")[2]+"\t"+ligand+"\t"+str(tmpdist1)+"\n")
                    tmpdist0 = float(copy.copy(tmpdist1))

def Centroids():
    '''
	@description
   		for loop centroid data, find distane b/w and if it < 8.5
    	not include 2 interacting HOH and same residues
    '''
    with open(re.sub(".pdb","",file1)+"_centroidNetSC",'w') as out:
        for residue1 in sorted(centroid.keys()):
            for residue2 in sorted(centroid.keys()):
                if "HOH" in residue1 and "HOH" in residue2:
                    continue
                if residue1==residue2:
                    continue
                x = MDdist3D(np.array(centroid[residue1]),np.array(centroid[residue2])) 
                if x < 8.5:
                    out.write(re.sub("-","",residue1)+"\t"+re.sub("-","",residue2)+"\t"+"SCSC"+"\t"+"1"+"\t"+"CENTROID\t"+"CENTROID\t"+"CENTROID\t"+str(x)+"\n")
    '''
    find the shortest distance b/w centroid and ligand
    '''
    with open(re.sub(".pdb","",file1)+"_centroidNetLigand",'w') as out:
        for residue1 in centroid.keys():
            if "HOH" in residue1:
                continue
            tmpdist0 = 150
            for ligand in sorted(ligandCentroids.keys()):
                tmpdist1 = MDdist3D(np.array(centroid[residue1]),ligandCentroids[ligand])
                if tmpdist1 > tmpdist0:
                    continue
                out.write(re.sub("-","",residue1)+"\t"+ligand+"\t"+str(tmpdist1)+"\n")
                tmpdist0 = float(copy.copy(tmpdist1))   



def pdb2peptide(file1):
    '''
    @description
    	Go through ATOM lines in pdb file and define 2 variables
    	if acids in one chain write it else not write Write data like 
    @input
    	input pdb file, ex. 4eiy.pdb
    @output
    	file, ex. 4eiy_polypeptidePairs
    ''' 
    total_res = []
    for res, id, seg in zip(u1.atoms.resnames, u1.atoms.resids, u1.atoms.segids):
        out = res + str(id) + seg
        total_res.append([res, id, seg])
        
    with open(file1.replace(".pdb","")+"_polypeptidePairs",'w') as out:
	    for i in range(len(total_res) - 1):
	        if total_res[i][1] == total_res[i+1][1] - 1:
	            out.write(''.join([str(elem) for elem in total_res[i]]) +"\t"+''.join([str(elem) for elem in total_res[i+1]]) +"\tSCORE\n")        

def Bfactor(file1):
	'''
	@description 
		select CA atoms and tempfactors
		write it to file, ex. 4eiy_Bfactor file 
	'''
	CA_atoms = u.select_atoms('name CA')	
	with open(file1.replace('.pdb', '_Bfactor'),'w') as bfact:
		for atom in CA_atoms.atoms:
		    bfact.write(atom.resname+str(atom.resid)+atom.segid+"\t"+str(atom.tempfactor)+"\n")


if __name__ == '__main__':
	t0 = time.time()
	file1 = sys.argv[1] # input file.pdb
	
	u = mda.Universe(file1)
	u1 = u.select_atoms('protein and not name OXT')

	acids_class = [] # for each residue create class AA1 with parameters
	for res in u1.residues:
	    acids_class.append(AminoAcid(res))
	print('Finish creation of acids class', time.time() - t0)

	resnames = u.atoms.resnames
	resids = u.atoms.resids
	segids = u.atoms.segids
	names = u.atoms.names
	coords = dict(zip(prs.getTotalResidue(u)+names, u.atoms.positions))
	chain = {} 
	for res, atom in zip(prs.getTotalResidue(u), names):
	    if atom in ["N","O","C","CA","HA2","HA3"]:
	        if res[0:3] == 'GLY' and atom in ["O","CA","HA2","HA3","N","NH"]:
	            chain[res+atom] = 'SC'
	        else:
	            chain[res+atom] = 'MC'
	    else:
	        chain[res+atom] = 'SC'

	h = hbonds(u, selection1='protein', selection2= 'protein', distance=2.5, distance_type='hydrogen', angle=120) 
	h.run()
	h.generate_table() 
	print('Generate hydrogen table + define chain list', time.time()-t0)

	# define some common selections 
	hoh_atoms =  u.select_atoms('resname HOH')
	metalls = u.select_atoms('resname {}'.format(' '.join(list(prs.metals))))
	protein = u.select_atoms('protein')
	allatoms = u.select_atoms('all')
	# allAtoms = list(prs.getTotalResidue(atoms) + atoms.atoms.names)
	# define ligands as rest of allatoms - water, Me, protein
	ligands = allatoms - protein - hoh_atoms - metalls
	ligandCentroids = ligands.center_of_geometry(compound='residues')
	ligandCentroids  = dict(zip(prs.getTotalResidue(ligands, flag=False), ligandCentroids))

	centroid_data = u.select_atoms('protein and not resname DG DC DT DA and not backbone or (resname GLY and not name N C O)') + hoh_atoms
	# calc weighted center of residues, where centroid_data.masses are weights
	center = centroid_data.center(centroid_data.masses, compound='residues')
	centroid = dict(zip(prs.getTotalResidue(centroid_data, flag=False), center))

	print('Create acid and ligand centroids', time.time()-t0)

	Bfactor(file1)
	print('Create bfactor file', time.time()-t0)

	SecondaryStructure(file1)
	print('End _secondaryStructure', time.time()-t0)

	RSA(file1)
	print('End RSA', time.time()-t0)

	PiPi_writing(file1)
	print('Found pipi interactions', time.time()-t0)

	PiCation_writing(file1)
	print('Found pication interactions', time.time()-t0)

	SB_writing(file1)
	print('Found salt bridges', time.time()-t0)

	Disulfide_writing(file1)
	print('Found disulfide bonds', time.time() -t0)

	MetalBonds_calculation1()
	print('Found metal bonds', time.time()-t0)

	VdW_calculation()
	print('Found Van der Waals bonds', time.time()-t0)

	Hydrogen_calculation()
	print('Found hydrogen bonds', time.time()-t0)

	DNA_bonds(file1)
	print('Found dna bonds', time.time()-t0)
	
	Ligand_bonds(file1)
	Centroids()
	print('Found ligand and centroid bonds', time.time()-t0)

'''Extra notes 

'''

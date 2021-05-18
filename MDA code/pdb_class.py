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
import MDAnalysis.lib.mdamath as MDAmath
import time
import parameters as prs 

# TODO: make preprocessing script + installation readme
'''
    Before need to preprocess pdb file
    1. Add hydrogens (phenix.reduce)
    2. Protonate water (phenix.ready_set)
    3. And create phi file (stride)
'''

def prepare_secondary_structure_file(PDB, phi_data):
    '''
    @description
        phi_file contains info about type of secStructures and phi-angles for amino acids. These data defined by 'ASG' lines.
        1. Iterate phi_file and define prevStruct as first ASG line
        2. If it is Coil or Bridge change to CoilA, Turn -> TurnA
        3. Iterate in phi file, define residue + id + segid, phi angle and name of secstructure
        4. If residue is PRO and structure is not aHelix and (next line not last and next structure is not aHelix), secStruct=new struct
        5. Rename structure
    @input
        phi_data - processed PHI file, list of ASG lines
    @output 
        secondary structure file (ex. 4eiy_secondaryStructure) in format 
            Acid_info | Structure_type | phi_angle
            HIS26A    | Coil           | 90
    '''
    # secStructure = open('SecondaryStructure', 'w')
    secStructure = open(PDB.replace('.pdb', '_secondaryStructure2'), 'w')

    phi = {} 
    counter = 0
    currentSecStructure = ''
    currentAngle = -1

    for i in range(len(phi_data)):
        #here point i == 0 to maintain prevStruct and not to recall for each i
        if i == 0:  prevStruct = phi_data[0][33:42].strip()
        if prevStruct == 'Coil' or prevStruct == 'Bridge':  prevStruct = 'CoilA'
        elif prevStruct == 'Turn':  prevStrcut = 'TurnA'

        line = phi_data[i]   
        res, resid, segid = line[5:8].strip(), line[11:15].strip(), line[9]
        phi_angle, structure = float(line[42:50].strip()), line[33:42].strip()
        phi[res+resid+segid] = phi_angle
        secStruct = structure

        if res == 'PRO' and structure != 'aHelix':
            if i+1 != len(phi_data) and phi_data[i+1][33:42].strip() != "aHelix":
                secStruct = phi_data[i+1][33:42].strip()

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
                  (phi_angle==prs.FULL_ANGLE and currentAngle > prs.ZERO_ANGLE)):

                if prevStruct == "TurnA":
                    secStruct="TurnB"; prevStruct = "TurnB"
                else:
                    secStruct="TurnA"; prevStruct = "TurnA"

        if ("Coil" in secStruct or "Turn" in secStruct) and (res=="GLY" and phi_angle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE):
            phiangle = phi_angle
            # secStructure.write(res+resid+segid+"\t"+secStruct+str(counter)+"\t"+str(phiangle)+"\n")
            secStructure.write(res+':'+resid+':'+segid+'\t'+secStruct+str(counter)+'\t'+str(phiangle)+'\n')
            currentSecStructure = secStruct                
            counter+=1
        elif secStruct != currentSecStructure:
            counter+=1
            phiangle = line[42:51].strip()
            if float(phiangle)==prs.FULL_ANGLE:
                phiangle=0
            # secStructure.write(res+resid+segid+"\t"+secStruct+str(counter)+"\t"+str(phiangle)+"\n")
            secStructure.write(res+':'+resid+':'+segid+'\t'+secStruct+str(counter)+'\t'+str(phiangle)+'\n')
            currentSecStructure = secStruct
        else:
            # secStructure.write(res+resid+segid+"\t"+secStruct+str(counter)+"\t"+str(phi_angle)+"\n")
            secStructure.write(res+':'+resid+':'+segid+'\t'+secStruct+str(counter)+'\t'+str(phiangle)+'\n')
        currentAngle = phi_angle
        if currentAngle==prs.FULL_ANGLE:
            curentAngle=-1
    secStructure.close()

def prepare_rsa_file(PDB, phi_data):
    '''
    @description
        Read phi file. Take solvent area values from phi file and divide by max solvent area of aminoacid
        Max solvent area taken from Rose 1985(https://sci-hub.se/10.1126/science.4023714) and defined in parametres.py
    @input
        phi_data - processed PHI file, list of ASG lines
    @output 
        rsa file (ex. 4eiy_rsa) in format
            Acid_info | RSA_value 
            HIS26A    | 0.47      
    '''
    with open(PDB.replace('.pdb','.rsa'), 'w') as rsafile:
        for line in phi_data:
            if line [5:8] in prs.area.keys():
                # rsafile.write(line[5:8]+line[11:15].strip()+line[9]+"\t"+str(float(line[64:69].strip())/prs.area[line[5:8]])+"\n") 
                rsafile.write(line[5:8]+':'+line[11:15].strip()+':'+line[9]+'\t'+str(float(line[64:69].strip())/prs.area[line[5:8]])+'\n') 

class AminoAcid:
    '''
    @description
        Define parameters of interactions for each aminoacid 
        - parameters, like coords of atoms involved in certain interaction, normalVecs, etc.
        - pipi, pication, disulf, salt bridge bonds options
    @input
        aminoacid - MDA residue object
    '''
    def __init__(self, aminoacid):        
        self.acid = aminoacid  # mda.residue object  
        self.resname = aminoacid.resname
        self.atoms = aminoacid.atoms
        self.res = aminoacid.resname+':'+str(aminoacid.resid)+':'+aminoacid.segid   #HIS:26:A
        self.centroid, self.normal = self.calc_centroid_normal_coords()
        self.pipi = self.PiPi()
        self.pication = self.PiCation()
        self.disulf = self.Disulf()
        self.sb = self.SaltBridge()

    def calc_centroid_normal_coords(self):
        '''
        @description
            Consider amino acids involved in pi interactions (HIS, TYR, PHE, TRP).
            Select acid's ring atoms (atom_groups) and for ring plane calculate coords of center of geometry (centroid) and normal vecs to ring plane of acid.

        '''
        centroid_dict = {} 
        normal_dict = {} 

        # TODO: be consistent with other functions, e.g. self.resname in prs.pication_aa
        if self.resname in prs.acids_set['HIS']:
            select = ['name CG ND1 CD2 CE1 NE2']
        elif self.resname in prs.acids_set['TYR']+prs.acids_set['PHE']:
            select = ['name CG CD1 CD2 CE1 CE2 CZ']
        elif self.resname in prs.acids_set['TRP']:
            select = ['name CD2 CE2 CE3 CZ2 CZ3 CH2', 'name CG CD1 CD2 NE1 CE2']
        else:
            return None, None

        atom_groups = [self.atoms.select_atoms(sel) for sel in select]
        centroids = [atom_group.center_of_geometry(compound='group') for atom_group in atom_groups]
        normals = [prs.calc_norm_vecs(atom_group, centroid) for atom_group, centroid in zip(atom_groups, centroids)] 

        if self.resname in prs.acids_set['TRP']:
            for i, name in enumerate([':NORMAL', ':PENT']):
                centroid_dict[self.res+name] = centroids[i]
                normal_dict[self.res+name] = normals[i]
        else:
            centroid_dict[self.res+':NORMAL'] = centroids[0]
            normal_dict[self.res+':NORMAL'] = normals[0]
            #{HIS:26:A:HEX : [x, y, z]}
        return centroid_dict, normal_dict


    def PiPi(self):
        '''
        @description
            Consider only amino acids involved in pi interactions (HIS, TYR, PHE, TRP).
            For acid return coords of self.centroid and self.normal  
        '''
        if self.resname in prs.pication_acids:
            resname_key = list(self.centroid.keys())    #['HIS:26:A:NORMAL']
            centroid = list(self.centroid.values())
            norm = list(self.normal.values())
            return resname_key, centroid, norm

    def PiCation(self):
        '''
        @description
            Consider amino acids having cations NZ, CZ (ARG, LYS).
            For acid return cation atom names and positions
        '''
        if self.resname in prs.cation_acids:
            cations = self.atoms.select_atoms('name NZ CZ')
            cation_names = self.res+':'+cations.names[0]
            # cation_names = cations.resnames[0]+':'+str(cations.resids[0])+':'+cations.segids[0]+':'+cations.names[0]
            cation_coord = cations.positions[0]
            return cation_names, cation_coord

    def Disulf(self):
        '''
        @description
            for CYS return coords of S* atoms
        @input
            - self.res, mda.residue object
        '''
        if self.resname in prs.disulf_aa:
            atom = self.atoms.select_atoms('name SG')
            atom_name = self.res+':'+atom.names[0]
            # atom_name = atom.resnames[0]+':'+str(atom.resids[0])+':'+atom.segids[0]+':'+atom.names[0]
            atom_coords = atom.positions[0]
            return atom_name, atom_coords

    def SaltBridge(self):
        '''
        @description
            selection of basic (ARG, LYS, HIS) and acidic (ASP, GLU) residues involved in salt bridges
            return mda.atoms and coords of atoms N* and O*
        @input
            - self.res, mda.residue object
        '''
        if self.resname in prs.basic_salt_aa:
            select = 'name NH* NZ* NE* ND*'
        elif self.resname in prs.acidic_salt_aa:
            select = 'name OE* OD*'
        else: 
            return None

        atoms = self.atoms.select_atoms(select)
        # use here prs.rplcAtom it doesn't matter which atom name OE1 OE2 use
        # atom_names = [self.res+':'+prs.rplcAtom(name) for name in atoms.names]
        atom_names = [self.res+':'+name for name in atoms.names] 
        atom_coords = atoms.positions
        return atom_names, atom_coords


# TODO: think about output format
def find_pipi_bonds(res1, res2):
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
    if res1.pipi is None or res2.pipi is None:
        return None
    res1_name, centroid1, norm1 = res1.pipi # ['TRP:28:A:NORMAL'], origin and normal coords
    res2_name, centroid2, norm2 = res2.pipi
    # acid1 = ':'.join(res1_name[0].split(':')[0:3])
    # acid2 = ':'.join(res2_name[0].split(':')[0:3])

    out = []
    for i in range(len(res1_name)):
        for j in range(len(res2_name)):
            pipiangle = np.degrees(MDAmath.angle(norm1[i][0], norm2[j][0]))
            dist = MDdist3D(centroid1[i], centroid2[j])
            # condition for angle and distance b/w ring planes
            if pipiangle > prs.RIGHT_ANGLE:
                pipiangle = prs.STRAIGHT_ANGLE - pipiangle
            if ((dist <= prs.PIPI_D1 and pipiangle < prs.PIPI_ANG1) or 
                (dist <=prs.PIPI_D2 and prs.PIPI_ANG2 < pipiangle < prs.PIPI_ANG3)):
                # out.append(acid1+'\t'+acid2+'\tSCSC'+prs.PIPI_CONST+'PIPI\t'+res1_name[i].split(':')[3]+'\t'+res2_name[j].split(':')[3]+'\n')
                out.append(res1_name[i]+'\t'+res2_name[j]+'\tSCSC'+prs.PIPI_CONST+'PIPI\t'+'\n')
    return out

def find_pication_bonds(res1, res2):
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
    if res1.pipi is not None and res2.pication is not None:
        res_name, centroid, norm = res1.pipi        
        cation, pos = res2.pication
    elif res1.pication is not None and res2.pipi is not None:
        res_name, centroid, norm = res2.pipi
        cation, pos = res1.pication
    else: 
        return None

    out = []
    for i in range(len(res_name)):
        # acid = ':'.join(res_name[i].split(':')[0:3])
        catvec = pos - centroid[i]
        catangle = np.degrees(MDAmath.angle(norm[i][0], catvec))

        if catangle > prs.RIGHT_ANGLE: 
            catangle = prs.STRAIGHT_ANGLE - catangle
        if MDdist3D(centroid[i], pos) < prs.PICAT_DIST and catangle < prs.PICAT_ANG:
            # bond = acid+'\t'+':'.join(cation.split(':')[:3])+'\tSCSC'+prs.PICATION_CONST+'PICAT\t'+res_name[i].split(':')[3]+'\t'+cation.split(':')[3]+'\n'
            # out.append(bond)
            out.append(res_name[i]+'\t'+cation+'\tSCSC'+prs.PICATION_CONST+'PICAT\t'+'\n')
    return out

def find_salt_bridges(res1, res2):
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
    # global saltBridges
    # saltBridges = []
    # check not None residues and call coords, atoms involved in salt bridges
    if res1.sb is None or res2.sb is None:
        return None
    else:
        res1_name, coords1 = res1.sb
        res2_name, coords2 = res2.sb
        # acid1 = res1.res
        # acid2 = res2.res

    out = []
    # check basic and acidic residues
    name1, name2 = res1.resname, res2.resname
    if all(name in prs.basic_salt_aa for name in [name1, name2]) or all(name in prs.acidic_salt_aa for name in [name1, name2]):
        return None
    if np.abs(res1.acid.resid - res2.acid.resid) == prs.seq_dist_cutoff:
        return None

    # pairs, distances = mda.lib.distances.capped_distance(coords1, coords2, max_cutoff=prs.SALT_DISTANCE, min_cutoff=0)

    for i in range(len(res1_name)):
        for j in range(len(res2_name)):
            dist = MDdist3D(coords1[i], coords2[j])
            if dist < prs.SALT_DISTANCE:
                saltBridges.append(res1_name[i]+'\t'+res2_name[j])
                out.append(res1_name[i]+'\t'+res2_name[j]+'\tSCSC'+prs.SALT_CONST+'SB\t'+'\n')
    return list(set(out))

    # if distances.size != 0:
    #   for k, [i, j] in enumerate(pairs):
    #       # use prs.rplsAtom, replace NH1 and NH2 to NH1/2, 
    #       # a1, a2 = atoms1[i].name, atoms2[j].name
    #       # a1, a2 = atoms1[i].split(':')[3], atoms2[j].split(':')[3]
    #       a1, a2 = res1_names[i], res2_names[j]
    #       # atom1 = acid1+':'+prs.rplcAtom(a1)
    #       # atom2 = acid2+':'+prs.rplcAtom(a2)
    #       # TODO: why to write both? > use saltBridges in hydrogen_bonds, so idk how bonds written in hb  
    #       saltBridges.append(a1+'\t'+a2)
    #       saltBridges.append(a2+'\t'+a1)
    #       # add to out list non repeating residues with same ids
    #       # ARG-1034A ASP-1073A must appear once
    #       # out.append(atom1+'\t'+atom2+prs.SALT_CONST+prs.rplcAtom(atoms1[i].name)+'\t'+prs.rplcAtom(atoms2[j].name)+'\t'+'\n')
    #       out.append(a1+'\t'+a2+'\t'+'\tSCSC'+prs.SALT_CONST+'SB\t'+'\n')
    #       # out.append(''.join(atom1.split(':')[0:3])+'\t'+''.join(atom2.split(':')[0:3])+
    #                  # '\tSCSC'+prs.SALT_CONST+'SB\t'+atom1.split(':')[3]+'\t'+atom2.split(':')[3]+'\n')
    #       # print(out)
    #   return out[0]
def find_disulfide_bonds(res1, res2):
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
        acid1_name, coords1 = res1.disulf
        acid2_name, coords2 = res2.disulf
        # coord1, coord2 = res1.disulf, res2.disulf
        
    # acid1 = ''.join(res1.res.split(':'))
    # acid2 = ''.join(res2.res.split(':'))
    dist = MDdist3D(coords1, coords2)
    if dist < prs.DISULF_D:
        bond = acid1_name+'\t'+acid2_name+'\tSCSC'+prs.DISULF_CONST+'SS\t'+'\n'
        # bond = ':'.join(acid1.split(':')[:3])+'\t'+':'.join(acid2.split(':')[:3])+'\tSCSC'+prs.DISULF_CONST+'SS\tSG\tSG\n'
        return bond

def find_hydrogen_bonds(PDB, u, segids, coords, chain, saltBridges):
    '''
    @description
        Generate table of hydrogen bonds (MDAnalysis.HydrogenBondAnalysis) with parameters: 
            selection1-['protein'], selection2-['protein'] - Find hydrogen bonds among protein amino acids
            distance_type='hydrogen' - Measures hydrogen bond lenghs b/w donor hydrogen and acceptor atom
            distance - Distance cutoff b/w donor hydrogen and acceptor atoms <= 2.5; 
            angle - Angle cutoff, D-H-A (donor-hydrogen-acceptor) angle >= 120
        Table contains info about donor hydrogen and acceptor atoms
            | time | donor_index | acceptor_index | donor_resnm | donor_resid | 
            | donor_atom | acceptor_resnm | acceptor_resid | acceptor_atom | distance | angle |
        Traverse table; define hydrogen, donor, acceptor atoms and find geometric parameters of hydrogen bonds, distances and angles
        Constraints: 
            Not include bonds b/w donor and acceptor that are involved in salt bridges.
            Hydrogen and acceptor must be from non-adjacent residues
            Confine distances b/w  hydrogen-acceptor and donor-acceptor

    @input
        u - MDAnalysis Universe used to call HydrogenBondAnalysis
        segids - list of atoms' segids defined as u.atoms.segids
        coords - dict of atom names and their positions
        chain - dict of atom names and chain type (MC or SC)
    @return 
        file with hydrogen bonds in format
        acid1       | acid2         | Chain_type | Distance | Bond_type
        GLU:28:A:H  | HIS:26:A:ND1  | MCSC       | 15.52    | HB    
    '''
    h = hbonds(u, selection1='protein', selection2= 'protein', distance_type='hydrogen', distance=2.5, angle=120) 
    h.run()
    h.generate_table()
    print(saltBridges)
    tmp = []
    hydrogenBonds = []
    acceptorBonds = {}
    for hbond in h.table:
        hydrogen = hbond[3]+':'+str(hbond[4])+':'+segids[hbond[1]]+':'+hbond[5] # HIS:28:A:H
        acceptor = hbond[6]+':'+str(hbond[7])+':'+segids[hbond[2]]+':'+hbond[8] # GLU:46:A:N
        dH = hbond[3]+':'+hbond[5]  # HIS:H, select resname and atom name to check whether it is in hydrogen dict
        dA = hbond[6]+':'+hbond[8]  # GLU:N, acceptoor atom, ancedent
        if dH not in prs.hydrogen:
            continue        
        donor = hbond[3]+':'+str(hbond[4])+':'+segids[hbond[1]]+':'+prs.hydrogen[dH] # HIS:28:A:N
        if not prs.isDonor(donor):
            continue
        if (donor+'\t'+acceptor) in saltBridges or (acceptor+'\t'+donor) in saltBridges: # no bonds from salt bridges
            continue
        if np.abs(hbond[4]-hbond[7]) > prs.seq_dist_cutoff:
        # if hydrogen.split(':')[:3] != acceptor.split(':')[:3] and (np.abs(hbond[4]-hbond[7]) != prs.seq_dist_cutoff):
            if (dA in prs.nextAcid) and prs.isAcceptor(acceptor):
                alpha = hbond[10] # D-H_A angle 
                d1 = hbond[9] # hydrogen-acceptor distance, cutoff=3.5
                d2 = MDdist3D(coords[donor], coords[acceptor]) # donor-acceptor distance, cutoff=4.5
                if d1 > prs.H_A_dist_cutoff or d2 > prs.D_A_dist_cutoff: 
                    continue

                # Define acceptor ancedent (AA) atom and its coords
                aNext = prs.nextAcid[dA]
                if len(aNext)==2: 
                    neighbor1 = ":".join(acceptor.split(":")[0:3])+":"+aNext[0]
                    neighbor2 = ":".join(acceptor.split(":")[0:3])+":"+aNext[1]
                    AA = (coords[neighbor1]+coords[neighbor2])/2
                else:
                    AA = coords[":".join(acceptor.split(":")[0:3])+":"+aNext[0]]

                # Define beta angle b/w AA-acceptor-hydrogen
                a = AA - coords[acceptor]
                b = coords[hydrogen] - coords[acceptor]
                beta = np.degrees(mda.lib.mdamath.angle(a, b))
                beta = prs.STRAIGHT_ANGLE-beta if beta>prs.RIGHT_ANGLE else beta

                # Define gamma angle b/w AA-acceptor-donor
                a = AA - coords[acceptor]
                b = coords[donor] - coords[acceptor]
                gamma = np.degrees(mda.lib.mdamath.angle(a, b)) 

                if [(prs.hydrogen[dH]!='S' and hbond[8]!='S' and d1 < prs.hb_d1 and d2 < prs.hb_d2 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]=='S' and hbond[8]!='S' and d1 < prs.hb_d1 and d2 < prs.hb_d2 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]!='S' and hbond[8]=='S' and d1 < prs.hb_d11 and d2 < prs.hb_d21 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]=='S' and hbond[8]=='S' and d1 < prs.hb_d12 and d2 < prs.hb_d22 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE)]:
                    
                    if acceptor not in acceptorBonds.keys():
                        acceptorBonds[acceptor] = []
                    if hbond[8] in prs.acceptor_atoms1:
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 2:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:2]
                    if hbond[8] in prs.acceptor_atoms2 and hbond[6]!="HOH": 
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 1:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:1] 
                    if hbond[8] == "O" and hbond[6]=="HOH":
                        acceptorBonds[acceptor].append(d1)
                    # beta = prs.STRAIGHT_ANGLE-beta if beta>prs.RIGHT_ANGLE else beta

                    if prs.SPHyb(donor)=="SP3" and prs.SPHyb(acceptor)=="SP3":
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP3', sphyb2='SP3', alpha=alpha, beta=beta) 

                    elif prs.SPHyb(donor)=="SP3" and prs.SPHyb(acceptor)=="SP2":
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP3', sphyb2='SP2', alpha=alpha, beta=beta)

                    elif prs.SPHyb(donor)=="SP2" and prs.SPHyb(acceptor)=="SP3":
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP2', sphyb2='SP3', alpha=alpha, beta=beta)

                    elif prs.SPHyb(donor)=="SP2" and prs.SPHyb(acceptor)=="SP2":
                        normalVecDonor = prs.normalDonorVecToPlane(donor, coords)
                        if normalVecDonor is None:
                            continue
                        normalVecAcceptor = prs.normalAcceptorVecToPlane(acceptor, coords)
                        if normalVecAcceptor is None:
                            continue
                        psi = np.degrees(mda.lib.mdamath.angle(normalVecDonor,normalVecAcceptor))
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP2', sphyb2='SP2', alpha=alpha, beta=beta, psi=psi)

                    # hydrogenBonds.append(hydrogen+"\t"+acceptor+"\t"+chain[donor]+chain[acceptor]+"\t"+str(d1)+"\t"+str(d2)+"\t"+str(alpha)+"\t"+str(beta)+"\t"+str(gamma)+"\t"+donor+"\t"+str(E)+"\t"+donor.split(":")[3]+"\t"+acceptor.split(":")[3]+"\n")
                    hydrogenBonds.append(hydrogen+'\t'+acceptor+'\t'+donor+'\t'+chain[donor]+chain[acceptor]+'\t'+str(d1)+'\t'+str(E))

    out = open(PDB.replace('.pdb','_hb'), 'w')
    finalHydrogenBonds = []
    donorBonds = {}
    for interaction in hydrogenBonds:
        hydrogen, acceptor, donor, chain, d1, E  = interaction.split("\t")[0:6] 
        bond = hydrogen+'\t'+acceptor+'\t'+chain+'\t'+E+'\tHB\t\n'

        if donor not in donorBonds:
            donorBonds[donor] = 0
        if donor.split(':')[0]=="HOH" or acceptor.split(':')[0]=="HOH":
            out.write(bond)

        if len(acceptorBonds[acceptor])==2 and (d1==str(acceptorBonds[acceptor][0]) or d1==str(acceptorBonds[acceptor][1])):
            if donor.split(':')[3] not in (prs.donor_atoms1+prs.donor_atoms2) and donorBonds[donor]==0:
                out.write(bond)
                donorBonds[donor]+=1
            elif donor.split(':')[3] in prs.donor_atoms1 and donorBonds[donor] < 2:
                out.write(bond)
                donorBonds[dA]+=1
            elif donor.split(':')[3] in prs.donor_atoms2 and donorBonds[donor] < 3:
                out.write(bond)
                donorBonds[donor]+=1
        elif d1 == str(acceptorBonds[acceptor][0]):
            if donor.split(':')[3] not in (prs.donor_atoms1+prs.donor_atoms2) and donorBonds[donor]==0:
                out.write(bond)
                donorBonds[donor]+=1
            elif donor.split(':')[3] in prs.donor_atoms1 and donorBonds[donor] < 2:
                out.write(bond)
                donorBonds[donor]+=1
            elif donor.split(':')[3] in prs.donor_atoms2 and donorBonds[donor] < 3:
                out.write(bond)
                donorBonds[donor]+=1
            elif donor.split(':')[3] in ["0"]:
                out.write(bond)
    out.close() 

def find_VanderWaals_bonds(PDB, protein, allatoms, chain):
    '''
    @description
        This is a method which uses mainly mda.objects. It takes moroe time then 2nd method 
        Calcualate capped distance b/w protein atoms with distance cutoff (4.25), which was selected as sum of max radii among atoms
        Modify and traverse pairs. Compare sum radii of two atoms (Rm) and distance b/w them (R)
        For atoms with the difference (R-Rm < 0.5) and non-adjacent residues calculate Energy value and write into file
    @input
        protein - MDA object, prot = u.select_atoms('protein').atoms
        allatoms - list of atoms in format ['HIS:26:A:N']
        chain - dict of atoms and chain types
    @return
        file with vanderwaals bonds in format, which then can be modified not to include bonds with negative energies
        Atom1      | Atom2       | Chain_type | Energy
        HIS:26:A:N | GLU:58:A:CD | MCSC       | 3.470
    '''
    pairs, distances = mda.lib.distances.capped_distance(protein.positions, protein.positions, max_cutoff=4.25, min_cutoff=0)
    # create array with sorted first numbers and other atoms' index from pairs
    van_atom = np.split(pairs[:, 1], np.cumsum(np.unique(pairs[:, 0], return_counts=True)[1])[:-1])
    vdw1, vdw2 = {}, {}

    for i, k in enumerate(van_atom):
        atom1 = protein[i]
        acid1 = allatoms[i]
        for j in k:
            atom2 = protein[j]
            acid2 = allatoms[j]
            if acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2] in vdw1 or acid2+'\t'+acid1+'\t'+chain[acid2]+chain[acid1] in vdw1:
                continue
            
            if atom1.resid!=atom2.resid and atom1.name in prs.radii and atom2.name in prs.radii:
                rm = prs.radii[atom1.name]+prs.radii[atom2.name]
                r = MDdist3D(atom1.position, atom2.position)

                if r-rm < 0.5 and not (np.abs(atom1.resid-atom2.resid)==1 and chain[acid1]=='MC' and chain[acid2]=='MC'):
                    if not acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2] in vdw1:
                        vdw1[acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2]] = []

                    E = prs.energy_vdw(rm , r)
                    vdw1[acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2]].append(E)

                    # if (('C' in atom1.name and 'C' in atom2.name) or 
                    #   ('C' in atom1.name and atom2.name in ['NE2','OE1','ND2','OD1'] and atom2.resname in ['GLN', 'ASN']) or 
                    #   ('C' in atom2.name and atom1.name in ['NE2','OE1','ND2','OD1'] and atom1.resname in ['GLN', 'ASN'])):
                    #   if not acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2] in vdw2:
                    #       vdw2[acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2]] = []
                    #   vdw2[acid1+'\t'+acid2+'\t'+chain[acid1]+chain[acid2]].append(E)
                    

    with open(PDB.replace(".pdb","_VanderWaals"),'w') as out:
        for contact in vdw1:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw1[contact]))+'\n')

    # with open(PDB.replace('.pdb', '_VanderWaals_new'), 'w') as out:
    #   for contact in vdw1:
    #       if not (sum(vdw1[contact]) < 0):
    #           out.write(contact+'\t'+str(sum(vdw1[contact]))+'\n')

    # with open(PDB.replace(".pdb","_VanderWaals1"),'w') as out:
    #   for contact in vdw2:
    #       if not (sum(vdw2[contact])<0 and abs(int(contact.split(':')[1]) - int(contact.split(':')[5]))==1):
    #           out.write(contact+'\t'+str(sum(vdw2[contact]))+'\n')

     
def find_vdw_bonds(PDB, protein, resnames, resids, segids, names, coords, chain):
    '''
    @description
        This is a method which uses mainly mda.objects. It takes moroe time then 2nd method 
        Calcualate capped distance b/w protein atoms with distance cutoff (4.25), which was selected as sum of max radii among atoms
        Modify and traverse pairs. Compare sum radii of two atoms (Rm) and distance b/w them (R)
        For atoms with the difference (R-Rm < 0.5) and non-adjacent residues calculate Energy value and write into file
    @input
        protein - MDA object, prot = u.select_atoms('protein').atoms
        allatoms - list of atoms in format ['HIS:26:A:N']
        chain - dict of atoms and chain types
    @return
        file with vanderwaals bonds in format, which then can be modified not to include bonds with negative energies
        Atom1      | Atom2       | Chain_type | Energy
        HIS:26:A:N | GLU:58:A:CD | MCSC       | 3.470
    @input
        protein - MDAnalysis object, u.select_atoms('protein')
        pairs - pairs of interacted atoms created by capped_distance in main
        file1 - pdb file used for generating HB file
        resnames - created in main array of resnames as u.atoms.resnames
        resids - created in main array of resids as u.atoms.resids
        segids - created in main array of segids as u.atoms.segids
        names - created in main array of names as u.atoms.names
        coords - dictionary of atoms and positions 
        chain - dictionary of atoms and their type of chain (MC or SC)

    '''

    pairs, distances = mda.lib.distances.capped_distance(protein.positions, protein.positions, max_cutoff=4.25, min_cutoff=0)
    # create array with sorted first numbers and other atoms' index from pairs
    van_atom = np.split(pairs[:, 1], np.cumsum(np.unique(pairs[:, 0], return_counts=True)[1])[:-1])
    vdw1, vdw2 = {}, {}

    for i, k in enumerate(van_atom):
        elem1 = resnames[i]+':'+str(resids[i])+':'+segids[i]+':'+names[i]
        res1 = resnames[i]+':'+str(resids[i])+':'+segids[i]
        atom1 = names[i]
        for j in k:
            elem2 = resnames[j]+':'+str(resids[j])+':'+segids[j]+':'+names[j]
            res2 = resnames[j]+':'+str(resids[j])+':' +str(segids[j])
            atom2 = names[j]
            if res2+'\t'+res1+'\t'+chain[elem2]+chain[elem1]+'\t'+atom2+'\t'+atom1 in vdw1:
                continue
            if (not res1==res2 and atom1 in prs.radii and atom2 in prs.radii):
                rm = prs.radii[atom1] + prs.radii[atom2] 
                r = MDdist3D(coords[elem1], coords[elem2]) 

                if r - rm < .5 and not (np.abs(int(elem1.split(':')[1]) - int(elem2.split(':')[1]))==1 and chain[elem1]=="MC" and chain[elem2]=="MC"):

                    if not res1+'\t'+res2+'\t'+chain[elem1]+chain[elem2]+'\t'+atom1+'\t'+atom2 in vdw1:
                        vdw1[res1+'\t'+res2+'\t'+chain[elem1]+chain[elem2]+'\t'+atom1+'\t'+atom2] = []
                    E = prs.energy_vdw(rm , r) 
                    vdw1[res1+'\t'+res2+'\t'+chain[elem1]+chain[elem2]+'\t'+atom1+'\t'+atom2].append(E)
                    if (("C" in atom1 and "C" in atom2) or 
                        ("C" in atom1 and atom2 in ["NE2","OE1","ND2","OD1"] and res2.split(" ")[0] in ["GLN","ASN"]) or 
                        (atom1 in ["NE2","OE1","ND2","OD1"] and res1.split(" ")[0] in ["GLN","ASN"] and "C" in atom2)):
                    
                        if not res1+'\t'+res2+'\t'+chain[elem1]+chain[elem2]+'\t'+atom1+'\t'+atom2 in vdw2:
                            vdw2[res1+'\t'+res2+'\t'+chain[elem1]+chain[elem2]+'\t'+atom1+'\t'+atom2] = []
                        vdw2[res1+'\t'+res2+'\t'+chain[elem1]+chain[elem2]+'\t'+atom1+'\t'+atom2].append(E)
   
    with open(PDB.replace(".pdb","_vdw"),'w') as out:
        for contact in vdw1:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw1[contact]))+'\n')

    with open(PDB.replace(".pdb","_vdw"),'w') as out:
        for contact in vdw2:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw1[contact]))+'\n')

def find_metal_bonds(PDB, metalls, acids_class):
    '''
    @description
        For each metal find atoms from acids_class that are bonded to that metal
        The distance b/w metal and atom must be < 3.5 (metal_dist_cutoff) and correspond to distance cutoffs b/w acid atoms and metals
        
        Write to list (metal2atom) appropriate bonds b/w metall and atoms
        Traverse 2 atoms in metal2atom, if angle b/w them [90, 180] write these atoms which are connected through metal
    @input
        metalls - mda selection of metalls 
        acids_class - created list of class of aminoacids
    @return
        file in format 
        Acid1      | Acid2      | Chain_type | Energy | Bond_type 
        HIS:26:A:C | GLU:58:A:N | MCSC       | 3      | METAL
    '''
    metal2atom=[]
    metalBonds = []
    for metal in metalls:
        # metall = met.resname+':'+str(met.resid)+':'+met.segid+':'+met.name
        met_pos = metal.position
        for i in range(len(acids_class)):
            acid = acids_class[i]
            for atom in acid.atoms:
                dist = MDdist3D(met_pos, atom.position)
                if dist < prs.metal_dist_cutoff and acid.resname in prs.METALLS_DICT:
                    met_dist, atoms = prs.METALLS_DICT[acid.resname]
                    if metal.resname in met_dist:
                        d = met_dist[metal.resname]
                        if acid.resname is 'ASP' and atom.name in atoms:
                            if dist < prs.GLUASPmonoDist[metal.resname]:
                                metal2atom.append(atom)
                        elif acid.resname is 'GLU' and atom.name in atoms:
                            if dist < prs.GLUASPmonoDist[metal.resname]:
                                metal2atom.append(atom)
                        elif dist <= d and atom.name in atoms:
                            metal2atom.append(atom)

        for i in range(len(metal2atom)):
            for j in range(i+1, len(metal2atom)):
                atom1, atom2 = metal2atom[i], metal2atom[j]
                atom_1 = atom1.resname+':'+str(atom1.resid)+':'+atom1.segid+':'+atom1.name
                atom_2 = atom2.resname+':'+str(atom2.resid)+':'+atom2.segid+':'+atom2.name
                atom1_pos = atom1.position
                atom2_pos = atom2.position
                coords1 = met_pos - atom1_pos
                coords2 = met_pos - atom2_pos

                angle = np.degrees(MDAmath.angle(coords1, coords2))
                if prs.RIGHT_ANGLE < angle < prs.STRAIGHT_ANGLE:
                    metalBonds.append(atom_1+'\t'+atom_2+'\t'+chain[atom_1]+chain[atom_2]+'\t3\tMETAL\n')
                    # metalBonds.append(''.join(atom_1.split(':')[:3])+'\t'+''.join(atom_2.split(':')[:3])+'\t'+
                    #               chain[atom_1]+chain[atom_2]+'\t3\tMETAL\t'+''.join(atom_1.split(':')[3])+'\t'
                    #               +''.join(atom_2.split(':')[3])+'\n')
                                    #   metall+'\t'+str(angle))
    
    with open(PDB.replace('.pdb', '_metal'), 'a') as out:
        for bond in metalBonds:
            out.write(bond)                
            

def find_dna_bonds(u, allatoms_data):
    '''
    @description
        Select DNA resnames (DA, DG, DC, DT).
        Iterate through 2 atoms around DNA resname, if angle b/w these 2 atoms suits to conditions.
        Append to list DNAbindingPairs and write output file
    @input
        u - MDAnalysis Universe object
        allatoms_data - dict of atoms and their chain types and coords
    @return 
        file in format
        DNA_atom | Protein_atom | Chaintype | Value | BondType |  
    '''
    # dna = u.select_atoms('resname DA DG DC DT') # select dna acids
    # allatoms_data - also contains information about DNA atoms (DA, DG, DC, DT)
    dna = u.select_atoms('resname DA DG DC DT')
    atoms_around_dna = u.select_atoms('protein and around 7.5 resname DA DG DC DT') # only consider protein atoms
    around_atoms = [i+':'+str(j)+':'+k+':'+l for i, j, k, l in zip(atoms_around_dna.resnames, atoms_around_dna.resids, atoms_around_dna.segids, atoms_around_dna.names)]
    
    dna_atoms = [i+':'+str(j)+':'+k+':'+l for i, j, k, l in zip(dna.resnames, dna.resids, dna.segids, dna.names)]
    dna_data = {atom : pos for atom, pos in zip(dna_atoms, dna.positions)}
    DNAbindingPairs = []

    for dna_atom, dna_coords in dna_data.items():  # DA:9:D:N1
        for atom1 in around_atoms:
            bound = False
            if (("N" in dna_atom.split(':')[3] and (prs.isDonor(atom1) or prs.isAcceptor(atom1))) or 
                ("O" in dna_atom.split(':')[3] and (prs.isDonor(atom1) or prs.isAcceptor(atom1))) or 
                (("N" in dna_atom.split(':')[3] or "O" in dna_atom.split(':')[3]) and "S" in atom1.split(':')[3]) or 
                ("C" in dna_atom.split(':')[3] and "O" in atom1.split(':')[3] and prs.isAcceptor(atom1)) or 
                ("O" in dna_atom.split(':')[3] and atom1.split(':')[3] in ["NZ","NH1","NH2","OD1","OD2","OE1","OE2"]) or 
                ("C" in dna_atom.split(':')[3] and "C" in atom1.split(':')[3])):

                for atom2 in around_atoms:
                    if atom2 == atom1 or atom2.split(':')[3] == 'H': # check isH function
                        continue
                    vector1 = dna_coords - allatoms_data[atom2]['coords']
                    vector2 = allatoms_data[atom2]['coords'] - allatoms_data[atom1]['coords']
                    omega = np.degrees(mda.lib.mdamath.angle(vector1, vector2))
                    if omega <= prs.RIGHT_ANGLE:
                        bound = True
                    else:
                        bound = False
                        break
            if bound:
                DNAbindingPairs.append([dna_atom, atom1])

    with open('DNAbonds', 'w') as out:
        for dna_atom, prot_atom in DNAbindingPairs:
            out.write(dna_atom+'\t'+prot_atom+'\tMC'+allatoms_data[prot_atom]['chain']+'\t10\tDNA\tNT\n')


def find_ligand_atom_bonds_new(PDB, ligand_centroids, atoms_dict):
    '''
    @description
        Consider sidechained atom, for it create 'distances' dict with 'dist' and 'ligand' keys
        'dist' key contains information about distance to most closed ligand 
        'ligand' key contains informtaion about that ligand
        Search for min distnace b/w atom and ligands, and write it to dictionary
    @input
        allatoms - list of atoms' names, ex. 'HIS:26:A:H'
        ligand_centroids - dict of ligand names and centroid coords
        chain - dict of atoms and their chain types
        coords - dict of atoms and their coordinates
    @return 
        file in format
        Atom       | Ligand     | Distance
        ALA:0:A:CB | CLR:2403:A | 30.8
    '''
    ligand_distances = {}
    with open('Ligands_new', 'w') as out:
        for atom, atom_data in atoms_dict.items():
            if atom_data['chain'] == 'SC':
                ligand_distances[atom] = {}
                ligand_distances[atom]['dist'] = prs.ligand_dist #tmp distance - 150
                ligand_distances[atom]['ligand'] = '' #tmp ligand name
                for ligand, ligand_coords in ligand_centroids.items():
                    dist = MDdist3D(atom_data['coords'], ligand_coords)
                    if distance > ligand_distances[atom]['dist']:
                        continue
                    ligand_distances[atom]['dist'] = distance
                    ligand_distances[atom]['ligand'] = ligand
                    out.write(atom+'\t'+ligand+'\t'+str(ligand_distances[atom]['dist'])+'\n')

        # for atom in allatoms:
        # for atom, atom_coords in coords.items():
        #     if chain[atom] == 'SC':
        #         ligand_distances[atom] = {}
        #         ligand_distances[atom]['dist'] = prs.ligand_dist #tmp distance - 150
        #         ligand_distances[atom]['ligand'] = '' #tmp ligand name
        #         for ligand, ligand_coords in ligand_centroids.items():
        #             distance = MDdist3D(atom_coords, ligand_coords)
        #             if distance > ligand_distances[atom]['dist']:
        #                 continue
        #             ligand_distances[atom]['dist'] = distance
        #             ligand_distances[atom]['ligand'] = ligand
        #             out.write(atom+'\t'+ligand+'\t'+str(ligand_distances[atom]['dist'])+'\n')
    

def find_ligand_atom_bonds_old(allatoms, chain, coords): # old version 
    '''
    @description
        iterate through protein atoms and ligands, finding distance b/w them and write to file
        replacing distance to less value try to find the less distance b/w atom and ligand
    '''
    with open('Ligands_old','w') as out:
        for atom in allatoms:
            if chain[atom] == 'SC':
                tmpdist0 = prs.ligand_dist
                for ligand in sorted(ligand_centroids.keys()):
                    tmpdist1 = MDdist3D(coords[atom], ligand_centroids[ligand])
                    if tmpdist1 > tmpdist0:
                        continue
                    out.write(atom+'\t'+ligand+'\t'+str(tmpdist1)+'\n')
                    tmpdist0 = float(copy.copy(tmpdist1))

# TODO: not clear
# def Centroids():
    # with open(re.sub(".pdb","",file1)+"_centroidNetSC",'w') as out:
    #     for residue1 in sorted(centroid.keys()):
    #         for residue2 in sorted(centroid.keys()):
    #             if "HOH" in residue1 and "HOH" in residue2:
    #                 continue
    #             if residue1==residue2:
    #                 continue
    #             x = MDdist3D(np.array(centroid[residue1]),np.array(centroid[residue2])) 
    #             if x < 8.5:
    #                 out.write(re.sub("-","",residue1)+"\t"+re.sub("-","",residue2)+"\t"+"SCSC"+"\t"+"1"+"\t"+"CENTROID\t"+"CENTROID\t"+"CENTROID\t"+str(x)+"\n")
    # '''
    # find the shortest distance b/w centroid and ligand
    # '''
    # with open(re.sub(".pdb","",file1)+"_centroidNetLigand",'w') as out:
    #     for residue1 in centroid.keys():
    #         if "HOH" in residue1:
    #             continue
    #         tmpdist0 = 150
    #         for ligand in sorted(ligandCentroids.keys()):
    #             tmpdist1 = MDdist3D(np.array(centroid[residue1]),ligandCentroids[ligand])
    #             if tmpdist1 > tmpdist0:
    #                 continue
    #             out.write(re.sub("-","",residue1)+"\t"+ligand+"\t"+str(tmpdist1)+"\n")
    #             tmpdist0 = float(copy.copy(tmpdist1))   

def find_centroid_centroid_bonds(protein_centroids):
    '''
    @description
        In main define centroids and their coordinates
        Here find distances b/w them and confine by centroid_dist (8.5 A)
    @input
        centroid_coords - dict of residues names and centroid coords
    @return
        file in format
        Centroid1 | Centroid2 | Chain_type | Energy | Distance
        ALA:0:A   | GLY:5:A   | SCSC       | 1      | 10 
    '''
    with open(PDB.replace('.pdb', '_centroidNetSC'), 'w') as out:
        for centroid1, coords1 in protein_centroids.items():
            for centroid2, coords2 in protein_centroids.items():
                dist = MDdist3D(coords1, coords2)
                if 0 < dist < prs.centroid_dist:
                    out.write(centroid1+'\t'+centroid2+'\tSCSC\t1\t'+str(dist)+'\n')
       
def find_centroid_ligand_bonds(protein_centroids, ligand_centroids):
    '''
    @description
        For each centroid key create 'distances' dict with 'dist' and 'ligand' keys
        'dist' key contains information about distance to most closed ligand 
        'ligand' key contains informtaion about that ligand
        Search for min distnace b/w residue centroid and ligand centroid, and write it to dictionary
    @input
        centroid_coords - dict of residues names and their centroid coords
        ligand_centroids - dict of ligands and their centroid coords
    '''            
    ligand_distances = {}
    with open(PDB.replace('.pdb', '_centroidNetLigand'), 'w') as out:
        for centroid, centroid_coords in protein_centroids.items():
            ligand_distances[centroid] = {}
            ligand_distances[centroid]['dist'] = prs.ligand_dist
            ligand_distances[centroid]['ligand'] = ''
            for ligand, ligand_coords in ligand_centroids.items():
                dist = MDdist3D(ligand_coords, centroid_coords)
                if dist > ligand_distances[centroid]['dist']:
                    continue
                ligand_distances[centroid]['dist'] = dist
                ligand_distances[centroid]['ligand'] = ligand
                out.write(centroid+'\t'+ligand+'\t'+str(distances[centroid]['dist'])+'\n')
                # out.write(''.join(centroid.split(':'))+'\t'+''.join((distances[centroid]['ligand']).split(':'))+'\t'+str(distances[centroid]['dist'])+'\n')

def find_peptide_pairs(protein):
	'''
 	@description
		Find peptide pairs b/w two successive residues in aminoacid sequence

 	@input
		protein - MDAnalysis object, select_atoms('protein')
	    PDB file, ex. 4eiy.pdb
	''' 
	residues = protein.residues
	with open('polypeptidePairs', 'w') as out:
	    for i in range(protein.n_residues - 1):
	    	if residues[i].resid == residues[i+1].resid - 1:
	    		out.write(residues[i].resname+':'+str(residues[i].resid)+':'+residues[i].segid+'\t'+
	    			residues[i+1].resname+':'+str(residues[i+1].resid)+':'+residues[i+1].segid+'\tMCMC\t10\tPP\n')

def prepare_bfactor_file(protein):
    '''
    @description
    	Select "CA" atoms of residues and extract B-factor(tempfactor)
    	And write to file ("Bfactor") in format
		Acid     | B-factor
		HIS:26:A | 14.48
    @input 
        protein - MDAnalysis object, select_atoms('protein')
    '''
    CA_atoms = protein.select_atoms('name CA')    
    with open('Bfactor','w') as bfact:
        for atom in CA_atoms.atoms:
            bfact.write(atom.resname+':'+str(atom.resid)+':'+atom.segid+"\t"+str(atom.tempfactor)+"\n")

# TODO: not clear
def select_vdw_bonds(vdw_file): #1BTL_vdw file 
    '''
    @description
        Read vdw_file and process VanderWaals bonds.
        Compute total sum of energies for same connected amino acids,
        if this total energy is negative -> don't include interactions b/w these amino acids
    @input
        vdw_file - created file with VanderWaals bonds, in format
        Atom1      | Atom2       | Chain_type | Energy
        HIS:26:A:N | GLU:58:A:CD | MCSC       | 3.470
    @return 
        processed file of VanderWaals bonds with positive total energy, in format
        Atom1      | Atom2       | Chain_type | Energy
        HIS:26:A:N | GLU:58:A:CD | MCSC       | 3.470
    '''
    with open(vdw_file, 'r') as vdw_out:
        lines = vdw_out.readlines()
        bond_energies = {} # dict of bonds and summed total energy
        
        for line in lines:
            acid1, acid2, chain, energy = line.split('\t')
            res1 = ':'.join(acid1.split(':')[0:3])
            res2 = ':'.join(acid2.split(':')[0:3])
            bond1 = res1+':'+res2
            bond2 = res2+':'+res1
            if bond1 not in s and bond2 not in s:
                bond_energies[bond1] = float(energy)
            elif bond1 in s and bond2 not in s:
                bond_energies[bond1] += float(energy)
            elif bond1 not in s and bond2 in s:
                bond_energies[bond2] += float(energy)
        # Define bonds with negative total energy
        neg_vdw_bonds = dict((bond, energy) for bond, energy in bond_energies.items() if energy <= 0)

    with open('VdW_noRepulse', 'w') as out:
        for line in lines:
            if any(bond.split('\t')[0] in line and bond.split('\t')[1] in line for bond in list(neg_vdw_bonds.keys())): # don't write acids with negative energies
                continue
            else:
                out.write(line.replace('\n', '\t')+'VDW\n')

# def main(pdb_file, phi_file): 
# def main(u, *args, **args):
# 	# IF WE INPUT UNIVERSE WITH PDB TOPOLOGY AND FRAMES_DICT OF COORDINATES
# 	# pdb_file - 1BTL.pdb
# 	# phi_file - 1BTL.phi
# 	phifile = open(phi_file, 'r').readlines()
# 	phi_data = [line for line in phifile if 'ASG' in line] 

# 	# DEFINE mda objects for protein, water, metall, ligand molecules
# 	u = mda.Universe(pdb_file)
# 	protein = u.select_atoms('protein')
# 	hoh = u.select_atoms('resname HOH')
# 	metalls = u.select_atoms('resname {}'.format(' '.join(list(prs.metals))))
# 	ligands = u.select_atoms('not protein and not resname HOH') - metalls
# 	# ligands = u - protein - hoh - metalls

# 	prot = protein.atoms # define protein atoms
#     prot_resnames, prot_resids, prot_segids, prot_atoms = prot.resnames, prot.resids, prot.segids, prot.names
#     allatoms = [(i+':'+str(j)+':'+k+':'+l) for i, j, k, l in zip(prot_resnames, prot_resids, prot_segids, prot_atoms)] # format 'HIS:26:A:N'

#     allatoms_data = {} # define dict of protein atoms' coords and chain types
#     for atom, position in zip(allatoms, prot.positions):
#         allatoms_data[atom] = {}
#         allatoms_data[atom]['coords'] = position
#         allatoms_data[atom]['chain'] = ''
#         if atom.split(':')[3] in ["N","O","C","CA","HA2","HA3"]:
#             if atom.split(':')[0] == 'GLY' and atom.split(':')[3] in ["O","CA","HA2","HA3","N","NH"]:
#                 allatoms_data[atom]['chain'] = 'SC'
#             else:
#                 allatoms_data[atom]['chain'] = 'MC'
#         else:
#             allatoms_data[atom]['chain'] = 'SC'

#     acids_class = [AminoAcid(res) for res in prot.residues]

#     with open(PDB.replace('.pdb', '_bonds'), 'w') as net:
#         for i in range(len(acids_class)):
#             for j in range(i+1, len(acids_class)):
#                 bonds = find_pipi_bonds(acids_class[i], acids_class[j])
#                 if bonds: 
#                     for bond in bonds:
#                         net.write(bond)
#                 bonds = find_pication_bonds(acids_class[i], acids_class[j])
#                 if bonds: 
#                     for bond in bonds:
#                         net.write(bond)
#                 bonds = find_salt_bridges(acids_class[i], acids_class[j])
#                 if bonds: 
#                     net.write(bonds[0]) #why we write only one connection b/w acids
#                     # for bond in bonds:
#                         # net.write(bond)
#                 bonds = find_disulfide_bonds(acids_class[i], acids_class[j])
#                 if bonds: 
#                     for bond in bonds:
#                         net.write(bond)

if __name__ == '__main__':
    t0 = time.time()
    # PDB = sys.argv[1] # input file.pdb
    # PHI = PDB.replace('.pdb', '.phi')

    # TEST definition
    PDB = '../Test_data/1BTL.pdb'   #sys.argv[1]
    PHI = '../Test_data/1BTL.phi'   #sys.argv[2]

    if os.path.exists(PDB) and os.path.exists(PHI):
        phifile = open(PHI, 'r').readlines()
        phi_data = [line for line in phifile if 'ASG' in line] 
        print('PDB and PHI files are ready')
    else:
        print('Point appropriate PDB and PHI files')
        exit()
    
    prepare_secondary_structure_file(PDB, phi_data)

    prepare_rsa_file(PDB, phi_data)
    

    # ### MAIN ###
    # Define mda.objects and create array of AminoAcid classes
    u = mda.Universe(PDB)
    protein = u.select_atoms('protein') # u1 = u.select_atoms('protein and not name OXT')
    hoh = u.select_atoms('resname HOH')
    metalls = u.select_atoms('resname {}'.format(' '.join(list(prs.metals))))
    ligands = u.select_atoms('not protein and not resname HOH') - metalls

    prot = protein.atoms #define for each atom of protein
    prot_resnames, prot_resids, prot_segids, prot_atoms = prot.resnames, prot.resids, prot.segids, prot.names

    allatoms = [(i+':'+str(j)+':'+k+':'+l) for i, j, k, l in zip(prot_resnames, prot_resids, prot_segids, prot_atoms)] # HIS:26:A:N

    residues = [(res+':'+str(rid)+':'+sid) for res, rid, sid in zip(prot_resnames, prot_resids, prot_segids)] # HIS:26:A
    # coords = {res+':'+atom : pos for res, atom, pos in zip(residues, prot_atoms, prot.positions)} # {HIS:26:A:N : array([2.61,1.454,10.018]}
    coords = {atom : pos for atom, pos in zip(allatoms, prot.positions)}
    chain = {} #like {'HIS:26:A:N': 'MC'}
    for res, atom in zip(residues, prot_atoms):
        if atom in ["N","O","C","CA","HA2","HA3"]:
            if res[0:3] == 'GLY' and atom in ["O","CA","HA2","HA3","N","NH"]:
                chain[res+':'+atom] = 'SC'
            else:
                chain[res+':'+atom] = 'MC'
        else:
            chain[res+':'+atom] = 'SC'

    allatoms_data = {}
    for atom, position in zip(allatoms, prot.positions):
        allatoms_data[atom] = {}
        allatoms_data[atom]['coords'] = position
        allatoms_data[atom]['chain'] = ''
        if atom.split(':')[3] in ["N","O","C","CA","HA2","HA3"]:
            if atom.split(':')[0] == 'GLY' and atom.split(':')[3] in ["O","CA","HA2","HA3","N","NH"]:
                allatoms_data[atom]['chain'] = 'SC'
            else:
                allatoms_data[atom]['chain'] = 'MC'
        else:
            allatoms_data[atom]['chain'] = 'SC'

    # for residues of protein create AminoAcid class
    acids_class = [AminoAcid(res) for res in prot.residues]

    # prepare_bfactor_file(PDB, prot)
    saltBridges = []
    # Search interactions
    with open(PDB.replace('.pdb', '_bonds'), 'w') as net:
        for i in range(len(acids_class)):
            for j in range(i+1, len(acids_class)):
                bonds = find_pipi_bonds(acids_class[i], acids_class[j])
                if bonds: 
                    for bond in bonds:
                        net.write(bond)
                bonds = find_pication_bonds(acids_class[i], acids_class[j])
                if bonds: 
                    for bond in bonds:
                        net.write(bond)
                bonds = find_salt_bridges(acids_class[i], acids_class[j])
                if bonds: 
                    net.write(bonds[0]) #why we write only one connection b/w acids
                    # for bond in bonds:
                        # net.write(bond)
                bonds = find_disulfide_bonds(acids_class[i], acids_class[j])
                if bonds: 
                    for bond in bonds:
                        net.write(bond)
    print(len(saltBridges))
    # protein = u1.residues
    # prot_resnames, prot_resids, prot_segids, prot_names = protein.resnames, protein.resids, protein.segids, protein.names

    # h = hbonds(u, selection1='protein', selection2= 'protein', distance=2.5, distance_type='hydrogen', angle=120) 
    # h.run()
    # h.generate_table() 
    find_hydrogen_bonds(PDB, u, prot_segids, coords, chain, saltBridges)

    find_vdw_bonds(PDB, prot, prot_resnames, prot_resids, prot_segids, prot_atoms, coords, chain) # old version

    find_VanderWaals_bonds(PDB, prot, allatoms, chain) # new version
    # # VandWaals_awk_replacement(PDB.replace('.pdb','')+'_vdw')
    # # VandWaals_awk_replacement(PDB.replace('.pdb','')+'_vdw2')

    find_metal_bonds(PDB, metalls, acids_class) 

    # # # Concatenate created files into one 
    # # pdb = PDB.replace('.pdb', '')
    # # print(pdb)
    # # cmd = f'cat {pdb}_bonds {pdb}_hb {pdb}_vdw_noRepulse {pdb}_vdw_noRepulse2 {pdb}_metal > {pdb}_net'
    # # os.system(cmd) 

    # # Delete unneccessary files
    # # cmd = f'rm {pdb}_bonds {pdb}_hb {pdb}_vdw_noRepulse {pdb}_vdw_noRepulse2 {pdb}_metal {pdb}_vdw {pdb}_vdw2'
    # # os.system(cmd)

    # LIGANDS
    ligand_names = [i+':'+str(j)+':'+k for i, j, k in zip(ligands.residues.resnames, ligands.residues.resids, ligands.residues.segids)]
    ligand_centroids = dict(zip(ligand_names, ligands.center_of_geometry(compound='group')))

    # CENTROIDS
    # Calc weighted center of residues, where centroid_data.masses are weights
    centroid_data = u.select_atoms('protein and not resname DG DC DT DA and not backbone or (resname GLY and not name N C O)') 
    center_coords = centroid_data.center(centroid_data.masses, compound='residues')
    centroid = centroid_data.residues
    centroid_names = [i+':'+str(j)+':'+k for i, j, k in zip(centroid.resnames, centroid.resids, centroid.segids)]
    protein_centroids = dict(zip(centroid_names, center_coords))
    
    find_ligand_atom_bonds_new(PDB, ligand_centroids, allatoms_data)
    find_ligand_atom_bonds_old(allatoms, chain, coords)

    # find_centroid_bonds(centroid_coords)
    find_centroid_centroid_bonds(protein_centroids)
    find_centroid_ligand_bonds(protein_centroids, ligand_centroids)
    # find_centroid_ligand_bonds(PDB, centroid_coords, ligand_centroids)

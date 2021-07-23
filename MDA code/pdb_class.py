import MDAnalysis as mda
import numpy as np
import string
import sys, os, re, copy, math
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.lib.distances import calc_bonds as MDdist3D
import MDAnalysis.lib.mdamath as MDAmath
import time
import parameters as prs 

# INCLUDE POLYPEPTIDE PAIRS, add PP1, PP2 for proper work of energetics, centroid scriptss

# TODO: make preprocessing script + installation readme
'''
    Before need to preprocess pdb file
    1. Add hydrogens (phenix.reduce)
    2. Protonate water (phenix.ready_set)
    3. And create phi file (stride)
'''

def prepare_secondary_structure_file(PDB, phi_file):
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
    print("Prepare secondary structure file...")
    phifile = [line for line in open(phi_file, 'r').readlines() if 'ASG' in line]  
    secStructure = open(PDB.replace('.pdb', '_secondaryStructure'), 'w')

    counter = 0
    currentSecStructure = ""
    currentAngle = -1

    for i in range(len(phifile)):
        #here point i == 0 to maintain prevStruct and not to recall for each i
        if i == 0:  prevStruct = phifile[0][33:42].strip()
        if prevStruct == 'Coil' or prevStruct == 'Bridge':  prevStruct = 'CoilA'
        elif prevStruct == 'Turn':  prevStrcut = 'TurnA'

        line = phifile[i]   
        resname, resid, segid = line[5:8].strip(), line[11:15].strip(), line[9]
        phi_angle = float(line[42:51].strip())
        secStruct = line[33:42].strip()

        if resname == "PRO" and line[33:42].strip() == "aHelix":
            if i+1 != len(phifile) and phifile[i+1][33:42].strip() != "aHelix":
                secStruct = phifile[i+1][33:42].strip()
        if secStruct=="Bridge" or secStruct=="Coil":
            if ((phi_angle > prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or 
                (phi_angle < prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE) or 
                (phi_angle==prs.FULL_ANGLE and currentAngle < prs.ZERO_ANGLE)):
                if prevStruct == "CoilA":
                    secStruct="CoilA"; prevStruct = "CoilA"
                else:
                    secStruct="CoilB"; prevStruct = "CoilB"
            elif ((phi_angle < prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE) or
                (phi_angle > prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or
                (phi_angle==prs.FULL_ANGLE and currentAngle > prs.ZERO_ANGLE)):
                if prevStruct == "CoilA":
                    secStruct="CoilB"; prevStruct = "CoilB"
                else:
                    secStruct="CoilA"; prevStruct = "CoilA"
        if secStruct=="Turn":
            if ((phi_angle > prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or 
                (phi_angle < prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE) or 
                (phi_angle== prs.FULL_ANGLE and currentAngle < prs.ZERO_ANGLE)):
                if prevStruct == "TurnA":
                    secStruct="TurnA"; prevStruct = "TurnA"
                else:
                    secStruct="TurnB"; prevStruct = "TurnB"
            elif ((phi_angle < prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE) or 
                (phi_angle > prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE) or 
                (phi_angle== prs.FULL_ANGLE and currentAngle > prs.ZERO_ANGLE)):
                if prevStruct == "TurnA":
                    secStruct="TurnB"; prevStruct = "TurnB"
                else:
                    secStruct="TurnA"; prevStruct = "TurnA"

        if (("Coil" in secStruct or "Turn" in secStruct) and 
            (resname=="GLY" and phi_angle > prs.ZERO_ANGLE and phi_angle!=prs.FULL_ANGLE)):
            phiangle = line[42:51].strip()
            if float(phiangle)==prs.FULL_ANGLE:
                phiangle=prs.ZERO_ANGLE
            secStructure.write(resname+':'+resid+':'+segid+'\t'+secStruct+str(counter)+'\t'+str(phiangle)+'\n')
            currentSecStructure = secStruct                
            counter+=1
        elif secStruct != currentSecStructure:
            counter+=1
            phiangle = line[42:51].strip()
            if float(phiangle)==prs.FULL_ANGLE:
                phiangle=prs.ZERO_ANGLE
            secStructure.write(resname+':'+resid+':'+segid+'\t'+secStruct+str(counter)+'\t'+str(phiangle)+'\n')
            currentSecStructure = secStruct
        else:
            secStructure.write(resname+':'+resid+':'+segid+'\t'+secStruct+str(counter)+'\t'+line[42:51].strip()+'\n')
        currentAngle = float(line[42:51])
        if currentAngle==prs.FULL_ANGLE:
            curentAngle=-1
    secStructure.close()

def prepare_rsa_file(PDB, phi_file):
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
    print('Prepare RSA file...')
    phi_data = [line for line in open(phi_file, 'r').readlines() if 'ASG' in line]  
    with open(PDB.replace('.pdb','.rsa'), 'w') as rsafile:
        for line in phi_data:
            if line [5:8] in prs.area.keys():
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
        # self.res = aminoacid.resname+':'+str(aminoacid.resid)+':'+aminoacid.segid   #HIS:26:A
        self.res = aminoacid.resname+':'+str(aminoacid.resid)+':'+self.atoms.chainIDs[0]
        self.centroid, self.normal = self.calc_centroid_normal_coords()
        self.pipi = self.PiPi()
        self.pication = self.PiCation()
        self.disulf = self.Disulf()
        self.sb = self.SaltBridge()

    def calc_centroid_normal_coords(self):
        '''
        @description
            Consider amino acids involved in pi interactions (HIS, TYR, PHE, TRP).
            Select acid's ring atoms (atom_groups) and for ring plane calculate coords of center of geometry (centroid) and normal vecs to ring plane of acid
        '''
        centroid_dict = {} 
        normal_dict = {} 

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
            for i, name in enumerate([':HEX', ':PENT']):
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
            resname_key = list(self.centroid.keys())    # ex. ['HIS:26:A:NORMAL']
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
            if cations.names.size != 0: # some residues cannot contain NZ CZ atoms, returns zero array([])
                cation_names = self.res+':'+cations.names[0]
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
        atom_names = [self.res+':'+name for name in atoms.names] 
        atom_coords = atoms.positions
        return atom_names, atom_coords

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
        - list of found pipi bonds to file, [TRP:229:A:HEX  TRP:290:A:HEX   SCSC    10  PIPI]
    '''
    if res1.pipi is None or res2.pipi is None:
        return None
    res1_name, centroid1, norm1 = res1.pipi # ['TRP:28:A:NORMAL'], origin and normal coords
    res2_name, centroid2, norm2 = res2.pipi

    out = []
    for i in range(len(res1_name)):
        for j in range(len(res2_name)):
            pipiangle = np.degrees(MDAmath.angle(norm1[i][0], norm2[j][0]))
            dist = MDdist3D(centroid1[i], centroid2[j])
            # condition for angle and distance b/w ring planes
            if pipiangle > prs.RIGHT_ANGLE:
                pipiangle = prs.STRAIGHT_ANGLE - pipiangle
            if ((dist <= prs.PIPI_D1 and pipiangle < prs.PIPI_ANG1) or (dist <=prs.PIPI_D2 and prs.PIPI_ANG2 < pipiangle < prs.PIPI_ANG3)):
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
        - list of found pication bonds to file, [TRP:290:A:HEX  ARG:259:A:CZ    SCSC    9.4 PICAT]
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
        catvec = pos - centroid[i]
        catangle = np.degrees(MDAmath.angle(norm[i][0], catvec))

        if catangle > prs.RIGHT_ANGLE: 
            catangle = prs.STRAIGHT_ANGLE - catangle
        
        distance = MDdist3D(centroid[i], pos) # distance b/w centroid positions and cation pos
        if distance < prs.PICAT_DIST and catangle < prs.PICAT_ANG:
            out.append(res_name[i]+'\t'+cation+'\tSCSC'+prs.PICATION_CONST+'PICAT\t'+'\n')
    return out

saltBridges = []
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
        - list of all salt bridges, [LYS:34:A:NZ   ASP:38:A:OD2]
        - list of found salt bridges to file,  [LYS:32:A:NZ    ASP:35:A:OD2    SCSC    20  SB]
    '''
    # check not None residues and call coords, atoms involved in salt bridges
    if res1.sb is None or res2.sb is None:
        return None
    else:
        res1_name, coords1 = res1.sb
        res2_name, coords2 = res2.sb

    out = []
    # check basic and acidic residues
    name1, name2 = res1.resname, res2.resname
    if all(name in prs.basic_salt_aa for name in [name1, name2]) or all(name in prs.acidic_salt_aa for name in [name1, name2]):
        return None
    if np.abs(res1.acid.resid - res2.acid.resid) == prs.seq_dist_cutoff:
        return None

    for i in range(len(res1_name)):
        for j in range(len(res2_name)):
            dist = MDdist3D(coords1[i], coords2[j])
            if dist < prs.SALT_DISTANCE:
                saltBridges.append(res1_name[i]+'\t'+res2_name[j])
                out.append(res1_name[i]+'\t'+res2_name[j]+'\tSCSC'+prs.SALT_CONST+'SB\t'+'\n')
    return out

def find_disulfide_bonds(res1, res2):
    '''
    @description
        for atoms of CYS find disulfide bonds based on distance cutoff
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @return 
        - list of found disulfide bonds, [CYS:77:A:SG    CYS:123:A:SG    SCSC    167 SS]
    '''
    # check if residues have disulfide params and return coords of SG atoms
    if res1.disulf is None or res2.disulf is None:  
        return None
    else:   
        acid1_name, coords1 = res1.disulf
        acid2_name, coords2 = res2.disulf
        
    dist = MDdist3D(coords1, coords2)
    if dist < prs.DISULF_D:
        bond = acid1_name+'\t'+acid2_name+'\tSCSC'+prs.DISULF_CONST+'SS\t'+'\n'
    return bond

def find_hydrogen_bonds_new(PDB, u, segids, allatoms_data, saltBridges):
    '''
    https://docs.mdanalysis.org/dev/documentation_pages/analysis/hydrogenbonds.html
    @description
        Use new version of HydrogenBondsAnalysis module
        Generate table of hydrogen bonds (MDAnalysis.HydrogenBondAnalysis) with parameters: 
            universe=u
            donors_sel, acceptors_sel - select N* and O* atoms for donors and acceptors and find hydrogen bonds
            d_h_a_angle_cutoff[90] - Angle cutoff, D-H-A (donor-hydrogen-acceptor), set to 90 to include water bonds
            d_h_cutoff[1.2] - Distance cutoff b/w Donor and Hydrogen atoms
            d_a_cutoff[3.5] - Distance cutoff b/w Donor and Acceptor atoms
        Table contains info about donor hydrogen and acceptor atoms
            | frame | donor_id | hydrogen_id | acceptor_id | distance | angle | 
        We traverse table; define hydrogen, donor, acceptor atoms and find geometric parameters of hydrogen bonds, distances and angles
        Constraints: 
            Not include bonds b/w donor and acceptor that are involved in salt bridges.
            Hydrogen and acceptor must be from non-adjacent residues
            Confine distances b/w  hydrogen-acceptor and donor-acceptor
    @input
        PDB - used for write hydrogen_bonds file
        u - MDAnalysis Universe used to call HydrogenBondAnalysis
        segids - list of atoms' segids defined as u.atoms.segids
        allatoms_data - dict of atom name and their positions and chain types (MC or SC)
        saltBdridges - list of salt bridges, hydrogen bonds should not be salt bridges
    @return 
        file with hydrogen bonds in format
        acid1       | acid2         | Chain_type | Distance | Bond_type
        GLU:28:A:H  | HIS:26:A:ND1  | MCSC       | 15.52    | HB    
    '''
    hbonds = HBA(universe=u, 
                donors_sel='all and name N* O*',
                hydrogens_sel='all and name H*',
                acceptors_sel='all and name N* O*',
                d_h_a_angle_cutoff=90, 
                d_h_cutoff=1.2,
                d_a_cutoff=3.5,
                )
    hbonds.run()
    tmp = []
    hydrogenBonds = []
    acceptorBonds = {}
    for hbond in hbonds.hbonds:
        d, h, a, d_a_dist, d_h_a_angle = hbond[1:] #d_h_a_angle, return indexes of atoms from u.atoms
        d, h, a = u.atoms[int(d)], u.atoms[int(h)], u.atoms[int(a)]
        # Whether use segid obj or segids dict based on indexes
        donor = ':'.join([d.resname, str(d.resid), d.segid, d.name])
        hydrogen = ':'.join([h.resname, str(h.resid), h.segid, h.name])
        acceptor = ':'.join([a.resname, str(a.resid), a.segid, a.name])
        acc_res, acc_name = a.resname, a.name
        dA = ':'.join([a.resname, a.name])
        dH = ':'.join([h.resname, h.name])
        if dH not in prs.hydrogen:
            continue
        if not prs.isDonor(donor):
            continue
        if (donor+'\t'+acceptor) in saltBridges or (acceptor+'\t'+donor) in saltBridges: # no bonds from salt bridges
            continue

        if np.abs(h.resid - a.resid) > prs.seq_dist_cutoff:
            if (dA in prs.nextAcid) and prs.isAcceptor(acceptor):
                alpha = d_h_a_angle                
                # d1 = MDdist3D(h.position, a.position) 
                # d2 = d_a_dist
                # here was bug in calculation of d1 and d2, for 1btl SER-82A-HG and ALA-79A-O distances werenot same as for MDdist3D function, 
                # while for other bonds it was OK
                d1 = MDdist3D(allatoms_data[hydrogen]['coords'], allatoms_data[acceptor]['coords'])
                d2 = MDdist3D(allatoms_data[donor]['coords'], allatoms_data[acceptor]['coords'])
                if d1 > prs.H_A_dist_cutoff or d2 > prs.D_A_dist_cutoff: 
                    continue

                # Define acceptor ancedent (AA) atom and its coords
                aNext = prs.nextAcid[dA]
                if len(aNext)==2: 
                    neighbor1 = ":".join(acceptor.split(":")[0:3])+":"+aNext[0]
                    neighbor2 = ":".join(acceptor.split(":")[0:3])+":"+aNext[1]
                    AA_coords = (allatoms_data[neighbor1]['coords']+allatoms_data[neighbor2]['coords'])/2
                else:
                    neighbor1 = ":".join(acceptor.split(":")[0:3])+":"+aNext[0]
                    AA_coords = allatoms_data[neighbor1]['coords']
                
                # Define beta angle b/w AA-acceptor-hydrogen
                a = AA_coords - allatoms_data[acceptor]['coords']
                b = allatoms_data[hydrogen]['coords'] - allatoms_data[acceptor]['coords']
                beta = np.degrees(mda.lib.mdamath.angle(a, b))

                # Define gamma angle b/w AA-acceptor-donor
                a = AA_coords - allatoms_data[acceptor]['coords']
                b = allatoms_data[donor]['coords'] - allatoms_data[acceptor]['coords']
                gamma = np.degrees(mda.lib.mdamath.angle(a, b)) 

                if ((prs.hydrogen[dH]!='S' and acc_name!='S' and d1 < prs.hb_d1 and d2 < prs.hb_d2 and alpha > 90 and beta > 90 and gamma > 90) or
                    (prs.hydrogen[dH]=='S' and acc_name!='S' and d1 < prs.hb_d1 and d2 < prs.hb_d2 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]!='S' and acc_name=='S' and d1 < prs.hb_d11 and d2 < prs.hb_d21 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                    (prs.hydrogen[dH]=='S' and acc_name=='S' and d1 < prs.hb_d12 and d2 < prs.hb_d22 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE)):
                    
                    if acceptor not in acceptorBonds.keys():
                        acceptorBonds[acceptor] = []
                    if acceptor.split(':')[3] in prs.acceptor_atoms1:
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 2:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:2]
                    if acceptor.split(':')[3] in prs.acceptor_atoms2 and acceptor.split(':')[0]!="HOH": 
                        acceptorBonds[acceptor].append(d1)
                        if len(acceptorBonds[acceptor]) > 1:
                            acceptorBonds[acceptor].sort()
                            acceptorBonds[acceptor] = acceptorBonds[acceptor][0:1] 
                    if acceptor.split(':')[3] == "O" and acceptor.split(':')[0]=="HOH":
                        acceptorBonds[acceptor].append(d1)
                    # beta = prs.STRAIGHT_ANGLE-beta if beta<prs.RIGHT_ANGLE else beta

                    if prs.SPHyb(donor)=="SP3" and prs.SPHyb(acceptor)=="SP3":
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP3', sphyb2='SP3', alpha=alpha, beta=beta) 

                    elif prs.SPHyb(donor)=="SP3" and prs.SPHyb(acceptor)=="SP2":
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP3', sphyb2='SP2', alpha=alpha, beta=beta)

                    elif prs.SPHyb(donor)=="SP2" and prs.SPHyb(acceptor)=="SP3":
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP2', sphyb2='SP3', alpha=alpha, beta=beta)

                    elif prs.SPHyb(donor)=="SP2" and prs.SPHyb(acceptor)=="SP2":
                        normalVecDonor = prs.normalDonorVecToPlane(donor, allatoms_data)
                        if normalVecDonor is None:
                            continue
                        normalVecAcceptor = prs.normalAcceptorVecToPlane(acceptor, allatoms_data)
                        if normalVecAcceptor is None:
                            continue
                        psi = np.degrees(mda.lib.mdamath.angle(normalVecDonor,normalVecAcceptor))
                        E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP2', sphyb2='SP2', alpha=alpha, beta=beta, psi=psi)

                    hydrogenBonds.append(hydrogen+'\t'+acceptor+'\t'+donor+'\t'+allatoms_data[donor]['chain']+allatoms_data[acceptor]['chain']+'\t'+str(d1)+'\t'+str(E))

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
                donorBonds[donor]+=1
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

def find_VanderWaals_bonds(PDB, protein, allatoms, allatoms_data):
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
            if (acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain'] in vdw1 or 
                acid2+'\t'+acid1+'\t'+allatoms_data[acid2]['chain']+allatoms_data[acid1]['chain'] in vdw1):
                continue
            
            if atom1.resid!=atom2.resid and atom1.name in prs.radii and atom2.name in prs.radii:
                rm = prs.radii[atom1.name]+prs.radii[atom2.name]
                r = MDdist3D(atom1.position, atom2.position)

                if r-rm < 0.5 and not (np.abs(atom1.resid-atom2.resid)==1 and allatoms_data[acid1]['chain']=='MC' and allatoms_data[acid2]['chain']=='MC'):
                    if not acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain'] in vdw1:
                        vdw1[acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain']] = []

                    E = prs.energy_vdw(rm , r)
                    vdw1[acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain']].append(E)

                    if (('C' in atom1.name and 'C' in atom2.name) or 
                        ('C' in atom1.name and atom2.name in ['NE2','OE1','ND2','OD1'] and atom2.resname in ['GLN', 'ASN']) or 
                        ('C' in atom2.name and atom1.name in ['NE2','OE1','ND2','OD1'] and atom1.resname in ['GLN', 'ASN'])):
                        if not acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain'] in vdw2:
                            vdw2[acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain']] = []
                        vdw2[acid1+'\t'+acid2+'\t'+allatoms_data[acid1]['chain']+allatoms_data[acid2]['chain']].append(E)
                    

    with open(PDB.replace(".pdb","_VanderWaals"),'w') as out:
        for contact in vdw1:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw1[contact]))+'\tVDW\n')

    # with open(PDB.replace('.pdb', '_VanderWaals_new'), 'w') as out:
    #   for contact in vdw1:
    #       if not (sum(vdw1[contact]) < 0):
    #           out.write(contact+'\t'+str(sum(vdw1[contact]))+'\n')

    with open(PDB.replace(".pdb","_VanderWaals2"),'w') as out:
        for contact in vdw2:
            # if not (sum(vdw2[contact])<0 and abs(int(contact.split(':')[1]) - int(contact.split(':')[5]))==1):
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw2[contact]))+'\tVDW2\n')

     
def find_vdw_bonds(PDB, protein, allatoms, allatoms_data):
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
        elem1 = allatoms[i]
        res1 = ':'.join(elem1.split(':')[:3])
        atom1 = elem1.split(':')[3]
        for j in k:
            elem2 = allatoms[j]
            res2 = ':'.join(elem2.split(':')[:3])
            atom2 = elem2.split(':')[3]
            if (elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain'] in vdw1 or
                elem2+'\t'+elem1+'\t'+allatoms_data[elem2]['chain']+allatoms_data[elem1]['chain'] in vdw1):
                continue
            if (not res1==res2 and atom1 in prs.radii and atom2 in prs.radii):
                rm = prs.radii[atom1] + prs.radii[atom2] 
                r = MDdist3D(allatoms_data[elem1]['coords'], allatoms_data[elem2]['coords'])

                if r - rm < .5 and not (np.abs(int(elem1.split(':')[1]) - int(elem2.split(':')[1]))==1 and allatoms_data[elem1]['chain']=="MC" and allatoms_data[elem2]['chain']=="MC"):
                    # if not res1+'\t'+res2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain']+'\t'+atom1+'\t'+atom2 in vdw1:
                        # vdw1[res1+'\t'+res2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain']+'\t'+atom1+'\t'+atom2] = []
                    if not elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain'] in vdw1:
                        vdw1[elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain']] = []

                    E = prs.energy_vdw(rm , r) 
                    vdw1[elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain']].append(E)

                    if (("C" in atom1 and "C" in atom2) or 
                        ("C" in atom1 and atom2 in ["NE2","OE1","ND2","OD1"] and res2.split(":")[0] in ["GLN","ASN"]) or 
                        ("C" in atom2 and atom1 in ["NE2","OE1","ND2","OD1"] and res1.split(":")[0] in ["GLN","ASN"])):
                    
                        if not elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain'] in vdw2:
                            vdw2[elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain']] = []

                        vdw2[elem1+'\t'+elem2+'\t'+allatoms_data[elem1]['chain']+allatoms_data[elem2]['chain']].append(E)
   
    with open(PDB.replace(".pdb","_vdw"),'w') as out:
        for contact in vdw1:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw1[contact]))+'\tVDW\n')

    with open(PDB.replace(".pdb","_vdw2"),'w') as out:
        for contact in vdw2:
            if not (sum(vdw1[contact])<0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1]))==1):
                out.write(contact+'\t'+str(sum(vdw1[contact]))+'\tVDW2\n')

def find_metal_bonds(metalls, acids_class, allatoms_data, PDB):
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
                atom_1 = atom1.resname+':'+str(atom1.resid)+':'+atom1.chainID+':'+atom1.name
                atom_2 = atom2.resname+':'+str(atom2.resid)+':'+atom2.chainID+':'+atom2.name
                atom1_pos = atom1.position
                atom2_pos = atom2.position
                coords1 = met_pos - atom1_pos
                coords2 = met_pos - atom2_pos

                angle = np.degrees(MDAmath.angle(coords1, coords2))
                if prs.RIGHT_ANGLE < angle < prs.STRAIGHT_ANGLE:
                    # metalBonds.append(atom_1+'\t'+atom_2+'\t'+chain[atom_1]+chain[atom_2]+'\t3\tMETAL\n')
                    metalBonds.append(atom_1+'\t'+atom_2+'\t'+allatoms_data[atom_1]['chain']+allatoms_data[atom_2]['chain']+'\t3\tMETAL\n')
                    # metalBonds.append(''.join(atom_1.split(':')[:3])+'\t'+''.join(atom_2.split(':')[:3])+'\t'+
                    #               chain[atom_1]+chain[atom_2]+'\t3\tMETAL\t'+''.join(atom_1.split(':')[3])+'\t'
                    #               +''.join(atom_2.split(':')[3])+'\n')
                                    #   metall+'\t'+str(angle))
    
    with open(PDB.replace('.pdb', '_metal'), 'a') as out:
        for bond in list(set(metalBonds)):
            out.write(bond)                
            

def find_dna_bonds(u, allatoms_data, PDB):
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

    with open(PDB.replace('.pdb', '_DNAbonds'), 'w') as out:
        for dna_atom, prot_atom in DNAbindingPairs:
            out.write(dna_atom+'\t'+prot_atom+'\tMC'+allatoms_data[prot_atom]['chain']+'\t10\tDNA\tNT\n')


def find_ligand_atom_bonds_new(ligand_centroids, allatoms_data, PDB):
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
    with open(PDB.replace('.pdb', '_ligands'), 'w') as out:
        for atom, atom_data in allatoms_data.items():
            if atom_data['chain'] == 'SC':
                ligand_distances[atom] = {}
                ligand_distances[atom]['dist'] = prs.ligand_dist #tmp distance - 150
                ligand_distances[atom]['ligand'] = '' #tmp ligand name
                for ligand, ligand_coords in ligand_centroids.items():
                    dist = MDdist3D(atom_data['coords'], ligand_coords)
                    if dist > ligand_distances[atom]['dist']:
                        continue
                    ligand_distances[atom]['dist'] = dist
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
    

def find_ligand_atom_bonds_old(allatoms, allatoms_data, ligand_centroids, PDB): # old version 
    '''
    @description
        iterate through protein atoms and ligands, finding distance b/w them and write to file
        replacing distance to less value try to find the less distance b/w atom and ligand
    '''
    with open(PDB.replace('.pdb', '_ligands_old'),'w') as out:
        for atom in allatoms:
            # if chain[atom] == 'SC':
            if allatoms_data[atom]['chain'] == 'SC':
                tmpdist0 = prs.ligand_dist
                for ligand in sorted(ligand_centroids.keys()):
                    tmpdist1 = MDdist3D(allatoms_data[atom]['coords'], ligand_centroids[ligand])
                    # tmpdist1 = MDdist3D(coords[atom], ligand_centroids[ligand])
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

def find_centroid_centroid_bonds(protein_centroids, PDB):
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
       
def find_centroid_ligand_bonds(protein_centroids, ligand_centroids, PDB):
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
                out.write(centroid+'\t'+ligand+'\t'+str(ligand_distances[centroid]['dist'])+'\n')
                # out.write(''.join(centroid.split(':'))+'\t'+''.join((distances[centroid]['ligand']).split(':'))+'\t'+str(distances[centroid]['dist'])+'\n')

def find_peptide_pairs(protein, PDB):
    '''
    @description
        Find peptide pairs b/w two successive residues in aminoacid sequence

    @input
        protein - MDAnalysis object, select_atoms('protein')
        PDB file, ex. 4eiy.pdb
    ''' 
    residues = protein.residues
    with open(PDB.replace('.pdb', '_polypeptidePairs'), 'w') as out:
        for i in range(protein.n_residues - 1):
            segid1 = residues[i].atoms.chainIDs[0]
            segid2 = residues[i+1].atoms.chainIDs[0]
            # when separate trajectory to pdb files, in new created pdbs Residue objects dont contain chainIDs, so select from atoms
            if residues[i].resid == residues[i+1].resid - 1:
                out.write(residues[i].resname+':'+str(residues[i].resid)+':'+segid1+':PP1'+'\t'+
                    residues[i+1].resname+':'+str(residues[i+1].resid)+':'+segid2+':PP2'+'\tMCMC\t10\tPP\n')

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
            acid1, acid2, chain, energy, _ = line.split('\t')
            res1 = ':'.join(acid1.split(':')[0:3])
            res2 = ':'.join(acid2.split(':')[0:3])
            bond1 = res1+'\t'+res2
            bond2 = res2+'\t'+res1
            if bond1 not in bond_energies and bond2 not in bond_energies:
                bond_energies[bond1] = float(energy)
            elif bond1 in bond_energies and bond2 not in bond_energies:
                bond_energies[bond1] += float(energy)
            elif bond1 not in bond_energies and bond2 in bond_energies:
                bond_energies[bond2] += float(energy)
        # Define bonds with negative total energy
        neg_vdw_bonds = dict((bond, energy) for bond, energy in bond_energies.items() if energy <= 0)

    with open(vdw_file+'_noRepulse', 'w') as out:
        for line in lines:
            if any(bond.split('\t')[0] in line and bond.split('\t')[1] in line for bond in list(neg_vdw_bonds.keys())): # don't write acids with negative energies
                continue
            else:
                out.write(line.replace('\n', '\t')+'\n')

def find_waterbonds(PDB, hb_file, allatoms_data):
    hb_file = open(hb_file, 'r').readlines()
    HOH = {}
    HOH2 = {}
    HOH_energy = {}

    for line in hb_file:
        line = line.split('\t')
        if 'HOH' in line[0] and 'HOH' not in line[1]:
            water = ':'.join(line[0].split(':')[:3])
            if water not in HOH:
                HOH[water] = []
            if (line[1], water) not in HOH_energy:
                HOH_energy[line[1], water] = []
            HOH[water].append(line[1])
            HOH[water] = list(set(HOH[water]))
            HOH_energy[line[1], water].append(float(line[3]))
            HOH_energy[line[1], water] = list(set(HOH_energy[line[1], water]))

        if 'HOH' in line[1] and 'HOH' not in line[0]:
            water = ':'.join(line[1].split(':')[:3])
            if not water in HOH:
                HOH[water] = []
            if (line[0], water) not in HOH_energy:
                HOH_energy[line[0], water] = []
            HOH[water].append(line[0])
            HOH[water] = list(set(HOH[water]))
            HOH_energy[line[0], water].append(float(line[3]))
            HOH_energy[line[0], water] = list(set(HOH_energy[line[0], water]))

        if 'HOH' in line[0] and 'HOH' in line[1]:
            water1 = ':'.join(line[0].split(':')[:3])
            water2 = ':'.join(line[1].split(':')[:3])
            if water1 not in HOH2:
                HOH2[water1] = []
            if water2 not in HOH2:
                HOH2[water2] = []
            # if (water1, water2) not in HOH2_energy:
                # HOH2_energy[water1, water2] = []

            HOH2[water1].append(water2)
            HOH2[water2].append(water1)
            HOH2[water1] = list(set(HOH2[water1]))
            HOH2[water2] = list(set(HOH2[water2]))

    import itertools
    with open(PDB.replace('.pdb', 'waterbonds'), 'w') as out:
        for water in HOH:
            if len(HOH[water]) > 1:
                x = list(itertools.combinations(HOH[water], 2))
                for item in x:
                    if item[0].split(':')[:3] != item[1].split(':')[:3]:
                        energy = sum([*HOH_energy[item[0], water], *HOH_energy[item[1], water]])/2
                        out.write(item[0]+'\t'+item[1]+'\t'+allatoms_data[item[0]]['chain']+allatoms_data[item[1]]['chain']+'\t'+
                            str(energy) + '\t'+water+'\n')
                    
    with open(PDB.replace('.pdb', 'waterbonds'), 'a') as out:
        for water1 in HOH:
            for atom1 in HOH[water1]:
                if water1 in HOH2:
                    for water2 in HOH2[water1]:
                        if water2 in HOH:
                            for atom2 in HOH[water2]:
                                if atom1 != atom2 and atom1.split(':')[:3] != atom2.split(':')[:3]:
                                    energy = sum([*HOH_energy[atom1, water1], *HOH_energy[atom2, water2]])/2
                                    out.write(atom1+'\t'+atom2+'\t'+allatoms_data[atom1]['chain']+
                                              allatoms_data[atom2]['chain']+'\t'+str(energy)+'\t'+water1+'\t'+water2+'\n')


def main(PDB, u):
    '''
    @description
        Define main variables: chain and coordinates of protein atoms; water, metall, ligand AtomGroups
        Write output files with found pipi, pication, hydrogen, VanderWaals, ligand bonds
    @input
        PDB - pdb file name, in case of trajectory file it would be pdb file name with number of trajectory frame
              ex. 1BTL.pdb or 1BTL_3.pdb
        u - MDAnalysis Universe with topology file and (trajectory frame's coordinates) 
    '''
    PHI = PDB.replace('pdb', 'phi')
    prepare_secondary_structure_file(PDB, phi_file=PHI)
    prepare_rsa_file(PDB, phi_file=PHI)
    
    protein = u.select_atoms('protein and not name OXT')
    hoh = u.select_atoms('resname HOH or resname SOL')
    metalls = u.select_atoms('resname {}'.format(' '.join(list(prs.metals))))
    ligands = u.select_atoms('not protein and not resname HOH and not resname SOL') - metalls

    prot = u.atoms # protein.atoms 
    prot_resnames, prot_resids, prot_segids, prot_atoms = prot.resnames, prot.resids, prot.chainIDs, prot.names # chainIDs vs segids
    allatoms = [(i+':'+str(j)+':'+k+':'+l) for i, j, k, l in zip(prot_resnames, prot_resids, prot_segids, prot_atoms)] # format 'HIS:26:A:N'

    allatoms_data = {} # define dict of protein atoms' coords and chain types
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

    print('Create for each acid `class AminoAcid` object...')
    acids_class = [AminoAcid(res) for res in prot.residues]

    print('Find pi-pi; pi-cation; disulfide bonds and salt bridges...')
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
                    net.write(bonds[0]) 
                bonds = find_disulfide_bonds(acids_class[i], acids_class[j])
                if bonds: 
                    for bond in bonds:
                        net.write(bond)

    # find_hydrogen_bonds(PDB, u, prot_segids, allatoms_data, saltBridges)
    find_hydrogen_bonds_new(PDB, u, prot_segids, allatoms_data, saltBridges)
    # find_waterbonds(PDB, PDB.replace('.pdb', '_hb'), allatoms_data)
    # cmd = '''awk '{split($1,x,":");split($2,y,":");if(x[1]x[2]!=y[1]y[2]){print}}' "waterbonds" > "1btl_waterbonds"'''

    # find_vdw_bonds(PDB, prot, prot_resnames, prot_resids, prot_segids, prot_atoms, coords, chain) # old version
    # find_vdw_bonds(PDB, prot, allatoms, allatoms_data)
    find_vdw_bonds(PDB, protein.atoms, allatoms, allatoms_data)
    # find_VanderWaals_bonds(PDB, prot, allatoms, allatoms_data) # new version
    # select_vdw_bonds(PDB.replace('.pdb', '_VanderWaals'))
    # select_vdw_bonds(PDB.replace('.pdb', '_VanderWaals2'))
    select_vdw_bonds(PDB.replace('.pdb', '_vdw'))
    select_vdw_bonds(PDB.replace('.pdb', '_vdw2'))

    find_metal_bonds(metalls, acids_class, allatoms_data, PDB) 
    # find_dna_bonds(u, allatoms_data)
    find_peptide_pairs(protein.atoms, PDB)

    print('FIND HB and VDW')
    pdb = PDB.replace('.pdb', '')
    print(pdb)
    cmd = f"cat {pdb}_bonds waterbonds {pdb}_hb {pdb}_vdw_noRepulse {pdb}_vdw2_noRepulse {pdb}_metal {pdb}_polypeptidePairs > {pdb}_net"
    os.system(cmd)
    cmd = f"cat {pdb}_hb1 | sed /HOH/d >>{pdb}_net"
    os.system(cmd)

    ligand_names = [i+':'+str(j)+':'+k for i, j, k in zip(ligands.residues.resnames, ligands.residues.resids, ligands.residues.chainIDs[0])] # chainIDs vs segids
    ligand_centroids = dict(zip(ligand_names, ligands.center_of_geometry(compound='residues')))

    # CENTROIDS
    # Calc weighted center of residues, where centroid_data.masses are weights
    centroid_data = u.select_atoms('protein and not resname DG DC DT DA and not backbone or (resname GLY and not name N C O)') 
    center_coords = centroid_data.center(centroid_data.masses, compound='residues')
    centroid = centroid_data.residues
    centroid_names = [i+':'+str(j)+':'+k[0] for i, j, k in zip(centroid.resnames, centroid.resids, centroid.chainIDs)] # chainIDs vs segids. here chainIDs for residues -> take k[0]
    protein_centroids = dict(zip(centroid_names, center_coords))
    
    find_ligand_atom_bonds_new(ligand_centroids, allatoms_data, PDB)
    # find_ligand_atom_bonds_old(allatoms, allatoms_data, ligand_centroids, PDB)

    # # find_centroid_bonds(centroid_coords)
    find_centroid_centroid_bonds(protein_centroids, PDB)
    find_centroid_ligand_bonds(protein_centroids, ligand_centroids, PDB)
    # find_peptide_pairs(protein)
    # prepare_bfactor_file(protein)
    
    # Delete Files
    # cmd = f"rm -rf {pdb}_bonds {pdb}_hb {pdb}_vdw {pdb}_vdw2 {pdb}_vdw_noRepulse {pdb}_vdw2_noRepulse {pdb}_metal {pdb}_polypeptidePairs"
    # os.system(cmd)

    # Choose positive energized bonds
    cmd = f"awk '$4>0' {pdb}_net > {pdb}tmp"
    os.system(cmd)
    cmd = f"mv {pdb}tmp {pdb}_net"
    os.system(cmd)

if __name__ == '__main__':
    t0 = time.time()
    PDB = sys.argv[1]
    PHI = PDB.replace('pdb', 'phi')    
    # os.system(f"stride {PDB} > {PHI}") # stride must be in usr/local/bin

    if os.path.exists(PDB) and os.path.exists(PHI):
        print('PDB and PHI files are ready')
    else:
        print('Point appropriate PDB and PHI files')
        exit()

    u = mda.Universe(PDB)
    main(PDB, u)

    print(time.time()-t0)



    ### OUTPUT FILES
    # 1BTL_secondaryStructure2      [HIS:26:A    CoilA1  0]
    # 1BTL.rsa                      [HIS:26:A   0.474]
    # 1BTL_Bonds with pipi, picat, saltbridge, disulfide bonds  [LYS:32:A:NZ    ASP:35:A:OD2    SCSC    20  SB ]
    # 1BTL_hb                       [GLU:28:A:H HIS:26:A:ND1    MCSC    15.524611902380917  HB  ]
    # 1BTL_vdw, 1BTL_vdw2 - prev vdw bonds      [HIS:26:A:N GLU:58:A:OE1    MCSC    2.561345646191593   VDW]
    # 1BTL_VanderWaals, 1BTL_VanderWaals2 - current vdw bonds
    # VDW + *_noRepulse files - deleted negative bonds
    # 1BTL_metal                    [ASP:52:A:OD1    SER:91:A:OG  SCSC    NA:2402:A:NA 125.47308026698812]
    # 1BTL_ligand_new, 1BTL_ligand_old  [ALA:0:A:CB   CLR:2403:A   30.84465258857814]
    # 1BTL_centroidNetSC            [HIS:26:A   PRO:27:A    SCSC    1   4.734220978924919]
    # 1BTL_centroidNetLigand        [ALA:0:A  CLR:2403:A   30.928275454825783]
    ### MISS DNA FILE
    ### 1BTL.phi processed by stiride 1BTL.pdb

    ### OPTIONAL FILES
    # Bfactor               [HIS:26:A  14.479999542236328]
    # polypeptidePairs      [HIS:26:A    PRO:27:A    MCMC    10  PP]

import MDAnalysis as mda
import numpy as np
import itertools
import sys
import os
import warnings
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.lib.distances import calc_bonds as mdadist
import MDAnalysis.lib.mdamath as mdamath
import time
import parameters as prs

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def prepare_secondary_structure_file(pdb, phi_file):
    """
    @description
        phi_file created by stride (http://webclu.bio.wzw.tum.de/stride/)
        phi_file has secondary structure and phi angles info for amino acids (defined by 'ASG' lines)
        Traverse phi_file lines and define residue's name, id, phi_angle and secondary structure type
        Based on phi_angle values and sec_structure type write them into output file (*_secondaryStructure)
    @input
        phi_file - file with secondary structure info
    @output
        secondary structure file (ex. 1btl_secondaryStructure) in format ...
            Acid    |   Structure_type | phi_angle
            HIS26A  |   Coil           | 90
    """
    print("Prepare secondary structure file.")

    phifile = [line for line in open(phi_file, 'r').readlines() if 'ASG' in line]
    with open(pdb.replace('.pdb', '_secondaryStructure'), 'w') as secStructure:
        counter = 0
        currentSecStructure = ""
        currentAngle = -1
        prevStruct = ''

        for i in range(len(phifile)):
            # define prevStruct var as first structure, if it's Coil or Bridge replace to CoilA
            # then define secStruct and compare with prevStruct, if they are not same add +1 to counter
            if i == 0:
                prevStruct = phifile[0][33:42].strip()
            if prevStruct == 'Coil' or prevStruct == 'Bridge':
                prevStruct = 'CoilA'
            elif prevStruct == 'Turn':
                prevStruct = 'TurnA'

            line = phifile[i]
            resname, resid, segid = line[5:8].strip(), line[11:15].strip(), line[9]
            phi_angle = float(line[42:51].strip())
            secStruct = line[33:42].strip()

            if resname == "PRO" and line[33:42].strip() == "aHelix":
                if i + 1 != len(phifile) and phifile[i + 1][33:42].strip() != "aHelix":
                    secStruct = phifile[i + 1][33:42].strip()
            if secStruct == "Bridge" or secStruct == "Coil":
                if ((phi_angle > prs.ZERO_ANGLE > currentAngle and phi_angle != prs.FULL_ANGLE) or
                        (phi_angle < prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE) or
                        (phi_angle == prs.FULL_ANGLE and currentAngle < prs.ZERO_ANGLE)):
                    if prevStruct == "CoilA":
                        secStruct = "CoilA"
                        prevStruct = "CoilA"
                    else:
                        secStruct = "CoilB"
                        prevStruct = "CoilB"
                elif ((phi_angle < prs.ZERO_ANGLE < currentAngle) or
                      (phi_angle > prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE and phi_angle != prs.FULL_ANGLE) or
                      (phi_angle == prs.FULL_ANGLE and currentAngle > prs.ZERO_ANGLE)):
                    if prevStruct == "CoilA":
                        secStruct = "CoilB"
                        prevStruct = "CoilB"
                    else:
                        secStruct = "CoilA"
                        prevStruct = "CoilA"
            if secStruct == "Turn":
                if ((phi_angle > prs.ZERO_ANGLE > currentAngle and phi_angle != prs.FULL_ANGLE) or
                        (phi_angle < prs.ZERO_ANGLE and currentAngle < prs.ZERO_ANGLE) or
                        (phi_angle == prs.FULL_ANGLE and currentAngle < prs.ZERO_ANGLE)):
                    if prevStruct == "TurnA":
                        secStruct = "TurnA"
                        prevStruct = "TurnA"
                    else:
                        secStruct = "TurnB"
                        prevStruct = "TurnB"
                elif ((phi_angle < prs.ZERO_ANGLE < currentAngle) or
                      (phi_angle > prs.ZERO_ANGLE and currentAngle > prs.ZERO_ANGLE and phi_angle != prs.FULL_ANGLE) or
                      (phi_angle == prs.FULL_ANGLE and currentAngle > prs.ZERO_ANGLE)):
                    if prevStruct == "TurnA":
                        secStruct = "TurnB"
                        prevStruct = "TurnB"
                    else:
                        secStruct = "TurnA"
                        prevStruct = "TurnA"

            if ("Coil" in secStruct or "Turn" in secStruct) and (resname == "GLY" and phi_angle > prs.ZERO_ANGLE and phi_angle != prs.FULL_ANGLE):
                phiangle = float(line[42:51].strip())
                if phiangle == prs.FULL_ANGLE:
                    phiangle = prs.ZERO_ANGLE
                secStructure.write(resname + ':' + resid + ':' + segid + '\t' + secStruct + str(counter) + '\t' + str(phiangle) + '\n')
                currentSecStructure = secStruct
                counter += 1
            elif secStruct != currentSecStructure:
                counter += 1
                phiangle = float(line[42:51].strip())
                if phiangle == prs.FULL_ANGLE:
                    phiangle = prs.ZERO_ANGLE
                secStructure.write(resname + ':' + resid + ':' + segid + '\t' + secStruct + str(counter) + '\t' + str(phiangle) + '\n')
                currentSecStructure = secStruct
            else:
                secStructure.write(resname + ':' + resid + ':' + segid + '\t' + secStruct + str(counter) + '\t' + line[42:51].strip() + '\n')

            currentAngle = float(line[42:51])
            if currentAngle == prs.FULL_ANGLE:
                currentAngle = -1


def prepare_rsa_file(pdb, phi_file):
    """
    @description
        phi_file created by stride (http://webclu.bio.wzw.tum.de/stride/)
        phi_file has secondary structure and phi angles info for amino acids (defined by 'ASG' lines)
        Take solvent area values from phi file and divide by max solvent area of amino acid, taken
        from Rose 1985 (https://sci-hub.se/10.1126/science.4023714) and defined in parameters.py
    @input
        phi_file - file with secondary structure info
    @output
        rsa file (ex. 1btl_rsa) in format
            Acid    |   RSA_value
            HIS26A  |   0.47
    """
    print('Prepare RSA file.')
    
    phi_data = [line for line in open(phi_file, 'r').readlines() if 'ASG' in line]
    with open(pdb.replace('.pdb', '.rsa'), 'w') as rsa_file:
        for line in phi_data:
            if line[5:8] in prs.area.keys():
                # solvent area - line[64:69] in phi_file
                rsa_file.write(line[5:8] + ':' + line[11:15].strip() + ':' + line[9] + '\t' + str(float(line[64:69].strip()) / prs.area[line[5:8]]) + '\n')


class AminoAcid:
    """
    @description
        For each input amino acid as mda residue object define
        coordinates of atoms involved in interactions, normal vectors, centroids
        Whether amino acid's atoms take part in PiPi, PiCation, Disulfide or Salt bridge interactions
        return corresponding parameters or return None
    @input
        aminoacid - MDA residue object
    """
    def __init__(self, amino):
        self.acid = amino
        self.resname = amino.resname
        self.atoms = amino.atoms
        self.res = amino.resname + ':' + str(amino.resid) + ':' + self.atoms.chainIDs[0]
        self.centroid, self.normal = self.calc_centroid_normal_coords()
        self.pipi = self.PiPi()
        self.pication = self.PiCation()
        self.disulf = self.Disulfide()
        self.sb = self.SaltBridge()

    def calc_centroid_normal_coords(self):
        """
        @description
            Consider amino acids involved in PiPi interactions (HIS, TYR, PHE, TRP).
            For each acid select corresponding ring plane atoms and
            Calculate center of geometry (centroid.group) and normal vector coordinates orthogonal to ring plane
            and write into dictionaries, for HIS, TYR, PHE add 'NORMAL' to residue, for TRP - 'HEX', 'PENT'
        """
        if self.resname in prs.acids_set['HIS']:
            select = ['name CG ND1 CD2 CE1 NE2']
        elif self.resname in prs.acids_set['TYR'] + prs.acids_set['PHE']:
            select = ['name CG CD1 CD2 CE1 CE2 CZ']
        elif self.resname in prs.acids_set['TRP']:
            select = ['name CD2 CE2 CE3 CZ2 CZ3 CH2', 'name CG CD1 CD2 NE1 CE2']
        else:
            return None, None

        # Calculate centroid and normal vector coordinates
        atom_groups = [self.atoms.select_atoms(sel) for sel in select]
        centroids = [atom_group.center_of_geometry(compound='group') for atom_group in atom_groups]
        normals = [prs.calc_norm_vecs(atom_group, centroid) for atom_group, centroid in zip(atom_groups, centroids)]

        centroid_dict = {}
        normal_dict = {}
        if self.resname in prs.acids_set['TRP']:
            for i, name in enumerate([':HEX', ':PENT']):
                centroid_dict[self.res + name] = centroids[i]
                normal_dict[self.res + name] = normals[i]
        else:
            centroid_dict[self.res + ':NORMAL'] = centroids[0]
            normal_dict[self.res + ':NORMAL'] = normals[0]  # ex. {HIS:26:A:NORMAL : [x, y, z]}

        return centroid_dict, normal_dict

    def PiPi(self):
        """
        @description
            Consider amino acids involved in PiPi interactions (HIS, TYR, PHE, TRP).
            For these acids return name of residue, coordinates of centroid (self.centroid) and normal vector (self.normal)
        """
        if self.resname in prs.pication_acids:
            resname_key = list(self.centroid.keys())  # ex. ['TRP:26:A:HEX', 'TRP:26:A:PENT]
            centroid = list(self.centroid.values())
            norm = list(self.normal.values())
            return resname_key, centroid, norm

    def PiCation(self):
        """
        @description
            Consider amino acids having cation atoms NZ (LYS), CZ (ARG)
            and return cation atom names and its position
        """
        if self.resname in prs.cation_acids:
            cations = self.atoms.select_atoms('name NZ CZ')
            if cations.names.size != 0:  # check, b/c some residues cannot contain NZ CZ atoms
                cation_name = self.res + ':' + cations.names[0]  # ex. ARG:232:A:CZ
                cation_coords = cations.positions[0]
                return cation_name, cation_coords

    def Disulfide(self):
        """
        @description
            Consider amino acid CYS and its S* atoms and
            return name and coordinates of S* atom
        """
        if self.resname in prs.disulf_acids:
            atom = self.atoms.select_atoms('name SG')
            atom_name = self.res + ':' + atom.names[0]
            atom_coords = atom.positions[0]
            return atom_name, atom_coords

    def SaltBridge(self):
        """
        @description
            Consider basic (ARG, LYS, HIS) and acidic (ASP, GLU) amino acids involved in salt bridges
            for basic select N* atoms, for acidic O* atoms and return atoms names and coordinates
        """
        if self.resname in prs.basic_salt_aa:
            select = 'name NH* NZ* NE* ND*'
        elif self.resname in prs.acidic_salt_aa:
            select = 'name OE* OD*'
        else:
            return None

        atoms = self.atoms.select_atoms(select)
        atom_names = [self.res + ':' + name for name in atoms.names]
        atom_coords = atoms.positions
        return atom_names, atom_coords


def find_pipi_bonds(res1, res2):
    """
    @description
        Find PiPi interactions b/w atoms with not None 'pipi' attribute
        Calculate angle b/w normal vectors and distance b/w centroids
        based on cutoff conditions return PiPi bonds with energy (10)
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @output
        - list of PiPi bonds in format: [TRP:229:A:HEX  TRP:290:A:HEX   SCSC    10  PIPI]
    """
    if res1.pipi is None or res2.pipi is None:
        return None
    res1_name, centroid1, norm1 = res1.pipi  # ['TRP:28:A:HEX'], origin and normal coords
    res2_name, centroid2, norm2 = res2.pipi

    pipi_bonds = []
    for i in range(len(res1_name)):
        for j in range(len(res2_name)):
            pipiangle = np.degrees(mdamath.angle(norm1[i][0], norm2[j][0]))
            dist = mdadist(centroid1[i], centroid2[j])

            if pipiangle > prs.RIGHT_ANGLE:
                pipiangle = prs.STRAIGHT_ANGLE - pipiangle
            if (dist <= prs.PIPI_D1 and pipiangle < prs.PIPI_ANG1) or (dist <= prs.PIPI_D2 and prs.PIPI_ANG2 < pipiangle < prs.PIPI_ANG3):
                pipi_bonds.append(res1_name[i] + '\t' + res2_name[j] + '\tSCSC\t' + prs.PIPI_CONST + '\tPIPI\t' + '\n')
    return pipi_bonds


def find_pication_bonds(res1, res2):
    """
    @description
        Find PiCation interactions b/w Pi (HIS, TYR, PHE, TRP) and Cation-contained (ARG, LYS) acids
        Calculate vector coordinates and angle b/w cation atom and centroid ring of HIS, TYR, PHE, TRP
        Based on distance (6.6A) and angle cutoff (30°) conditions return PiCation bonds with energy (9.4)
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @return
        - list of PiCation bonds in format: [TRP:290:A:HEX  ARG:259:A:CZ    SCSC    9.4 PICAT]
    """
    if res1.pipi is not None and res2.pication is not None:
        res_name, centroid, norm = res1.pipi
        cation, pos = res2.pication
    elif res1.pication is not None and res2.pipi is not None:
        res_name, centroid, norm = res2.pipi
        cation, pos = res1.pication
    else:
        return None

    pication_bonds = []
    for i in range(len(res_name)):
        cat_vec = pos - centroid[i]
        cat_angle = np.degrees(mdamath.angle(norm[i][0], cat_vec))
        distance = mdadist(centroid[i], pos)  # distance b/w centroid positions and cation pos

        if cat_angle > prs.RIGHT_ANGLE:
            cat_angle = prs.STRAIGHT_ANGLE - cat_angle
        if distance < prs.PICAT_D and cat_angle < prs.PICAT_ANG:
            pication_bonds.append(res_name[i] + '\t' + cation + '\tSCSC\t' + prs.PICATION_CONST + '\tPICAT\t' + '\n')
    return pication_bonds


saltBridges = []


def find_salt_bridges(res1, res2):
    """
    @description
        Find SaltBridges b/w atoms with not None 'sb' attribute
        Consider only basic and acidic not neighbored residues
        Calculate distance b/w atoms and based on cutoff (4A) return 'salt_bridges' list
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @output
        - list of salt_bridges in format: [LYS:32:A:NZ    ASP:35:A:OD2    SCSC    20  SB]
    """
    if res1.sb is None or res2.sb is None:
        return None
    else:
        res1_name, coords1 = res1.sb
        res2_name, coords2 = res2.sb

    # check whether residues are basic or acidic
    name1, name2 = res1.resname, res2.resname
    if all(name in prs.basic_salt_aa for name in [name1, name2]) or all(name in prs.acidic_salt_aa for name in [name1, name2]):
        return None
    if np.abs(res1.acid.resid - res2.acid.resid) == 1:  # prs.seq_dist_cutoff
        return None

    salt_bridges = []
    for i in range(len(res1_name)):
        for j in range(len(res2_name)):
            dist = mdadist(coords1[i], coords2[j])
            if dist < prs.SALT_DISTANCE:
                saltBridges.append(res1_name[i] + '\t' + res2_name[j])
                salt_bridges.append(res1_name[i] + '\t' + res2_name[j] + '\tSCSC\t' + prs.SALT_CONST + '\tSB\t' + '\n')
    return salt_bridges


def find_disulfide_bonds(res1, res2):
    """
    @description
        Find disulfide bonds b/w SG atoms of CYS with not None 'disulf' attribute
        Calculate distance b/w SG atoms and based on cutoff (3A) return bonds
    @input
        - res1, AminoAcid(mda.residue) class
        - res2, AminoAcid(mda.residue) class
    @return
        - list of disulfide bonds in format: [CYS:77:A:SG    CYS:123:A:SG    SCSC    167 SS]
    """
    if res1.disulf is None or res2.disulf is None:
        return None
    else:
        acid1_name, coords1 = res1.disulf
        acid2_name, coords2 = res2.disulf

    dist = mdadist(coords1, coords2)
    if dist < prs.DISULF_D:
        bond = acid1_name + '\t' + acid2_name + '\tSCSC\t' + prs.DISULF_CONST + '\tSS\t' + '\n'
        return bond


def find_hydrogen_bonds(pdb, uni, allatoms_data, salt_bridges):
    """https://docs.mdanalysis.org/dev/documentation_pages/analysis/hydrogenbonds.html
    http://pubs.sciepub.com/ajme/3/2/3/index.html
    @description
        Use new version of HydrogenBondsAnalysis module
        Generate table of hydrogen bonds (MDAnalysis.HydrogenBondAnalysis) with parameters:
            universe=uni
            donors_sel, acceptors_sel - select N* and O* atoms for donors and acceptors
            hydrogens_sel  - select H* atoms for hydrogens
            d_h_a_angle_cutoff[90] - Angle cutoff, D-H-A (donor-hydrogen-acceptor), set to 90 to include water bonds
            d_h_cutoff[1.2] - Distance cutoff b/w Donor and Hydrogen atoms
            d_a_cutoff[3.5] - Distance cutoff b/w Donor and Acceptor atoms
        Table contains info about donor hydrogen and acceptor atoms
            | frame | donor_id | hydrogen_id | acceptor_id | distance | angle |
        Traverse table; define hydrogen, donor, acceptor atoms and find geometric parameters of hydrogen bonds, distances and angles
        Constraints:
            Not include bonds b/w donor and acceptor that are involved in salt_bridges.
            Hydrogen and acceptor must be from non-adjacent residues
            Confine distances b/w  hydrogen-acceptor and donor-acceptor
    @input
        pdb - used for write hydrogen_bonds file
        uni - MDAnalysis Universe used to call HydrogenBondAnalysis
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
        saltBridges - list of salt bridges, hydrogen bonds should not be salt bridges
    @return
        file with hydrogen bonds in format
            acid1       | acid2         | Chain_type | Distance | Bond_type
            GLU:28:A:H  | HIS:26:A:ND1  | MCSC       | 15.52    | HB
    """
    print('Search hydrogen bonds in protein and among water molecules')

    hbonds = HBA(universe=uni,
                 donors_sel='all and name N* O*',
                 hydrogens_sel='all and name H*',
                 acceptors_sel='all and name N* O*',
                 d_h_a_angle_cutoff=90,
                 d_h_cutoff=1.2,
                 d_a_cutoff=3.5)
    hbonds.run()

    hydrogenBonds = []
    acceptorBonds = {}
    for hbond in hbonds.hbonds:
        # d-donor, h-hydrogen, a-acceptor
        d, h, a, d_a_dist, d_h_a_angle = hbond[1:]  # return indexes of uni.atoms
        d, h, a = uni.atoms[int(d)], uni.atoms[int(h)], uni.atoms[int(a)]
        donor = ':'.join([d.resname, str(d.resid), d.chainID, d.name])
        hydrogen = ':'.join([h.resname, str(h.resid), h.chainID, h.name])
        acceptor = ':'.join([a.resname, str(a.resid), a.chainID, a.name])

        acc_res, acc_name = a.resname, a.name  # acceptor residue name and atom name
        dA = ':'.join([a.resname, a.name])
        dH = ':'.join([h.resname, h.name])
        if dH not in prs.hydrogen:  # hydrogen atom must be in prs.hydrogen
            continue
        if not prs.isDonor(donor):  # check whether is donor or not
            continue
        if (donor + '\t' + acceptor) in salt_bridges or (acceptor + '\t' + donor) in salt_bridges:  # not include bonds from salt bridges
            continue

        if np.abs(h.resid - a.resid) > 1:  # prs.seq_dist_cutoff, not include non-adjacent acids
            if (dA in prs.nextAcid) and prs.isAcceptor(acceptor):
                alpha = d_h_a_angle
                d1 = mdadist(allatoms_data[hydrogen]['coords'], allatoms_data[acceptor]['coords'])
                d2 = mdadist(allatoms_data[donor]['coords'], allatoms_data[acceptor]['coords'])
                if d1 > prs.H_A_dist_cutoff or d2 > prs.D_A_dist_cutoff:
                    continue

                # Define acceptor antecedent (AA) atom and its coords
                AA_coords = None
                aNext = prs.nextAcid[dA]
                if len(aNext) == 2:
                    neighbor1 = ":".join(acceptor.split(":")[0:3]) + ":" + aNext[0]
                    neighbor2 = ":".join(acceptor.split(":")[0:3]) + ":" + aNext[1]
                    if neighbor1 in allatoms_data and neighbor2 in allatoms_data:
                        AA_coords = (allatoms_data[neighbor1]['coords'] + allatoms_data[neighbor2]['coords']) / 2
                else:
                    neighbor1 = ":".join(acceptor.split(":")[0:3]) + ":" + aNext[0]
                    AA_coords = allatoms_data[neighbor1]['coords']
                # here problem with sol in frames when add them based on distance
                # not all atoms of one water mol are added, that's why it is not
                # read aNext in neighbor1

                if AA_coords is not None:
                    # Define beta angle b/w AA--acceptor--hydrogen
                    a = AA_coords - allatoms_data[acceptor]['coords']
                    b = allatoms_data[hydrogen]['coords'] - allatoms_data[acceptor]['coords']
                    beta = np.degrees(mda.lib.mdamath.angle(a, b))

                    # Define gamma angle b/w AA--acceptor--donor
                    a = AA_coords - allatoms_data[acceptor]['coords']
                    b = allatoms_data[donor]['coords'] - allatoms_data[acceptor]['coords']
                    gamma = np.degrees(mda.lib.mdamath.angle(a, b))

                    if ((prs.hydrogen[dH] != 'S' and acc_name != 'S' and d1 < prs.hb_d1 and d2 < prs.hb_d2 and alpha > 90 and beta > 90 and gamma > 90) or
                            (prs.hydrogen[dH] == 'S' and acc_name != 'S' and d1 < prs.hb_d1 and d2 < prs.hb_d2 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                            (prs.hydrogen[dH] != 'S' and acc_name == 'S' and d1 < prs.hb_d11 and d2 < prs.hb_d21 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE) or
                            (prs.hydrogen[dH] == 'S' and acc_name == 'S' and d1 < prs.hb_d12 and d2 < prs.hb_d22 and alpha > prs.RIGHT_ANGLE and beta > prs.RIGHT_ANGLE and gamma > prs.RIGHT_ANGLE)):

                        if acceptor not in acceptorBonds.keys():
                            acceptorBonds[acceptor] = []
                        if acceptor.split(':')[3] in prs.acceptor_atoms1:
                            acceptorBonds[acceptor].append(d1)
                            if len(acceptorBonds[acceptor]) > 2:
                                acceptorBonds[acceptor].sort()
                                acceptorBonds[acceptor] = acceptorBonds[acceptor][0:2]
                        if acceptor.split(':')[3] in prs.acceptor_atoms2 and (acceptor.split(':')[0] != 'SOL' or acceptor.split(':')[0] != 'HOH'):
                            acceptorBonds[acceptor].append(d1)
                            if len(acceptorBonds[acceptor]) > 1:
                                acceptorBonds[acceptor].sort()
                                acceptorBonds[acceptor] = acceptorBonds[acceptor][0:1]
                        if acceptor.split(':')[3] == "O" and acceptor.split(':')[0] == "HOH":
                            acceptorBonds[acceptor].append(d1)
                        if acceptor.split(':')[3] == 'OW' and acceptor.split(':')[0] == 'SOL':  # proceed frames
                            acceptorBonds[acceptor].append(d1)
                        # beta = prs.STRAIGHT_ANGLE-beta if beta<prs.RIGHT_ANGLE else beta

                        E = None
                        if prs.SPHyb(donor) == "SP3" and prs.SPHyb(acceptor) == "SP3":
                            E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP3', sphyb2='SP3', alpha=alpha, beta=beta)

                        elif prs.SPHyb(donor) == "SP3" and prs.SPHyb(acceptor) == "SP2":
                            E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP3', sphyb2='SP2', alpha=alpha, beta=beta)

                        elif prs.SPHyb(donor) == "SP2" and prs.SPHyb(acceptor) == "SP3":
                            E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP2', sphyb2='SP3', alpha=alpha, beta=beta)

                        elif prs.SPHyb(donor) == "SP2" and prs.SPHyb(acceptor) == "SP2":
                            normalVecDonor = prs.normalDonorVecToPlane(donor, allatoms_data)
                            if normalVecDonor is None:
                                continue
                            normalVecAcceptor = prs.normalAcceptorVecToPlane(acceptor, allatoms_data)
                            if normalVecAcceptor is None:
                                continue
                            psi = np.degrees(mda.lib.mdamath.angle(normalVecDonor, normalVecAcceptor))
                            E = prs.HydrogenBondEnergy(dist=d2, sphyb1='SP2', sphyb2='SP2', alpha=alpha, beta=beta, psi=psi)
                        if E is not None:
                            hydrogenBonds.append(hydrogen + '\t' + acceptor + '\t' + donor + '\t' + allatoms_data[donor]['chain'] + allatoms_data[acceptor]['chain'] + '\t' + str(d1) + '\t' + str(E))

    with open(pdb.replace('.pdb', '_hb'), 'w') as out:
        donorBonds = {}
        for interaction in hydrogenBonds:
            hydrogen, acceptor, donor, chain, d1, E = interaction.split("\t")[0:6]
            bond = hydrogen + '\t' + acceptor + '\t' + chain + '\t' + E + '\tHB\t\n'

            if donor not in donorBonds:
                donorBonds[donor] = 0
            if donor.split(':')[0] == "HOH" or acceptor.split(':')[0] == "HOH":
                out.write(bond)
            if donor.split(':')[0] == 'SOL' or acceptor.split(':')[0] == 'SOL':
                out.write(bond)

            if len(acceptorBonds[acceptor]) == 2 and (d1 == str(acceptorBonds[acceptor][0]) or d1 == str(acceptorBonds[acceptor][1])):
                if donor.split(':')[3] not in (prs.donor_atoms1 + prs.donor_atoms2) and donorBonds[donor] == 0:
                    out.write(bond)
                    donorBonds[donor] += 1
                elif donor.split(':')[3] in prs.donor_atoms1 and donorBonds[donor] < 2:
                    out.write(bond)
                    donorBonds[donor] += 1
                elif donor.split(':')[3] in prs.donor_atoms2 and donorBonds[donor] < 3:
                    out.write(bond)
                    donorBonds[donor] += 1
            elif d1 == str(acceptorBonds[acceptor][0]):
                if donor.split(':')[3] not in (prs.donor_atoms1 + prs.donor_atoms2) and donorBonds[donor] == 0:
                    out.write(bond)
                    donorBonds[donor] += 1
                elif donor.split(':')[3] in prs.donor_atoms1 and donorBonds[donor] < 2:
                    out.write(bond)
                    donorBonds[donor] += 1
                elif donor.split(':')[3] in prs.donor_atoms2 and donorBonds[donor] < 3:
                    out.write(bond)
                    donorBonds[donor] += 1
                elif donor.split(':')[3] in ["0", 'OW']:
                    out.write(bond)


def find_vdw_bonds(pdb, protein, allatoms, allatoms_data):
    """
    @description
        Calculate capped distance b/w protein atoms with distance cutoff (4.25), which was selected as sum of max radii among atoms
        Traverse pairs and calculate summed radii of two considered atoms (rm) and distance b/w them (r)
        For interacted atoms with the difference (r-rm < 0.5) and non-adjacent residues calculate Energy value and write into file
    @input
        protein - mda object, protein = uni.select_atoms('protein').atoms
        allatoms - list of atoms in format ['HIS:26:A:N']
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    @return
        file with VanderWaals bonds in format, which then can be modified not to include bonds with negative energies
            Atom1      | Atom2       | Chain_type | Energy
            HIS:26:A:N | GLU:58:A:CD | MCSC       | 3.470
    """
    print('Search VanderWaals bonds')
    pairs, distances = mda.lib.distances.capped_distance(protein.positions, protein.positions, max_cutoff=4.25, min_cutoff=0)
    # create array with sorted first numbers and other atoms' index from pairs
    van_atom = np.split(pairs[:, 1], np.cumsum(np.unique(pairs[:, 0], return_counts=True)[1])[:-1])
    vdw1, vdw2 = {}, {}

    for i, k in enumerate(van_atom):
        elem1 = allatoms[i]  # ex. HIS:26:A:CB
        res1 = ':'.join(elem1.split(':')[:3])  # ex. HIS:26:A
        atom1 = elem1.split(':')[3]  # ex. CB
        for j in k:
            elem2 = allatoms[j]
            res2 = ':'.join(elem2.split(':')[:3])
            atom2 = elem2.split(':')[3]
            if (elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain'] in vdw1 or
                    elem2 + '\t' + elem1 + '\t' + allatoms_data[elem2]['chain'] + allatoms_data[elem1]['chain'] in vdw1):
                continue
            if not res1 == res2 and atom1 in prs.radii and atom2 in prs.radii:
                rm = prs.radii[atom1] + prs.radii[atom2]
                r = mdadist(allatoms_data[elem1]['coords'], allatoms_data[elem2]['coords'])
                # not include adjacent residues and atoms with MC type
                if r - rm < .5 and not (np.abs(int(elem1.split(':')[1]) - int(elem2.split(':')[1])) == 1 and allatoms_data[elem1]['chain'] == "MC" and allatoms_data[elem2]['chain'] == "MC"):
                    if not elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain'] in vdw1:
                        vdw1[elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain']] = []

                    E = prs.energy_vdw(rm, r)
                    vdw1[elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain']].append(E)

                    if (("C" in atom1 and "C" in atom2) or
                            ("C" in atom1 and atom2 in ["NE2", "OE1", "ND2", "OD1"] and res2.split(":")[0] in ["GLN", "ASN"]) or
                            ("C" in atom2 and atom1 in ["NE2", "OE1", "ND2", "OD1"] and res1.split(":")[0] in ["GLN", "ASN"])):

                        if not elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain'] in vdw2:
                            vdw2[elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain']] = []

                        vdw2[elem1 + '\t' + elem2 + '\t' + allatoms_data[elem1]['chain'] + allatoms_data[elem2]['chain']].append(E)

    with open(pdb.replace(".pdb", "_vdw"), 'w') as out:
        for contact in vdw1:
            if not (sum(vdw1[contact]) < 0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1])) == 1):
                out.write(contact + '\t' + str(sum(vdw1[contact])) + '\tVDW\n')

    with open(pdb.replace(".pdb", "_vdw2"), 'w') as out:
        for contact in vdw2:
            if not (sum(vdw1[contact]) < 0 and abs(int(contact.split('\t')[0].split(':')[1]) - int(contact.split('\t')[1].split(':')[1])) == 1):
                out.write(contact + '\t' + str(sum(vdw1[contact])) + '\tVDW2\n')


def find_metal_bonds(metals, acids_class, allatoms_data, pdb):
    """
    @description
        For each metal find atoms from acids_class that are bonded to that metal
        The distance b/w metal and atom must be < 3.5A (metal_dist_cutoff) and
        correspond to distance cutoffs b/w acid atoms and metals (prs.METALS_DICT)
        Write to list (metal2atom) appropriate bonds b/w metal and atoms
        Traverse 2 atoms in metal2atom, if angle b/w them is in [90, 180] write these atoms which are connected through metal
    @input
        metals - mda selection of metals
        acids_class - created list of class of amino acids
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    @return
        file in format
            Acid1      | Acid2      | Chain_type | Energy | Bond_type
            HIS:26:A:C | GLU:58:A:N | MCSC       | 3      | METAL
    """
    print('Search bonds b/w metals and atoms')
    metal2atom = []
    metalBonds = []
    for metal in metals:
        met_pos = metal.position
        for i in range(len(acids_class)):
            acid = acids_class[i]
            for atom in acid.atoms:
                dist = mdadist(met_pos, atom.position)
                if dist < prs.metal_dist_cutoff and acid.resname in prs.METALLS_DICT:
                    met_dist, atoms = prs.METALLS_DICT[acid.resname]
                    if metal.resname in met_dist:
                        d = met_dist[metal.resname]
                        if acid.resname == 'ASP' and atom.name in atoms:
                            if dist < prs.GLUASPmonoDist[metal.resname]:
                                metal2atom.append(atom)  # append atom mda object
                        elif acid.resname == 'GLU' and atom.name in atoms:
                            if dist < prs.GLUASPmonoDist[metal.resname]:
                                metal2atom.append(atom)
                        elif dist <= d and atom.name in atoms:
                            metal2atom.append(atom)

        for i in range(len(metal2atom)):
            for j in range(i + 1, len(metal2atom)):
                atom1, atom2 = metal2atom[i], metal2atom[j]
                atom_1 = atom1.resname + ':' + str(atom1.resid) + ':' + atom1.chainID + ':' + atom1.name
                atom_2 = atom2.resname + ':' + str(atom2.resid) + ':' + atom2.chainID + ':' + atom2.name
                atom1_pos = atom1.position
                atom2_pos = atom2.position
                coords1 = met_pos - atom1_pos  # distance b/w metal and atom1
                coords2 = met_pos - atom2_pos  # distance b/w metal and atom2
                angle = np.degrees(mdamath.angle(coords1, coords2))

                if prs.RIGHT_ANGLE < angle < prs.STRAIGHT_ANGLE:
                    metalBonds.append(atom_1 + '\t' + atom_2 + '\t' + allatoms_data[atom_1]['chain'] + allatoms_data[atom_2]['chain'] + '\t3\tMETAL\n')

    with open(pdb.replace('.pdb', '_metal'), 'a') as out:
        for bond in list(set(metalBonds)):
            out.write(bond)


def find_dna_bonds(uni, allatoms_data, pdb):
    """
    @description
        Select DNA resnames (DA, DG, DC, DT).
        Iterate through 2 atoms around DNA resname, if angle b/w these 2 atoms suits to conditions.
        Append to list DNAbindingPairs and write output file
    @input
        uni - MDAnalysis Universe object
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    @return
        file in format
        DNA_atom | Protein_atom | Chaintype | Value | BondType |
    """
    print('Search DNA bonds.')

    dna = uni.select_atoms('resname DA DG DC DT')
    atoms_around_dna = uni.select_atoms('protein and around 7.5 resname DA DG DC DT')  # only consider protein atoms
    around_atoms = [i + ':' + str(j) + ':' + k + ':' + l for i, j, k, l in zip(atoms_around_dna.resnames, atoms_around_dna.resids, atoms_around_dna.segids, atoms_around_dna.names)]

    dna_atoms = [i + ':' + str(j) + ':' + k + ':' + l for i, j, k, l in zip(dna.resnames, dna.resids, dna.segids, dna.names)]
    dna_data = {atom: pos for atom, pos in zip(dna_atoms, dna.positions)}
    DNAbindingPairs = []

    for dna_atom, dna_coords in dna_data.items():  # DA:9:D:N1
        for atom1 in around_atoms:
            bound = False
            if (("N" in dna_atom.split(':')[3] and (prs.isDonor(atom1) or prs.isAcceptor(atom1))) or
                    ("O" in dna_atom.split(':')[3] and (prs.isDonor(atom1) or prs.isAcceptor(atom1))) or
                    (("N" in dna_atom.split(':')[3] or "O" in dna_atom.split(':')[3]) and "S" in atom1.split(':')[3]) or
                    ("C" in dna_atom.split(':')[3] and "O" in atom1.split(':')[3] and prs.isAcceptor(atom1)) or
                    ("O" in dna_atom.split(':')[3] and atom1.split(':')[3] in ["NZ", "NH1", "NH2", "OD1", "OD2", "OE1", "OE2"]) or
                    ("C" in dna_atom.split(':')[3] and "C" in atom1.split(':')[3])):

                for atom2 in around_atoms:
                    if atom2 == atom1 or atom2.split(':')[3] == 'H':
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

    with open(pdb.replace('.pdb', '_DNAbonds'), 'w') as out:
        for dna_atom, prot_atom in DNAbindingPairs:
            out.write(dna_atom + '\t' + prot_atom + '\tMC' + allatoms_data[prot_atom]['chain'] + '\t10\tDNA\tNT\n')


def find_ligand_atom_bonds(ligand_centroids, allatoms_data, pdb):  # arg - prot_atoms to consider only protein residues
    """
    @description
        For each protein 'SC' atom find closest ligand centroid
        Create 'ligand_distances' dict and write tmp info about min distance to ligand and ligand
        Write ligand-atom bonds into file
    @input
        ligand_centroids - dict of ligand names and centroid coords
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    @return
        file in format:
            Atom       | Ligand     | Distance
            ALA:0:A:CB | CLR:2403:A | 30.8
    """
    print('Search bonds b/w protein atoms and closest ligands.')

    ligand_distances = {}
    with open(pdb.replace('.pdb', '_ligands'), 'w') as out:
        for atom, atom_data in allatoms_data.items():
            if atom_data['chain'] == 'SC':
                ligand_distances[atom] = {}
                ligand_distances[atom]['dist'] = prs.ligand_dist  # tmp distance - 150
                ligand_distances[atom]['ligand'] = ''  # tmp ligand name
                for ligand, ligand_coords in ligand_centroids.items():
                    dist = mdadist(atom_data['coords'], ligand_coords)
                    if dist > ligand_distances[atom]['dist']:
                        continue
                    ligand_distances[atom]['dist'] = dist
                    ligand_distances[atom]['ligand'] = ligand
                    out.write(atom + '\t' + ligand + '\t' + str(ligand_distances[atom]['dist']) + '\n')


def find_centroid_centroid_bonds(protein_centroids, pdb):
    """
    @description
        Traverse protein centroids, calculate distance b/w them
        based on distance cutoff (8.5A) find bonds and write into file
    @input
        protein_centroids - dict of residues names and centroid coords
    @return
        file in format:
            Residue1    | Residue2  | Chain_type | Energy | Distance
            ALA:0:A     | GLY:5:A   | SCSC       | 1      | 10
    """
    print('Search protein centroid-centroid interactions')
    with open(pdb.replace('.pdb', '_centroidNetSC'), 'w') as out:
        for centroid1, coords1 in protein_centroids.items():
            for centroid2, coords2 in protein_centroids.items():
                dist = mdadist(coords1, coords2)
                if 0 < dist < prs.centroid_dist:
                    out.write(centroid1 + '\t' + centroid2 + '\tSCSC\t1\tCENTROID\tCENTROID\tCENTROID\t' + str(dist) + '\n')


def find_centroid_ligand_bonds(protein_centroids, ligand_centroids, pdb):
    """
    @description
        For each protein centroid find closest ligand centroid
        Create 'ligand_distances' dict and write tmp info about min distance to ligand and ligand
        Search for min distance b/w residue centroid and ligand centroid, and write it into file
    @input
        centroid_coords - dict of residues names and their centroid coords
        ligand_centroids - dict of ligands and their centroid coords
    @output
        file in format:
            Residue |   Ligand  |   Distance
    """
    print('Search protein centroid-closest ligand interactions')

    ligand_distances = {}
    with open(pdb.replace('.pdb', '_centroidNetLigand'), 'w') as out:
        for centroid, centroid_coords in protein_centroids.items():
            ligand_distances[centroid] = {}
            ligand_distances[centroid]['dist'] = prs.ligand_dist  # tmp distance 150
            ligand_distances[centroid]['ligand'] = ''  # tmp ligand
            for ligand, ligand_coords in ligand_centroids.items():
                dist = mdadist(ligand_coords, centroid_coords)
                if dist > ligand_distances[centroid]['dist']:
                    continue
                ligand_distances[centroid]['dist'] = dist
                ligand_distances[centroid]['ligand'] = ligand
                out.write(centroid + '\t' + ligand + '\t' + str(ligand_distances[centroid]['dist']) + '\n')


def find_peptide_pairs(protein, pdb):
    """
    @description
        Find peptide pairs b/w two successive residues in amino acid sequence
    @input
        protein - MDAnalysis object, u.select_atoms('protein')
    @output
        file with peptide pairs in format: [HIS:26:A:PP1	PRO:27:A:PP2	MCMC	10	PP]
    """
    residues = protein.residues
    with open(pdb.replace('.pdb', '_polypeptidePairs'), 'w') as out:
        for i in range(protein.n_residues - 1):
            segid1 = residues[i].atoms.chainIDs[0]
            segid2 = residues[i + 1].atoms.chainIDs[0]
            # when separate trajectory to pdb files, in new created pdbs Residue objects dont contain chainIDs, so select from atoms
            if residues[i].resid == residues[i + 1].resid - 1:
                out.write(residues[i].resname + ':' + str(residues[i].resid) + ':' + segid1 + ':PP1' + '\t' +
                          residues[i + 1].resname + ':' + str(residues[i + 1].resid) + ':' + segid2 + ':PP2' + '\tMCMC\t10\tPP\n')


# def prepare_bfactor_file(protein):
#     """
#     @description
#         Select "CA" atoms of residues and extract B-factor(tempfactor)
#     @input
#         protein - MDAnalysis object, u.select_atoms('protein')
#     @output
#         file in format:
#             Acid     | B-factor
#             HIS:26:A | 14.48
#     """
#     CA_atoms = protein.select_atoms('name CA')
#     with open('Bfactor', 'w') as bfact:
#         for atom in CA_atoms.atoms:
#             bfact.write(atom.resname + ':' + str(atom.resid) + ':' + atom.segid + "\t" + str(atom.tempfactor) + "\n")


def filter_neg_vdw_bonds(vdw_file):
    """
    @description
        Read vdw_file and process VanderWaals bonds.
        Compute total sum of energies for same connected amino acids,
        if this total energy is negative -> don't include interactions b/w these amino acids
    @input
        vdw_file - a created file by function 'find_vds_bonds'
    @output
        processed file of VanderWaals bonds with positive total energy, in format
            Atom1      | Atom2       | Chain_type | Energy
            HIS:26:A:N | GLU:58:A:CD | MCSC       | 3.470
    """
    with open(vdw_file, 'r') as vdw_out:
        lines = vdw_out.readlines()
        bond_energies = {}  # dict of bonds and summed total energy
        for line in lines:
            acid1, acid2, chain, energy, _ = line.split('\t')
            res1 = ':'.join(acid1.split(':')[0:3])
            res2 = ':'.join(acid2.split(':')[0:3])
            bond1 = res1 + '\t' + res2
            bond2 = res2 + '\t' + res1
            if bond1 not in bond_energies and bond2 not in bond_energies:
                bond_energies[bond1] = float(energy)
            elif bond1 in bond_energies and bond2 not in bond_energies:
                bond_energies[bond1] += float(energy)
            elif bond1 not in bond_energies and bond2 in bond_energies:
                bond_energies[bond2] += float(energy)
        # Define bonds with negative total energy
        neg_vdw_bonds = dict((bond, energy) for bond, energy in bond_energies.items() if energy <= 0)

    with open(vdw_file + '_noRepulse', 'w') as out:
        for line in lines:  # travers vdw_out.readlines
            if any(bond.split('\t')[0] in line and bond.split('\t')[1] in line for bond in list(neg_vdw_bonds.keys())):
                continue
            else:
                out.write(line.replace('\n', '\t') + '\n')


def find_waterbonds(pdb, hb_file, allatoms_data):
    """
    @description
        Traverse lines of hb_file, create dicts 'HOH' of water and 'HOH_energy' of energies
        Search bonds b/w atoms connected through one or two atoms of water
    @input
        hb_file - file with hydrogen bonds contained bonds with water molecules (HOH, SOL)
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    @output
        file in format:
            HIS:26:A    SO4:291:A   30.52
    """
    print('Find hydrogen bonds connected via water atoms')

    hb_file = open(hb_file, 'r').readlines()
    HOH = {}  # key - water, value - atom
    HOH2 = {}  # key - water, value - atom
    HOH_energy = {}  # key - atom and water, value - energy from hb_file

    for line in hb_file:
        line = line.split('\t')
        if ('HOH' in line[0] and 'HOH' not in line[1]) or ('SOL' in line[0] and 'SOL' not in line[1]):
            water = ':'.join(line[0].split(':')[:3])
            if water not in HOH:
                HOH[water] = []
            if (line[1], water) not in HOH_energy:
                HOH_energy[line[1], water] = []
            HOH[water].append(line[1])
            HOH[water] = list(set(HOH[water]))
            HOH_energy[line[1], water].append(float(line[3]))
            HOH_energy[line[1], water] = list(set(HOH_energy[line[1], water]))

        if ('HOH' in line[1] and 'HOH' not in line[0]) or ('SOL' in line[1] and 'SOL' not in line[0]):
            water = ':'.join(line[1].split(':')[:3])
            if water not in HOH:
                HOH[water] = []
            if (line[0], water) not in HOH_energy:
                HOH_energy[line[0], water] = []
            HOH[water].append(line[0])
            HOH[water] = list(set(HOH[water]))
            HOH_energy[line[0], water].append(float(line[3]))
            HOH_energy[line[0], water] = list(set(HOH_energy[line[0], water]))

        if ('HOH' in line[0] and 'HOH' in line[1]) or ('SOL' in line[0] and 'SOL' in line[1]):
            water1 = ':'.join(line[0].split(':')[:3])
            water2 = ':'.join(line[1].split(':')[:3])
            if water1 not in HOH2:
                HOH2[water1] = []
            if water2 not in HOH2:
                HOH2[water2] = []
            HOH2[water1].append(water2)
            HOH2[water2].append(water1)
            HOH2[water1] = list(set(HOH2[water1]))
            HOH2[water2] = list(set(HOH2[water2]))

    # handle single water interactions b/w atoms
    with open(pdb.replace('.pdb', '_waterbonds'), 'w') as out:
        for water in HOH:
            if len(HOH[water]) > 1:
                x = list(itertools.combinations(HOH[water], 2))
                for item in x:
                    if item[0].split(':')[:3] != item[1].split(':')[:3]:
                        energy = sum([*HOH_energy[item[0], water], *HOH_energy[item[1], water]]) / 2
                        out.write(item[0] + '\t' + item[1] + '\t' + allatoms_data[item[0]]['chain'] + allatoms_data[item[1]]['chain'] + '\t' + str(energy) + '\tWATER\n')
    # handle double water interactions b/w atoms
    with open(pdb.replace('.pdb', '_waterbonds'), 'a') as out:
        for water1 in HOH:
            for atom1 in HOH[water1]:
                if water1 in HOH2:
                    for water2 in HOH2[water1]:
                        if water2 in HOH:
                            for atom2 in HOH[water2]:
                                if atom1 != atom2 and atom1.split(':')[:3] != atom2.split(':')[:3]:
                                    energy = sum([*HOH_energy[atom1, water1], *HOH_energy[atom2, water2]]) / 2
                                    out.write(atom1 + '\t' + atom2 + '\t' + allatoms_data[atom1]['chain'] + allatoms_data[atom2]['chain'] + '\t' + str(energy) + '\tWATER\n')


def create_atoms_data(u):
    print('Create dict of atoms coords and chain types')
    prot = u.atoms  # protein.atoms
    prot_resnames, prot_resids, prot_segids, prot_atoms = prot.resnames, prot.resids, prot.chainIDs, prot.names
    allatoms = [(i + ':' + str(j) + ':' + k + ':' + l) for i, j, k, l in zip(prot_resnames, prot_resids, prot_segids, prot_atoms)]  # format 'HIS:26:A:N'

    allatoms_data = {}  # define dict of atoms' coords and chain types
    for atom, position in zip(allatoms, prot.positions):
        allatoms_data[atom] = {}
        allatoms_data[atom]['coords'] = position
        allatoms_data[atom]['chain'] = ''
        if atom.split(':')[3] in ["N", "O", "C", "CA", "HA2", "HA3"]:
            if atom.split(':')[0] == 'GLY' and atom.split(':')[3] in ["O", "CA", "HA2", "HA3", "N", "NH"]:
                allatoms_data[atom]['chain'] = 'SC'
            else:
                allatoms_data[atom]['chain'] = 'MC'
        else:
            allatoms_data[atom]['chain'] = 'SC'
    return allatoms, allatoms_data


def main(pdb):
    """
    @description
        Create phi_file using stride program. Define Universe, protein mda objects
        Define atoms positions and side chain types. Search PiPi, PiCation, Hydrogen, VanderWaals, SaltBridges bonds
        Find centroid interactions
    @input
        pdb - pdb file
    """
    phi = pdb.replace('pdb', 'phi')
    os.system(f"stride {pdb} > {phi}")  # stride must be located in usr/local/bin
    prepare_secondary_structure_file(pdb, phi)
    prepare_rsa_file(pdb, phi)

    u = mda.Universe(pdb)
    protein = u.select_atoms('protein and not name OXT')
    # hoh = u.select_atoms('resname HOH or resname SOL')
    metals = u.select_atoms('resname {}'.format(' '.join(list(prs.metals))))
    ligands = u.select_atoms('not protein and not resname HOH and not resname SOL') - metals
    # protein = u.select_atoms('protein and resid 0:72 and not name OXT')  # for 2oob ubiquitin
    # ligands = u.select_atoms('resid 929:972')  # for 2oob ubiquitin

    allatoms, allatoms_data = create_atoms_data(u)

    print('Create for each residue `AminoAcid` class object')
    acids_class = [AminoAcid(res) for res in protein.residues]

    print('Find pi-pi; pi-cation; disulfide bonds and salt bridges')
    with open(pdb.replace('.pdb', '_bonds'), 'w') as net:
        for i in range(len(acids_class)):
            for j in range(i + 1, len(acids_class)):
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

    find_hydrogen_bonds(pdb, u, allatoms_data, saltBridges)
    find_waterbonds(pdb, pdb.replace('.pdb', '_hb'), allatoms_data)
    find_vdw_bonds(pdb, protein.atoms, allatoms, allatoms_data)
    filter_neg_vdw_bonds(pdb.replace('.pdb', '_vdw'))
    filter_neg_vdw_bonds(pdb.replace('.pdb', '_vdw2'))
    # find_metal_bonds(metals, acids_class, allatoms_data, pdb)
    find_peptide_pairs(protein.atoms, pdb)
    # find_dna_bonds(u, allatoms_data, pdb)

    # Calculate weighted center of residues, where centroid_data.masses are weights
    centroid_data = protein.select_atoms('(resname GLY and not name N C O) or (protein and not backbone)')
    center_coords = centroid_data.center(centroid_data.masses, compound='residues')
    centroid = centroid_data.residues
    centroid_names = [i + ':' + str(j) + ':' + k[0] for i, j, k in zip(centroid.resnames, centroid.resids, centroid.chainIDs)]  # chainIDs vs segids. here chainIDs for residues -> take k[0]
    protein_centroids = dict(zip(centroid_names, center_coords))

    find_centroid_centroid_bonds(protein_centroids, pdb)

    if ligands.n_residues:  # for existed ligands
        ligand_names = [i + ':' + str(j) + ':' + k for i, j, k in zip(ligands.residues.resnames, ligands.residues.resids, ligands.residues.chainIDs[0])]  # chainIDs vs segids
        ligand_centroids = dict(zip(ligand_names, ligands.center_of_geometry(compound='residues')))
        find_ligand_atom_bonds(ligand_centroids, allatoms_data, pdb)
        find_centroid_ligand_bonds(protein_centroids, ligand_centroids, pdb)
    else:
        open(pdb.replace('.pdb', '_ligands'), 'w').close()
        open(pdb.replace('.pdb', '_centroidNetLigand'), 'w').close()

    pdb = pdb.replace('.pdb', '')  # path/to/pdb
    # Merge all found bonds to net_file; select bonds with E>0
    cmd = f"cat {pdb}_polypeptidePairs {pdb}_vdw_noRepulse {pdb}_bonds {pdb}_waterbonds {pdb}_vdw2_noRepulse {pdb}_metal  > {pdb}_net; \
            cat {pdb}_hb | sed '/HOH/d; /SOL/d' >> {pdb}_net; \
            awk '$4>0' {pdb}_net > {pdb}tmp; \
            mv {pdb}tmp {pdb}_net; \
            rm -rf {pdb}_bonds {pdb}_vdw {pdb}_vdw2 {pdb}_vdw_noRepulse {pdb}_vdw2_noRepulse {pdb}_metal {pdb}_polypeptidePairs {pdb}_hb {pdb}_waterbonds"
    os.system(cmd)


if __name__ == '__main__':
    t0 = time.time()

    PDB = sys.argv[1]  # input PDB structure
    main(PDB)

    print(time.time() - t0)

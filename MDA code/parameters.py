import numpy as np
import math

# constant values of angles and distances for salt bridges, pipi, pication
# VanderWaals and hydrogen bonds, for ligands and centroid interactions
SALT_CONST = '20'
SALT_DISTANCE = 4

PIPI_CONST = '10'
PIPI_D1 = 4.4
PIPI_ANG1 = 30
PIPI_D2 = 5.5
PIPI_ANG2 = 60
PIPI_ANG3 = 120

PICATION_CONST = '9.4'
PICAT_D = 6.6
PICAT_ANG = 30

DISULF_CONST = '167'
DISULF_D = 3

VDW_D = .5
FULL_ANGLE = 360
STRAIGHT_ANGLE = 180
RIGHT_ANGLE = 90
ZERO_ANGLE = 0

seq_dist_cutoff = 1
H_A_dist_cutoff = 3.5  # used in hydrogen bonds as hydrogen-acceptor distance cutoff
D_A_dist_cutoff = 4.5  # used in hydrogen bonds as donor-acceptor distance cutoff
hb_d1 = 2.5
hb_d11 = 3.2
hb_d12 = 3.3
hb_d2 = 3.9
hb_d21 = 4.1
hb_d22 = 4.5

ligand_dist = 150
centroid_dist = 8.5

acceptor_atoms1 = ["OD1", "OD2", "OE1", "OE2", "OG", "OG1", "SG"]
acceptor_atoms2 = ["O", "ND1", "NE2", "OH", 'OW']
donor_atoms1 = ["NH1", "NH2", "ND2", "NE2"]
donor_atoms2 = ["NZ"]


def energy_vdw(rm, r):
    """
    energy formula for VanderWaals bonds
    r - summed radii of 2 atoms
    rm - distance b/w 2 atoms
    """
    E = (-.997 * ((rm / r) ** 12 - 2 * (rm / r) ** 6)) * 4
    return E


def HydrogenBondEnergy(dist, sphyb1, sphyb2, alpha, beta, psi=None):
    """
    Energy for various hybridization
    Ef = V*[5*(do/d)^12 - 6*(do/d)^10]*F where do is 2.8, d is the distance between donor and acceptor atom (ie d2), V is 8
    sp3 donor - sp3 acceptor F=cos^2(alpha)cos^2(beta-109.5)
    sp3 donor - sp2 acceptor F=cos^2(alpha)cos^2(beta)
    sp2 donor - sp3 acceptor F=cos^4(alpha)
    sp2 donor - sp2 acceptor F=cos^2(alpha)cos^2(max[alpha,beta])

    dist - distance b/w donor and acceptor atoms
    sphyb1 - type of donor hybridization
    sphyb2 - type of acceptor hybridization
    alpha, beta, psi - parameters of hydrogen bonding
    """
    E = -33.47 * (5 * (2.8 / dist) ** 12 - 6 * (2.8 / dist) ** 10)
    if sphyb1 == 'SP3' and sphyb2 == 'SP3':
        F = ((math.cos(math.radians(alpha))) ** 2) * (math.cos(math.radians(beta - 109.5))) ** 2
    elif sphyb1 == 'SP3' and sphyb2 == 'SP2':
        F = ((math.cos(math.radians(alpha))) ** 2) * (math.cos(math.radians(beta))) ** 2
    elif sphyb1 == 'SP2' and sphyb2 == 'SP3':
        F = ((math.cos(math.radians(alpha))) ** 4)
    # elif sphyb1 == 'SP2' and sphyb2 == 'SP2':
    else:
        F = ((math.cos(math.radians(alpha))) ** 2) * (math.cos(math.radians(max([beta, psi])))) ** 2
    E_total = E * F
    return E_total


# max solvent surface area taken from Rose 1985 to calculate RSA values
area = {'ALA': 118.1, 'ARG': 256.0, 'ASN': 165.5, 'ASP': 158.7, 'CYS': 146.1, 'GLN': 193.2,
        'GLU': 186.2, 'GLY': 88.1, 'HIS': 202.5, 'ILE': 181.0, 'LEU': 193.1, 'LYS': 225.8,
        'MET': 203.4, 'PHE': 222.8, 'PRO': 146.8, 'SER': 129.8, 'THR': 152.5, 'TRP': 266.3,
        'TYR': 236.8, 'VAL': 164.5}


def calc_norm_vecs(res_com, origin):
    """
    @description
        Function used in calculations of normal vectors to ring planes for pipi, pication bonds
        Traverse atoms, find cross product, vector perpendicular to vectors a, b
    @input
        res_com - mda.selection of ring atoms of HIS, TYR, PHE, TRP
        origin - centroid, center of geometry

    """

    def normalize(pos):
        return [pos[i] / np.linalg.norm(pos[i]) for i in range(len(pos))]

    num = res_com.n_atoms
    vec = 0
    for i in range(0, num):
        a = normalize(res_com.positions[i::num] - res_com.positions[(i + 2) % num::num])
        b = normalize(res_com.positions[(i + 1) % num::num] - res_com.positions[(i + 2) % num::num])
        normalVec = np.cross(a, b)
        for j in range(len([origin])):
            if normalVec[j][0] < 0:
                normalVec[j] = normalVec[j] * -1
        vec += normalVec / num
    return vec


# amino acids variations
acids_set = {'HIS': ['HIS', 'HISD', 'HISE', 'HISH', 'HIS1'],
             'ARG': ['ARG', 'ARGN'],
             'LYS': ['LYS', 'LYSN'],
             'TYR': ['TYR'],
             'PHE': ['PHE'],
             'TRP': ['TRP']}

pication_acids = acids_set['HIS'] + acids_set['TYR'] + acids_set['PHE'] + acids_set['TRP']
cation_acids = acids_set['ARG'] + acids_set['LYS']
disulfide_acids = ['CYS', 'CYS2']
basic_salt_acids = ['ARG', 'ARGN', 'LYS', 'HIS', 'HISD', 'HISE', 'HISH', 'HIS1']
acidic_salt_acids = ['ASP', 'ASPH', 'GLU', 'GLUH']
# salt_acids = ['ARG', 'ARGN', 'LYS', 'LYSN', 'HIS', 'HISD', 'HISE', 'HISH', 'HIS1', 'ASP', 'ASPH', 'GLU', 'GLUH']

# dictionaries for donors and acceptors in hydrogen bonds
donor_dict = {
    'HIS': {'ND1': ['CG', 'CE1'], 'NE2': ['CD2', 'CE1']},
    'ARG': {'NE': ['CD', 'CZ'], 'NH1': ['HH11', 'CZ'], 'NH2': ['HH21', 'CZ']},
    'TRP': {'NE1': ['CD1', 'CE2']},
    'TYR': {'OH': ['CZ', 'CE2']},
    'ASN': {'ND2': ['HD21', 'CG']},
    'GLN': {'NE2': ['HE21', 'CD']},
    'HOH': {'O': ['H1', 'H2']},
    'SOL': {'OW': ['HW1', 'HW2']},
    'N': ['CA', 'H']
}

acceptor_dict = {
    'ASN': {'OD1': ['CG', 'CB']},
    'ASP': {'OD1': ['CG', 'OD2'], 'OD2': ['CG', 'OD1']},
    'GLN': {'OE1': ['CD', 'CG']},
    'GLU': {'OE1': ['CD', 'OE2'], 'OE2': ['CD', 'OE1']},
    'HIS': {'ND1': ['CG', 'CE1'], 'ND2': ['CD2', 'CE1']},
    'TYR': {'OH': ['CZ', 'CE1']},
    'HOH': {'O': ['H1', 'H2']},
    'SOL': {'OW': ['HW1', 'HW2']},
    'O': ['C', 'CA']
}


def normalDonorVecToPlane(donor, allatoms_data):
    """
    @description
        Calculate coordinates of normal vector to donor atom
    @input
        atom - donor atom in format HIS:26:A:CB
        allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    """
    res = ':'.join(donor.split(':')[0:3]) + ':'  # HIS:26:A:
    acid = donor.split(':')[0]  # HIS
    atom = donor.split(':')[3]  # CB
    if acid in donor_dict and atom in donor_dict[acid]:
        a1, a2 = donor_dict[acid][atom]
        a = allatoms_data[res + a1]['coords'] - allatoms_data[res + atom]['coords']
        b = allatoms_data[res + a2]['coords'] - allatoms_data[res + atom]['coords']
        normal = np.cross(a, b)
    elif atom == 'N':
        a1, a2 = donor_dict[atom]
        a = allatoms_data[res + a1]['coords'] - allatoms_data[res + atom]['coords']
        b = allatoms_data[res + a2]['coords'] - allatoms_data[res + atom]['coords']
        normal = np.cross(a, b)
    else:
        return None
    return normal


def normalAcceptorVecToPlane(acceptor, allatoms_data):
    """
    @description
      Calculate coordinates of normal vector to acceptor atom
    @input
      atom - donor atom in format HIS:26:A:CB
      allatoms_data - dict of atoms names, their positions and chain types (MC or SC)
    """
    res = ':'.join(acceptor.split(':')[0:3]) + ':'
    acid = acceptor.split(':')[0]
    atom = acceptor.split(':')[3]
    if acid in acceptor_dict and atom in acceptor_dict[acid]:
        a1, a2 = acceptor_dict[acid][atom]
        a = allatoms_data[res + a1]['coords'] - allatoms_data[res + atom]['coords']
        b = allatoms_data[res + a2]['coords'] - allatoms_data[res + atom]['coords']
        normal = np.cross(a, b)
    elif atom == 'O':
        a1, a2 = acceptor_dict[atom]
        a = allatoms_data[res + a1]['coords'] - allatoms_data[res + atom]['coords']
        b = allatoms_data[res + a2]['coords'] - allatoms_data[res + atom]['coords']
        normal = np.cross(a, b)
    else:
        return None
    return normal


# Define dictionaries for hydrogen bonds search
hydrogen = {
    "ARG:H": "N", "ARG:HE": "NE", "ARG:HH11": "NH1", "ARG:HH12": "NH1", "ARG:HH21": "NH2", "ARG:HH22": "NH2",
    "HIS:H": "N", "HIS:HE1": "NE2",
    "LYS:H": "N", "LYS:HZ1": "NZ", "LYS:HZ2": "NZ", "LYS:HZ3": "NZ",
    "ASP:H": "N",
    "GLU:H": "N",
    "SER:H": "N", "SER:HG": "OG",
    "THR:H": "N", "THR:HG1": "OG1",
    "ASN:H": "N", "ASN:HD21": "ND2", "ASN:HD22": "ND2",
    "GLN:H": "N", "GLN:HE21": "NE2", "GLN:HE22": "NE2",
    "CYS:H": "N", "CYS:HG": "SG",
    "SEC:H": "N",
    "GLY:H": "N",
    "ALA:H": "N",
    "VAL:H": "N",
    "ILE:H": "N",
    "LEU:H": "N",
    "MET:H": "N",
    "PHE:H": "N",
    "TYR:H": "N", "TYR:HH": "OH",
    "TRP:H": "N", "TRP:HE1": "NE1",
    "HOH:H1": "O", "HOH:H2": "O",
    'SOL:HW1': 'OW', 'SOL:HW2': 'OW'}

nextAcid = {
    "ARG:N": ["CA"], "ARG:NE": ["CD", "CZ"], "ARG:NH1": ["CZ"], "ARG:NH2": ["CZ"], "ARG:O": ["C"],
    "HIS:N": ["CA"], "HIS:NE2": ["CG", "CE1"], "HIS:ND1": ["CD2", "CE1"], "HIS:O": ["C"],
    "LYS:N": ["CA"], "LYS:NZ": ["CE"], "LYS:O": ["C"],
    "ASP:N": ["CA"], "ASP:OD1": ["CG"], "ASP:OD2": ["CG"], "ASP:O": ["C"],
    "GLU:N": ["CA"], "GLU:OE1": ["CD"], "GLU:OE2": ["CD"], "GLU:O": ["C"],
    "SER:N": ["CA"], "SER:OG": ["CB"], "SER:O": ["C"],
    "THR:N": ["CA"], "THR:OG1": ["CB"], "THR:O": ["C"],
    "ASN:N": ["CA"], "ASN:OD1": ["CG"], "ASN:ND2": ["CG"], "ASN:O": ["C"],
    "GLN:N": ["CA"], "GLN:OE1": ["CD"], "GLN:NE2": ["CD"], "GLN:O": ["C"],
    "CYS:N": ["CA"], "CYS:O": ["C"], "CYS:SG": ["CB"],
    "SEC:N": ["CA"], "SEC:O": ["C"],
    "GLY:N": ["CA"], "GLY:O": ["C"],
    "PRO:N": ["CA", "CB"], "PRO:O": ["C"],
    "ALA:N": ["CA"], "ALA:O": ["C"],
    "VAL:N": ["CA"], "VAL:O": ["C"],
    "ILE:N": ["CA"], "ILE:O": ["C"],
    "LEU:N": ["CA"], "LEU:O": ["C"],
    "MET:N": ["CA"], "MET:O": ["C"], "MET:SD": ["CG", "CE"],
    "PHE:N": ["CA"], "PHE:O": ["C"],
    "TYR:N": ["CA"], "TYR:OH": ["CZ"], "TYR:O": ["C"],
    "TRP:N": ["CA"], "TRP:NE1": ["CD1", "CE2"], "TRP:O": ["C"],
    "HOH:O": ["H1", "H2"],
    "SOL:OW": ['HW1', 'HW2'],
}


def isAcceptor(atom):
    """ Check whether atom is acceptor or not """
    acid = atom.split(':')[0] + ':' + atom.split(':')[3]
    SC_acceptors = ["ASN:OD1", "ASP:OD1", "ASP:OD2", "GLN:OE1", "GLU:OE1", "GLU:OE2", "HIS:ND1",
                    "HIS:NE2", "SER:OG", "THR:OG1", "TYR:OH", "CYS:SG", "HOH:O", 'SOL:OW']
    MC_acceptors = ["ARG:O", "HIS:O", "LYS:O", "ASP:O", "GLU:O", "SER:O", "THR:O", "ASN:O", "GLN:O", "CYS:O",
                    "SEC:O", "GLY:O", "PRO:O", "ALA:O", "VAL:O", "ILE:O", "LEU:O", "MET:O", "PHE:O", "TYR:O", "TRP:O"]
    if acid in SC_acceptors or acid in MC_acceptors:
        return True
    else:
        return False


def isDonor(atom):
    """ Check whether atom is donor or not """
    acid = atom.split(':')[0] + ':' + atom.split(':')[3]  # 3
    SC_donors = ["ARG:NE", "ARG:NH1", "ARG:NH2", "ASN:ND2", "GLN:NE2", "HIS:ND1", "HIS:NE2",
                 "LYS:NZ", "SER:OG", "THR:OG1", "TRP:NE1", "TYR:OH", "CYS:SG", "HOH:O", 'SOL:OW']
    MC_donors = ["ARG:N", "HIS:N", "LYS:N", "ASP:N", "GLU:N", "SER:N", "THR:N", "ASN:N", "GLN:N", "CYS:N",
                 "SEC:N", "GLY:N", "PRO:N", "ALA:N", "VAL:N", "ILE:N", "LEU:N", "MET:N", "PHE:N", "TYR:N", "TRP:N"]
    if acid in SC_donors or acid in MC_donors:
        return True
    else:
        return False


def SPHyb(elem):
    """ Return type of hybridization, SP3 or SP2 """
    acid = elem.split(':')[0]
    atom = elem.split(':')[3]
    if ((acid == "SER" and atom == "OG") or (acid == "THR" and atom == "OG1") or
            (acid == "CYS" and atom == "SG") or (acid == "LYS" and atom == "NZ") or
            (acid == 'SOL' and atom == 'OW') or (acid == "HOH" and atom == "O")):
        return "SP3"
    else:
        return "SP2"


# van der waals radii of atoms
radii = {'C': 1.76, 'CA': 1.87, 'CB': 1.87, 'CD': 1.815, 'CD1': 1.815,
         'CD2': 1.815, 'CE': 1.87, 'CE1': 1.76, 'CE2': 1.76, 'CE3': 1.76,
         'CG': 1.815, 'CG1': 1.87, 'CG2': 1.87, 'CH2': 1.76, 'CZ': 1.76,
         'CZ2': 1.76, 'CZ3': 1.76, 'D11': 1.8, 'D12': 1.8, 'D13': 1.8,
         'D21': 1.8, 'D22': 1.8, 'D23': 1.8, 'E21': 1.8, 'E22': 1.8,
         'G11': 1.8, 'G12': 1.8, 'G13': 1.8, 'G21': 1.8, 'G22': 1.8,
         'G23': 1.8, 'N': 1.65, 'ND1': 1.5, 'ND2': 1.65, 'NE': 1.65,
         'NE1': 1.65, 'NE2': 1.65, 'NH1': 1.65, 'NH2': 1.65, 'NZ': 1.5,
         'O': 1.4, 'OD1': 1.4, 'OD2': 1.4, 'OE1': 1.4, 'OE2': 1.4,
         'OG': 1.4, 'OG1': 1.4, 'OH': 1.4, 'SD': 1.85, 'SG': 1.85}

# distances b/w metals and amino were taken from TABLE.4
# https://sci-hub.se/https://doi.org/10.1107/S0907444906014594#
# For SER and THR, expect M—O distances between those for water and for monodentate carboxylate
# For TYR, expect M—O distances that are significantly shorter (by 0.1 A˚ ) than for monodentate carboxylate
# one problem why we add for all distances + 0.5A
MCOdist = {"NA": 2.88, "MG": 2.76, "K": 3.24, "CA": 2.86, "MN": 2.69, "FE": 2.54, "CO": 2.58, "CU": 2.54, "ZN": 2.57}
GLUASPmonoDist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 3.29, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
GLUASPbiDist = {"NA": 3, "MG": 2.7, "K": 3.4, "CA": 3.4, "MN": 2.7, "FE": 2.7, "CO": 2.7, "CU": 2.7, "ZN": 2.7}  # table 3
ASNdist = {"NA": 2.88, "MG": 2.76, "K": 3.24, "CA": 2.86, "MN": 2.69, "FE": 2.54, "CO": 2.58, "CU": 2.54, "ZN": 2.57}
GLNdist = {"NA": 2.88, "MG": 2.76, "K": 3.24, "CA": 2.86, "MN": 2.69, "FE": 2.54, "CO": 2.58, "CU": 2.54, "ZN": 2.57}
SERdist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 2.89, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
THRdist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 2.89, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
TYRdist = {"NA": 2.91, "MG": 2.57, "K": 3.31, "CA": 2.89, "MN": 2.69, "FE": 2.59, "CO": 2.59, "CU": 2.63, "ZN": 2.59}
HISdist = {"MN": 2.71, "FE": 2.66, "CO": 2.64, "CU": 2.52, "ZN": 2.53}
CYSdist = {"MN": 2.85, "FE": 2.8, "CO": 2.75, "CU": 2.65, "ZN": 2.81}
METdist = {"FE": 3, "CU": 3.3}
LYSdist = {"ZN": 2.8}

metals = ["NA", "MG", "K", "CA", "MN", "FE", "CO", "CU", "ZN", "S", "CL", 'NA+']
metal_dist_cutoff = 3.5

METALS_DICT = {
    'O': MCOdist,
    'GLU': [GLUASPmonoDist, ['OE1', 'OE2']],
    'ASP': [GLUASPbiDist, ['OD1', 'OD2']],
    'ASN': [ASNdist, 'OD1'],
    'GLN': [GLNdist, 'OE1'],
    'SER': [SERdist, 'OG'],
    'THR': [THRdist, 'OG1'],
    'TYR': [TYRdist, 'OH'],
    'HIS': [HISdist, ['NE2', 'ND1']],
    'CYS': [CYSdist, 'SG'],
    'MET': [METdist, 'SD'],
    'LYS': [LYSdist, 'NZ']
}

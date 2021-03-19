import numpy as np
import math

# metalls
metals = ["NA","MG","K","CA","MN","FE","CO","CU","ZN", "S"]

# constant values of angles and distances
SALT_CONST = '\t20\t' 
SALT_DISTANCE = 4

PIPI_CONST = '\t10\t'
PIPI_D1 = 4.4
PIPI_ANG1 = 30
PIPI_D2 = 5.5
PIPI_ANG2 = 60
PIPI_ANG3 = 120

PICATION_CONST = '\t9.4\t'
PICAT_DIST = 6.6 #change from PICAT_D
PICAT_ANG = 30

DISULF_CONST = '\t167\t'
DISULF_D = 3

VDW_D = .5
FULL_ANGLE = 360
STRAIGHT_ANGLE = 180
RIGHT_ANGLE = 90
ZERO_ANGLE = 0


def energy_vdw(rm, r): # energy formula for vdw bonds
    # rm - distance b/w 2 atoms
    # r- sum vdw radii of 2 atoms
    E = (-.997*((rm/r)**12-2*(rm/r)**6))*4
    return E

def energy_hb(d2, sphyb1, sphyb2, alpha, beta, psi): # energy of hybridization for hydrogen bonds
    '''
    d2 - distance b/w donor and acceptor
    sphyb1 - type of donor hybridization 
    sphyb2 - type of acceptor hybridization
    alpha, beta, psi - parameters of hydrogen bonding
    '''
    E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)
    if sphyb1=='SP3' and sphyb2=='SP3':
        F = ((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(beta-109.5)))**2
    elif sphyb1=='SP3' and sphyb2=='SP2':
        F = ((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(beta)))**2
    elif sphyb1=='SP2' and sphyb2=='SP3':
        F = ((math.cos(math.radians(alpha)))**4) 
    elif sphyb1=='SP2' and sphyb2=='SP2':    
        F = ((math.cos(math.radians(alpha)))**2)*(math.cos(math.radians(max([beta,psi]))))**2                    
    E_total = E*F       
    return E_total

# max solvent surface area taken from Rose 1985 to calculate RSA values
area = {'ALA': 118.1, 'ARG': 256.0, 'ASN': 165.5, 'ASP': 158.7, 'CYS': 146.1, 'GLN': 193.2,
        'GLU': 186.2, 'GLY': 88.1, 'HIS': 202.5, 'ILE': 181.0, 'LEU': 193.1, 'LYS': 225.8,
        'MET': 203.4, 'PHE': 222.8, 'PRO': 146.8, 'SER': 129.8, 'THR': 152.5, 'TRP': 266.3,
        'TYR': 236.8, 'VAL': 164.5}

# hydrogen constants
# class hydrogen_vars:
#     d1 = 2.5; d11 = 3.2; d12 = 3.3
#     d2 = 3.9; d21 = 4.1; d22 = 4.5
#     def energy(d2):
#         E = -33.47*(5*(2.8/d2)**12 - 6*(2.8/d2)**10)
#         return E

def calc_norm_vecs(res_com, origin):
    '''Function used in calc origin, normals for pipi, pication
    acid - mda.residue object
    res_com - mda.selection
    origin - center of geometry
    
    '''
    num = res_com.n_atoms
    vec = 0
    for i in range(0, num):
        a = normalize(res_com.positions[i::num] - res_com.positions[(i+2)%num::num])
        b = normalize(res_com.positions[(i+1)%num::num] - res_com.positions[(i+2)%num::num])
        normalVec = np.cross(a, b)
        for j in range(len([origin])):
            if normalVec[j][0]<0:
                normalVec[j] = normalVec[j]*-1
        vec+=normalVec/num
    return vec


def isNeighbor(atom1,atom2):
  # figure out whether 2 atoms are neighbors if aminoacid ids differ by 1
    if len(atom1.split(':')) == 4 :
        # index1 = -int(atom1.split(':')[2][0:-1])
        index1 = -int(atom1.split(':')[1])
    else: index1 = atom1.split(':')[1][0:-1]
    if len(atom2.split(':')) == 4:
        # index2 = -int(atom2.split(':')[2][0:-1])
        index2 = -int(atom2.split(':')[1])
    else:
        index2 = atom2.split(':')[1][0:-1]
    if abs(int(index1) - int(index2)) == 1:
        return True
    else: return False
    
# aminoacids variations
# normVec4acids = {} # maybe create dict of acids and atoms for normalVecPIPI
his = ['HIS', 'HISD', 'HISE', 'HISH', 'HIS1']
salt_aa = ['ARG', 'ARGN', 'LYS', 'LYSN', 'HIS', 'HISD', 'HISE', 'HISH', 'HIS1', 'ASP', 'ASPH', 'GLU', 'GLUH']
pication_aa = his + ['TYR', 'TRP', 'PHE']
cation_aa = ['ARG', 'ARGN', 'LYS', 'LYSN']
disulf_aa = ['CYS', 'CYS2']
basic_salt_aa = ['ARG', 'ARGN', 'LYS', 'HIS', 'HISD', 'HISE', 'HISH', 'HIS1']
acidic_salt_aa = ['ASP', 'ASPH', 'GLU', 'GLUH']

def normalize(a):
    return [a[i] / np.linalg.norm(a[i]) for i in range(len(a))]

# HYDROGEN BONDS
donor_dict = {'HIS' : {'ND1': ['CG', 'CE1'], 
                    'NE2': ['CD2', 'CE1']},
            'ARG' : {'NE' : ['CD', 'CZ'],
                        'NH1': ['HH11', 'CZ'],
                        'NH2': ['HH21', 'CZ']},
            'TRP' : {'NE1': ['CD1', 'CE2']}, 
            'TYR' : {'OH' : ['CZ', 'CE2']},
            'ASN' : {'ND2': ['HD21', 'CG']}, 
            'GLN' : {'NE2': ['HE21', 'CD']}, 
            'HOH' : {'O'  : ['H1', 'H2']}, 
            'N' : ['CA', 'H']}

acceptor_dict = {'ASN' : {'OD1' : ['CG', 'CB']}, 
                 'ASP' : {'OD1' : ['CG', 'OD2'],
                          'OD2' : ['CG', 'OD1']},
                 'GLN' : {'OE1' : ['CD', 'CG']}, 
                 'GLU' : {'OE1' : ['CD', 'OE2'],
                          'OE2' : ['CD', 'OE1']},
                 'HIS' : {'ND1' : ['CG', 'CE1'], 
                          'ND2' : ['CD2', 'CE1']},
                 'TYR' : {'OH' : ['CZ', 'CE1']},
                 'HOH' : {'O' : ['H1', 'H2']}, 
                 'O' : ['C', 'CA']}


def normalDonorVecToPlane1(A, coords):
    '''
    calculate coords of normal to donor acid
    Take as input name of Donor (A) and coord dictionary
    '''
    res = ':'.join(A.split(':')[0:3]) + ':'
    r = A.split(':')[0]
    atom = A.split(':')[3]
    if (r in donor_dict and atom in donor_dict[r]):
        a1, a2 =donor_dict[r][atom]
        a = coords[res+a1] - coords[res+atom]
        b = coords[res+a2] - coords[res+atom]
        normal = np.cross(a, b)
    elif atom == 'N':
        a1, a2 = donor_dict[atom]
        a = coords[res+a1] - coords[res+atom]
        b = coords[res+a2] - coords[res+atom]
        normal = np.cross(a, b)
    else:
        return None
    return normal

def normalAcceptorVecToPlane1(A, coords):
    '''
    calculate coords of normal to Acceptor acid
    Take as input name of Acceptor (A) and coord dictionary
    '''
    res = ':'.join(A.split(':')[0:3]) + ':'
    r = A.split(':')[0]
    atom = A.split(':')[3]
    if (r in acceptor_dict and atom in acceptor_dict[r]):
        a1, a2 = acceptor_dict[r][atom]
        a = coords[res+a1] - coords[res+atom]
        b = coords[res+a2] - coords[res+atom]
        normal = np.cross(a, b)
    elif atom == 'O':
        a1, a2 = acceptor_dict[atom]
        a = coords[res+a1] - coords[res+atom]
        b = coords[res+a2] - coords[res+atom]
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
"HOH:H1": "O", "HOH:H2": "O"}

nextAcid = {
"ARG:N": ["CA"],"ARG:NE": ["CD","CZ"], "ARG:NH1": ["CZ"], "ARG:NH2": ["CZ"], "ARG:O": ["C"],
"HIS:N": ["CA"], "HIS:NE2": ["CG","CE1"], "HIS:ND1": ["CD2","CE1"], "HIS:O": ["C"],
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
"PRO:N": ["CA","CB"],"PRO:O": ["C"],
"ALA:N": ["CA"], "ALA:O": ["C"],
"VAL:N": ["CA"], "VAL:O": ["C"],
"ILE:N": ["CA"], "ILE:O": ["C"],
"LEU:N": ["CA"], "LEU:O": ["C"],
"MET:N": ["CA"], "MET:O": ["C"], "MET:SD": ["CG","CE"],
"PHE:N": ["CA"], "PHE:O": ["C"],
"TYR:N": ["CA"], "TYR:OH": ["CZ"], "TYR:O": ["C"],
"TRP:N": ["CA"], "TRP:NE1": ["CD1","CE2"], "TRP:O": ["C"],
"HOH:O": ["H1","H2"]
}

# check if acid+atom is Acceptor 
def isAcceptor(elem):
    acid = elem.split(':')[0]
    atom = elem.split(':')[3]
    SCacceptors = ["ASN:OD1","ASP:OD1","ASP:OD2","GLN:OE1","GLU:OE1","GLU:OE2","HIS:ND1","HIS:NE2","SER:OG","THR:OG1","TYR:OH","CYS:SG","HOH:O"]
    MCacceptors = ["ARG:O","HIS:O","LYS:O","ASP:O","GLU:O","SER:O","THR:O","ASN:O","GLN:O","CYS:O","SEC:O","GLY:O","PRO:O","ALA:O","VAL:O","ILE:O","LEU:O","MET:O","PHE:O","TYR:O","TRP:O"]
    # if atom.split(':')[0]+':'+atom.split(':')[2] in SCacceptors or atom.split(':')[0]+':'+atom.split(':')[2] in MCacceptors:
    #     return(True)
    # else:
    #     return(False)
    if acid+':'+atom in SCacceptors or acid+':'+atom in MCacceptors:
        return(True)
    else:
        return(False)


# check if acid+atom is Donor
def isDonor(elem):
    acid = elem.split(':')[0]
    atom = elem.split(':')[3]
    SCdonors = ["ARG:NE","ARG:NH1","ARG:NH2","ASN:ND2","GLN:NE2","HIS:ND1","HIS:NE2","LYS:NZ","SER:OG","THR:OG1","TRP:NE1","TYR:OH","CYS:SG","HOH:O"]
    MCdonors = ["ARG:N","HIS:N","LYS:N","ASP:N","GLU:N","SER:N","THR:N","ASN:N","GLN:N","CYS:N","SEC:N","GLY:N","PRO:N","ALA:N","VAL:N","ILE:N","LEU:N","MET:N","PHE:N","TYR:N","TRP:N"]
    # if atom.split(':')[0]+':'+atom.split(':')[2] in SCdonors or atom.split(':')[0]+':'+atom.split(':')[2] in MCdonors:
    #     return(True)
    # else:
    #     return(False)
    if acid+':'+atom in SCdonors or acid+':'+atom in MCdonors:
        return(True)
    else:
        return(False)



# return type of hybridization
def SPHyb(elem):
    acid = elem.split(':')[0]
    atom = elem.split(':')[3]
    if ((acid=="SER" and atom=="OG") or 
        (acid=="THR" and atom=="OG1") or 
        (acid=="CYS" and atom=="SG") or 
        (acid=="LYS" and atom=="NZ") or 
        (acid=="HOH" and atom=="O")):
        return("SP3")
    else:
        return("SP2")

# van der waals radii of atoms
radii= {'C': 1.76, 'CA': 1.87, 'CB': 1.87, 'CD': 1.815, 'CD1': 1.815,
       'CD2': 1.815, 'CE': 1.87, 'CE1': 1.76, 'CE2': 1.76, 'CE3': 1.76,
       'CG': 1.815, 'CG1': 1.87, 'CG2': 1.87, 'CH2': 1.76, 'CZ': 1.76,
       'CZ2': 1.76, 'CZ3': 1.76, 'D11': 1.8, 'D12': 1.8, 'D13': 1.8,
       'D21': 1.8, 'D22': 1.8, 'D23': 1.8, 'E21': 1.8, 'E22': 1.8,
       'G11': 1.8, 'G12': 1.8, 'G13': 1.8, 'G21': 1.8, 'G22': 1.8,
       'G23': 1.8, 'N': 1.65, 'ND1': 1.5, 'ND2': 1.65, 'NE': 1.65,
       'NE1': 1.65, 'NE2': 1.65, 'NH1': 1.65, 'NH2': 1.65, 'NZ': 1.5,
       'O': 1.4, 'OD1': 1.4, 'OD2': 1.4, 'OE1': 1.4, 'OE2': 1.4,
       'OG': 1.4, 'OG1': 1.4, 'OH': 1.4, 'SD': 1.85, 'SG': 1.85}

# distances b/w metals and aminos
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

METALLS_DICT = {
                'O' : MCOdist,
                'GLU' : [GLUASPmonoDist, ['OE1', 'OE2']],
                'ASP' : [GLUASPbiDist, ['OD1', 'OD2']], 
                'ASN' : [ASNdist, 'OD1'],
                'GLN' : [GLNdist, 'OE1'], 
                'SER' : [SERdist, 'OG'],
                'THR' : [THRdist, 'OG1'], 
                'TYR' : [TYRdist, 'OH'], 
                'HIS' : [HISdist, ['NE2', 'ND1']], 
                'CYS' : [CYSdist, 'SG'],
                'MET' : [METdist, 'SD'], 
                'LYS' : [LYSdist, 'NZ']
                }

# replace atom 
def rplcAtom(atom):
    if atom == 'NH1' or atom == 'NH2':   item = 'NH1/2'
    elif atom == 'OD1' or atom == 'OD2': item = 'OD1/2'
    elif atom == 'OE1' or atom == 'OE2': item = 'OE1/2'
    else: item = atom
    return item

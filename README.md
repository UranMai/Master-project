# Master-project
    Before using this code preprocess pdb file
    1. Add hydrogens (phenix.reduce)
    2. Protonate water (phenix.ready_set)
    3. And create phi file (stride)
    Then run this script >>> python pdb2edge_class.py PDB.pdb




> pdb2edge_class.py \
# Functions documentation

### SecondaryStructure Function
- **Description**\
Read PDB and Phi file \
Extract phi angles from phi file for each acid on conditions
- **Parameteres**\
Input PDB file 
- **Returns**\
''pdb_secondaryStructure'' file in format\
''pdb_secondaryStructure2'' file in format
  acid + resid + segid | secStructure | phi_angle
  ------------- | ------------- | ------------- 
  HIS26A	| CoilA1	| 0
  PRO27A	| aHelix2	| -57.19
  GLU28A	| aHelix2 | -65.68 

### RSA Function
- **Description**\
Read Phi file \
Extract phi areas from phi file and divided by area of aminoacids
- **Parameteres**\
Input PDB file 
area of acid variable from params1.py
- **Returns**\
"rsa" file. 
  acid + resid + segid | area 
  ------------- | -------------
  HIS26A | 0.4740740740740741
  PRO27A | 0.5933242506811989

### Aminoacid class (AA class)
- **Desciption**\
`def __init__(self, aminoacid)`\
Define bonds parameteres\
`def Calc_ORIGIN_NORM(self, acid)`\
return coordinates dict of norm vectors and origin
### `def PiPi(self, res, origin_dict, vec_dict)`

    - description
    return name of acid and coords of acids involved in PIPI interactions (HIS, TYR, PHE, TRP)
    check if acid is right 
    - parameteres
    res - acid name from __init__
    origin_dict, vec_dict - coords from def Calc_ORIGIN_NORM()
    
### `def PiCation(self, res, origin_dict, vec_dict)` 

    - description
    return name of acid and coords of acids involved in PICation interactions (HIS, TYR, PHE, TRP)
    and cation atoms NZ, CZ of acids ARG and LYS
    check if acid is right 
    - parameteres
    res - acid name from __init__
    origin_dict, vec_dict - coords from def Calc_ORIGIN_NORM()
    
### `def Disulf(self, res)`

    - description
    check if residues is in disulfide bonds involved
    - return 
    coordinates of atoms SG of CYS
    
### `def Salt_Bridge(self, res)`

    - description
    Define basic (ARG, LYS, HIS) and acidic (ASP, GLU) acids
    - return
    atoms and coords of these acids
    
- **Parameteres**\
mda.residues

### `def PIPI(res1, res2)`
- description
check if pipi exist and call pipi parametre of class residue\
iterate through pipi acids\
find angle between norm vectores of pipi and distance between origins
if angle and distance are under conditions write it to array
- parameteres 
residue1 class and residue2 class
- return
return array of acids that get along with conditions

### `def PiPi_writing(res1, res2)`
- description
take data from `def PIPI(res1, res2)` function and write it
- paramteres
input pdb file 
- return
write pdb_pipi2 file

### `PiCation(res1, res2)`
- description
check and call options of pication class\
find angle and coords of vector between cation and ring\
if angle and distance are under conditions write it to array
- parameteres
residue1 class and residue2 class
- return 
return array of acids that get along with conditions

### `def PiCation_writing(res1, res2)`
- description
take data from `def PiCation(res1, res2)` function and write it
- paramteres
input pdb file 
- return
write pdb_pication2 file

### `SaltBridge(res1, res2)`
- description
    - take coords of involved acids\
    - find pairs and distances via mda.capped_distance
    - if no distance return None
    - else iterate through pairs and condition
- parameteres
residue1 class and residue2 class
- return
return list of data and saltBridges

### `SB_writing(file1)`
- description
take data from `def SaltBridge(res1, res2)` function and write it
- parameteres
input pdb file
- return 
pdb_disulf file
- return 
write pdb_salrBridges_Barlow file\
list of saltBridges that will be used in 

### `def Disulfide(res1, res2)`
- description
    - Delete the check on res1==res2\
    - Take coords of SG atoms of CYS acid
    - find distance between SG atoms of 2 residues
    - if distance less than CONST_distance
- parameters
residue1 class and residue2 class
- return 
res1+res2

### `def Disulfide_writing(file1)`
- descirption
take data from `def Disulfide(res1, res2)` function and write it
- parameters
input pdb file
- return 
pdb_disulf file

### `def Hydrogen_calculation()`
- desciption
    - take saltBridges from function
    - iterate through created hydrogen table h.table, which contains donor and acceptor
    - define parameteres of hydrogen bonding and write it to hydrogenBonds list
    - write pdb_hb file
-return 

### `def VdW_calculation()`
- description
    - Select atoms in mda.Universe
    - capped_distance
    - Iterate through atoms, find distance between atoms and write it
- return 
pdb_vdw file and pdb_vdw2 file with own condtions

### `def Calc_MetalBonds(met, AA1)`
- desciption
Iterate through AA1.acid atoms\
Conditions 
- parameteres
met - metal atom\
AA1 - class of aminoacids
- return
metal-atom bonds

### `def write_metalbonds(met, metal2atom)`
- description

- parameteres
met - metal atom\
metal2atom - output of function `def Calc_MetalBonds`

### `def MetalBonds_calculation()`
- return 
pdb_metal file 

### `def ligand_bonds()`


### `def Centroid_calculation`

### `def pdb2peptide`

### `def Bfactor()`
- desciption



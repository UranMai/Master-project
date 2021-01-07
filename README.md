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
`__init__(self, aminoacid)`\
Define bonds parameteres\
`Calc_ORIGIN_NORM(self, acid)`\
return coordinates dict of norm vectors and origin
`def PiPi(self, res, origin_dict, vec_dict)`
    - description
    return name of acid and coords of acids involved in PIPI interactions (HIS, TYR, PHE, TRP)
    - parameteres
    
    
- **Parameteres**\
mda.residues

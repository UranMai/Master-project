### Run
```sh
python cluster_frames.py PDB XTC
```

This code allow to center protein removing artifacts from periodic boundary conditions and cluster trajectory frames by RMSD value.\
For each cluster centroid run ```pdb_class.py``` to find aminoacid interactions in a protein

### **Input**

PDB and XTC are topology and trajectory files used from
1. MDAnalysisTests or
2. ```run_1.tpr``` and ```run_1.xtc``` from google.drive after MD production

\
Now it runs on Test_data/1BTL.pdb, creating ```'1BTL_net'``` in format ```[LYS:32:A:NZ	ASP:35:A:OD2	SCSC	20	SB]```\
and ligand files
```sh
python pdb_class.py
```

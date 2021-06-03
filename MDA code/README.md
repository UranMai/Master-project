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
python pdb_class.py ../Test_data/1btl_clean.pdb
```

# Run
Consider topology (**run_1.tpr**) and trajectory files (**run_1.xtc**). Scripts process: clusters centroid frames, each k-th frame and all frames.
After run main function from ```pdb_class.py``` on each frame and compute protein bonds
```sh
python cluster_frames ../Test_data/1btl_clean.pdb ../Test_data/run_1.tpr ../Test_data/run_1.xtc
```
where ```1btl_clean.pdb``` - pdb file after **PDBFixer**, ```run_1.tpr``` - topology file, ```run_1.xtc``` - trajectory file after Gromacs MD production on 1byl_clean.pdb

To run energetics, centroid scripts take 2 inputs: PDB file - ```../Test_data/1btl_clean.pdb```; and ```../Test_data/1btl_500.pdb```
finalscore.py takes input: ```../Test_data/1btl_500.pdb```

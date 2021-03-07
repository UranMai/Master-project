https://github.com/UranMai/Master-project/wiki

### Preprocessing PDB file
1. Install **Stride** from [here](http://webclu.bio.wzw.tum.de/stride/install.html) and copy binary file stride to `usr/local/bin`
2. Install **Phenix** from [here](https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/phenix/install-setup-run.html)
3. Run ```bash preprocess.sh PDBfile```


[phenix.reduce](https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/phenix/reference/hydrogens.html)

### Run code
```
pip install MDAnalysis
pip install igraph
```

1. Find protein interactions in pdb file, (ex. 1BTL.pdb)
```
python pdb_class.py file.pdb
```
2. Calculate energetics and centroid scores
```
python energetics.py file.pdb
python centroid.py file.pdb
```

Before using code, pdb file need in preprocessing.
* Add hydrogens using **phenix.reduce**
* Protonate water using **phenix.ready_set**
* Make phi file using **Stride**

| Total time | pdb2edge.py | MDA |
| ---------- | ----------- | --- |
| 6vxx.pdb | ~3200s | ~600s | 
| 1BTL.pdb | ~40s | ~16s |
| 4eiy.pdb | ~100s | ~30s |

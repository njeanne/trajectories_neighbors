# Molecular Dynamics Trajectory atoms neighboring analysis

From molecular dynamics trajectory files (`*.nc`), the script performs a trajectory analysis to search contacts. It 
looks for the atoms neighboring of two different residues. 

The neighbors' atoms are written to a CSV file with the residues they belong to.

## Conda environment

A [conda](https://docs.conda.io/projects/conda/en/latest/index.html) YAML environment file is provided: 
`conda_env/contacts_env.yml`. The file contains all the dependencies to run the script.

The conda environment is generated using the command:
```shell script
# create the environment
conda env create -f conda_env/contacts_env.yml

# activate the environment
conda activate contacts
```

## Usage

The script can be tested with the test data provided in the `data` directory, which contains the trajectory 
file `test_data_20-frames.nc` and the topology associated file `test_data.parm`. The commands are:

```shell script
conda activate contacts

# run with 4 processors in parallel
mpirun -np 4 python ./trajectories_neighbors.py --out results/test --sample test --topology data/test_data.parm \
--frames test_data_20-frames.nc:5-20 --nanoseconds 1 --distance-contacts 3.0 --proportion-contacts 50.0 --roi 7 12 \
--padding 0 data/test_data_20-frames.nc

conda deactivate
```

The optional parameters used are:
- `--frames test_data_20-frames.nc:5-20`: frames selection from 5 to 20 from the file `test_data_20-frames.nc`.
- `--proportion-contacts 50.0`: a contact is validated only if it is at least present in 50% of the frames 500 to 2000.
- `--distance-contacts 3.0`: maximal distance in Angstroms between 2 atoms of different residues.
- `--roi 7 12`: the region of interest by residues index (1-index.)
- `--padding 0`: the number of residues to ignore on each side of the region of interest.
- `--nanoseconds`: the molecular dynamics simulation duration.

If the analysis resumes a previous analysis, use the `--resume` parameter with the path of the YAML file of the 
previous analysis as argument:

```shell script
# run with 4 processors in parallel
mpirun -np 4 python ./trajectories_neighbors.py --out results/test --sample test \
--resume tests/expected/analysis_resumed.yaml --topology data/test_data.parm --frames test_data_20-frames.nc:5-20 \
--nanoseconds 1 --distance-contacts 3.0 --proportion-contacts 50.0 --roi 7 12 --padding 0 data/test_data_20-frames.nc
```

## Outputs

The script outputs is a CSV file of the validated contacts:

|neighbors          |residue 1 position|residue 1|atom 1|residue 2 position|residue 2|atom 2|frames with contact|total frames|proportion frames (%)|
|-------------------|------------------|---------|------|------------------|---------|------|-------------------|------------|---------------------|
|PRO7_N-ALA6_N      |7                 |PRO      |N     |6                 |ALA      |N     |17                 |20          |85.0                 |
|PRO7_N-ALA6_CA     |7                 |PRO      |N     |6                 |ALA      |CA    |20                 |20          |100.0                |
|PRO7_N-ALA6_C      |7                 |PRO      |N     |6                 |ALA      |C     |20                 |20          |100.0                |
|PRO7_N-ALA6_O      |7                 |PRO      |N     |6                 |ALA      |O     |20                 |20          |100.0                |
|PRO7_CD-ALA6_CA    |7                 |PRO      |CD    |6                 |ALA      |CA    |20                 |20          |100.0                |
|PRO7_CD-ALA6_HA    |7                 |PRO      |CD    |6                 |ALA      |HA    |20                 |20          |100.0                |
|PRO7_CD-ALA6_C     |7                 |PRO      |CD    |6                 |ALA      |C     |20                 |20          |100.0                |
|PRO7_CD-ALA6_O     |7                 |PRO      |CD    |6                 |ALA      |O     |20                 |20          |100.0                |
|PRO7_HD2-ALA6_CA   |7                 |PRO      |HD2   |6                 |ALA      |CA    |14                 |20          |70.0                 |
|PRO7_HD2-ALA6_HA   |7                 |PRO      |HD2   |6                 |ALA      |HA    |20                 |20          |100.0                |
|PRO7_HD2-ALA6_C    |7                 |PRO      |HD2   |6                 |ALA      |C     |20                 |20          |100.0                |
|PRO7_HD3-ALA6_CA   |7                 |PRO      |HD3   |6                 |ALA      |CA    |11                 |20          |55.0                 |
|PRO7_HD3-ALA6_HA   |7                 |PRO      |HD3   |6                 |ALA      |HA    |19                 |20          |95.0                 |
|PRO7_HD3-ALA6_HB1  |7                 |PRO      |HD3   |6                 |ALA      |HB1   |11                 |20          |55.0                 |
|PRO7_HD3-ALA6_C    |7                 |PRO      |HD3   |6                 |ALA      |C     |19                 |20          |95.0                 |
|PRO7_CG-ALA6_HA    |7                 |PRO      |CG    |6                 |ALA      |HA    |20                 |20          |100.0                |
|PRO7_CG-ALA6_C     |7                 |PRO      |CG    |6                 |ALA      |C     |18                 |20          |90.0                 |
|PRO7_HA-ALA6_C     |7                 |PRO      |HA    |6                 |ALA      |C     |20                 |20          |100.0                |
|PRO7_HA-ALA6_O     |7                 |PRO      |HA    |6                 |ALA      |O     |19                 |20          |95.0                 |
|PRO7_C-ALA6_C      |7                 |PRO      |C     |6                 |ALA      |C     |20                 |20          |100.0                |
|PRO7_C-ALA6_O      |7                 |PRO      |C     |6                 |ALA      |O     |19                 |20          |95.0                 |
|ALA12_HA-PRO13_N   |12                |ALA      |HA    |13                |PRO      |N     |20                 |20          |100.0                |
|ALA12_HA-PRO13_CD  |12                |ALA      |HA    |13                |PRO      |CD    |10                 |20          |50.0                 |
|ALA12_HA-PRO13_HD3 |12                |ALA      |HA    |13                |PRO      |HD3   |13                 |20          |65.0                 |
|ALA12_CB-PRO13_N   |12                |ALA      |CB    |13                |PRO      |N     |20                 |20          |100.0                |
|ALA12_CB-PRO13_CD  |12                |ALA      |CB    |13                |PRO      |CD    |20                 |20          |100.0                |
|ALA12_CB-PRO13_HD2 |12                |ALA      |CB    |13                |PRO      |HD2   |17                 |20          |85.0                 |
|ALA12_CB-PRO13_HD3 |12                |ALA      |CB    |13                |PRO      |HD3   |20                 |20          |100.0                |
|ALA12_HB2-PRO13_HD2|12                |ALA      |HB2   |13                |PRO      |HD2   |16                 |20          |80.0                 |
|ALA12_O-PRO13_N    |12                |ALA      |O     |13                |PRO      |N     |20                 |20          |100.0                |
|ALA12_O-PRO13_CD   |12                |ALA      |O     |13                |PRO      |CD    |20                 |20          |100.0                |
|ALA12_O-PRO13_HD2  |12                |ALA      |O     |13                |PRO      |HD2   |16                 |20          |80.0                 |
|ALA12_O-PRO13_HD3  |12                |ALA      |O     |13                |PRO      |HD3   |19                 |20          |95.0                 |
|ALA12_O-PRO13_CA   |12                |ALA      |O     |13                |PRO      |CA    |20                 |20          |100.0                |
|ALA12_O-PRO13_HA   |12                |ALA      |O     |13                |PRO      |HA    |20                 |20          |100.0                |

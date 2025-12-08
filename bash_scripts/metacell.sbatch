#!/bin/bash
#SBATCH --job-name="metacell_fast2b_2x_scANVI_basis"
#SBATCH --output="metacell_fast2b_2x_scANVI_basis.log"
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -t 01:00:00

# Set working directory where the python script is located
WORK_DIR="/ocean/projects/cis240075p/asachan/datasets/TA_muscle/ERCC1_KO_mice/integrated_samples"
cd $WORK_DIR || { echo "Error: Could not change directory to $WORK_DIR"; exit 1; }

# Load necessary modules and activate the environment
# Load GCC first to provide newer libstdc++ (GLIBCXX_3.4.29)
module load gcc/13.3.1-p20240614

module load anaconda3/2022.10
eval "$(conda shell.bash hook)"  # Proper conda initialization
conda activate metasheller-py311

# Define file paths as bash variables
adata_file = "/ocean/projects/cis240075p/asachan/datasets/TA_muscle/ERCC1_KO_mice/integrated_samples/seacells_F_FastIIB_FastIIX.h5ad"

# Run metacell
echo "Running metacell for batch fast2b_2x_scANVI_basis"

########################## Python script ##########################
python -c "
import scanpy as sc 
import metashells as ms
from metashells.seashells import SeaShells 

import os 

results_dir = './meta_cell_outs_fast2b_2x_scANVI_basis'
os.makedirs(results_dir, exist_ok=True)

adata = sc.read_h5ad('/ocean/projects/cis240075p/asachan/datasets/TA_muscle/ERCC1_KO_mice/integrated_samples/adata_harmony_v2.h5ad')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

conch = SeaShells(
    adata_full = adata, 
    results_dir = results_dir,
    cells_per_metacell = 75, # recommended is 75 but can go as low as ~20 if really needed
    sample_col = 'sample_id'
)
conch.run()
"

# Log completion of the batch
echo "Done with metacell for batch fast2b_2x_scANVI_basis"
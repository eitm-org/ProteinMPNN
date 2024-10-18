import os
import glob
import json
from tqdm import tqdm

def map_motif_sequence(motif_pdbs_dir, pdbs_dir, output_dir, verbose=False):
    """
    Map motif sequence information into PDB files of generated structures 
    in preparation for later conditional inverse folding.

    Args:
        motif_pdbs_dir:
            Directory containing motif structures, where each PDB file (corresponding to 
            the same filename in the pdbs directory) contains the motif structure, aligned 
            in residue indices with the generated structure.
        pdbs_dir:
            Directory containing generated structures in the PDB format.
        output_dir:
            Base output directory.
        verbose:
            Whether to print detailed progress information.

    Returns:
        processed_pdb_dir:
            Output directory (specified as [output_dir]/processed_pdbs), where each file 
            contains the generated structure in the PDB format, with mapped motif sequence 
            information.
    """

    # Process
    for pdb_filepath in tqdm(
        glob.glob(os.path.join(pdbs_dir, '*.pdb')),
        desc='Mapping motif sequence', disable=not verbose
    ):

        # Parse
        domain_name = pdb_filepath.split('/')[-1].split('.')[0]

        # Create residue index to name mapping
        motif_pdb_filepath = os.path.join(motif_pdbs_dir, f'{domain_name}.pdb')
        with open(motif_pdb_filepath) as file:
            residue_name_dict = dict([
                (int(line[22:26]), line[17:20]) for line in file
                if line.startswith('ATOM') and line[12:16].strip() == 'CA'
            ])

        # Update
        lines = []
        with open(pdb_filepath) as file:
            for line in file:
                assert line.startswith('ATOM') and line[21] == 'A'
                residue_index = int(line[22:26])
                residue_name = line[17:20]
                if residue_index in residue_name_dict:
                    residue_name = residue_name_dict[residue_index]
                lines.append(line[:17] + residue_name + line[20:])

        # Save
        processed_pdb_filepath = os.path.join(output_dir, f'{domain_name}.pdb')
        with open(processed_pdb_filepath, 'w') as file:
            file.write(''.join(lines))

def create_fixed_positions_dict(motif_pdbs_dir, processed_pdbs_dir, output_path, verbose=False):
    """
    Run conditional inverse folding to obtain sequences.

    Args:
        motif_pdbs_dir:
            Directory containing motif structures, where each PDB file (corresponding to 
            the same filename in the pdbs directory) contains the motif structure, aligned 
            in residue indices with the generated structure.
        processed_pdb_dir:
            Directory containing processed PDB files, where each file contains the generated 
            structure in the PDB format, with mapped motif sequence information.
        output_dir:
            Base output directory.
        verbose:
            Whether to print detailed progress information.

    Returns:
        fixed_positions_dict:
            Output file (specified as [output_dir]/fixed_positions.jsonl) of fixed residue indices to 
            pass into proteinMPNN for motif scaffolding sequence prediction.
    """
    # Process
    full_fixed_positions_dict = {}
    for processed_pdb_filepath in tqdm(
        glob.glob(os.path.join(processed_pdbs_dir, '*.pdb')),
        desc='Inverse folding', disable=not verbose
    ):
        domain_name = processed_pdb_filepath.split('/')[-1].split('.')[0]
        # Create fixed positions dictionary
        with open(os.path.join(motif_pdbs_dir, f'{domain_name}.pdb')) as file:
            fixed_residue_indices = [
                int(line[22:26]) for line in file 
                if line.startswith('ATOM') and line[12:16].strip() == 'CA'
            ]
        fixed_positions_dict = {}
        fixed_positions_dict[domain_name] = {
            'A': fixed_residue_indices
        }
        full_fixed_positions_dict = full_fixed_positions_dict | fixed_positions_dict
        # Predict
    with open(output_path,  'w') as file:
        file.write(json.dumps(full_fixed_positions_dict) + '\n')

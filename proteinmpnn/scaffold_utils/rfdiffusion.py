import os
import glob
import json
from tqdm import tqdm


def save_motif_pdbs(pdbs_dir, verbose=False):
    """
    Parse rfdiffusion generated pdbs and save pdbs of only the motif into a new directroy called motif_pdbs
    """
    for pdb_filepath in tqdm(
        glob.glob(os.path.join(pdbs_dir, '*.pdb')),
        desc='Save motif pdbs', disable=not verbose
    ):
        domain_name = pdb_filepath.split('/')[-1].split('.')[0]
        lines = []
        with open(pdb_filepath) as file:
            for line in file:
                if line.startswith('ATOM') and line[60:66].strip() == '1.00':
                    lines.append(line)
        # Save
        output_dir = os.path.join(pdbs_dir, 'motif_pdbs')
        if not os.path.isdir(output_dir): os.makedirs(output_dir)
        motif_pdb_filepath = os.path.join(output_dir, f'{domain_name}.pdb')
        with open(motif_pdb_filepath, 'w+') as file:
            file.write(''.join(lines))


def create_fixed_positions_dict(pdbs_dir, output_path, verbose=False):
    """
    Run conditional inverse folding to obtain sequences.

    Args:
       pdbs_dir:
            Directory containing rfdiffusion generated PDB files, where each file contains the generated 
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
        glob.glob(os.path.join(pdbs_dir, '*.pdb')),
        desc='Create fixed positions dictionary', disable=not verbose
    ):
        domain_name = processed_pdb_filepath.split('/')[-1].split('.')[0]
        # Create fixed positions dictionary
        with open(processed_pdb_filepath) as file:
            fixed_residue_indices = [
                int(line[22:26]) for line in file 
                if line.startswith('ATOM')  and line[60:66].strip() == '1.00'
            ]
        fixed_positions_dict = {}
        fixed_positions_dict[domain_name] = {
            'A': fixed_residue_indices
        }
        full_fixed_positions_dict = full_fixed_positions_dict | fixed_positions_dict
        # Predict
    with open(output_path,  'w') as file:
        file.write(json.dumps(full_fixed_positions_dict) + '\n')

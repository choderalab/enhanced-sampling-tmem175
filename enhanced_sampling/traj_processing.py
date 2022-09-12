import mdtraj as md
import os


def process_traj(traj_path, top_path, output_dir, stride=100):
    print(f"Loading traj from {traj_path}, {top_path}")
    traj = md.load_xtc(traj_path,
                       top=top_path)
    print("Removing solvent")
    ref = traj.remove_solvent()
    print("Superposing")
    protein_idx = ref.topology.select('protein')
    ref.superpose(reference=ref,
                  frame=0,
                  atom_indices=protein_idx)

    print("Saving start and end pdbs")
    ref[0].save_pdb(os.path.join(output_dir, "start.pdb"))
    ref[-1].save_pdb(os.path.join(output_dir, "end.pdb"))
    print("Saving strided trajectory")
    strided = ref[::stride]
    strided.save_xtc(os.path.join(output_dir, f"trajectory_{stride:03d}.xtc"))
    return True

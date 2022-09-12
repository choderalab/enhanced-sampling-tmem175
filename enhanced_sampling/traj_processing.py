import mdtraj as md
import os

def process_traj(traj_path, top_path, output_dir, stride=100):
    traj = md.load_xtc(traj_path,
                       top=top_path)
    ref = traj.remove_solvent()
    ref[0].save_pdb(os.path.join(output_dir, "start.pdb"))
    ref[-1].save_pdb(os.path.join(output_dir, "end.pdb"))
    strided = ref[::stride]
    strided.save_xtc(os.path.join(output_dir, f"trajectory_{stride:03d}.xtc"))
    return True


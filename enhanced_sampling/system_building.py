import openmm, json, os, yaml, mdtraj, bz2
from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PDBxFile
from simtk import unit


def load_input_dir(input_dir,
                   charmm_param_dir=None,
                   load_psf=True,
                   positions="state"
                   ):
    assert os.path.exists(input_dir)

    ## load psf
    if load_psf:
        psf_path = os.path.join(input_dir, "step5_input.psf")
        print(f"Loading psf from {psf_path}")
        psf = CharmmPsfFile(psf_path)
    else:
        psf = None

    if charmm_param_dir:
        print(f"Loading Charmm params from {charmm_param_dir}")
        param_paths = [os.path.join(charmm_param_dir, path) for path in os.listdir(charmm_param_dir)]
        params = CharmmParameterSet(*param_paths)
    else:
        params = None

    if positions == "state":
        input_state_path = os.path.join(input_dir, "state.xml.bz2")
        with bz2.open(input_state_path, 'rb') as infile:
            state = openmm.XmlSerializer.deserialize(infile.read().decode())

        if load_psf:
            x, y, z = state.getPeriodicBoxVectors()
            psf.setBox(x[0], y[1], z[2])

        positions = state.getPositions()
    elif positions == "pdb":
        state = None
        input_pdb_path = os.path.join(input_dir, "step5_input.pdb")
        print(f"loading positions from {input_pdb_path}")
        pdb = PDBFile(input_pdb_path)
        positions = pdb.positions

        if load_psf:
            print(f"loading box vectors")
            ## get box size from sysinfo.dat
            x, y, z = load_xyz_from_datfile(os.path.join(input_dir, "sysinfo.dat"))
            psf.setBox(x, y, z)

    elif positions == "pdbx":
        state = None
        input_pdb_path = os.path.join(input_dir, "final_frame.cif")
        print(f"loading positions from {input_pdb_path}")
        pdb = PDBxFile(input_pdb_path)
        positions = pdb.positions

    else:
        state = None
        positions = None
        pdb = None

    return {"psf": psf, "params": params, "state": state, "positions": positions, "pdb": pdb}


def load_xyz_from_datfile(dat_file_path):
    with open(dat_file_path) as file:
        data = json.load(file)
        x, y, z = map(float, data['dimensions'][:3]) * unit.angstroms
    print(f'Setting box size x: {x}, y: {y}, z:{z}')
    return (x, y, z)

def build_virtual_bond(psf, params):
    vbond_selections = params.virtual_bond
    print(vbond_selections)
    topology = mdtraj.Topology.from_openmm(psf.topology)
    print(topology)
    atom1 = topology.select(vbond_selections[0])[0]
    atom2 = topology.select(vbond_selections[1])[0]

    force = openmm.CustomBondForce("0*r")
    force.addBond(int(atom1), int(atom2))

    return force


def get_platform_from_params(params):
    platform_name = params.platform
    platform = openmm.Platform.getPlatformByName(platform_name)

    if not platform_name == 'CPU':
        print(f"Setting precision to {params.precision}")
        platform.setPropertyDefaultValue('Precision', params.precision)
        print(f"Setting DeterministicForces to 'true`")
        platform.setPropertyDefaultValue('DeterministicForces', 'true')

    return platform
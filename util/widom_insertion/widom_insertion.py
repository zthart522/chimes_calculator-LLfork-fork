import os

# Each worker process is single-threaded; parallelism comes from processes.
# Must be set before numpy is imported to take effect.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
import time
import multiprocessing as mp
from multiprocessing import util as mp_util

import numpy as np

sys.path.append(os.path.normpath(os.getcwd()))

if not os.path.exists("config.py"):
    print("Error: Cannot find config file")
    print("       Searched for: config.py in the current directory")
    exit(0)

import config

# The instance-based wrapper lives in serial_interface/api (not chimesFF/api)
api_dir = config.CHMS_REPO + "/serial_interface/api"
wrapper_file = api_dir + "/chimescalc_serial_py.py"

if not os.path.exists(wrapper_file):
    print("Error: Cannot find chimescalc_serial_py.py")
    print("       Searched for: ", wrapper_file)
    print("       Check CHMS_REPO path provided in config.py. Full paths must be used.")
    exit(0)

sys.path.insert(0, api_dir)
import chimescalc_serial_py as wrapper_py

# The instance API uses the _dl library
if hasattr(config, 'CHMS_LIB'):
    LIB_PATH = config.CHMS_LIB + "/libchimescalc_dl.so"
else:
    LIB_PATH = config.CHMS_REPO + "/build/libchimescalc_dl.so"

if not os.path.exists(LIB_PATH):
    print("Error: Cannot find chimes_calculator library file")
    print("       Searched for: ", LIB_PATH)
    print("       Ensure the library has been compiled (e.g. via ./install from " + config.CHMS_REPO + ")")
    exit(0)

# ---------------------------------------------------------
# Settings (module level so worker processes see them)
# ---------------------------------------------------------
PARAM_FILE = config.PARAM_FILE
TEMPERATURE = config.TEMP
KB = 8.617333262e-5 * 23.06            # kcal/mol/K; must match CHIMES energy units
BETA = 1.0 / (KB * TEMPERATURE)

INSERTED_TYPE = "C"
SEED = config.SEED
N_WIDOM = config.N_WIDOM                # insertions PER FRAME
SIGMA = 0.45                            # used for BOTH Z and the sampler
N_Z_SAMPLES = 100000
Z_CHUNK = 10000                         # chunk size to limit memory in compute_Z


def _prepare_frame(x, y, z, cell_a, cell_b, cell_c):
    """Prepare atom positions, orthogonal box lengths, and volume."""
    atom_pos = np.column_stack((x, y, z))

    box = np.array([
        np.linalg.norm(cell_a),
        np.linalg.norm(cell_b),
        np.linalg.norm(cell_c),
    ], dtype=float)

    V = np.prod(box)
    return atom_pos, box, V


def _density_factor(points, atom_pos, box, sigma):
    """
    Evaluate f(r) = max(0, 1 - sum_i exp(-r_i^2 / (2*sigma^2))).

    points can be a single position with shape (3,) or an array
    of positions with shape (n_points, 3).
    """
    points = np.asarray(points, dtype=float)
    dr = points[..., None, :] - atom_pos
    dr -= box * np.round(dr / box)
    r2 = np.sum(dr**2, axis=-1)

    return np.maximum(
        0.0,
        1.0 - np.sum(np.exp(-r2 / (2.0 * sigma**2)), axis=-1)
    )


def compute_Z(x, y, z, cell_a, cell_b, cell_c, rng,
              sigma=0.8, n_samples=100000, chunk=10000):
    """Estimate dimensionless Z = average of f(r) over the cell (chunked)."""
    atom_pos, box, V = _prepare_frame(
        x, y, z, cell_a, cell_b, cell_c
    )

    total = 0.0
    done = 0
    while done < n_samples:
        n = min(chunk, n_samples - done)
        points = rng.random((n, 3)) * box
        total += np.sum(_density_factor(points, atom_pos, box, sigma))
        done += n
    return total / n_samples


def evaluate_p_raw(x_new, y_new, z_new, x, y, z,
                   cell_a, cell_b, cell_c, sigma=0.8):
    """Evaluate unnormalized density p_raw(r) = f(r) / V."""
    atom_pos, box, V = _prepare_frame(
        x, y, z, cell_a, cell_b, cell_c
    )

    r = np.array([x_new, y_new, z_new], dtype=float)
    return _density_factor(r, atom_pos, box, sigma) / V


def sample_uniform_excluding_atoms(
    x, y, z, cell_a, cell_b, cell_c, rng, sigma=0.8
):
    """
    Rejection-sample positions with density proportional to f(r) / V.
    Assumes an orthogonal, axis-aligned periodic cell with origin (0,0,0).
    """
    atom_pos, box, V = _prepare_frame(
        x, y, z, cell_a, cell_b, cell_c
    )

    while True:
        r = rng.random(3) * box
        f = _density_factor(r, atom_pos, box, sigma)

        # Uniform proposal density is 1/V; acceptance probability is f.
        if rng.random() < f:
            p_raw = f / V
            return r, p_raw


def read_xyz(path, frame=None, skip=1):
    """
    Read every skip-th frame from a multi-frame XYZ file.

    Parameters
    ----------
    path : str
        Input XYZ file.
    frame : int or None
        Maximum number of frames to read after applying skip.
        None reads all selected frames.
    skip : int
        Read every skip-th frame (1 = every frame).

    Returns
    -------
    cell_as, cell_bs, cell_cs : lists
        Cell vectors for each selected frame.
    xs, ys, zs : lists of lists
        Cartesian coordinates for each frame.
    atom_typess : list of lists
        Element symbols for each frame.
    natomss : list
        Number of atoms in each frame.
    """

    if skip < 1:
        raise ValueError("skip must be >= 1")
    if frame is not None and frame < 0:
        raise ValueError("frame must be >= 0 or None")

    # Fixed orthogonal periodic cell, Angstrom.
    cell_a = [7.906793224, 0.0, 0.0]
    cell_b = [0.0, 12.259055612, 0.0]
    cell_c = [0.0, 0.0, 12.259055612]

    cell_as, cell_bs, cell_cs = [], [], []
    xs, ys, zs = [], [], []
    atom_typess, natomss = [], []

    with open(path, "r") as fin:
        iframe = 0
        nselected = 0

        while True:
            line = fin.readline()
            if not line:
                break

            # First line of an XYZ frame: number of atoms.
            try:
                natoms = int(line.strip())
            except ValueError:
                raise ValueError(
                    "Expected atom count at XYZ frame beginning "
                    "(frame index {}).".format(iframe)
                )

            # Second line is the comment/timestep.
            comment = fin.readline()
            if not comment:
                raise ValueError("Unexpected EOF after atom count.")

            # Select frames 0, skip, 2*skip, ...
            selected = (iframe % skip == 0)

            if selected and frame is not None and nselected >= frame:
                break

            x, y, z = [], [], []
            atom_types = []

            for _ in range(natoms):
                line = fin.readline()
                if not line:
                    raise ValueError("Unexpected EOF in atom coordinates.")

                if not selected:
                    continue

                fields = line.split()
                if len(fields) < 4:
                    raise ValueError("Invalid XYZ atom line: " + line.strip())

                atom_types.append(fields[0])
                x.append(float(fields[1]))
                y.append(float(fields[2]))
                z.append(float(fields[3]))

            if selected:
                cell_as.append(cell_a.copy())
                cell_bs.append(cell_b.copy())
                cell_cs.append(cell_c.copy())

                xs.append(x)
                ys.append(y)
                zs.append(z)
                atom_typess.append(atom_types)
                natomss.append(natoms)

                nselected += 1

            iframe += 1

    print("Read {} frames; selected {} frames.".format(
        iframe, len(natomss)
    ))

    return (
        cell_as, cell_bs, cell_cs,
        xs, ys, zs, atom_typess, natomss
    )


# ---------------------------------------------------------
# Worker-side code
# ---------------------------------------------------------
# One ChIMES instance per worker process. The ctypes library handle and
# instance pointer cannot be pickled, so they are created inside each
# worker by the pool initializer and never cross process boundaries.
_chimes_ptr = None


def _close_worker_chimes():
    global _chimes_ptr
    if _chimes_ptr is not None:
        wrapper_py.chimes_close_instance(_chimes_ptr)
        _chimes_ptr = None


def _init_worker(lib_path, param_file):
    global _chimes_ptr
    wrapper_py.chimes_wrapper = wrapper_py.init_chimes_wrapper(lib_path)
    _chimes_ptr = wrapper_py.chimes_open_instance()
    wrapper_py.set_chimes_instance(_chimes_ptr, small=True)
    wrapper_py.init_chimes_instance(_chimes_ptr, param_file, rank=0)

    # Pool workers are usually terminated rather than exited normally,
    # so atexit is unreliable; Finalize is the multiprocessing-safe hook.
    mp_util.Finalize(None, _close_worker_chimes, exitpriority=10)


def process_frame(task):
    """
    Run all Widom insertions for one frame.
    Returns (iframe, weights, text_block, Z, secs), where text_block holds
    the lines for widom_insertion.type_X.dat for this frame.
    """
    (iframe, seed_seq, cell_a, cell_b, cell_c,
     x, y, z, atom_types, natoms) = task

    start = time.time()
    rng = np.random.default_rng(seed_seq)   # independent stream per frame

    cell_a = np.asarray(cell_a, dtype=float)
    cell_b = np.asarray(cell_b, dtype=float)
    cell_c = np.asarray(cell_c, dtype=float)

    Z = compute_Z(x, y, z, cell_a, cell_b, cell_c, rng,
                  sigma=SIGMA, n_samples=N_Z_SAMPLES, chunk=Z_CHUNK)
    V = abs(np.dot(cell_a, np.cross(cell_b, cell_c)))

    weights = np.empty(N_WIDOM, dtype=float)
    lines = []

    for i in range(N_WIDOM):
        r_new, q_sampled = sample_uniform_excluding_atoms(
            x, y, z, cell_a, cell_b, cell_c, rng, sigma=SIGMA
        )
        x_new, y_new, z_new = r_new

        result = wrapper_py.calculate_chimes_insertion_instance(
            _chimes_ptr,
            natoms,
            x, y, z,
            atom_types,
            INSERTED_TYPE,
            float(x_new), float(y_new), float(z_new),
            cell_a, cell_b, cell_c
        )

        dE = result[4]
        boltzmann_factor = np.exp(-BETA * dE)
        weighted_factor = boltzmann_factor * Z / (V * q_sampled)
        weights[i] = weighted_factor

        # frame and iter are 1-based
        lines.append(
            f"{iframe + 1} {i + 1} "
            f"{x_new:.8f} {y_new:.8f} {z_new:.8f} "
            f"{boltzmann_factor:.16e} {weighted_factor:.16e}"
        )

    return iframe, weights, "\n".join(lines) + "\n", Z, time.time() - start


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def default_nproc():
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        if var in os.environ:
            return int(os.environ[var])
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def main():
    # -----------------------------------------------------
    # Convert trajectory and read frames
    # -----------------------------------------------------
    os.system(f"python3 {config.CHMS_REPO}/util/widom_insertion/src/lmp2xyz.py {config.TRAJ_FILE} widom_traj.xyz")

    max_frames = getattr(config, "MAX_FRAMES", None)
    frames = read_xyz('widom_traj.xyz', frame=max_frames, skip=config.FRAME_SKIP)
    nframes = len(frames[0])

    max_workers = getattr(config, "MAX_WORKERS", None)
    if max_workers is None:
        max_workers = default_nproc()
    nproc = max(1, min(int(max_workers), nframes))
    print(f"Running {nframes} frames on up to {nproc} worker processes "
          f"(MAX_WORKERS = {max_workers}).", flush=True)

    # One independent, reproducible RNG stream per frame. Results do not
    # depend on the number of workers or on scheduling order.
    seeds = np.random.SeedSequence(SEED).spawn(nframes)

    # tasks[i] = (i, seed, cell_a, cell_b, cell_c, x, y, z, atom_types, natoms)
    tasks = [
        (i, seeds[i], *(frames[k][i] for k in range(8)))
        for i in range(nframes)
    ]

    per_frame_averages = np.zeros(nframes)
    per_frame_mu = np.full(nframes, np.nan)
    per_frame_weights = [None] * nframes

    outfile_name = f"widom_insertion.type_{INSERTED_TYPE}.dat"

    pending = {}          # finished frames waiting to be written in order
    next_to_write = 0
    ndone = 0
    t0 = time.time()

    with open(outfile_name, "w") as fout:

        fout.write("# frame iter x_new y_new z_new "
                   "boltzmann_factor weighted_boltzmann_factor\n")

        # Do NOT load the ChIMES library in the parent before forking.
        with mp.Pool(processes=nproc,
                     initializer=_init_worker,
                     initargs=(LIB_PATH, PARAM_FILE)) as pool:

            for iframe, weights, text, Z, secs in pool.imap_unordered(
                    process_frame, tasks, chunksize=1):

                per_frame_weights[iframe] = weights
                avg = np.mean(weights)
                per_frame_averages[iframe] = avg
                if avg > 0.0:
                    per_frame_mu[iframe] = -(1.0 / BETA) * np.log(avg)

                # Write frames to the file in frame order, as soon as the
                # next expected frame is available.
                pending[iframe] = text
                while next_to_write in pending:
                    fout.write(pending.pop(next_to_write))
                    next_to_write += 1
                fout.flush()

                ndone += 1
                print(
                    f"[{ndone}/{nframes}] Frame {iframe + 1}: "
                    f"Z = {Z:.4f}, Avg = {avg:.8e}, "
                    f"mu = {per_frame_mu[iframe]:.8f}, "
                    f"Time = {secs:.2f} s, "
                    f"Elapsed = {time.time() - t0:.1f} s",
                    flush=True
                )

    # -----------------------------------------------------
    # Aggregate
    # -----------------------------------------------------
    all_weights = np.concatenate(per_frame_weights)

    overall_avg = np.mean(all_weights)
    overall_mu = -(1.0 / BETA) * np.log(overall_avg)
    ess = np.sum(all_weights) ** 2 / np.sum(all_weights ** 2)

    print("\nOverall results")
    print(f"Number of frames: {nframes}")
    print(f"Insertions per frame: {N_WIDOM}")
    print(f"Overall average: {overall_avg:.8e}")
    print(f"Chemical potential: {overall_mu:.8f}")
    print(f"Effective sample size: {ess:.1f} of {len(all_weights)}")

    np.savetxt(
        "per_frame_averages.dat",
        np.column_stack((
            np.arange(1, nframes + 1),
            per_frame_averages,
            per_frame_mu
        )),
        header="frame average_boltzmann_factor chemical_potential"
    )

    print(f"\nPer-insertion data written to: {outfile_name}")
    print("Per-frame summary written to: per_frame_averages.dat")


if __name__ == "__main__":
    main()
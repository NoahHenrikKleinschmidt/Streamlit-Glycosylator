import glycosylator as gl
import buildamol as bam
import streamlit as st
import os, zipfile
from pathlib import Path
import stmol
from tempfile import NamedTemporaryFile
from time import time

# INITIALIZE SESSION STATE
st.session_state.setdefault("glycans", {})
st.session_state.setdefault("glycan_optimized", False)
st.session_state.setdefault("glycan_conformers", [])
st.session_state.setdefault("glycan", None)
st.session_state.setdefault("scaffold", None)
st.session_state.setdefault("glyco_sites", [])
st.session_state.setdefault("scaffold_type", "scaffold")
st.session_state.setdefault("glycans_to_extend", [])
st.session_state.setdefault("available_sequons", gl.available_sequons())
st.session_state.setdefault("scaffold_conformers", [])
st.session_state.setdefault("scaffold_graph_edges", None)
st.session_state.setdefault("shield", None)
st.session_state.setdefault("env_hyperparameters", {})


def working_glycan(g=None):
    if g is not None:
        st.session_state["glycan"] = g
    return st.session_state["glycan"]


def scaffold(s=None):
    if s is not None:
        st.session_state["scaffold"] = s
        st.session_state["glyco_sites"].clear()
        st.session_state["glycans_to_extend"].clear()
        st.session_state["scaffold_conformers"].clear()
        reset_scaffold_graph()
        st.session_state["shield"] = None
        st.session_state["glycans_to_extend"].clear()
        st.session_state["scaffold_type"] = "scaffold"
    return st.session_state["scaffold"]


def reset_scaffold_graph():
    st.session_state["scaffold_graph_edges"] = None


def scaffold_type(s=None):
    if s is not None:
        st.session_state["scaffold_type"] = s
    return st.session_state["scaffold_type"]


_SUPPORTED_RESIDUE_DICT = {
    "Asparagine": "ASN",
    "Serine": "SER",
    "Threonine": "THR",
    "Ceramide": "CER",
}

SUPPORTED_RESIDUES = set(_SUPPORTED_RESIDUE_DICT.values())


def make_glycan(glycan_name, glycan_input):
    try:
        G = gl.glycan(glycan_input.strip(), glycan_name.strip())
        st.session_state["glycan"] = G
        st.session_state["glycans"][glycan_name] = G
        st.session_state["glycan_optimized"] = False
        return G
    except Exception as e:
        st.error(f"Error: {e}")
        return None


def load_glycan(glycan_file, infer_bonds):
    try:
        f = NamedTemporaryFile(suffix=".pdb", delete=False)
        with open(f, "w") as f:
            f.write(glycan_file.getvalue().decode("utf-8"))

        G = gl.glycan(f, glycan_file.name)
        if infer_bonds or G.count_bonds() == 0:
            G.reindex()
            G.infer_bonds(restrict_residues=False)
        G.infer_glycan_tree()
        st.session_state["glycan"] = G
        st.session_state["glycans"][G.id] = G
        st.session_state["glycan_optimized"] = False
        f.close()

        return G
    except Exception as e:
        st.error(f"Error: {e}")
        return None


def load_scaffold(scaffold_file, infer_existing_glycans, type: str = "scaffold"):
    try:
        pdb = scaffold_file.getvalue().decode("utf-8")

        if type == "protein":
            G = gl.Protein._from_pdb_string(pdb)
        elif type == "membrane":
            G = gl.Membrane._from_pdb_string(pdb)
        else:
            G = gl.Scaffold._from_pdb_string(pdb)

        st.session_state["scaffold_type"] = type

        if infer_existing_glycans:
            G.find_glycans()
            for glycan in G.glycans:
                st.session_state["glycans"][glycan.id] = glycan
        scaffold(G)
        return G
    except Exception as e:
        st.error(f"Error: {e}")
        return None


def extend_glycan(glycan_input):
    G = working_glycan()
    try:
        G.extend_to(glycan_input)
        return G
    except Exception as e:
        st.error(f"Error: {e}")
        return None


def extend_glycans_on_scaffold(glycan_input):

    S = scaffold()
    for glycan in st.session_state["glycans_to_extend"]:
        try:
            S.extend_glycan(glycan, glycan_input.strip())
        except Exception as e:
            st.error(f"Error on glycan {glycan}: {e}")
            raise e
    reset_scaffold_graph()


def show_glycan_2d3d(container=st):
    G = working_glycan()

    e1, e2, e3 = container.columns((1, 2, 1))
    e1.markdown("### 2D SNFG Representation")
    e1.pyplot(G.draw2d().figure)

    e2.markdown("### 3D Model")
    with e2:
        stmol.showmol(G.py3dmol().view)

    e3.markdown("### Glycan Histogram")
    e3.bar_chart(G.hist(), x="residue", y="count")


def show_scaffold_3d(type: str = None, S=None, container=st):
    S = S or scaffold()
    if S is None:
        return
    if type == "protein":
        v = S.py3dmol("cartoon", "spectrum")
    elif type == "membrane":
        v = S.py3dmol("sphere", "beige")
    else:
        v = S.py3dmol()
    stmol.showmol(v.view)


def show_protein_snfg(container=st):
    S = scaffold()
    if S is None:
        return
    container.pyplot(S.snfg().fig)


def optimize_glycan(container=st):
    columns = container.columns((1, 3))
    methods = {
        "atom distances": gl.optimizers.DistanceRotatron,
        "total energy": gl.optimizers.ForceFieldRotatron,
        "residue overlap": gl.optimizers.OverlapRotatron,
    }
    method = columns[0].selectbox(
        "Optimization metric",
        list(methods.keys()),
        help="Select the optimization metric to use.",
    )

    algorithms = {
        "particle swarm": "swarm",
        "genetic algorithm": "genetic",
        "simulated annealing": "anneal",
        "gradient descent": "scipy",
    }
    algorithm = columns[0].selectbox(
        "Optimization algorithm",
        list(algorithms.keys()),
        help="Select the optimization algorithm to use.",
    )

    n_conformers = columns[0].number_input(
        "Number of conformers",
        min_value=1,
        max_value=10,
        value=1,
        help="Select the number of conformers to generate.",
    )

    n_parallel = columns[0].number_input(
        "Number of parallel optimizations",
        min_value=1,
        max_value=10,
        value=1,
        help="Select the number of parallel optimizations to run.",
    )

    show_overlays = columns[0].checkbox(
        "Show overlays", value=False, help="Show overlays of the conformers."
    )

    optimize_button = columns[0].button("Optimize Model")

    if optimize_button:
        with st.spinner("Optimizing glycan model..."):

            G = working_glycan()

            graph = G.get_atom_graph()
            edges = G.get_residue_connections()
            edges = graph.direct_edges(G.get_root(), edges)

            env = methods[method](graph, edges)

            if n_parallel > 1:
                from concurrent.futures import ThreadPoolExecutor, as_completed

                with ThreadPoolExecutor(max_workers=n_parallel) as executor:
                    conformers = list(
                        executor.map(
                            lambda i: gl.optimizers.optimize(
                                G.copy(), env, algorithm=algorithms[algorithm]
                            ),
                            range(n_conformers),
                        )
                    )
                # conformers = gl.optimizers.parallel_optimize(
                #     G, [env] * n_conformers, n_processes=n_parallel, unify_final=False
                # )
            elif n_conformers > 1:
                conformers = [
                    gl.optimizers.optimize(
                        G.copy(), env, algorithm=algorithms[algorithm]
                    )
                    for i in range(n_conformers)
                ]
            else:
                conformers = [
                    gl.optimizers.optimize(
                        G.copy(), env, algorithm=algorithms[algorithm]
                    )
                    for i in range(n_conformers)
                ]
            st.session_state["glycan_conformers"] = conformers
            st.session_state["glycan_optimized"] = True

        with st.spinner("Drawing conformers..."):
            if show_overlays:
                with columns[1]:
                    v = G.py3dmol(color="gray")
                    s = dict(v.style)
                    s["stick"]["color"] = "red"
                    s["stick"]["opacity"] = 0.8
                    for conformer in conformers:
                        v.add(conformer, style=s)

                    stmol.showmol(v.view)

        return conformers


def export(prefix, filename, format, container=st):
    if prefix == "glycan":
        obj = working_glycan()
        conformers = st.session_state["glycan_conformers"]
    elif prefix == "scaffold":
        obj = scaffold()
        conformers = st.session_state.get("scaffold_conformers", None)

    if obj is None:
        container.error(f"No {prefix} found.")
        return

    format = format.lower().strip()
    if len(filename) == 0:
        container.error(
            f"Please select a format and provide a filename to export the {prefix}."
        )
        return

    filename = Path(filename)

    zip_fname = filename.with_suffix(".zip")

    can_make_zip = False
    if format == "pdb":

        data = bam.utils.pdb.encode_pdb(obj)
        main_fname = filename.with_suffix(".pdb")

        if conformers:
            can_make_zip = True
            files = []
            for i, conformer in enumerate(conformers):
                fname = filename.with_suffix(f".conf{i}.pdb")
                conformer.to_pdb(fname)
                files.append(fname)

            create_zip_file(files, zip_fname)
            for f in files:
                os.remove(f)

    elif format == "molfile":
        obj.to_molfile("tmp.mol")
        data = open("tmp.mol", "r").read()
        os.remove("tmp.mol")
        main_fname = filename.with_suffix(".mol")

        if conformers:
            can_make_zip = True
            files = []

            for i, conformer in enumerate(conformers):
                fname = filename.with_suffix(f".conf{i}.mol")
                conformer.to_molfile(fname)
                files.append(fname)

            create_zip_file(files, zip_fname)
            for f in files:
                os.remove(f)

    elif format == "pickle":
        obj.save("tmp.pkl")
        data = open("tmp.pkl", "rb").read()
        os.remove("tmp.pkl")
        main_fname = filename.with_suffix(".pkl")

        if conformers:
            can_make_zip = True
            files = []

            for i, conformer in enumerate(conformers):
                fname = filename.with_suffix(f".conf{i}.pkl")
                conformer.save(fname)
                files.append(fname)

            create_zip_file(files, zip_fname)
            for f in files:
                os.remove(f)

    elif format == "smiles":
        data = obj.to_smiles()
        main_fname = filename.with_suffix(".smi")

        if conformers:
            container.warning("SMILES format will not be exported for conformers.")

    elif format == "iupac":
        data = obj.to_iupac()
        main_fname = filename.with_suffix(".txt")

        if conformers:
            container.warning("IUPAC format will not be exported for conformers.")

    container.download_button(
        "Download Main File",
        data=data,
        file_name=(main_fname.name),
    )

    if can_make_zip:
        container.download_button(
            "Download ZIP of all conformers",
            data=open(zip_fname, "rb").read(),
            file_name=(zip_fname.name),
        )
        os.remove(zip_fname)


def create_zip_file(filenames, zip_filename):
    with zipfile.ZipFile(zip_filename, "w") as zip_file:
        for filename in filenames:
            zip_file.write(filename)


def attach_glycan():

    G = working_glycan()
    S = scaffold()

    sites = st.session_state["glyco_sites"]
    if not sites:
        st.error(
            "No glycosylation sites found. Please identify glycosylation sites first."
        )
        return

    _sites = {}
    for s in sites:
        if any(s is root.parent for root in S.get_glycans().keys()):
            st.info(f"Skipping site {s} as it is already glycosylated.")
            continue

        if s.resname not in _sites:
            _sites[s.resname] = []
        _sites[s.resname].append(s)

    hydrogenator = gl.structural.Hydrogenator()
    for name, residues in _sites.items():
        if name not in SUPPORTED_RESIDUES:
            st.error(f"Residue {name} is not supported for glycosylation.")
            continue

        # first make sure that all residues are not yet occupied
        residues = [
            i
            for i in residues
            if not any(i is root.parent for root in S.get_glycans().keys())
        ]

        # now we must make sure that all the scaffold residues are also
        # filled with bonds and the target atoms have hydrogen neighbors
        link = gl.get_linkage(name + "-glyco")

        to_delete = []
        for residue in residues:
            anchor = S.get_atom(link._stitch_ref_atoms[0], residue=residue)
            if len(S.get_neighbors(anchor)) == 0 or not link.can_be_target(S, residue):
                S.infer_bonds_for(residue)
                if not S.get_hydrogen(anchor):
                    hydrogenator.add_hydrogens(anchor, S)
                    if (
                        len(S.get_neighbors(anchor, filter=lambda x: x.element == "H"))
                        == 2
                    ):
                        S.remove_atoms(S.get_hydrogen(anchor))
                H = S.get_hydrogen(anchor)
                if H is None:
                    st.error(
                        f"Could not find or make a suitable hydrogen for atom {anchor} for residue {residue} to attach glycan."
                    )
                    to_delete.append(residue)
                    continue
                H.id = link.deletes[0][0]

        residues = [i for i in residues if i not in to_delete]

        if len(residues) == 0:
            continue

        S.attach(mol=G, residues=residues, inplace=True, chain="new")

    # make sure to reset since now the system has changed
    reset_scaffold_graph()


def show_glycosylated_3d(type: str, container=st):
    S = scaffold()
    if S is None:
        return
    if type == "protein":
        v = S.py3dmol("cartoon", "spectrum", glycans=False)
    elif type == "membrane":
        v = S.py3dmol("sphere", glycans=False)
    else:
        v = S.py3dmol()
    for glycan in S.glycans:
        v.add(glycan, style={"stick": {}})
    stmol.showmol(v.view)


def optimize_scaffold(container=st):
    columns = container.columns((1, 3))
    methods = {
        "atom distances": gl.optimizers.DistanceRotatron,
        "total energy": gl.optimizers.ForceFieldRotatron,
        "residue overlap": gl.optimizers.OverlapRotatron,
    }
    method = columns[0].selectbox(
        "Optimization metric",
        list(methods.keys()),
        help="Select the optimization metric to use.",
    )

    algorithms = {
        "particle swarm": "swarm",
        "genetic algorithm": "genetic",
        "simulated annealing": "anneal",
        "gradient descent": "scipy",
    }
    _algorithms = {
        "particle swarm": gl.optimizers.swarm_optimize,
        "genetic algorithm": gl.optimizers.genetic_optimize,
        "simulated annealing": gl.optimizers.anneal_optimize,
        "gradient descent": gl.optimizers.scipy_optimize,
    }
    algorithm = columns[0].selectbox(
        "Optimization algorithm",
        list(algorithms.keys()),
        help="Select the optimization algorithm to use.",
    )

    n_conformers = columns[0].number_input(
        "Number of conformers",
        min_value=1,
        max_value=10,
        value=1,
        help="Select the number of conformers to generate.",
    )

    n_parallel = columns[0].number_input(
        "Number of parallel optimizations",
        min_value=1,
        max_value=10,
        value=1,
        help="Select the number of parallel optimizations to run. Parallel optimization can make the process faster but because the optimizations will be run independently and later combined it will be less accurate than optimizing the entire system at once (if only optimizing for one single conformation; this does not apply if using parallel optimization to generate multiple conformers).",
    )

    only_clashing = columns[0].checkbox(
        "Only clashing glycans",
        value=True,
        help="Only optimize glycans that are clashing with the scaffold or themselves.",
    )

    use_numba = columns[0].checkbox(
        "Use Numba",
        value=False,
        help="Use Numba for faster optimization. Requires Numba to be installed. It may actually slow down the optimization depending on the size of the input and algorithm used. Use only with large systems.",
    )

    hyperparam_ext = columns[1].expander("Hyperparameters", expanded=False)

    if method == "atom distances":
        _dist_rot_hyperparams(hyperparam_ext)
    elif method == "total energy":
        _forcefield_rot_hyperparams(hyperparam_ext)
    elif method == "residue overlap":
        _overlap_rot_hyperparams(hyperparam_ext)
    if algorithm == "genetic algorithm":
        _genetic_hyperparams(hyperparam_ext)
    elif algorithm == "particle swarm":
        _swarm_hyperparams(hyperparam_ext)
    elif algorithm == "simulated annealing":
        _anneal_hyperparams(hyperparam_ext)
    elif algorithm == "gradient descent":
        _scipy_hyperparams(hyperparam_ext)

    optimize_button = columns[0].button("Optimize Model")
    clear_cache = columns[0].button(
        "Clear Cache",
        help="Clear the cache of the scaffold graph (if you make changes to the scaffold and the graph is not adjusted automatically).",
    )
    highlight_clashing = columns[0].button(
        "Show Clashing Glycans",
        help="Show a 3D view of all glycans that are clashing with the scaffold (red).",
    )

    if clear_cache:
        st.session_state["scaffold_graph_edges"] = None

    if highlight_clashing:
        with columns[1]:
            highlight_clashing_glycans()

    if optimize_button:

        if use_numba:
            gl.utils.use_numba()

        with st.spinner(
            "Optimizing glycosylated model (this is your chance to get coffee or go for lunch)..."
        ):
            t1 = time()
            st.write("Preparing graph...")
            current_graph = st.session_state["scaffold_graph_edges"]
            S = scaffold()
            if not current_graph or current_graph[0] != S.id:
                graph, edges = gl.optimizers.make_scaffold_graph(
                    S,
                    only_clashing_glycans=only_clashing,
                    include_root=True,
                    slice=4,
                )
                st.session_state["scaffold_graph_edges"] = (
                    graph._molecule.id,
                    graph,
                    edges,
                )
            elif current_graph[0] == S.id:
                graph, edges = current_graph[1:]

            st.write(
                f"Graph has {len(graph.nodes)} nodes and {len(edges)} edges to optimize."
            )
            if len(edges) == 0:
                st.error(
                    "No edges to optimize! Try setting the 'only clashing glycans' option to False to optimize all glycans."
                )
                return

            st.write("Preparing environment...")
            env = methods[method](
                graph, edges, **st.session_state["env_hyperparameters"]
            )
            st.write("Optimizing conformers...")
            if n_parallel > 1:
                from concurrent.futures import ThreadPoolExecutor, as_completed

                if n_conformers > 1:
                    st.write(
                        f"Running {n_parallel} independent parallel optimizations for different conformations."
                    )
                    env = gl.optimizers.ScaffoldRotatron(env)
                    with ThreadPoolExecutor(max_workers=n_parallel) as executor:
                        conformers = list(
                            executor.map(
                                lambda i: gl.optimizers.optimize(
                                    S.copy(),
                                    env,
                                    algorithm=algorithms[algorithm],
                                    **st.session_state.get("algo_hyperparameters", {}),
                                ),
                                range(n_conformers),
                            )
                        )
                else:

                    envs = gl.optimizers.split_environment(env, n_parallel)
                    envs = [gl.optimizers.ScaffoldRotatron(e) for e in envs]

                    algo = _algorithms[algorithm]
                    st.write(
                        f"Split the optimization into {n_parallel} independent parallel optimizations."
                    )
                    with ThreadPoolExecutor(max_workers=n_parallel) as executor:
                        solutions = list(
                            executor.map(
                                lambda e: algo(
                                    e,
                                    **st.session_state.get("algo_hyperparameters", {}),
                                ),
                                envs,
                            )
                        )
                    out = S.copy()
                    for (sol, _), e in zip(solutions, envs):
                        out = gl.optimizers.apply_rotatron_solution(sol, e, out)
                    conformers = [out]

            else:
                env = gl.optimizers.ScaffoldRotatron(env)
                conformers = [
                    gl.optimizers.optimize(
                        S.copy(),
                        env,
                        algorithm=algorithms[algorithm],
                        **st.session_state.get("algo_hyperparameters", {}),
                    )
                    for i in range(n_conformers)
                ]

        t2 = time()
        st.write(f"Optimization took {t2 - t1:.2f} seconds.")
        gl.utils.dont_use_numba()
        st.session_state["scaffold_conformers"] = conformers


def shield_scaffold(repeats, edge_samples, angle_step, visualize, container=st):
    S = scaffold()
    if edge_samples == 0:
        edge_samples = "all"

    shield = gl.GlycoShield(S)
    shield.simulate(
        repeats=repeats,
        angle_step=angle_step,
        edge_samples=edge_samples,
        visualize_conformations=visualize,
    )
    st.session_state["shield"] = shield


def highlight_clashing_glycans(container=st):
    S = scaffold()
    if S is None:
        return
    if st.session_state["scaffold_type"] == "protein":
        v = S.py3dmol("cartoon", "gray", glycans=False)
    elif st.session_state["scaffold_type"] == "membrane":
        v = S.py3dmol("sphere", color="gray", glycans=False)
    else:
        v = S.py3dmol(color="gray", glycans=False)

    for glycan in S.glycans:
        if glycan.clashes_with_scaffold():
            color = "red"
        else:
            color = "limegreen"
        v.add(glycan, style={"stick": {"color": color}})
    stmol.showmol(v.view)
    container.info("Red glycans are clashing with the scaffold. Green glycans are not.")


def _dist_rot_hyperparams(container=st):
    container.write("### Distance Rotatron Hyperparameters")
    container.write(
        "If the default values do not work, you can try playing around with the env_hyperparameters. Especially in cases where multiple glycans are close-by it may be necessary to adjust these values."
    )
    pushback = container.number_input(
        "Pushback",
        min_value=0.1,
        max_value=10.0,
        value=3.0,
        help="The amount to pushback for close-by atoms. This will disfavour clashes but makes the optimization more erratic and less suited for gradient descent.",
    )
    unfold = container.number_input(
        "Unfold",
        min_value=0.1,
        max_value=10.0,
        value=2.0,
        help="The amount to unfold to apply to the entire structure. This will favour extending the structure as a whole and makes the optimization smoother. But it may incentivize leaving some clashes 'for the greater good'.",
    )
    st.session_state["env_hyperparameters"] = {"pushback": pushback, "unfold": unfold}


def _overlap_rot_hyperparams(container=st):
    container.write("### Overlap Rotatron Hyperparameters")
    container.write(
        "If the default values do not work, you can try playing around with the env_hyperparameters."
    )

    artificial_spread = container.number_input(
        "Artificial Spread",
        min_value=1.0,
        max_value=20.0,
        value=4.0,
        help="The amount to artificially spread for the Gaussians of rotational units. Especially for for larger glycans that are 'medium distanced' it may be that the rotational units are too far away to compute an overlap without artificially extending the glycans.",
    )
    st.session_state["env_hyperparameters"] = {"artificial_spread": artificial_spread}


def _forcefield_rot_hyperparams(container=st):
    container.write("### Force Field Rotatron Hyperparameters")
    container.info(
        "There are no env_hyperparameters to adjust for the force field rotatron."
    )


def _genetic_hyperparams(container=st):
    container.write("### Genetic Algorithm Hyperparameters")
    max_generations = container.number_input(
        "Max Generations",
        min_value=1,
        max_value=1000000,
        value=100,
        help="The maximum number of generations to run the genetic algorithm.",
    )
    population_size = container.number_input(
        "Population Size",
        min_value=1,
        max_value=1000,
        value=50,
        help="The size of the population to run the genetic algorithm.",
    )
    parents = container.number_input(
        "Number of Parents",
        min_value=1,
        max_value=1000,
        value=10,
        help="The number of parents to select from the population.",
    )
    children = container.number_input(
        "Number of Children",
        min_value=1,
        max_value=1000,
        value=10,
        help="The number of children to generate from the parents.",
    )
    mutants = container.number_input(
        "Number of Mutants",
        min_value=1,
        max_value=1000,
        value=3,
        help="The number of random mutants to generate from the population.",
    )
    st.session_state["algo_hyperparameters"] = {
        "max_generations": max_generations,
        "population_size": population_size,
        "parents": parents,
        "children": children,
        "mutants": mutants,
    }


def _swarm_hyperparams(container=st):
    container.write("### Particle Swarm Hyperparameters")
    n_particles = container.number_input(
        "Number of Particles",
        min_value=1,
        max_value=1000,
        value=50,
        help="The number of particles to run the particle swarm optimization.",
    )
    max_steps = container.number_input(
        "Max Iterations",
        min_value=1,
        max_value=1000,
        value=30,
        help="The maximum number of iterations to run the particle swarm optimization.",
    )
    c1 = container.number_input(
        "C1",
        min_value=0.0,
        max_value=1.0,
        value=0.5,
        help="The cognitive parameter for the particle swarm optimization.",
    )
    c2 = container.number_input(
        "C2",
        min_value=0.0,
        max_value=1.0,
        value=0.3,
        help="The social parameter for the particle swarm optimization.",
    )
    w = container.number_input(
        "Inertia",
        min_value=0.0,
        max_value=1.0,
        value=0.9,
        help="The inertia parameter for the particle swarm optimization.",
    )
    st.session_state["algo_hyperparameters"] = {
        "n_particles": n_particles,
        "max_steps": max_steps,
        "c1": c1,
        "c2": c2,
        "w": w,
    }


def _anneal_hyperparams(container=st):
    container.write("### Simulated Annealing Hyperparameters")
    n_particles = container.number_input(
        "Number of Particles",
        min_value=1,
        max_value=1000,
        value=50,
        help="The number of particles to run the simulated annealing optimization.",
    )
    max_steps = container.number_input(
        "Max Iterations",
        min_value=1,
        max_value=1000,
        value=100,
        help="The maximum number of iterations to run the simulated annealing optimization.",
    )
    variance = container.number_input(
        "Variance",
        min_value=0.0,
        max_value=1.0,
        value=0.2,
        help="The variance for the simulated annealing optimization.",
    )
    cooldown_rate = container.number_input(
        "Cooling Rate",
        min_value=0.0,
        max_value=1.0,
        value=0.003,
        help="The cooling rate for the simulated annealing optimization.",
    )
    st.session_state["algo_hyperparameters"] = {
        "n_particles": n_particles,
        "max_steps": max_steps,
        "variance": variance,
        "cooldown_rate": cooldown_rate,
    }


def _scipy_hyperparams(container=st):
    method = container.selectbox(
        "Method",
        ["L-BFGS-B", "Powell", "Nelder-Mead", "CG", "BFGS", "TNC", "COBYLA", "SLSQP"],
        help="The optimization method to use for the scipy optimization.",
    )
    steps = container.number_input(
        "Max Iterations",
        min_value=1,
        max_value=10000000,
        value=100000,
        help="The maximum number of iterations to run the scipy optimization.",
    )
    st.session_state["algo_hyperparameters"] = {
        "method": method,
        "steps": steps,
    }

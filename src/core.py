import glycosylator as gl
import buildamol as bam
import streamlit as st
import stmol
from tempfile import NamedTemporaryFile
from time import time

st.session_state.setdefault("glycans", {})
st.session_state.setdefault("glycan_optimized", False)
st.session_state.setdefault("conformers", [])
st.session_state.setdefault("glycan", None)
st.session_state.setdefault("scaffold", None)
st.session_state.setdefault("scaffold_glycosylated", None)
st.session_state.setdefault("glyco_sites", [])

SUPPORTED_RESIDUES = set(
    (
        "ASN",
        "SER",
        "THR",
    )
)


def make_glycan(glycan_name, glycan_input):
    try:
        G = gl.glycan(glycan_input, glycan_name)
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


def load_scaffold(scaffold_file, infer_bonds, infer_existing_glycans):
    try:
        f = NamedTemporaryFile(suffix=".pdb", delete=False)
        with open(f, "w") as f:
            f.write(scaffold_file.getvalue().decode("utf-8"))

        G = gl.Scaffold.from_pdb(f)

        if G.count_bonds() == 0 or infer_bonds:
            G.reindex()
            G.infer_bonds(restrict_residues=False)

        if infer_existing_glycans:
            G.find_glycans()

        st.session_state["scaffold"] = G
        st.session_state["scaffold_glycosylated"] = None
        st.session_state["glyco_sites"] = []
        st.session_state["scaffold_conformers"] = []
        st.session_state["scaffold_graph_edges"] = None
        f.close()
        return G
    except Exception as e:
        st.error(f"Error: {e}")
        return None


def extend_glycan(glycan_input):
    G = st.session_state["glycan"]
    try:
        G.extend_to(glycan_input)
        return G
    except Exception as e:
        st.error(f"Error: {e}")
        return None


def show_glycan_2d3d(container=st):
    G = st.session_state["glycan"]

    e1, e2, e3 = container.columns((1, 2, 1))
    # e1 = container.expander("2D SNFG Representation", expanded=False)
    e1.markdown("### 2D SNFG Representation")
    e1.pyplot(G.draw2d().figure)

    # e2 = container.expander("3D Model", expanded=True)
    e2.markdown("### 3D Model")
    with e2:
        stmol.showmol(G.py3dmol().view)

    # e3 = container.expander("Glycan Histogram", expanded=False)
    e3.markdown("### Glycan Histogram")
    e3.bar_chart(G.hist(), x="residue", y="count")

    # columns[0].markdown("### 2D SNFG Representation")
    # columns[1].markdown("### 3D Model")
    # columns[1].plotly_chart(G.draw(atoms=False).figure, use_container_width=True)


def show_scaffold_3d(type: str = None, container=st):
    S = st.session_state["scaffold"]
    if S is None:
        return
    if type == "protein":
        v = S.py3dmol("cartoon", "spectrum")
    elif type == "membrane":
        v = S.py3dmol("sphere", "beige")
    else:
        v = S.py3dmol()
    for glycan in S.glycans.values():
        v.add(glycan, style={"stick": {}})
    stmol.showmol(v.view)

    # v = gl.utils.visual.MoleculeViewer3D()
    # v.draw_edges(*G.get_residue_connections(), showlegend=False)
    # container.plotly_chart(v.figure, use_container_width=True)


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

            G = st.session_state["glycan"]

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
            st.session_state["conformers"] = conformers
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
                # v = G.draw(atoms=False)
                # for conformer in conformers:
                #     v += conformer.draw(atoms=False, line_color="red")
                # columns[1].plotly_chart(v.figure, use_container_width=True)
        return conformers


def export_glycan(filename, export_conformers, container=st):
    G = st.session_state["glycan"]
    conformers = st.session_state.get("conformers", None)

    G.to_pdb(filename)
    container.success(f"Model saved to {filename}.pdb")

    if export_conformers and conformers:
        fname = ".".join(filename.split(".")[:-1])
        for i, conformer in enumerate(conformers):
            conformer.to_pdb(f"{fname}_conf{i}.pdb")
            container.success(f"Model saved to {filename}_conf{i}.pdb")


glycosite_find_modes = ["Existing Sequon", "Custom Sequon", "Residue Numbers"]


def find_glyco_sites(container=st):

    c1, c2 = container.columns((1, 2))
    method = c1.selectbox(
        "Method",
        glycosite_find_modes,
        help="Select the method to use for finding glycosylation sites.",
    )
    st.session_state["glycosite_find_method"] = method

    if method == glycosite_find_modes[0]:
        existing_sequons = gl.available_sequons()
        sequon = c1.selectbox(
            "Select sequon",
            existing_sequons,
            help="Select a pre-existing sequon to use for finding glycosylation sites.",
        )
        st.session_state["glycosite_find_sequon"] = sequon

    elif method == glycosite_find_modes[1]:
        sequon = c1.text_input(
            "Custom sequon",
            placeholder="(N)(?=[A-OQ-Z][ST])",
            help="Input a custom sequon to use for finding glycosylation sites. This needs to have one capturing group that matches the glycosylated residue.",
        )
        st.session_state["glycosite_find_sequon"] = sequon

    elif method == glycosite_find_modes[2]:
        residues = c1.text_area(
            "Residue numbers",
            placeholder="18, 75, 342, ...",
            help="Enter a list of residue numbers to identify glycosylation sites. Separate the numbers by commas or newlines.",
        )
        if residues:
            residues = [int(i) for i in residues.replace(",", "\n").split("\n")]
            st.session_state["glycosite_find_residues"] = residues

    show_sites = c1.checkbox(
        "Show glycosylation sites",
        value=False,
        help="Show the identified glycosylation sites.",
    )

    find_button = c1.button("Find Glycosylation Sites")
    S = st.session_state["scaffold"]
    if find_button:

        if method == glycosite_find_modes[0] or method == glycosite_find_modes[1]:
            sites = S.find_glycosylation_sites(
                sequon=st.session_state["glycosite_find_sequon"]
            )
            s = []
            for chain in sites:
                s.extend(sites[chain])
            st.session_state["glyco_sites"] = s
        elif method == glycosite_find_modes[2]:
            sites = S.get_residues(*st.session_state["glycosite_find_residues"])
            st.session_state["glyco_sites"] = sites

    if find_button and show_sites:
        sites = st.session_state.get("glyco_sites", [])
        if sites:

            v = S.py3dmol("cartoon", "spectrum")
            sites_mol = gl.Molecule.empty()
            sites_mol.add_residues(*(i.copy() for i in sites))
            sites_mol.infer_bonds()
            v.view.addModel(bam.utils.pdb.encode_pdb(sites_mol), "pdb")
            v.view.setStyle(
                {"model": 1},
                {
                    "sphere": {"color": "orange", "opacity": 0.8},
                },
            )
            # v.add(_v)
            # v.add(sites_mol, style={"sphere": {"color": "orange"}})

            with c2:
                stmol.showmol(v.view)
        else:
            c2.warning("No glycosylation sites found (yet).")

    e2 = c2.expander("Glycosylation Sites", expanded=not show_sites)
    sites = st.session_state.get("glyco_sites", {})
    if sites:

        e2.dataframe(
            {
                "Chain": [i.parent.id for i in sites],
                "Residue": [i.resname for i in sites],
                "Number": [i.serial_number for i in sites],
            }
        )

        for s in sites:
            if s.resname not in SUPPORTED_RESIDUES:
                e2.error(
                    f"Residue {s.resname} at position {s.serial_number} is not supported for glycosylation."
                )


def attach_glycan():

    G = st.session_state["glycan"]
    S = st.session_state["scaffold"]
    S_glyco = st.session_state.get("scaffold_glycosylated", None)
    if S_glyco is None or S_glyco.id != S.id:
        S_glyco = S.copy()
        st.session_state["scaffold_glycosylated"] = S_glyco

    sites = st.session_state.get("glyco_sites", [])
    if not sites:
        st.error(
            "No glycosylation sites found. Please identify glycosylation sites first."
        )
        return

    _sites = {}
    for s in sites:
        if s.resname not in _sites:
            _sites[s.resname] = []
        _sites[s.resname].append(s)

    for name, residues in _sites.items():
        if name not in SUPPORTED_RESIDUES:
            st.error(f"Residue {name} is not supported for glycosylation.")
            continue

        S_glyco.attach(mol=G, residues=residues, inplace=True, chain="new")


def show_glycosylated_3d(type: str, container=st):
    S = st.session_state["scaffold_glycosylated"]
    if S is None:
        return
    if type == "protein":
        v = S.py3dmol("cartoon", "spectrum")
    elif type == "membrane":
        v = S.py3dmol("sphere")
    else:
        v = S.py3dmol()
    for glycan in S.glycans.values():
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

    use_numba = columns[0].checkbox(
        "Use Numba",
        value=False,
        help="Use Numba for faster optimization. Requires Numba to be installed. It may actually slow down the optimization depending on the size of the input and algorithm used. Use only with large systems.",
    )

    hyperparam_ext = columns[1].expander("Hyperparameters", expanded=False)
    st.session_state.setdefault("env_hyperparameters", {})
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
    if optimize_button:

        if use_numba:
            gl.utils.use_numba()

        with st.spinner(
            "Optimizing glycosylated model (this is your chance to get coffee or go for lunch)..."
        ):
            t1 = time()
            st.write("Preparing graph...")
            current_graph = st.session_state.get("scaffold_graph_edges", None)
            S = st.session_state["scaffold_glycosylated"]
            if not current_graph or current_graph[0] != S.id:
                if len(S._internal_residues) == 0:
                    S.hollow_out()
                graph, edges = gl.optimizers.make_scaffold_graph(
                    S,
                    only_clashing_glycans=False,
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


def export_scaffold(filename, export_conformers, container=st):
    S = st.session_state["scaffold_glycosylated"]
    conformers = st.session_state.get("scaffold_conformers", None)

    S.to_pdb(filename)
    container.success(f"Model saved to {filename}.pdb")

    if export_conformers and conformers:
        fname = ".".join(filename.split(".")[:-1])
        for i, conformer in enumerate(conformers):
            conformer.to_pdb(f"{fname}_conf{i}.pdb")
            container.success(f"Model saved to {filename}_conf{i}.pdb")


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

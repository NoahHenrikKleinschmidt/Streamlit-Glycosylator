import streamlit as st
import core


def page_new_glycan():
    st.title("Generate Glycan Model")
    st.write(
        "Generate a 3D model for a glycan. You may either provide an IUPAC glycan sequence or use the GlyTouCan ID to fetch the glycan sequence. Alternatively, if you have a non-standard glycan you may provide the SMILES descriptor of the entire structure."
    )
    glycan_name = st.text_input(
        "Glycan Name",
        f"glycan-{len(st.session_state['glycans']) + 1}",
        help="Enter the name of the glycan for reference purposes.",
    )
    glycan_input = st.text_area(
        "Glycan Input",
        "",
        help="Enter the glycan sequence here in IUPAC format, as SMILES, or a GlyTouCan ID.",
    )
    show_glycan = st.checkbox("Show Glycan after generation", value=True)
    generate_button = st.button("Generate Model")

    if generate_button:
        G = core.make_glycan(glycan_name, glycan_input)
        st.success("Glycan model generated successfully.")

        if show_glycan:
            core.show_glycan_2d3d(st)

    if st.session_state["glycan"]:
        _optimize_glycan_subpage()
        _export_glycan_subpage()


def page_load_glycan():

    st.title("Load a Glycan Model from a file")
    st.write("Load a glycan model from a PDB file.")
    file_uploader = st.file_uploader(
        "Upload a PDB file",
        type=["pdb"],
        help="Upload a PDB file containing the glycan model.",
    )

    load_button = st.button("Load Model")
    infer_bonds = st.checkbox(
        "Infer bonds missing bonds",
        value=False,
        help="Infer bonds between atoms in the glycan if the PDB file does not have 'CONECT' lines.",
    )
    show_glycan = st.checkbox("Show Glycan after loading", value=False)

    if load_button and file_uploader is not None:
        G = core.load_glycan(file_uploader, infer_bonds)
        if G:
            st.success("Glycan model loaded successfully.")

            if show_glycan:
                core.show_glycan_2d3d(st)

    if st.session_state["glycan"]:

        st.markdown("---")
        st.markdown("## Extend Glycan Model")
        st.write(
            "If the glycan model is not complete, you can extend it by adding more residues. This can be done by providing an IUPAC sequence for the 'full' glycan that you want to build. The glycan model will be extended by adding the residues from the IUPAC sequence to the end of the current glycan model."
        )

        glycan_input = st.text_area(
            "Glycan Input",
            "",
            help="Enter the glycan sequence here in IUPAC format to extend the current glycan model.",
        )
        show_glycan = st.checkbox("Show Glycan after extension", value=True)
        extend_button = st.button("Extend Model")
        if extend_button:
            core.extend_glycan(glycan_input)
            st.success("Glycan model extended successfully.")

            if show_glycan:
                core.show_glycan_2d3d(st)

    if st.session_state["glycan"]:
        _optimize_glycan_subpage()
        _export_glycan_subpage()


def page_glycosylate_protein():
    st.title("Glycosylate a protein")
    st.write(
        "Glycosylate a protein by adding a glycan to one or more residues. You must provide a PDB file for the protein you want to glycosylate and have generated or loaded a glycan model."
    )

    col1, col2 = st.columns((1, 2.5))
    file_uploader = col1.file_uploader(
        "Upload a PDB file",
        type=["pdb"],
        help="Upload a PDB file containing the protein structure.",
    )

    infer_existing_glycans = col1.checkbox(
        "Infer existing glycans",
        value=False,
        help="Infer any already existing glycans in the protein structure.",
    )
    load_button = col1.button("Load Protein")
    show_protein = col1.button("Show Protein")

    if load_button:
        S = core.load_scaffold(file_uploader, infer_existing_glycans, "protein")

        if S:
            st.success("Protein loaded successfully.")

    if show_protein:
        with st.spinner("Rendering protein..."):
            with col2:
                core.show_scaffold_3d("protein")
                core.show_protein_snfg()

    st.markdown("---")
    _extend_glycans_subpage()

    st.markdown("---")
    st.markdown("## Identify Glycosylation sites")

    st.write(
        "To glycosylate the protein, first, select which residues should be modified. This can be done either by using a sequon pattern (pre-existing or user-defined regex) or by providing a list of residue numbers. Note that only one of the two methods can be used."
    )
    st.info(
        "Currently the app only supports N-linked (at Asparagines) and O-linked (at Serines and Threonines) glycosylation. If you need some non-canonical form of glycosylation, you will have to check out the Glycosylator library directly. It's not difficult to work with, so don't be afraid to try it out!"
    )

    _find_protein_glyco_sites_subpage()

    _attach_glycan_subpage()

    st.markdown("---")
    st.markdown("## [optional] Optimize Glycans")
    st.write(
        f"It is possible that the raw glycosylation output is not optimal and contains clashes with the {st.session_state['scaffold_type']} or between glycans. You can optimize the {st.session_state['scaffold_type']} structure to resolve these issues. Depending on the size of the {st.session_state['scaffold_type']} and the number of glycans attached this may take a while..."
    )

    core.optimize_scaffold(st)

    conformers = st.session_state.get("scaffold_conformers", None)
    if conformers:
        st.markdown("---")
        st.markdown("## Inspect Conformers")
        st.write(
            f"You can inspect the conformers generated for the glycosylated {st.session_state['scaffold_type']} structures. If you like you can also select one of the conformers to be the main one."
        )

        columns = st.columns((1, 5))
        conformer = columns[0].selectbox(
            "Select conformer",
            range(len(conformers)),
            help="Select the conformer to be displayed.",
        )

        select_main = columns[0].button(
            "Select as Main Conformer",
            help="Select the current conformer as the main conformer. This will be used from now on for any further glycosylations or optimizations.",
        )

        if select_main:
            st.session_state["scaffold"] = conformers[conformer]
            columns[0].success("Main conformer selected.")

        with columns[1]:
            core.show_scaffold_3d_compare_optimized(
                core.scaffold_type(), S=conformers[conformer]
            )

    _export_scaffold_subpage()


def _extend_glycans_subpage():

    st.markdown("## Extend existing glycans")
    st.write(
        "Extend any existing glycans by adding new residues to the end of the glycans. This can be done by providing a IUPAC sequence for the 'full' glycan that you want to build. The glycan model will be extended by adding the residues from the IUPAC sequence to the end of the current glycan model."
    )
    if not core.scaffold():
        st.error(
            f"No {st.session_state['scaffold_type']} loaded. Please load a {st.session_state['scaffold_type']} first."
        )
        return

    if len(core.scaffold().glycans) == 0:
        st.error(
            f"No glycans found in the {st.session_state['scaffold_type']}. If you are sure there are glycans already, be sure to select the 'Infer existing glycans' option when loading the {st.session_state['scaffold_type']}."
        )
        return

    S = core.scaffold()
    _glycans = {
        f"{glycan.id} @ {root.parent}": glycan
        for root, glycan in S.get_glycans().items()
    }
    glycans_to_modify = st.multiselect(
        "Select Glycans to Extend",
        list(_glycans.keys()),
        help="Select the glycans that should be extended.",
    )
    st.table(
        {
            "Glycan": list(_glycans.keys()),
            "Current IUPAC String": [i.to_iupac() for i in _glycans.values()],
        }
    )

    glycan_input = st.text_area(
        "Glycan Input",
        "",
        help="Enter the glycan sequence here in IUPAC format to extend the current glycans to.",
    )

    glycans_to_modify = [_glycans[g] for g in glycans_to_modify]
    st.session_state["glycans_to_extend"] = glycans_to_modify

    extend_glycans = st.button("Extend Glycans")
    if extend_glycans:
        core.extend_glycans_on_scaffold(glycan_input)
        st.success("Glycans extended successfully.")

    show_glycans = st.button("Show Glycans after extension")
    if show_glycans:
        with st.spinner("Rendering..."):
            core.show_scaffold_3d(
                st.session_state["scaffold_type"],
                S=st.session_state["scaffold_glycosylated"],
            )


def page_glycosylate_membrane():
    st.title("Glycosylate a Membrane")
    st.write(
        "Glycosylate a membrane by adding a glycan to one or more residues. You must provide a PDB file for the membrane you want to glycosylate and have generated or loaded a glycan model."
    )

    col1, col2 = st.columns((1, 2.5))
    file_uploader = col1.file_uploader(
        "Upload a PDB file",
        type=["pdb"],
        help="Upload a PDB file containing the membrane structure.",
    )

    infer_existing_glycans = col1.checkbox(
        "Infer existing glycans",
        value=False,
        help="Infer any already existing glycans on the membrane surface.",
    )
    load_button = col1.button("Load Membrane")
    show = col1.button("Show Membrane")

    if load_button:
        S = core.load_scaffold(file_uploader, infer_existing_glycans, "membrane")

        if S:
            st.success("Membrane loaded successfully.")

    if show:
        with st.spinner("Rendering membrane..."):
            with col2:
                core.show_scaffold_3d("membrane")

    st.markdown("---")
    _extend_glycans_subpage()

    st.markdown("---")
    st.markdown("## Identify Glycosylation sites")

    st.write(
        "To glycosylate the protein, first, select which residues should be modified. This can be done either by using a sequon pattern (pre-existing or user-defined regex) or by providing a list of residue numbers. Note that only one of the two methods can be used."
    )
    st.info(
        "Currently the app only supports N-linked (at Asparagines) and O-linked (at Serines and Threonines) glycosylation. If you need some non-canonical form of glycosylation, you will have to check out the Glycosylator library directly. It's not difficult to work with, so don't be afraid to try it out!"
    )

    _find_membrane_glyco_sites_subpage()

    _attach_glycan_subpage()

    st.markdown("---")
    st.markdown("## [optional] Optimize Glycans")
    st.write(
        "It is possible that the raw glycosylation output is not optimal and contains clashes with the protein or between glycans. You can optimize the protein structure to resolve these issues. Depending on the size of the protein and the number of glycans attached this may take a while..."
    )

    core.optimize_scaffold(st)

    conformers = st.session_state.get("scaffold_conformers", None)
    if conformers:
        st.markdown("---")
        st.markdown("## Inspect Conformers")
        st.write(
            "You can inspect the conformers generated for the glycosylated protein structures. If you like you can also select one of the conformers to be the main one."
        )

        columns = st.columns((1, 5))
        conformer = columns[0].selectbox(
            "Select conformer",
            range(len(conformers)),
            help="Select the conformer to be displayed.",
        )

        select_main = columns[0].button(
            "Select as Main Conformer",
            help="Select the current conformer as the main conformer. This will be used from now on for any further glycosylations or optimizations.",
        )

        if select_main:
            st.session_state["scaffold"] = conformers[conformer]
            columns[0].success("Main conformer selected.")

        C = conformers[conformer]
        v = C.py3dmol()
        with columns[1]:
            core.stmol.showmol(v.view)

    _export_scaffold_subpage()


def page_glycoshield():
    st.title("Glycan Quickshield")
    st.write(
        f"Simulate the shielding effect of the attached glycans on the {st.session_state['scaffold_type']} structure. This will compute an 'exposure' value for all {st.session_state['scaffold_type']} residues based on the shielding effect of the glycans. The exposure value is a measure of how much a residue is exposed to the solvent. The higher the value, the more exposed the residue is. The shielding effect is calculated by simulating the glycan conformations and measuring distances to the glycan residues. This will generate also a 3D visualization of the exposure values on the {st.session_state['scaffold_type']} structure."
    )

    if core.scaffold() is None or len(core.scaffold().glycans) == 0:
        st.error(
            f"No glycosylated {st.session_state['scaffold_type']} loaded. Please glycosylate a {st.session_state['scaffold_type']} first or load one with attached glycans."
        )
        return

    c1, c2 = st.columns((1, 2))
    edge_samples = c1.number_input(
        "Edge samples",
        value=3,
        min_value=0,
        max_value=10,
        help="The number of edges to subsample for each glycan while simulating. Set this to 0 to use *all* available edges.",
    )
    repeats = c1.number_input(
        "Repeats",
        value=3,
        min_value=1,
        max_value=100,
        help="The number of repeats to simulate the shielding effect. The more repeats, the more accurate the shielding effect will be. Subsampling will be performed for every repeat. If all edges are used only one repeat is performed and this input is ignored!",
    )

    angle_step = c1.number_input(
        "Angle step",
        value=50,
        min_value=10,
        max_value=180,
        help="The angle step to use when rotating the glycan conformations. The smaller the angle step, the more conformations will be generated and the more accurate the shielding effect will be. However, this will also increase the computation time (and memory load if conformations are also captures).",
    )

    visualize = c1.checkbox(
        "Capture Conformation Overlay",
        value=True,
        help="Visualize the shielding effect in 3D by overlaying all accepted conformations. This may take some time and require considerable amounts of memory depending on the size of the system and the number of glycans attached.",
    )

    shield_button = c1.button("Simulate Shielding")
    display_results = c1.button("Display Results")
    if shield_button:
        with st.spinner("Simulating..."):
            core.shield_scaffold(repeats, edge_samples, angle_step, visualize)

        st.success("Shielding simulation successful.")

    if shield_button or display_results:
        shield = st.session_state["shield"]
        if shield is None:
            return

        with c2:
            v = shield.py3dmol()
            core.stmol.showmol(v.view)

        st.plotly_chart(shield.plot_exposure("plotly"))
        e = st.expander("Exposure Values", expanded=False)
        e.markdown(
            "The exposure values for each residue in the protein. You can download the data as a CSV by clicking on the small icon in the top right corner of the table."
        )
        e.dataframe(shield.df)


def _export_scaffold_subpage(container=st, prefix=None):
    prefix = prefix or st.session_state["scaffold_type"]

    container.markdown("---")
    container.markdown(f"## Export Glycosylated {prefix.title()}")
    container.write(f"You can export the glycosylated {prefix} model(s) to PDB.")

    format = st.selectbox(
        "Format",
        ["PDB", "MOLFILE", "PICKLE"],
        help="Select the format to export the glycan model to.",
    )

    filename = container.text_input(
        "Filename",
        placeholder=f"my-glycosylated-{prefix}.pdb",
        help=f"Enter the filename for the PDB file to save the glycosylated {prefix} to.",
    )

    core.export(prefix, filename, format)


def _attach_glycan_subpage(container=st):

    container.markdown("---")
    container.markdown("## Glycosylate new residues")
    container.write(
        "Select the glycan to use when glycosylating. This requires that you have generated or loaded a glycan model previously. This glycan will be attached to all the residues selected in the previous step. If you need heterogeneous glycosylation, be sure to select only the residues that should have the same glycan attached at one time. The process can be repeated for different glycans."
    )

    glycans = st.session_state["glycans"]
    if not glycans:
        st.error("No glycans available. Please generate or load a glycan model first.")
        return

    c1, c2 = st.columns((1, 2))
    glycan_names = list(glycans.keys())
    glycan = c1.selectbox(
        "Select Glycan",
        glycan_names,
        help="Select the glycan to attach to the protein at the given residues.",
    )
    st.session_state["glycan"] = glycans[glycan]

    if not st.session_state["glyco_sites"]:
        c1.error(
            "No glycosylation sites found. Please identify glycosylation sites first."
        )
        return

    c1.info(
        f"The glycan '{glycan}' will be attached to residues: {st.session_state['glyco_sites']}."
    )

    attach_button = c1.button("Glycosylate")
    show_glycosylated = c1.button("Show the Structure")

    if attach_button:

        with st.spinner("Glycosylating..."):
            core.attach_glycan()

        st.success("Glycosylation successful.")

    if show_glycosylated:
        with st.spinner("Rendering..."):
            with c2:
                core.show_glycosylated_3d(st.session_state["scaffold_type"], c2)


def _export_glycan_subpage():
    st.markdown("---")
    st.markdown("## Export Glycan Model")
    st.write("You can export the glycan model(s).")

    format = st.selectbox(
        "Format",
        ["PDB", "MOLFILE", "IUPAC", "SMILES", "PICKLE"],
        help="Select the format to export the glycan model to.",
    )
    filename = st.text_input(
        "Filename",
        placeholder="my-new-glycan.pdb",
        help="Enter the filename for the glycan(s) to be saved to.",
    )

    core.export("glycan", filename, format)


def _optimize_glycan_subpage():
    st.markdown("---")
    st.markdown("## Optimize Glycan Model")
    st.write(
        "If there is something wrong with the 3D model, we can optimize it here. This may take some time depending on the size of the glycan."
    )
    conformers = core.optimize_glycan(st)

    if st.session_state["glycan_optimized"]:
        st.markdown("---")
        st.markdown("## Inspect Conformers")
        st.write(
            "You can inspect the conformers generated for the glycan model. If you like you can also select one of the conformers to be the main one."
        )

        columns = st.columns((1, 5))
        conformers = st.session_state["glycan_conformers"]
        conformer = columns[0].selectbox(
            "Select conformer",
            range(len(conformers)),
            help="Select the conformer to be displayed (red). The main conformer will be used for further analysis (black).",
        )

        select_main = columns[0].button(
            "Select as Main Conformer",
            help="Select the current conformer (red) as the main conformer.",
        )
        if select_main:
            st.session_state["glycan"] = conformers[conformer]
            st.session_state["glycans"][conformers[conformer].id] = conformers[
                conformer
            ]
            columns[0].success("Main conformer selected.")

        G = st.session_state["glycan"]

        with columns[1]:
            v = G.py3dmol(color="gray")
            s = dict(v.style)
            s["stick"]["color"] = "red"
            s["stick"]["opacity"] = 0.9
            v.add(conformers[conformer], style=s)
            core.stmol.showmol(v.view)

        # v = G.draw(atoms=False)
        # v += conformers[conformer].draw(atoms=False, line_color="red")
        # columns[1].plotly_chart(v.figure, use_container_width=True)


def _find_protein_glyco_sites_subpage(container=st):
    glycosite_find_modes = ["Existing Sequon", "Custom Sequon", "Residue Numbers"]

    c1, c2 = container.columns((1, 2))
    method = c1.selectbox(
        "Method",
        glycosite_find_modes,
        help="Select the method to use for finding glycosylation sites.",
    )
    st.session_state["glycosite_find_method"] = method

    if method == glycosite_find_modes[0]:
        existing_sequons = st.session_state["available_sequons"]
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

    find_button = c1.button("Find Glycosylation Sites")
    show_sites = c1.button(
        "Show glycosylation sites",
        help="Show the identified glycosylation sites.",
    )
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

    if show_sites:
        sites = st.session_state.get("glyco_sites", [])
        if sites:

            v = S.py3dmol("cartoon", "spectrum")
            for site in sites:
                v.add(site, style={"sphere": {"color": "orange", "opacity": 0.8}})

            with c2:
                e = c2.expander("3D View", expanded=True)
                with e:
                    core.stmol.showmol(v.view)

                if st.session_state["scaffold_type"] == "protein":
                    e2 = c2.expander("2D View", expanded=False)
                    e2.pyplot(S.snfg().fig)
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
            if s.resname not in core.SUPPORTED_RESIDUES:
                e2.error(
                    f"Residue {s.resname} at position {s.serial_number} is not supported for glycosylation."
                )


def _find_membrane_glyco_sites_subpage(container=st):
    glycosite_find_modes = ["Residue Names", "Residue Numbers"]
    c1, c2 = container.columns((1, 2))
    method = c1.selectbox(
        "Method",
        glycosite_find_modes,
        help="Select the method to use for finding glycosylation sites.",
    )
    st.session_state["glycosite_find_method"] = method

    residues = None
    if method == glycosite_find_modes[0]:
        residues = c1.multiselect(
            "Residue names",
            list(set(i.resname for i in st.session_state["scaffold"].get_residues())),
            help="Select one or more residues to use as glycosylation sites.",
        )

    elif method == glycosite_find_modes[1]:
        residues = c1.text_area(
            "Residue numbers",
            placeholder="18, 75, 342, ...",
            help="Enter a list of residue numbers to identify glycosylation sites. Separate the numbers by commas or newlines.",
        )
        if residues:
            residues = [int(i) for i in residues.replace(",", "\n").split("\n")]

    if residues:
        st.session_state["glycosite_find_residues"] = residues

    find_button = c1.button("Find Glycosylation Sites")
    show_sites = c1.button(
        "Show glycosylation sites",
        help="Show the identified glycosylation sites.",
    )
    S = st.session_state["scaffold"]
    if find_button:

        sites = S.get_residues(*st.session_state["glycosite_find_residues"])
        st.session_state["glyco_sites"] = sites

    if show_sites:
        sites = st.session_state.get("glyco_sites", [])
        if sites:

            v = S.py3dmol()
            for site in sites:
                v.add(site, style={"sphere": {"color": "orange", "opacity": 0.8}})

            with c2:
                e = c2.expander("3D View", expanded=True)
                with e:
                    core.stmol.showmol(v.view)
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
            if s.resname not in core.SUPPORTED_RESIDUES:
                e2.error(
                    f"Residue {s.resname} at position {s.serial_number} is not supported for glycosylation."
                )


app_pages = {
    "Generate New Glycan": page_new_glycan,
    "Load/Modify Existing Glycan": page_load_glycan,
    "Glycosylate Protein": page_glycosylate_protein,
    "Glycosylate Membrane": page_glycosylate_membrane,
    "Glycan Quickshield": page_glycoshield,
}


def build_ui():
    st.set_page_config(
        page_title="Glycosylator",
        page_icon="./src/light_logo.png",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    if st.get_option("theme.primaryColor") == "#000000":
        st.sidebar.image(
            "./src/dark_logo.png",
            use_column_width=True,
        )
    else:
        st.sidebar.image(
            "./src/light_logo.png",
            use_column_width=True,
        )

    st.sidebar.title("Navigation")
    st.sidebar.write("Select here the application you want to use")
    page = st.sidebar.selectbox("Application", list(app_pages.keys()))

    app_pages[page]()

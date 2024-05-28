import streamlit as st
import core


def page_generate():
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
    show_glycan = st.checkbox("Show Glycan after generation", value=False)
    generate_button = st.button("Generate Model")

    if generate_button:
        G = core.make_glycan(glycan_name, glycan_input)
        st.success("Glycan model generated successfully.")

        if show_glycan:
            core.show_glycan_2d3d(st)

    _optimize_glycan_subpage()

    _export_glycan_subpage()


def page_load():

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

    if st.session_state.get("glycan", None):

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
        show_glycan = st.checkbox("Show Glycan after extension", value=False)
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

    show_protein = col1.checkbox("Show Protein after loading", value=False)
    infer_bonds = col1.checkbox(
        "Infer bonds missing bonds",
        value=False,
        help="Infer bonds between atoms in the scaffold if the PDB file does not have 'CONECT' lines for amino acid residues.",
    )
    infer_existing_glycans = col1.checkbox(
        "Infer existing glycans",
        value=False,
        help="Infer any already existing glycans in the protein structure.",
    )
    load_button = col1.button("Load Protein")

    if load_button:
        S = core.load_scaffold(file_uploader, infer_bonds, infer_existing_glycans)

        if S:
            st.success("Protein loaded successfully.")

            if show_protein:
                with st.spinner("Rendering protein..."):
                    with col2:
                        core.show_scaffold_3d("protein")

    st.markdown("---")
    st.markdown("## Identify Glycosylation sites")

    st.write(
        "To glycosylate the protein, first, select which residues should be modified. This can be done either by using a sequon pattern (pre-existing or user-defined regex) or by providing a list of residue numbers. Note that only one of the two methods can be used."
    )
    st.info(
        "Currently the app only supports N-linked (at Asparagines) and O-linked (at Serines and Threonines) glycosylation. If you need some non-canonical form of glycosylation, you will have to check out the Glycosylator library directly. It's not difficult to work with, so don't be afraid to try it out!"
    )

    core.find_glyco_sites()

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
        v = C.py3dmol("cartoon", "spectrum")
        for glycan in C.glycans.values():
            v.add(glycan, style={"stick": {}})
        with columns[1]:
            core.stmol.showmol(v.view)

    _export_scaffold_subpage()


def _export_scaffold_subpage(container=st, prefix="protein"):

    container.markdown("---")
    container.markdown(f"## Export Glycosylated {prefix.title()}")
    container.write(f"You can export the glycosylated {prefix} model(s) to PDB.")

    filename = container.text_input(
        "Filename",
        value=f"my-glycosylated-{prefix}.pdb",
        help="Enter the filename for the PDB file to save the glycosylated protein model to.",
    )

    export_conformers = container.checkbox(
        "Export all conformers",
        value=True,
        help="Export all conformers to separate PDB files (they will inherit the same filename but add '_conf<N>' at the end).",
    )

    export_button = container.button("Export Model")

    if export_button:
        core.export_scaffold(filename, export_conformers)


def _attach_glycan_subpage(container=st):

    container.markdown("---")
    container.markdown("## Select Glycan to attach")
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

    c1.info(
        f"The glycan '{glycan}' will be attached to residues: {st.session_state['glyco_sites']}."
    )

    show_glycosylated = c1.checkbox(
        "Show the structure after glycosylation", value=False
    )
    attach_button = c1.button("Glycosylate")

    if attach_button:

        with st.spinner("Glycosylating..."):
            core.attach_glycan()

        st.success("Glycosylation successful.")

    if show_glycosylated:
        with st.spinner("Rendering..."):
            with c2:
                core.show_glycosylated_3d("protein", c2)


def show_page3():
    st.title("Page 3")
    # Add your content for Page 3 here


def _export_glycan_subpage():
    st.markdown("---")
    st.markdown("## Export Glycan Model")
    st.write("You can export the glycan model(s) to PDB.")

    filename = st.text_input(
        "Filename",
        value="my-glycan.pdb",
        help="Enter the filename for the PDB file to save the glycan model to.",
    )

    export_conformers = st.checkbox(
        "Export all conformers",
        value=True,
        help="Export all conformers to separate PDB files (they will inherit the same filename but add '_conf<N>' at the end).",
    )

    export_button = st.button("Export Model")

    if export_button:
        core.export_glycan(filename, export_conformers)


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
        conformers = st.session_state["conformers"]
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


app_pages = {
    "Generate New Glycan": page_generate,
    "Load Existing Glycan": page_load,
    "Glycosylate Protein": page_glycosylate_protein,
    "Glycosylate Membrane": show_page3,
}


def build_ui():
    st.set_page_config(
        page_title="Glycosylator",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.sidebar.title("Navigation")
    st.sidebar.write("Select here the application you want to use")
    page = st.sidebar.selectbox("Application", list(app_pages.keys()))

    app_pages[page]()

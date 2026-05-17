import marimo

__generated_with = "0.23.5"
app = marimo.App(width="full")


@app.cell
def _(mo):
    # DB_BASE_URLS = "https://bigsdb.pasteur.fr/api"
    BASE_API = {
        "PubMLST": "https://rest.pubmlst.org",
        "Pasteur": "https://bigsdb.pasteur.fr/api",
    }
    select_mlst_api = mo.ui.dropdown(options=BASE_API, value="PubMLST")
    select_mlst_api
    return (select_mlst_api,)


@app.cell
def _(get_available_dbs, mo, params_doras_int):
    available_dbs = get_available_dbs(params_doras_int)
    selected_db_key = mo.ui.dropdown(
        label="Choose a DB",
        options=available_dbs,
        # value="escherichia",
    ).form()
    return (selected_db_key,)


@app.cell
def _(selected_db_key):
    selected_db_key
    return


@app.cell
def _(mo, toml_mode):
    mo.stop(
        toml_mode is True,
    )

    # selected_db_key
    isolates_list_dirs = mo.ui.file_browser(
        label="Select Target folders", multiple=True, selection_mode="directory"
    )
    # create_toml.py
    exp_name_text = mo.ui.text(
        label="Name of the experiment", max_length=7, placeholder="Optional"
    )
    BASE_FOLDER = "/test_data/"
    return (exp_name_text,)


@app.cell
def _(db_manual_input, mo, requests, selected_db_key, toml_mode):
    mo.stop(toml_mode is True, mo.md("TOML mode"))
    mo.stop(
        selected_db_key.value is None and db_manual_input.value is None,
        mo.md("Select a DB"),
    )
    if db_manual_input.value:
        db_to_use = db_manual_input.value
    else:
        db_to_use = selected_db_key.value
    schemes_available = requests.get(url=f"{db_to_use}/schemes")
    dict_scheme_available = {
        i["description"]: i["scheme"] for i in schemes_available.json()["schemes"]
    }

    scheme_select = mo.ui.dropdown(
        options=dict_scheme_available, label="Choose a scheme:"
    )
    scheme_select
    # names["diphtheria"]
    return (scheme_select,)


@app.cell
def _(mo):
    toml_mode = mo.ui.switch(label="Load DORAS parameters from TOML")
    toml_mode
    return (toml_mode,)


@app.cell
def _():
    # Example usage

    # barcodes = [

    #     f"barcode{i:02d}_{j}h_medaka"
    #     for i in range(1, 21)
    #     # for j in [6]
    #     # f"bc{i:02d}_{j}h_medaka" for i in range(1, 21) for j in range(10, 11)
    # ]
    # mo.stop(bar)

    barcodes = ["barcode05"]
    # test_folder = "20250321_TimetoST_analysis/time_to_st_vres/WGS/"
    test_folder = "testdoras/"
    request_params = {"details": "true", "base64": "true"}
    return


@app.cell
def _():
    return


@app.cell
def _(load_params_from_toml, mo, toml_mode):
    mo.md("##Load toml or use UI to define params")
    if toml_mode.value:
        try:
            params_doras_toml = load_params_from_toml("./_doras_config.toml")
        except Exception as e:
            mo.md(f"Failed to load TOML: {e}")
    else:
        mo.md("Using manually defined parameters via UI.")
    return (params_doras_toml,)


@app.cell
def _(scheme_select):
    url = scheme_select.value
    return (url,)


@app.cell
def _(mo, scheme_select, select_mlst_api, selected_db_key):
    mo.md(f"""
    The URL for the selected db ({selected_db_key.value}) from {select_mlst_api.value} is: {scheme_select.value}
    """)
    return


@app.cell
def _(mo):
    mo.md("## Identity of the sequence")
    # "https://rest.pubmlst.org/db/pubmlst_efaecium_seqdef"

    db_manual_input = mo.ui.text(
        label="Manual input of the DB to use",
        kind="url",
        placeholder="https://rest.pubmlst.org/db/pubmlst_efaecium_seqdef",
    ).form()
    # mo.stop(selected_db_key.value is not None, mo.md("DB already selected"))
    db_manual_input
    return (db_manual_input,)


@app.cell
def _(mo, params_doras_manual, params_doras_toml, select_mlst_api, toml_mode):
    mo.stop(select_mlst_api.value is None, mo.md("Choose the base API"))
    # Use manual params if TOML mode is off
    if toml_mode.value:
        params_doras_int = params_doras_toml
    else:
        params_doras_int = params_doras_manual
    params_doras_int
    return (params_doras_int,)


@app.cell
def _(Path, config, exp_name_text, mo, request_form, toml):
    # Define the parameters as they are in extension_system_Test.py
    # FetchConsensusDataParams = {
    #     "mlst_genes_path": "./datasets/ecoli_mlst_shortnames.fasta",
    #     "fastq_files": "./datasets/exp12_wgs_barcode03_trim.fastq.gz",
    #     "test_folder": "test_data",
    #     "barcodes": ["medaka_ratio0.7"],
    #     "postprocessing": False
    # }
    # Path to the TOML file you want to create
    toml_file_path = Path(f"{exp_name_text.value}_doras_config.toml")

    mo.stop(request_form.value is None, "Please confirm fetching data")
    # Write the configuration dictionary to a TOML file
    with open(toml_file_path, "w") as toml_file:
        toml.dump(config, toml_file)

    print(f"TOML file created at {toml_file_path}")
    return


@app.cell
def _(mo):
    st_sample = mo.ui.number(start=1, stop=100, value=1)
    return (st_sample,)


@app.cell
def _(mo, st_sample):
    stop_sample = mo.ui.number(start=st_sample.value + 1, stop=100)
    range_samples = mo.md("First sample: {st_sample} Last Sample {stop_sample}").batch(
        st_sample=st_sample, stop_sample=stop_sample
    )
    return (range_samples,)


@app.cell
def _(mo, range_samples):
    refresh_query = mo.ui.refresh(
        options=["1m", "30s", "5s", "2m"],
        default_interval="1m",
    )
    mo.stop(
        range_samples.elements["st_sample"].value
        >= range_samples.elements["stop_sample"].value,
        mo.md("cannot have first sample higher than the last one"),
    )
    # range_samples
    return (refresh_query,)


@app.cell
def _():
    import altair as alt

    return (alt,)


@app.cell
def _(mo):
    details_switch = mo.ui.switch(label=r"##Alignment Details", value=True)
    return (details_switch,)


@app.cell
def _():
    return


@app.cell
def _(details_switch, mo, request_form):
    mo.hstack([details_switch, request_form])
    return


@app.cell
async def _(
    details_switch,
    mo,
    params_doras,
    query_consensus,
    refresh_query,
    request_form,
):
    refresh_query
    mo.stop(request_form.value is None, mo.md("Submit the form"))
    full_results = await query_consensus(params_doras, details=details_switch.value)
    return (full_results,)


@app.cell
def _():
    return


@app.cell
def _(create_heatmap_plot, create_heatmap_plot_ST, full_results, mo):
    if full_results.get("details") is not None:
        plots_details = create_heatmap_plot(full_results["details"])
        plots_sts = create_heatmap_plot_ST(full_results["st"])
        Full_plot = plots_details | plots_sts
    else:
        plots_sts = create_heatmap_plot_ST(full_results["st"])
        mo.ui.altair_chart(plots_sts)
    return (Full_plot,)


@app.cell(disabled=True)
def _(full_results, mo):
    mo.ui.table(full_results["st"], page_size=20)
    return


@app.cell
def _(full_results):
    full_results
    return


@app.cell
def _(Full_plot, mo):
    mo.ui.altair_chart(Full_plot)
    return


@app.cell(hide_code=True)
def code_queries():
    # import asyncio
    # from pydantic import BaseModel
    # import aiohttp
    # import logging
    # from functools import cache

    # # Set up logging with debug mode included
    # logger = logging.getLogger(__name__)
    # logging.basicConfig(
    #     level=logging.DEBUG,  # Changed to DEBUG for detailed logs
    #     format="%(asctime)s - %(levelname)s - %(message)s",
    # )

    # def generate_dataframe(response_list):
    #     data = []
    #     for response in response_list:
    #         exact_matches = response.get("exact_matches", [])
    #         best_match = response.get("best_match", {})
    #         if exact_matches:
    #             for match in exact_matches:
    #                 data.append(
    #                     {
    #                         "Isolate": response.get("Isolate", "None"),
    #                         "Loci": response.get("Loci", "None"),
    #                         "complete": True,
    #                         "allele_id": match.get("allele_id"),
    #                         "href": match.get("href"),
    #                     }
    #                 )
    #         elif best_match:
    #             data.append(
    #                 {
    #                     "Isolate": response.get("Isolate", "None"),
    #                     "Loci": response.get("Loci", "None"),
    #                     "complete": False,
    #                     "allele_id": best_match.get("allele_id"),
    #                     "identity": best_match.get("identity"),
    #                     "alignment": best_match.get("alignment"),
    #                 }
    #             )
    #     df = pd.DataFrame(data)
    #     return df

    # # @cache
    # def get_loci_status(params: FetchConsensusDataParams):
    #     response_list = []
    #     for bc in mo.status.progress_bar(params.sample_names):
    #         logger.info(f"Fetching alignment details from {params.selected_db}")
    #         logger.debug(bc)  # Changed to debug
    #         consensus = (
    #             params.path
    #             / str(bc)
    #             # / f"{bc}_extended_ref_postprocessing"
    #             / "consensus.fasta"
    #         )
    #         logging.debug(f"Fetching data for {consensus}")
    #         fa = SeqIO.parse(consensus, "fasta")

    #         for seq in fa:
    #             _sequence_data = {"sequence": f"{seq.seq}"}
    #             loci = seq.id.split("_")[0]
    #             _api_url = f"{selected_db}/loci/{loci}/sequence"
    #             try:
    #                 _response = requests.post(
    #                     _api_url,
    #                     data=json.dumps(_sequence_data),
    #                     params=json.dumps(params.request_params),
    #                 )
    #                 _text = _response.json()
    #                 logger.debug(_text)  # Changed to debug
    #                 _text["Isolate"] = bc
    #                 _text["Loci"] = loci
    #                 response_list.append(_text)
    #             except Exception as e:
    #                 logging.error(
    #                     f"Alignment details failed for {bc}"
    #                 )  # Fixed variable name

    #     return response_list

    # def get_loci_status_single(params: FetchConsensusDataParams):
    #     fa = SeqIO.parse(consensus, "fasta")
    #     params = {"details": "true", "base64": "true"}
    #     list_loci = []
    #     list_isolates = []
    #     for seq in fa:
    #         list_isolates.append(sample_name)
    #         _sequence_data = {"sequence": f"{seq.seq}"}
    #         loci = seq.id.split("_")[0]
    #         list_loci.append(loci)
    #         _api_url = f"{selected_db}/loci/{loci}/sequence"

    #         _response = requests.post(
    #             _api_url,
    #             data=json.dumps(_sequence_data),
    #             params=params,  # Corrected to use params
    #         )
    #         _text = _response.json()
    #         logger.debug(_text)  # Changed to debug
    #         logger.info(f"Processed loci: {_loci}")

    #     return list_loci, list_isolates

    # # @cache
    # async def fetch_consensus_data(params: FetchConsensusDataParams):
    #     st_status = {"Isolate": [], "ST": []}
    #     for bc in mo.status.progress_bar(params.sample_names):
    #         consensus_path = (
    #             params.path
    #             / bc
    #             # / f"{bc}_extended_ref_postprocessing"
    #             / "consensus.fasta"
    #         )
    #         if not consensus_path.exists():
    #             logger.error(f"File {consensus_path} does not exist.")
    #             continue

    #         with open(consensus_path, "r") as x:
    #             _fasta = x.read()

    #         _sequence_data = {
    #             "base64": "true",
    #             "sequence": base64.b64encode(_fasta.encode()).decode(),
    #         }
    #         st_status["ST"].append("Unknown")
    #         st_status["Isolate"].append(bc)
    #         try:
    #             async with aiohttp.ClientSession() as session:
    #                 async with (
    #                     session.post(
    #                         # f"{params.selected_db}/sequence",
    #                         f"{params.selected_db}/schemes/3/sequence",  # TODO add params for schemes High priority
    #                         data=json.dumps(_sequence_data),
    #                     ) as response
    #                 ):
    #                     _matches = await response.json()
    #                     print(_matches)
    #                     logger.debug(
    #                         f"ST status for {bc}: {st_status['Isolate']}"
    #                     )  # Changed to debug
    #                     logger.debug(f"{_matches}")
    #                     if _matches.get("fields"):
    #                         fields = _matches.get("fields")
    #                         st_status[
    #                             "ST"
    #                         ].pop()  # remove the default "Unknow if the ST is found"
    #                         st_status["ST"].append(fields.get("ST", "None"))
    #                         logger.info(
    #                             f"ST profile solved for barcode{bc}: ST{_matches.get('fields')['ST']}"
    #                         )
    #                     if _matches.get("exact_matches"):
    #                         logger.info(
    #                             f"Number of fully identified alleles: {len(_matches['exact_matches'])}"
    #                         )
    #                     logging.debug(_matches)  # Changed to debug

    #         except Exception as r:
    #             logger.error(f"Request problem for sample {bc}: {r}")

    #     return st_status

    # def verify_connection(BASEURL):
    #     ok = requests.get(url=BASEURL)
    #     if ok.status_code == 200:
    #         return ok.json()
    #     else:
    #         return False

    # def get_available_dbs(BASEURL):
    #     "Return a list of all available DBs (e.g E.Coli, Dipht, Clostridium...) available"
    #     responses = verify_connection(BASEURL)
    #     if responses:
    #         try:
    #             return {
    #                 response["name"]: response["databases"][1]["href"]
    #                 for response in responses
    #                 if len(response["databases"]) > 1
    #             }
    #         except Exception as e:
    #             logging.error(f"Some field missing in reponse {e}")

    # class FetchConsensusDataParams(BaseModel):
    #     selected_db: str
    #     path: Path
    #     sample_names: list[str]
    #     request_params: Dict
    #     mlst_genes_path: str
    #     fastq_files: str
    #     postprocessing: bool
    return


@app.cell
def _():
    import marimo as mo
    import requests
    import pandas as pd
    from pathlib import Path

    return Path, mo, pd, requests


@app.cell
def _():
    from src.config import load_params_from_toml, DorasParams
    from src.bigsdb_tools import query_consensus, get_available_dbs

    return (
        DorasParams,
        get_available_dbs,
        load_params_from_toml,
        query_consensus,
    )


@app.cell
def _(mo, ui_experiment_name, ui_genome_size):
    mo.hstack([ui_experiment_name, ui_genome_size])
    return


@app.cell
def _(barcodes_ui, mo):
    # UI elements for manual DORAS parameters
    ui_genome_size = mo.ui.number(
        label="Genome Size (estimated size of the genome of the bacterial species)",
        value=5_000_000,
        start=0,
        stop=20 * 10**6,
        step=1 * 10**5,
    )
    # Reuse existing base API selector

    # Available DBs based on selected base API
    # Reuse existing manual DB input
    # Reuse existing scheme selector

    # Path to
    ui_experiment_name = mo.ui.text(value="DORAS").form(
        label="experiment name (used for the output dir of DORAS)",
        submit_button_label="Confirm experiment name",
    )

    # Sample names from directory browser
    def get_sample_names():
        if barcodes_ui.value:
            return [barcodes_ui.name(i) for i in range(len(barcodes_ui.value))]
        return []

    # Request params UI
    ui_details_switch = mo.ui.switch(label="Include details", value=True)
    ui_base64_switch = mo.ui.switch(label="Include base64", value=True)

    # FASTQ files subfolder
    ui_fastq_path = mo.ui.file_browser(
        label="##FASTQ subfolder (location of the raw fastqs)",
        selection_mode="directory",
    )

    # Postprocessing switch
    ui_postprocess = mo.ui.switch(label="Postprocessing", value=False)

    # MLST genes path
    ui_mlst_path = mo.ui.file_browser(
        label="##MLST genes path (initial reference that will used to grow the extended reference)",
        initial_path="base_refs/",
        selection_mode="file",
    )
    return (
        ui_base64_switch,
        ui_details_switch,
        ui_experiment_name,
        ui_fastq_path,
        ui_genome_size,
        ui_mlst_path,
        ui_postprocess,
    )


@app.cell
def _(Path, mo, ui_experiment_name):
    ui_output_dir = mo.ui.file_browser(
        label=f"##Output directory (where the output dir will be created) (default: currrent dir, name {ui_experiment_name.value}",
        selection_mode="directory",
        initial_path=Path.cwd(),
    ).form(submit_button_label="Confirm output dir")
    return (ui_output_dir,)


@app.cell
def _(ui_output_dir):
    ui_output_dir.element.path()
    return


@app.cell
def _(mo, ui_output_dir):
    mo.md(f"""
    Output directory is {ui_output_dir.element.path()}
    """)
    return


@app.cell
def _(mo, ui_experiment_name, ui_fastq_path, ui_mlst_path, ui_output_dir):
    # ui_pubmlst = mo.vstack([select_mlst_api,ui_genome_size,ui_db_dropdown,ui_scheme_select])
    mo.stop(
        ui_experiment_name.value is None,
        output=mo.md("Confirm the experiment name"),
    )
    ui_paths = mo.vstack([ui_fastq_path, ui_mlst_path, ui_output_dir], gap=2)

    mo.hstack([ui_paths], align="start")
    return


@app.cell
def _(mo, ui_fastq_path):
    mo.stop(
        ui_fastq_path.path() is None,
        output=mo.md("##Enter a valid raw fastq path"),
    )
    for i, k in enumerate(ui_fastq_path.path().walk()):
        if i == 0:
            sample_names_list_options = k[1]
    return (sample_names_list_options,)


@app.cell
def _(mo, sample_names_list_options):
    _s = mo.ui.multiselect(
        label="## Select which barcodes contain the relevant samples",
        options=sample_names_list_options,
    )
    selected_samples_ui = _s.form(
        label="Confirm the barcodes",
        submit_button_label="Confirm",
        clear_on_submit=False,
        clear_button_label="OK",
    )
    selected_samples_ui
    return (selected_samples_ui,)


@app.cell
def _(
    DorasParams,
    Path,
    db_manual_input,
    mo,
    select_mlst_api,
    selected_db_key,
    selected_samples_ui,
    ui_base64_switch,
    ui_details_switch,
    ui_experiment_name,
    ui_fastq_path,
    ui_genome_size,
    ui_mlst_path,
    ui_output_dir,
    ui_postprocess,
    url,
):
    # Selected DB: manual if provided else dropdown
    mo.stop(selected_samples_ui.value is None, output=mo.md("Confirm the barcodes"))

    if db_manual_input.value:
        manual_db = db_manual_input.value.strip()
        selected_db_val = manual_db
    else:
        selected_db_val = selected_db_key.value  # ui_db_dropdown defined in UI cell

    # Sample names: parse from text area
    # Request params
    request_params_val = {
        "details": "true" if ui_details_switch.value else "false",
        "base64": "true" if ui_base64_switch.value else "false",
    }

    # Build DorasParams
    params_doras_manual = DorasParams(
        experiment_name=ui_experiment_name.value,
        genome_size=ui_genome_size.value,
        selected_db=selected_db_val,
        base_api=select_mlst_api.value,
        url=url,
        scheme=url,
        sample_names=selected_samples_ui.value,
        request_params=request_params_val,
        fastq_files_path=str(ui_fastq_path.path()),
        output_dir_path=str(Path(ui_output_dir.element.path())),
        postprocessing=ui_postprocess.value,
        mlst_genes_path=str(ui_mlst_path.path()),
    )
    # display the object

    # List of selected directories from FASTQ path browser
    return (params_doras_manual,)


@app.cell
def _(params_doras_manual):
    params_doras_manual.dict(exclude_defaults=False)
    return


@app.cell
def _(mo):
    ui_show_default_toml = mo.ui.switch(label="Show default parameters of the toml")
    return


@app.cell
def _():
    from src.config import save_params_to_toml

    return (save_params_to_toml,)


@app.cell
def _(mo):
    ui_save_toml = mo.ui.run_button(kind="warn", label="Save toml file")
    ui_save_toml
    return (ui_save_toml,)


@app.cell
def _(
    Path,
    mo,
    params_doras_manual,
    save_params_to_toml,
    ui_experiment_name,
    ui_save_toml,
):
    mo.stop(not ui_save_toml.value, output=mo.md("#Click save toml file first"))
    save_params_to_toml(
        params=params_doras_manual,
        file_path=f"{Path(params_doras_manual.output_dir_path) / ui_experiment_name.value}.toml",
    )
    return


@app.cell
def _(mo):
    ##RUn
    from src.config import Phase

    phase = mo.ui.dropdown(
        options=[Phase.EXTENSION, Phase.QUERY_ST], value=Phase.EXTENSION
    )
    phase
    return (phase,)


@app.cell
def _(Path, mo, params_doras_manual, phase):
    from src.doras_manager import DorasManager
    import asyncio

    mlst_genes_path = params_doras_manual.mlst_genes_path
    _test_folder = params_doras_manual.output_dir_path
    if params_doras_manual.test_mode:
        _test_folder = Path(params_doras_manual.output_dir_path)
    _size = 3000
    params_doras_manual.force_finalize = params_doras_manual.force_finalize

    manager = DorasManager(
        mlst_ref_genes_path=mlst_genes_path,
        barcodes=params_doras_manual.sample_names,
        folder_path=params_doras_manual.output_dir_path,
        mapq=params_doras_manual.min_map_quality,
        overwrite=False,
        target_extension_size=_size,
        params=params_doras_manual,
        clean_up=False,
        test_mode=params_doras_manual.test_mode,  # TODO add params for TEST_MODE
        phase=phase.value,
    )

    ui_run_doras = mo.ui.run_button(label="Run DORAS")
    return asyncio, manager, ui_run_doras


@app.cell
def _(ui_run_doras):
    ui_run_doras
    return


@app.cell
async def _(asyncio, manager, mo, ui_run_doras):
    mo.stop(not ui_run_doras.value, output="Press Run DORAS")

    processor_task = asyncio.create_task(manager.main())
    await processor_task
    return


@app.cell
def _(manager):
    print(manager.sts)
    return


@app.cell
def _(mo):
    mo.md("""
    Functions
    """)
    return


@app.cell
def code_plots(alt, pd):
    def create_heatmap_plot_ST(
        df: pd.DataFrame, x_column: str = "ST", y_column: str = "Isolate"
    ) -> alt.Chart:
        # Configure common options. We specify the aggregation
        # as a transform here so we can reuse it in both layers.
        base = alt.Chart(df).encode(
            alt.X(f"{x_column}:O"),
            alt.Y(f"{y_column}:O"),
        )

        color = (
            alt.when(alt.datum.ST != "Unknown")
            .then(alt.value("green"))
            .otherwise(alt.value("red"))
        )

        # Configure heatmap
        heatmap = base.mark_rect().encode(color=color)

        # Configure text
        _text = base.mark_text(baseline="middle").encode(
            alt.Text(f"{x_column}"), color=alt.value("white")
        )

        # Draw the chart
        return heatmap + _text

    def create_heatmap_plot(
        df: pd.DataFrame, x_column: str = "Loci", y_column: str = "Isolate"
    ) -> alt.Chart:
        # Configure common options. We specify the aggregation
        # as a transform here so we can reuse it in both layers.
        base = alt.Chart(df).encode(
            alt.X(f"{x_column}:O"),
            alt.Y(f"{y_column}:O"),
        )

        # Configure heatmap
        heatmap = base.mark_rect().encode(
            alt.Color("complete").scale(scheme="viridis").title("Complete or not")
        )
        color = (
            alt.when(alt.datum.complete == True)
            .then(alt.value("green"))
            .otherwise(alt.value("red"))
        )

        # Configure text
        if "identity" in df.columns:
            _text = (
                base.mark_text(baseline="middle")
                .encode(alt.Text("identity", format=".0f"), color=color)
                .transform_calculate(
                    identity="datum.identity != null ? datum.identity : ''"
                )
            )
        else:
            _text = base.mark_text(baseline="middle").encode(
                alt.Text("complete:N"), color=color
            )

        # Draw the chart
        return heatmap + _text

    # wgs_plot = create_heatmap_plot(df_loci_wgs)
    # as_plot = create_heatmap_plot(df_loci_as)
    # merged_wgs_as = create_heatmap_plot(df_loci_merged_wgs_as)
    return create_heatmap_plot, create_heatmap_plot_ST


if __name__ == "__main__":
    app.run()

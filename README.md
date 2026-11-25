# eWHALE bioinformatics workflows

Bioinformatics workflows and analysis scripts used by partners in the multinational [eWHALE project](https://ewhale.eu/) to process and analyze environmental DNA (eDNA) metabarcoding data.

This repository brings together workflows applied to samples from the Azores, Iceland, the French coast, and the Norwegian coast. These studies use complementary genetic markers to investigate marine biodiversity, including biodiversity associated with different whale species.

> **Project status:** Associated publications are currently in preparation. The code, parameters, and documentation presented here are subject to change as analyses and manuscripts are finalized. For reproducibility, record the exact commit used for each analysis.

## Workflows and study regions

| Folder or linked repository | Study region and application | Markers |
| --- | --- | --- |
| [`metabarcoding-pipeline`](metabarcoding-pipeline) | Azores and Iceland: analysis of samples collected around various whale species | COI, 18S rRNA, 16S rRNA |
| [`INRAE_pipeline_metabarcoding`](INRAE_pipeline_metabarcoding) | French coastal biodiversity monitoring in the English Channel | 12S rRNA |
| [`IMR_pipeline_metabarcoding`](IMR_pipeline_metabarcoding) | Norwegian coastal biodiversity monitoring | 12S rRNA, COI |

The workflows reflect the assays, sequencing data, and analytical requirements of the respective studies. Processing steps, filtering thresholds, reference databases, and output formats may therefore differ among workflows. Consult the relevant scripts and documentation before applying a workflow to a new dataset.

### Azores and Iceland: linked metabarcoding pipeline

The [`metabarcoding-pipeline`](https://github.com/laurenrodriguez/metabarcoding-pipeline) repository is maintained separately and included here as a **Git submodule**. It documents marker-specific sequence processing, taxonomic assignment, control-based filtering, replicate and sampling-event aggregation, and downstream ecological analyses.

See its [README](metabarcoding-pipeline/README.md) and [run-order documentation](metabarcoding-pipeline/docs/RUN_ORDER.md) for requirements, input formats, and analysis steps. Study-specific settings should be checked against the relevant dataset and publication.

### INRAE: French coastal metabarcoding

The [`INRAE_pipeline_metabarcoding`](INRAE_pipeline_metabarcoding) folder contains the numbered `Part0`–`Part4` shell scripts used for the 12S workflow. Review the scripts in numerical order and adapt paths, input files, reference databases, and software settings to the local computing environment.

### IMR: Norwegian coastal metabarcoding

The [`IMR_pipeline_metabarcoding`](IMR_pipeline_metabarcoding) folder contains:

- [`metabarcoding_pipeline.sh`](IMR_pipeline_metabarcoding/metabarcoding_pipeline.sh): sequence filtering, dereplication, chimera removal, clustering, curation, and taxonomic-assignment steps.
- [`config.ini`](IMR_pipeline_metabarcoding/config.ini): input files, working directories, marker selection, thread count, clustering settings, and reference-database paths.

The IMR configuration uses `CO1` to denote the COI marker. The current script contains an explicit stopping point before taxonomic assignment with `ecotag`; that stage requires review, including configuration of its input and output variables, before execution.

## Getting started

To download this repository **including the linked metabarcoding pipeline**, clone it with:

```bash
git clone --recurse-submodules https://github.com/eWHALE-DNA/bioinformatics.git
cd bioinformatics
```

If you have already cloned the repository, retrieve the submodule with:

```bash
git submodule update --init --recursive
```

These commands retrieve the submodule version recorded by this repository. The separately maintained pipeline may have newer commits; retain the recorded version when reproducing a specific analysis.

Before running a workflow:

1. Review its documentation, scripts, and required software.
2. Adapt input/output paths, sample metadata, primer and marker settings, and computing resources.
3. Obtain the required reference databases and record their versions or download dates.
4. Check that filtering and taxonomic-assignment settings are appropriate for the assay and study.
5. Retain the code version, configuration files, software versions, and execution logs used for the analysis.

These are research workflows with study-specific settings. Some stages require manual preparation or review, and the repository does not provide a single command that runs all partner workflows.

## Data availability

Raw sequence data generated from samples collected in the **Azores in 2024** are archived at NCBI under BioProject accession **[PRJNA1522437](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1522437)**.

Sequence data from other study regions will be deposited as the associated studies are finalized. Accession numbers and links to associated publications will be added here as they become available.

Raw sequence data and reference databases should be obtained separately from the appropriate archives or database providers.

## Citation and acknowledgment

If you use or adapt these workflows, please acknowledge the eWHALE project and the contributing workflow authors, and cite the associated study publications when available.

Please also cite the software used in your analysis. References are collected in [`Software Citations.md`](Software%20Citations.md), with additional information in the individual workflows. Record the repository URL and exact commit used; for the linked pipeline, record its own commit as well.

## Questions and contributions

For questions about the workflows, documentation, or reproducibility, please [open an issue](https://github.com/eWHALE-DNA/bioinformatics/issues) and identify the relevant folder, marker, and code version. Include relevant commands and error messages when reporting a problem.

## License

Code in this repository is distributed under the [MIT License](LICENSE). Consult the linked pipeline and third-party software for their respective licenses. Sequence datasets and reference databases remain subject to the terms specified by their providers.

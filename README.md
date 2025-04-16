<p align="center">
  <img width="auto" height="auto" src="https://github.com/Integrative-Transcriptomics/PTMVision/blob/main/app/ptmvision/static/resources/logo.png">
</p>

# About

<p align="justify">
PTMVision (Interactive Visualization of Post Translational Modifications) is an interactive web platform for effortless exploratory visual analysis of post translational modifications (PTMs) obtained from mass spectrometry (MS) data. It enables researchers to comprehend the intricate landscape of PTMs in order to unravel possible mechanisms related to cellular processes and diseases.
</p>

The webserver is available at [https://ptmvision-tuevis.cs.uni-tuebingen.de/](https://ptmvision-tuevis.cs.uni-tuebingen.de/).

# Key Features

- Variety of overview and summary visualizations that help identify patterns and trends in PTM data.
- Visualizations are dynamically linked and support features to highlight PTMs of interest.
- Integration of protein structure and annotation data from UniProt provide multidimensional view. This improves the interpretability of the data is significantly.
- Unique emphasis on interactions of residues in close contact within a single protein and associated modifications. This allows to study the dynamic interplay between PTMs and their putative effects on protein folding and function.
- Responsive and intuitive interface that allows users to effortlessly navigate the website and through their data.
- Supports high degree of data provenance by allowing various input formats and an easy export of session data to continue analysis without repeated processing input data.

# Supported File Formats

Can't find your favorite search engine? [Tell us!](https://github.com/Integrative-Transcriptomics/PTMVision/issues/new)

| **Search Engine**                         | **Version Tested**                              | **Columns Used**                                                                                                                              | **Postprocessing Applied by PTMVision**                                                                                                                                                                                                                |
| ----------------------------------------- | ----------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| ionbot                                    | v0.10.0, v0.11.0                                | uniprot_id, unexpected_modification, position                                                                                                 |                                                                                                                                                                                                                                                        |
| MSFragger (via FragPipe with PTMShepherd) | <br> MSFragger v4.0,<br> PTMShepherd v2.0.6<br> | <br> Modified Peptide, Peptide, Protein ID, Assigned<br> Modifications, Observed Modifications, MSFragger<br> Localization, Protein Start<br> | <br> Filter ambiguous localisations; Filter ambiguous<br> modifications; Filter combinations of modifications; Map<br> mass shifts and UniMod descriptions to UniMod IDs; Map PTM<br> position onto protein; Retrieve UniMod classification; \*\*.<br> |
| Sage                                      | v0.13.1                                         | \*                                                                                                                                            | <br> 1% FDR filtering (PSM level); Remove decoys; Remove peptides<br> matching to more than one protein; Map mass shifts to UniMod<br> IDs; Map PTM position onto protein; Mass shifts are not<br> localised, only variable mods are used.<br>         |
| MaxQuant                                  | <br> v2.4.13.0<br>                              | \*                                                                                                                                            | <br> Remove peptides matching to more than one protein; Map<br> MaxQuant modification names to UniMod names; Map PTM<br> position onto protein.<br>                                                                                                    |
| Spectronaut (PTM Site Report)             | v9                                              | <br> PTM.ProteinId, PTM.SiteLocation, PTM.ModificationTitle,<br> PTM.SiteAA<br>                                                               | <br> Modification mass shift and UniMod classification inferred<br> from the UniMod ID (mapped from PTM.ModificationTitle) and<br> the protein sequence.<br>                                                                                           |

<sup>\* See https://psm-utils.readthedocs.io/en/v0.6.0/ for details.</sup>

<sup>\*\* For FragPipe, localisations for unexpected modifications are derived from the "MSFragger Localization" column. This localisation might be only marginally better than the next best site, so caution is advised. As soon as the next MSFragger version is released, which reports localisation scores as well, PTMVision will include a way to filter modifications based on the scores.</sup>

## Custom Plain CSV

<p align="justify">
We support a plain CSV format to provide users the possibility to transfer individually compiled PTM data to PTMVision if no other format can be used. The following table provides an overview of the required columns and their content. In essence, the format provides the description of a modification at one position of a protein per line as comma separated values - thus, one file can contain information on several proteins. If <em>mass_shift</em> or <em>classification</em> are not provided, they are tried to be inferred from the UniMod ID and the protein sequence.
</p>

| **Column Name**         | **uniprot_id**                        | **position**                                                                           | **modification_unimod_name**                                                         | **modification_unimod_id**                                                         | **mass_shift (optional)**                                                                                            | **classification (optional)**                                                                                                                  |
| ----------------------- | ------------------------------------- | -------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| **Content Description** | The UniProt Accession of the protein. | <br> The position of the modification in the protein sequence<br> (1-based index).<br> | <br> The UniMod name of the modification observed on the<br> specified position.<br> | <br> The UniMod ID of the modification observed on the specified<br> position.<br> | <br> The mass shift of the modification in Dalton. This will be<br> inferred from the UniMod ID if not provided.<br> | <br> The UniMod classification of the modification. This will be<br> inferred from the UniMod ID and modified residue if not<br> provided.<br> |
| **Example Entry**       | P04075                                | 5                                                                                      | carbamidomethyl                                                                      | 4                                                                                  | 57.021464                                                                                                            | Artefact                                                                                                                                       |

<p align="justify">
Important note: Since UniMod classifications are amino acid specific, inferring the classification requires querying the protein sequence from UniProt. If the positions in the csv come from searches using old(er) sequences/entries, it is possible that the amino acid sequence changed in the meantime, which could lead to modifications being misplaced.
</p>

An example file for the plain csv format can be found [here](https://github.com/Integrative-Transcriptomics/PTMVision/blob/2b8dfcfb281150d4c634427ffc1d03d5873925fd/app/ptmvision/static/resources/ptmvision.plain-format.template.csv).

## Run PTMVision Locally

If required, PTMVision can be executed locally as a Docker container. After the repository has been cloned, simply run `docker build --tag ptmvision app` and `docker run -dp 127.0.0.1:5001:5001 ptmvision` in the root directory of the project. We currently do not provide a Docker image.

<p align="center">
  <img width="256" height="128" src="https://github.com/Integrative-Transcriptomics/PTMVision/blob/main/ptmvis/ptmvis/static/resources/logo.png">
</p>

# About

PTMVision is a webserver for analysis and exploration of PTM data derived from common proteomics search engines.
The webserver is available at [https://ptmvision-tuevis.cs.uni-tuebingen.de/](https://ptmvision-tuevis.cs.uni-tuebingen.de/).

# Supported file formats

Can't find your favorite search engine? [Tell us!](https://github.com/Integrative-Transcriptomics/PTMVision/issues/new)

|                                                                                                                   | file                     | version tested                      | Columns used                                                                                                                                                                  | postprocessing                                                                                                                                 | comment                                                     |
|-------------------------------------------------------------------------------------------------------------------|--------------------------|-------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------|
| ionbot                                                                                                            | modifications.csv | v0.10.0                             | "uniprot_id", "unexpected_modification",  "position"                                                                                                                          | -                                                                                                                                              |                                                             |
| Sage                                                                                                              | <run>.sage.tsv           | v0.13.1                             | *                                                                                                                                                                             | <ul><li>1% FDR Filtering (PSM level)</li><li>Remove decoys</li><li>Remove peptides matching to >1 protein</li>  <li>Map mass shifts to UniMod IDs Map PTM position onto protein</li></ul>   | mass shifts are not localised,  only variable mods are used |
| MSFragger                                                                                                         | psm.tsv                  | MSFragger v4.0 & PTMShepherd v2.0.6 | "Modified Peptide",       "Peptide",         "Protein ID", "Assigned Modifications", "Observed Modifications",            "MSFragger Localization",           "Protein Start" | <ul><li>Map mass shifts to UniMod IDs</li><li> Map PTM position onto protein </li></ul>                                                                                    | currently only variable mods are used                       |
| MaxQuant                                                                                                          | msms.txt                 | v2.4.13.0                           | *                                                                                                                                                                             | <ul><li>Remove peptides matching to >1 protein </li><li>Map MaxQuant modification names to UniMod names</li> <li>Map PTM position onto protein </li></ul>                          |                                                             |
| csv                                                                                                               | <name>.csv               | -                                   | uniprot_id, position, modification_unimod_name, modification_unimod_id, [optional: classification], [optional: annotation]                                                    | -                                                                                                                                              |                                                             |
| additional psm-utils  supported engines  (e.g. OpenMS,  MS Amanda,  Percolator,  X!Tandem,  Proteome Discoverer)  | *                        | -                                   | *                                                                                                                                                                             | <ul><li>1% FDR Filtering (PSM level)</li> <li>Removing decoys</li> <li>Removing peptides matching to >1 protein</li></ul>                                                           |                                                             |

\* see https://psm-utils.readthedocs.io/en/v0.6.0/ for details

An example file for the plain csv format can be found [here](https://github.com/Integrative-Transcriptomics/PTMVision/tree/psm_utils_parsing/ptmvis/ptmvis/backend/example_data).

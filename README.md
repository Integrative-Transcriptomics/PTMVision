<p align="center">
  <img width="auto" height="auto" src="https://github.com/Integrative-Transcriptomics/PTMVision/blob/main/app/ptmvision/static/resources/logo.png">
</p>

# PTMVision

**PTMVision** is an interactive web-based platform for the visual exploration of post-translational modifications (PTMs) identified in mass spectrometry-based proteomics data. It enables researchers to interpret the complex PTM landscape, supporting insights into molecular mechanisms of biological processes and disease.

🔗 **Web server**: [https://ptmvision-tuevis.cs.uni-tuebingen.de](https://ptmvision-tuevis.cs.uni-tuebingen.de)

---

## ✨ Features

- Rich **interactive visualizations** for sample-level and protein-level PTM data
- Linked plots for seamless navigation and exploration
- **3D structure integration** and contact map analysis for structural context
- Highlights residue-residue PTM interactions in close spatial proximity
- **Annotation integration** from UniProt and UniMod for interpretability
- Support for multiple **search engine output formats** and a custom CSV format
- **Session export/import** for reproducibility and collaboration
- Intuitive UI accessible to users without programming expertise

---

## 📁 Supported Input Formats

PTMVision accepts output from several widely-used search engines:

| **Search Engine**                         | **Tested Versions**                          | **Key Columns**                                                                                           | **Processing Steps**                                                                                                                                                                                                            |
|------------------------------------------|----------------------------------------------|------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **ionbot**                                | v0.10.0, v0.11.0                             | `uniprot_id`, `unexpected_modification`, `position`                                                        | Direct parsing of modification information                                                                                                                                                |
| **MSFragger (via FragPipe + PTMShepherd)** | MSFragger v4.0, PTMShepherd v2.0.6           | `Modified Peptide`, `Assigned Modifications`, `Localization`, `Protein ID`, `Protein Start`               | Filtering ambiguous or multi-localized sites; mass shift mapping to UniMod; UniProt sequence alignment; UniMod classification                                                          |
| **Sage**                                   | v0.13.1                                      | *See [psm-utils](https://psm-utils.readthedocs.io/en/v0.6.0/)*                                             | FDR filtering; decoy removal; UniMod mapping; PTM localization from variable mods only                                                                                                   |
| **MaxQuant**                               | v2.4.13.0                                    | *Search result columns used for mapping*                                                                   | Mapping MaxQuant-specific names to UniMod; resolving multi-mapping peptides; UniProt mapping                                                                                             |
| **Spectronaut (PTM Site Report)**          | v9                                           | `PTM.ProteinId`, `PTM.SiteLocation`, `PTM.ModificationTitle`, `PTM.SiteAA`                                 | Infers mass shifts and UniMod classifications from reported data and UniProt sequence                                                                                                   |

> **Note**: For FragPipe-derived data, localization confidence is inferred from the `MSFragger Localization` field. Use caution until future versions provide dedicated localization scores.

---

## 📄 Custom Plain CSV Format

Users can submit manually compiled PTM data using a structured CSV format when other formats are unavailable.

### Required Columns

| Column Name                 | Description                                                                                      |
|----------------------------|--------------------------------------------------------------------------------------------------|
| `uniprot_id`               | UniProt accession of the protein                                                                 |
| `position`                 | 1-based index of the modified residue                                                            |
| `modification_unimod_name`| Name of the UniMod modification                                                                  |
| `modification_unimod_id`  | UniMod ID                                                                                         |
| `mass_shift` *(optional)* | Monoisotopic mass shift in Daltons (inferred from UniMod ID if omitted)                         |
| `classification` *(optional)* | UniMod classification (inferred from ID and sequence context if omitted)                |

**Example entry**:
```
P04075,5,carbamidomethyl,4,57.021464,Artefact
```

> ⚠ **Note**: Classification relies on up-to-date UniProt sequences. Sequence version mismatches may lead to misaligned positions.

📄 [Download CSV template](https://github.com/Integrative-Transcriptomics/PTMVision/blob/main/app/ptmvision/static/resources/ptmvision.plain-format.template.csv)

---

## 🐳 Running Locally via Docker

To deploy PTMVision locally:

```bash
# Clone repository
git clone https://github.com/Integrative-Transcriptomics/PTMVision.git
cd PTMVision

# Build Docker image
docker build --tag ptmvision app

# Run locally
docker run -dp 127.0.0.1:5001:5001 ptmvision
```

> Note: Prebuilt Docker images are not currently provided.

---

## 🙋 Need Help?

- 📖 Read the [Usage Guide](https://ptmvision-tuevis.cs.uni-tuebingen.de/static/usage.html)
- 🐛 Found an issue or missing format? [Open a GitHub issue](https://github.com/Integrative-Transcriptomics/PTMVision/issues/new)

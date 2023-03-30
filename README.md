# Import Pan UKBB GSR data into AnVIL from Broad Institute

Visit the Pan-UK Biobank phenotype manifest [Google Sheets webpage](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288) to identify one or more phenotypes of interest for your analysis. Enter at least one phenocode from the Google Sheets webpage (e.g., "30600" is the phenocode corresponding to Albumin) in the JSON template. Then use the "download_pan_ukbb" workflow in AnVIL to import the data into the syntax of the PRIMED consortium. Finally, perform the standard GSR data validation workflows to import the data tables and confirm that requirements are met.

<br/>

In the JSON input to the AnVIL workflow, you can set <code>"download_pan_ukbb.population": ["all_available"]</code> to download the GSR data for all populations having data for that phenotype. Alternatively, you can choose specific populations of interest using a subset of the list <code>"download_pan_ukbb.population": ["meta", "metaHQ", "AFR", "AMR", "CSA", "EAS", "EUR", "MID"]</code>. If you choose a specific population that is unavailable for a particular phenotype, then no data will download for that population. If you specify multiple phenotypes, then the populations you specify will apply for each phenotype.

<br/>

Note that the population abbreviations from the [Pan UKBB website](https://pan.ukbb.broadinstitute.org/docs/technical-overview) are described as follows:
| Abbreviation | Description                   |
| ------------ | ----------------------------  |
| meta         | meta-analysis of any quality  |
| metaHQ       | meta-analysis of high quality |
| AFR          | African ancestry              |
| AMR          | Admixed American ancestry     |
| CSA          | Central/South Asian ancestry  |
| EAS          | East Asian ancestry           |
| EUR          | European ancestry             |
| MID          | Middle Eastern ancestry       |

<br/>

Once the data is downloaded to your AnVIL workspace, you can begin the next steps of performing the data validation workflows. As of March 2023, the data will pass the GSR data validation workflows, however this should be confirmed manually; the Pan UKBB GSR raw data syntax may change over time, as it is stored in Amazon Web Service.

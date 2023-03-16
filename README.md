# Import Pan UKBB GSR data into AnVIL from Broad Institute

Visit the Pan-UK Biobank phenotype manifest [Google Sheets webpage](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288) to identify one or more phenotypes of interest for your analysis. Enter at least one phenocode from the Google Sheets webpage (e.g., "30600" is the phenocode corresponding to Albumin) in the JSON template. Then use the "download_pan_ukbb" workflow in AnVIL to import the data into the syntax of the PRIMED consortium. Finally, perform the standard GSR data validation workflows to import the data tables and confirm that requirements are met.

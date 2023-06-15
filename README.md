# Import Pan UKBB GSR data into AnVIL from Broad Institute

Visit the Pan-UK Biobank phenotype manifest [Google Sheets webpage](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=1450719288) to identify one or more phenotypes of interest for your analysis. PRIMED consortium members can access an internal spreadsheet, which contains information about which traits were already downloaded. Then use the "download_pan_ukbb" workflow in AnVIL to import the data into the syntax of the PRIMED consortium.

<br/>

The "download_pan_ukbb" AnVIL workflow has two required inputs and seven optional inputs:

### (Required)

#### 1. bucket_name

In your AnVIL workspace, go to Dashboard &rarr; Cloud Information &rarr; Bucket Name. Copy this name as a string into the JSON template.

#### 2. phenocode

In the phenotype manifest, most traits can be uniquely identified using the "phenocode" column alone. If a trait you wish to download has a unique phenocode, then you can leave the "coding" and "modifier" inputs empty, which the WDL script will label as, "NO", by default.

### (Optional)

#### 4. coding

If a trait you wish to download does not have a unique phenocode, then you may need "coding", "modifier", or both to uniquely identify a trait. For example, to download one of the colorectal cancer traits, you need to set both: <br/>
<code>"download_pan_ukbb.phenocode": ["20001"]</code> <br/>
<code>"download_pan_ukbb.coding": ["1020"]</code> <br/>
because the particular trait is uniquely identified in the phenotype manifest by this phenocode and coding. Note that we can still leave "modifier" empty because it is not needed to identify this trait.

#### 5. modifier

Similar as above, to download one of the blood pressure traits, you need to set both: <br/>
<code>"download_pan_ukbb.phenocode": ["DBP"]</code> <br/>
<code>"download_pan_ukbb.modifier": ["combined_medadj_raw"]</code> <br/>
because the particular trait is uniquely identified in the phenotype manifest by this phenocode and modifier. Note that we can still leave "coding" empty because it is not needed to identify this trait.

#### 6. population

By default, all available data is downloaded. If you wish to download data for only a particular population, or only a meta analysis or high-quality meta analysis, you can specify which populations should be downloaded using a subset of the following list: <br/>
<code>"download_pan_ukbb.population": ["meta", "metaHQ", "AFR", "AMR", "CSA", "EAS", "EUR", "MID"]</code>.

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

#### 7. conceptID

If you know the concept ID of the trait you are downloading, please specify it here. If you are unsure which concept ID the trait corresponds to, this field can be left blank and it will be labeled "TBD".

#### 8. disk_gb

If you encounter an error regarding disk space, i.e. the instance does not have enough memory to save the raw and processed data, then you can manually increase disk space here. For most purposes, the default will suffice.

#### 9. mem_gb

If you encounter an error regarding memory, i.e. the R script is timing out or crashing without a clear coding error, then you can manually increase memory here to determine if that resolves the issue. For most purposes, the default will suffice. Increasing memory may be especially useful for traits with large sample sizes.

<br/>

### Note about syntax
Currently, the workflow only allows you to download one trait at a time. Therefore, if an input expects an array of strings, please only enter a single string within the array. For example, although <code>"download_pan_ukbb.phenocode": ["30000", "30010"]</code> <br/> is valid syntax for JSON array, you should only enter <code>"download_pan_ukbb.phenocode": ["30000"]</code> <br/> because the R script is currently only written to process one phenotype at a time. We made this decision to keep the workflow job history and output well-organized within AnVIL.

<br/>

### Next steps

Once the data are downloaded to your AnVIL workspace, you can begin the next steps of performing the data validation workflows. As of June 2023, the data will pass the GSR data validation workflows, however this should be confirmed manually. The Pan UKBB GSR raw data syntax may change over time, as it is stored in Amazon Web Service, so the R script supporting the docker image may require updates over time.

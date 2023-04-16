version 1.0

workflow download_pan_ukbb {
    input {
        Array[String]+ phenocode
        Array[String] population = ["all_available"]
        Array[String] conceptID = ["TBD"]
        Int disk_gb = 25
        Int mem_gb = 50
    }
    
    call create {
        input: phenocode = phenocode,
               population = population,
               conceptID = conceptID,
               disk_gb = disk_gb,
               mem_gb = mem_gb
    }
    
    call move {
        input: analysis_table = create.analysis_table,
               data_table = create.data_table,
               file_table = create.file_table,
               phenocode = phenocode,
               disk_gb = disk_gb,
               mem_gb = mem_gb
    }
    
    output {
        Array[File] analysis_table = move.analysis_table
        Array[File] data_table = move.data_table
        Array[File] file_table = move.file_table
    }
    
     meta {
          author: "UW Coordinating Center"
          email: "sdmorris@uw.edu"
    }
}

task create {
    input {
        Array[String] phenocode
        Array[String] population
        Array[String] conceptID
        Int disk_gb
        Int mem_gb
    }
    
    command <<<
        Rscript /usr/local/primed_ukbb_gsr/download_pan_ukbb.R \
            --phenocode ${sep=" " phenocode} \
            --population ${sep=" " population} \
            --conceptID ${sep=" " conceptID}
    >>>
    
    output {
        Array[File] analysis_table = glob("*_analysis.tsv")
        Array[File] file_table = glob("*_file.tsv")
        Array[File] data_table = glob("*_data.tsv.gz")
    }
    
    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.0"
        disks: "local-disk ${disk_gb} SSD"
        memory: "${mem_gb} GB"
    }
}

task move {
    input {
        Array[File] analysis_table
        Array[File] data_table
        Array[File] file_table
        Array[String] phenocode
        Int disk_gb
        Int mem_gb
    }
    
    # The command chunk below was adapted from:
    # https://support.terra.bio/hc/en-us/community/
    # posts/360068067031-Write-cromwell-output-to-
    # its-own-folder-instead-of-the-root-directory
    
    command <<<
        #!/bin/bash
        bucket = fc-bb562a6c-b341-4f67-8016-c36ffd74b988
        while read analysis_table
        do
          x = ${analysis_table%/}
          basename = $(basename $x)
          gsutil -m mv $x gs://${bucket}/UKBB-Data/${sep="" phenocode}/${basename}
        done < <(
        )
    >>>
    
    # gsutil ls -d gs://${bucket}
    
    output {
        Array[File] analysis_table = glob("*_analysis.tsv")
        Array[File] file_table = glob("*_file.tsv")
        Array[File] data_table = glob("*_data.tsv.gz")
    }
    
    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.0"
        disks: "local-disk ${disk_gb} SSD"
        memory: "${mem_gb} GB"
    }
}

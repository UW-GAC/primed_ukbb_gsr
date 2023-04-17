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
        input: analysis_table_in = create.analysis_table,
               data_table_in = create.data_table,
               file_table_in = create.file_table,
               phenocode = phenocode
    }
    
    output {
        Array[File] analysis_table = create.analysis_table
        Array[File] data_table = create.data_table
        Array[File] file_table = create.file_table
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
        Array[String] analysis_table_in
        Array[String] data_table_in
        Array[String] file_table_in
        Array[String] phenocode
    }
    
    # The command chunk below was adapted from:
    #
    # https://support.terra.bio/hc/en-us/community/
    # posts/360068067031-Write-cromwell-output-to-
    # its-own-folder-instead-of-the-root-directory
    #                   and
    # https://bioinformatics.stackexchange.com/
    # questions/18827/how-do-i-use-bash-for-loops
    # -to-loop-through-a-wdl-array
    
    command <<<
        #!/bin/bash
        files=('~{sep="' '" analysis_table_in}' '~{sep="' '" data_table_in}' '~{sep="' '" file_table_in}')
        line=1
        bucket="fc-bb562a6c-b341-4f67-8016-c36ffd74b988"
        for x in ${files[@]}; do
            echo "File"${line}
            ((line+=1))
            fname=$(basename ${x})
            echo ${fname}
            oldpath=${x}
            echo ${oldpath}
            newpath="gs://${bucket}/UKBB-Data/~{sep="AND" phenocode}/${fname}"
            echo ${newpath}
            gsutil -m mv ${oldpath} ${newpath}
        done;
    >>>
        
    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.0"
    }
}

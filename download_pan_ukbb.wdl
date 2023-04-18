version 1.0

workflow download_pan_ukbb {
    input {
        Array[String]+ phenocode
        String+ bucket_name
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
               file_table_in = create.file_table,
               data_table_in = create.data_table,
               bucket_name = bucket_name,
               phenocode = phenocode
    }
    
    output {
        Array[String] analysis_table = move.analysis_table
        Array[String] file_table = move.file_table
        Array[String] data_table = move.data_table
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
        Array[String] file_table_in
        Array[String] data_table_in
        Array[String] phenocode
        String bucket_name
    }
    
    command <<<
        #!/bin/bash
        echo "Beginning bash script..."
        bucket=~{bucket_name}
        files=('~{sep="' '" analysis_table_in}')
        for x in ${files[@]}; do
            fname=$(basename ${x})
            echo ${fname}
            oldpath=${x}
            echo ${oldpath}
            newpath="gs://${bucket}/UKBB-Data/~{sep="AND" phenocode}/${fname}"
            echo ${newpath} >> $analysis_table_in.txt
            gsutil -m mv ${oldpath} ${newpath}
        done;
        files=('~{sep="' '" file_table_in}')
        for x in ${files[@]}; do
            fname=$(basename ${x})
            echo ${fname}
            oldpath=${x}
            echo ${oldpath}
            newpath="gs://${bucket}/UKBB-Data/~{sep="AND" phenocode}/${fname}"
            echo ${newpath} >> $file_table_in.txt
            gsutil -m mv ${oldpath} ${newpath}
        done;
        files=('~{sep="' '" data_table_in}')
        for x in ${files[@]}; do
            fname=$(basename ${x})
            echo ${fname}
            oldpath=${x}
            echo ${oldpath}
            newpath="gs://${bucket}/UKBB-Data/~{sep="AND" phenocode}/${fname}"
            echo ${newpath} >> $data_table_in.txt
            gsutil -m mv ${oldpath} ${newpath}
        done;
    >>>
    
    output {
        Array[String] analysis_table = read_lines("analysis_table_in.txt")
        Array[String] file_table = read_lines("file_table_in.txt")
        Array[String] data_table = read_lines("data_table_in.txt")
    }
    
    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.0"
    }
}

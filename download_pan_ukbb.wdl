version 1.0

workflow download_pan_ukbb {
    input {
        Array[String]+ phenocode
        Array[String] population = ["all_available"]
        Array[String] conceptID = ["TBD"]
        Int disk_gb = 25
        Int mem_gb = 50
    }
    
    call folder {
        input: save = "NULL"
    }
    
    call results {
        input: phenocode = phenocode,
               population = population,
               conceptID = conceptID,
               disk_gb = disk_gb,
               mem_gb = mem_gb
    }
    
    output {
        String file_path = folder.file_path
        Array[File] analysis_table = results.analysis_table
        Array[File] data_table = results.data_table
        Array[File] file_table = results.file_table
    }
    
     meta {
          author: "UW Coordinating Center"
          email: "sdmorris@uw.edu"
    }
}

task folder {
    input {
        String save
    }
    
    command <<<
        Rscript; \
        write.table("", file = "get_filepath.tsv")
    >>>
    
    output {
        String file_path = glob("*get_filepath.tsv")
    }
    
    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.0"
        disks: "local-disk 1 SSD"
        memory: "1 GB"
    }
}

task results {
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
        memory: "~{mem_gb}GB"
    }
}

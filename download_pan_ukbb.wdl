version 1.0

workflow download_pan_ukbb {
    input {
        Array[String]+ phenocode
        Array[String] population = "all_available"
        Int disk_gb = 25
        Int mem_gb = 50
    }

    call results {
        input: phenocode = phenocode,
               population = population,
               disk_gb = disk_gb,
               mem_gb = mem_gb
    }

    output {
        Array[File] analysis_table = results.analysis_table
        Array[File] file_tables = results.file_tables
        Array[File] data_tables = results.data_tables
    }

     meta {
          author: "UW Coordinating Center"
          email: "sdmorris@uw.edu"
    }
}

task results {
    input {
        Array[String] phenocode
        Array[String] population
        Int disk_gb
        Int mem_gb
    }

    command {
        Rscript /usr/local/primed_ukbb_gsr/get_AWS_data.R \
            --phenocode ${sep=" " phenocode} \
            --population ${sep=" " population}
    }

    output {
        Array[File] analysis_table = glob("*_analysis.tsv")
        Array[File] file_tables = glob("*_file.tsv")
        Array[File] data_tables = glob("*_data.tsv.gz")
    }

    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.3"
        disks: "local-disk ${disk_gb} SSD"
        memory: "~{mem_gb}GB"
    }
}

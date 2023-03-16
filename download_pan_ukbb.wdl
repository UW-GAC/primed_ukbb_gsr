version 1.0

workflow download_pan_ukbb {
    input {
        Array[String]+ phenocode
    }

    call results {
        input: phenocode = phenocode
    }

    output {
        Array[File] analysis_table = results.analysis_table
        Array[File] file_tables    = results.file_tables
        Array[File] data_tables    = results.data_tables
    }

     meta {
          author: "UW Coordinating Center"
          email: "sdmorris@uw.edu"
    }
}

task results {
    input {
        Array[String] phenocode
    }

    command {
        Rscript /usr/local/primed_ukbb_gsr/get_AWS_data.R \
            --phenocode ${sep=" " phenocode}
    }

    output {
        Array[File] analysis_table = glob("*_analysis.tsv")
        Array[File] file_tables    = glob("*_file.tsv")
        Array[File] data_tables    = glob("*_data.tsv.gz")
    }

    runtime {
        docker: "uwgac/primed-pan-ukbb:0.1.0"
    }
}

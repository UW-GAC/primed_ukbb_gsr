version 1.0

workflow download_pan_ukbb {
    input {
        Float phenocode
    }

    call results {
        input: phenocode = phenocode
    }

     meta {
          author: "Grant Hopkins"
          email: "TBD"
    }
}

task results{
    input {
        Float phenocode
    }

    command {
        Rscript /usr/local/primed_ukbb_gsr/get_AWS_data.R \
            --phenocode ${phenocode}
    }

    runtime {
        docker: "uwgac/primed-file-checks:0.2.5"
    }
}
rule all:
    input:
        "results/uncertain_merged_all_cnv_annotations_filtered_per_ind_with_clinical.csv",
        "results/merged_per_individual_annotations.csv",
        "results/merged_annotations.csv"


rule download_classifycnv:
    output:
        directory("tools/ClassifyCNV")
    shell:
        """
        mkdir -p tools
        wget -O tools/ClassifyCNV.zip https://github.com/Genotek/ClassifyCNV/archive/refs/heads/master.zip
        unzip tools/ClassifyCNV.zip -d tools
        mv tools/ClassifyCNV-master tools/ClassifyCNV
        rm tools/ClassifyCNV.zip
        """


rule download_cnv_database:
    output:
        "databases/dgvMerged.txt"
    shell:
        """
        wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/dgvMerged.txt.gz
        gunzip -c dgvMerged.txt.gz > {output}
        rm dgvMerged.txt.gz
        """


rule download_brain_database:
    output:
        "databases/brain_database.tsv"
    shell:
        """
        wget -O {output} "https://www.proteinatlas.org/search/tissue_category_rna%3Abrain%3BTissue+enriched%2CGroup+enriched%2CTissue+enhanced+AND+sort_by%3Atissue+specific+score?format=tsv&download=yes"
        """


rule download_cytobands:
    output:
        "databases/cytoBand.txt"
    shell:
        """
        wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
        gunzip -c cytoBand.txt.gz > {output}
        rm cytoBand.txt.gz
        """


rule annotate_vcf:
    input:
        config="configs/config.yaml"
    output:
        annotated_vcf="results/merged_annotations.csv"
    shell:
        """
        python scripts/annotate_vcf.py --config {input.config}
        """


rule map_to_cytobands:
    input:
        annotated_vcf="results/merged_annotations.csv",
        cytoband_file="databases/cytoBand.txt"
    output:
        "results/cytoband_merged_annotations.csv"
    shell:
        """
        python scripts/map_to_cytobands.py {input.annotated_vcf} {input.cytoband_file} {output}
        """


rule define_brain_genes:
    input:
        vcf="results/cytoband_merged_annotations.csv",
        db_brain="databases/brain_database.tsv"
    output:
        "results/brain_cytoband_merged_annotations.csv"
    shell:
        """
        python scripts/define_brain_genes.py --input {input.vcf} --db {input.db_brain} --output {output}
        """


rule identify_rare_cnv:
    input:
        cnv_file="results/brain_cytoband_merged_annotations.csv",
        db_file="databases/dgvMerged.txt"
    output:
        "results/rare_brain_cytoband_merged_annotations.csv"
    shell:
        """
        python scripts/identify_rare_cnv.py --input {input.cnv_file} --db {input.db_file} --output {output}
        """


rule filter_cnv:
    input:
        csv_file="results/rare_brain_cytoband_merged_annotations.csv"
    output:
        "results/filtered_annotations.csv"
    shell:
        """
        python scripts/filter_cnv.py --csv_file {input.csv_file} --output_file {output} --cnvLength --cnvQual --cnvBinSupportRatio --cnvCopyRatio --Chromosome --Classification
        """


rule create_stat_per_individual:
    input:
        csv_file="results/filtered_annotations.csv",
        no_patho_file="results/uncertain_merged_all_cnv_annotations_filtered.csv"
    output:
        all_cnv="results/per_individual_annotations.csv",
        no_patho_cnv="results/uncertain_merged_all_cnv_annotations_filtered_per_ind.csv"
    shell:
        """
        python scripts/create_individual_stat.py {input.csv_file} {output.all_cnv}
        python scripts/create_individual_stat.py {input.no_patho_file} {output.no_patho_cnv}
        """


rule merge_genetic_and_clinical:
    input: 
        indiv_file="results/per_individual_annotations.csv",
        cnv_file="results/filtered_annotations.csv",
        no_patho_file="results/uncertain_merged_all_cnv_annotations_filtered_per_ind.csv"
    output:
        indiv_file="results/merged_per_individual_annotations.csv",
        cnv_file="results/merged_all_cnv_filtered_annotations.csv",
        no_patho_file="results/uncertain_merged_all_cnv_annotations_filtered_per_ind_with_clinical.csv"
    shell:
        """
        python scripts/merge_genetic_and_clinical.py {input.indiv_file} {output.indiv_file}
        python scripts/merge_genetic_and_clinical.py {input.cnv_file} {output.cnv_file}
        python scripts/merge_genetic_and_clinical.py {input.no_patho_file} {output.no_patho_file}
        """


rule exclude_pathogenic:
    input: 
        cnv_file="results/merged_all_cnv_annotations.csv",
        patho_file="results/pathogenic_cnv_with_marshal.csv"
    output: 
        no_patho="results/uncertain_merged_all_cnv_annotations_filtered.csv"
    shell:
        """
        python scripts/remove_pathogenic_cnv.py {input.cnv_file} {input.patho_file} {output.no_patho}
        """




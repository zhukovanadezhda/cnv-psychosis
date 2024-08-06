rule all:
    input:
        "results/filtered_annotations.csv"

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

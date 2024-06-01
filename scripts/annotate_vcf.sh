#!/bin/bash

DATA_PATH="/mnt/c/Users/nadia/projects/schizophrenia/data"
CLASSIFYCNV_PATH="/mnt/c/Users/nadia/projects/schizophrenia/ClassifyCNV"

# Check if CLASSIFYCNV_PATH is valid
if [ ! -d "$CLASSIFYCNV_PATH" ]; then
  echo "CLASSIFYCNV_PATH is not a valid directory."
  exit 1
fi

# Check if VCFS_PATH is valid
if [ ! -d ${DATA_PATH}/GZ ]; then
  echo "DATA_PATH is not a valid directory."
  exit 1
fi

# Count the number of .vcf.gz files in the directory
vcf_files_count=$(find ${DATA_PATH}/GZ -maxdepth 1 -type f \( -name "*.vcf.gz" -o -name "*.vcf" \) | wc -l)
echo "Found $vcf_files_count .vcf.gz and .vcf files in the directory."

processed_files=0

for vcf_gz_file in ${DATA_PATH}/GZ/*.vcf.gz; do
  # Check if there are any .vcf.gz files in the directory
  if [ -e "$vcf_gz_file" ]; then
    gzip -dk "$vcf_gz_file" && mv "${vcf_gz_file%.gz}" ./data/VCF
  fi
done

for vcf_file in ${DATA_PATH}/VCF/*.vcf; do
  # Check if there are any .vcf files in the directory
  if [ -e "$vcf_file" ]; then
    echo "Processing $(basename "$vcf_file")."
    ((processed_files++))
  fi

    # Filter variants to leave only FILTER=PASS
    filtered_vcf="./data/FILTERED_VCF/filtered_$(basename "$vcf_file")"
    bcftools view -i 'FILTER="PASS"' -O z -o "${filtered_vcf}" "${vcf_file}"

    # Convert the filtered VCF to ClassifyCNV input format
    bed_file="./data/BED/$(basename "$vcf_file").bed"
    python3 scripts/vcf2bed.py "${filtered_vcf}" "${bed_file}"

    # Annotate with ClassifyCNV
    result_dir=$(python3 "${CLASSIFYCNV_PATH}/ClassifyCNV.py" --infile "${bed_file}" --GenomeBuild hg38 2>&1 | grep -oP 'Results saved to \K.*')

    if [ -n "$result_dir" ]; then
      mv "${result_dir}/Scoresheet.txt" ./annotations/annotated_$(basename "$vcf_file")
    else
      echo "Error: No results found."
    fi
    echo "Processed $processed_files of $vcf_files_count files."
done

if [ -d "/mnt/c/Users/nadia/projects/schizophrenia/ClassifyCNV_results" ]; then
  rm -r "/mnt/c/Users/nadia/projects/schizophrenia/ClassifyCNV_results"
fi

echo "Annotation completed for all VCF files."
# The impact of copy number variants in individuals in early phases of psychosis
---  
  
#### **Nadezhda ZHUKOVA**<sup>1, 2</sup>  
*nadezhda.zhukova@inserm.fr*  

1 Université Paris Cité, 75013 Paris, France.  
2 Inserm U1266, Institut de Psychiatrie et Neurosciences de Paris, 75014 Paris, France.  

---

## Abstract

Psychotic disorders typically progress through adolescence in multiple stages: ultra-high-risk (UHR) state, first episode of psychosis (FEP), and chronic presentation. About 20% of UHR individuals transition to psychosis within 2 years (converters), while others either maintain attenuated symptoms or recover (non-converters). FEP can signal the onset of schizophrenia or bipolar disorder, and early intervention is beneficial. However, antipsychotic treatment may be harmful due to side effects, making accurate prediction of outcomes vital for targeted intervention. Given the high heritability of psychotic disorders (60%-80%), we hypothesized that genetics could help to identify prognostic biomarkers to predict psychotic transition. Therefore, this study investigates the link between copy number variants (CNVs) and psychosis in UHR individuals. 

A cohort of 387 UHR and FEP individuals was followed up for 12 months to determine their clinical outcomes. They underwent whole-genome sequencing followed by CNV detection, annotation, and statistical analysis. The analysis revealed the presence of eight known pathogenic or likely pathogenic CNVs in this population, as well as nine CNVs associated with schizophrenia in a large-scale pangenomic study by Marshall et al. However, our genome-wide analysis did not reveal significant CNVs, likely due to the small sample size and the low frequency of CNVs. Functional enrichment analysis of rare deletions indicated the involvement of neurotransmitter transport and membrane vesicle processes in psychosis, highlighted the role of neuron projections in high-risk individuals, and suggested a protective role of dendrite regulation in non-converters. Neither significant difference was found between groups, nor significant correlations between clinical scales and CNV parameters. The predictive model was developed to discriminate converters and non-converters, however, the performances were limited. Overall, CNVs have demonstrated limited predictive power, encouraging future integration of multiple data sources.

## Installation

### Clone the repository

```bash
git clone git@github.com:zhukovanadezhda/cnv-psychosis.git
cd cnv-psychosis
```
### Setup the conda environment

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) and [mamba](https://github.com/mamba-org/mamba). Create the `cnv-psychosis` conda environment:

```bash
mamba env create -f binder/environment.yml
```

### Load the environment

```bash
conda activate cnv-psychosis
```

Remark: to deactivate an active environment, use:

```bash
conda deactivate
```

## Usage

### Running Analysis

   To run the analysis workflows, use Snakemake:

   ```bash
   snakemake --snakefile Snakefile
   ```

## Contact

For questions or issues, please open an issue on GitHub or contact [nadiajuckova@gmail.com](mailto:nadiajuckova@gmail.com).

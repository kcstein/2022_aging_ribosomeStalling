# 2022_aging_ribosomeStalling

This repository contains the code and annotation files used to examine how aging impacts ribosome stalling and co-translational proteostasis as published in:

Stein KC, Morales-Polanco F, van der Lienden J, Rainbolt TK, Frydman J. Ageing exacerbates ribosome pausing to disrupt cotranslational proteostasis. Nature. 2022 Jan;601(7894):637-642. doi: 10.1038/s41586-021-04295-4. Epub 2022 Jan 19. [PMID: 35046576](https://pubmed.ncbi.nlm.nih.gov/35046576/)

Raw data and pre-processed codon counts tables can be downloaded from GEO: [GSE152850](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152850)


## Worm and Yeast directories
Code for processing and analyzing data from chronologically aged *S. cerevisiae* and *C. elegans*

### post_processing:
Import codon counts tables, calculate relevant metrics by position and gene, and identify age-dependent stall sites.

### gene_analysis:
**DESeq.R** Uses DESeq2 to determine significantly enriched genes.<br>
**Ubiquitination.R** Examines the association of age-dependent ribosome stalling with co-translation ubiquitination using data from [Duttler, et al. Mol Cell 2013](https://pubmed.ncbi.nlm.nih.gov/23583075/).

### ribosome_pausing:
Examine pausing at features of interest, including at residues and codons of the translatome, at polybasic regions, and at motifs previously associated with ribosome collisions.

### sequence_analysis:
Examine enrichment of sequences at age-dependent ribosome pause sites.

### reporter:
Examine how aging impacts translation elongation along a stalling reporter.

## Published directory

### Meydan_Guydosh_2020_disome
Code used to isolate disome positions originally identified in [Meydan and Guydosh, Mol Cell 2020](https://pubmed.ncbi.nlm.nih.gov/32615089/) to establish expectations of ribosome occupancy around polybasic regions and analyze association with age-dependent stalling.

### worm_ageDependent_aggregration_massSpec
Code used to analyze the relationship between age-dependent ribosome stalling and protein aggregation in *C. elegans*. Protein aggregation data was obtained from [Walther, et al. Cell 2015](https://pubmed.ncbi.nlm.nih.gov/25957690/) and [David, et al. PLoS Biol 2010](https://pubmed.ncbi.nlm.nih.gov/20711477/).

### yeast_lifespan_rqc
Code used to analyze the relationship between lifespan and RQC flux. Chronological lifespan data was obtained from [Powers, et al. Genes Dev 2006](https://pubmed.ncbi.nlm.nih.gov/16418483/) and RQC flux data was obtained from [Brandman, et al. Cell 2012](https://pubmed.ncbi.nlm.nih.gov/23178123/).

### Young_Green_2015_3AT
Code used to validate statistical methodology of identifying ribosome pause sites using data from [Young, et al. Cell 2015](https://pubmed.ncbi.nlm.nih.gov/26276635/).



## Figures directory

Specific visualization code for generating plots included in publication.



## Data_tables directory

Annotation of polybasic sites within the yeast and worm proteomes and summary tables used for analysis.


MovingAverage.R: function to average ribosome occupancy over user-defined window.


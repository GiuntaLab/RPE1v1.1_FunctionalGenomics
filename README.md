# RPE1 genome for Functional Genomics
This folder contains documentation, scripts and files relative to the main figures of the following manuscript: 

 Corda, Volpe et al. **Cell line matched reference enables high-precision functional genomics**. Nature Communications, 2025 (*in press*)

In this new work, we demonstrate the improvements when using the diploid reference genome of the human laboratory cell  line RPE-1, RPE1v1.1 (Volpe, Colantoni et al. [The reference genome of the human diploid cell line RPE-1](https://www.nature.com/articles/s41467-025-62428-z). *Nature Communications*, 2025; https://github.com/GiuntaLab/RPE1), including the two fully phased haplotypes, for sequencing and omics data analyses. \
RPE1v1.1 has been used as matched reference genomes to analyze sequencing data generated from the same cell line, an approach we refer to as *isogenomic reference genome*. The improvement in alignment quality using matched reads-reference enables high-precision mapping for profiling phased epigenome and methylome. The two publications were originally shared as preprints as one manuscript: 

**Preprint:** Volpe et al., [The complete diploid reference genome of RPE-1 identifies human phased epigenetic landscapes](https://www.biorxiv.org/content/10.1101/2023.11.01.565049v1), 2023. \
Documentation, scripts, files and figures relative to the main figures are available in the https://github.com/GiuntaLab/RPE1/tree/main/RPE1v1.1_FunctionalGenomics folder.

Altogether, the improvements in haplotype-resolved mapping, guide RNA (gRNA) design for genome engineering, transcriptome analysis and alignment of sequencing reads serve as a proof-of-concept calling for a comprehensive catalog of complete assemblies for commonly used cells for a widespread application of isogenomic reference genomes to enable high-precision multi-omics analyses.

# RPE1 Genome Assembly
## RPE1v1.1
Volpe, Colantoni et al. [The reference genome of the human diploid cell line RPE-1](https://www.nature.com/articles/s41467-025-62428-z). *Nature Communications*, 2025

Docmentation, scripts, files and figures relative to the main figures are available in the https://github.com/GiuntaLab/RPE1/tree/main/RPE1v1.1 folder.
Supplemental information is available in the **Zenodo repository**: https://doi.org/10.5281/zenodo.15789913.

## Genome availability

The RPE1v1.1 genome has been deposited in NCBI GenBank under the accession numbers [JBJQNK000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBJQNK000000000) (Hap1) and [JBJQNL000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBJQNL000000000) (Hap2), with links to BioProject accession numbers [PRJNA1193286](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1193286) (Hap1) and [PRJNA1193302](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1193302) (Hap2), both under the umbrella BioProject accession number [PRJNA1195024](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1195024).


To download the diploid RPE1v1.1 genome with the chromosome names used in the study, please use the following link, fill in the form, access the link that will be emailed to you and download the *GIUNTAlab_RPE1V1.1.gz* file: https://forms.gle/7TPggHvT1Na2G3oX7

RPE1v1.1 genome is also available in the **UCSC Genome Browser**:
- **Hap1**: https://genome.ucsc.edu/h/GCA_050656345.1
- **Hap2**: https://genome.ucsc.edu/h/GCA_050656315.1

## Sequencing data availability
PacBio HiFi and Illumina raw reads generated and used for this study are available under the BioProject accession number [PRJNA1359433](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1359433).
- **PacBio HiFi** raw reads used for genome assembly and the NucPlot analysis:  [SRR33464826](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464826), [SRR33464827](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464827), [SRR33464828](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464828), and [SRR33464829](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464829)
- **ONT** raw reads used for genome assembly: [SRR33464817](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464817), [SRR33464818](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464818), [SRR33464819](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464819), [SRR33464820](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464820), [SRR33464821](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464821), [SRR33464822](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464822), [SRR33464823](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464823), [SRR33464824](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464824), [SRR33464830](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464830), and [SRR33464831](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464831)
- **Hi-C** raw reads used for genome assembly: [SRR33464825](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464825)
- **PacBio HiFi** raw reads from RPE-1 cells expressing an inducible Cas9 and a PAN-centromeric sgRNA used for the NM and MAPQ score analysis: [SRR35988140](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR35988140)
- **Illumina** raw reads from RPE-1 cells expressing an inducible Cas9 and a non-targeting sgRNA used for the NM and MAPQ score analysis: [SRR35988139](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR35988139)


All the analyses were performed in the Giunta Lab HPC or using the University server at [HPC TeraStat2](https://www.dss.uniroma1.it/it/HPCTerastat2).


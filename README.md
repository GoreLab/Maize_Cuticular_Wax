# Maize_Cuticular_Wax
**PREPRINT**
This repository and website documents all analyses, summary, tables and figures associated with the following PREPRINT: Integrative multi-omic analysis identifies genes associated with cuticular wax biogenesis in adult maize leaves (a place holder for the link)

**Abstract**

Studying the genetic basis of leaf wax composition and its correlation with leaf cuticular conductance (gc) is crucial for improving crop water-use efficiency. The leaf cuticle, which comprises a cutin matrix and various waxes, functions as an extracellular hydrophobic layer, protecting against water loss upon stomatal closure. To address the limited understanding of genes associated with the natural variation of leaf cuticular waxes and their connection to gc, we conducted statistical genetic analyses using leaf transcriptomic, metabolomic, and physiological data sets collected from a maize (Zea mays L.) panel of ~300 inbred lines. Through a random forest analysis with 60 cuticular wax traits, it was shown that high molecular weight wax esters play an important role in predicting gc. Integrating results from genome-wide and transcriptome-wide studies (GWAS and TWAS) via a Fisher’s combined test revealed 231 candidate genes detected by all three association tests. Among these, 11 genes exhibit known or predicted roles in cuticle-related processes. Throughout the genome, multiple hotspots consisting of GWAS signals for several traits from one or more wax classes were discovered, identifying four additional plausible candidate genes and providing insights into the genetic basis of correlated wax traits. Establishing a partially shared genetic architecture, we identified 35 genes for both gc and at least one wax trait, with four considered plausible candidates. Our study uncovered the genetic control of maize leaf waxes, establishing a link between wax composition and gc, with implications for potentially breeding more water-use efficient maize.

**Data availability**

Much of the supporting data and output from the analyses documented here are too large for GitHub.

**Scripts 1 to 3 are related to the cuticular wax data processing, BLUP calculation, outlier removal, and random forest regression.**

1.process_wax_raw.R performs data cleaning, imputed zero values that were under detection limit, and incorporate the experimental design information.

2.0.wrapper_boxcox_modelSel.R calculates BLUPs for each of the 60 cuticular wax traits to combine phenotype from two environments and remove special variation.

2.1.organize_wax_blup.R formats BLUPs for the 60 cuticular wax traits by assigning proper accession IDs.

2.2.Transform_BLUP.R performs box-cox transformation for BLUPs for each of the 60 cuticular wax traits.

2.3.OutlierRM_for_Transformed_Wax.R performs outlier removal for the transformed BLUPs for each of the 60 cuticular wax traits.

3.RF_ce_pred_by_wax_updated.R performs random forest regression for cuticular conductance using 60 cuticular waxes.

**Script 4 to 10 perform GWAS, TWAS, Fisher’s combined test and candidate gene selection.**

4.model_selection_forGWAS_TWAS.R performs model selection for GWAS and TWAS analysis.

5.GAPIT_GWAS_wax_transBLUP_HMP3_v4.R performs GWAS using the GAPIT R package.

6.define_GWAS_peaks_20200925.R defines peak SNPs using the top 0.002% associated SNPs in GWAS for each of the 60 cuticular wax traits.

7.TWAS_wax.R performs TWAS using the rrBLUP R package for each of the 60 cuticular wax traits.

8.Fisher's_comb_wax.R performs Fisher’s combined test (FCT) for each of the 60 cuticular wax traits.

9.wax_cand_list_supp_tables.R selects candidate genes identified in GWAS, TWAS and FCT using the criteria described in the manuscript.

10.intersection_GWAS_TWAS_FCT.R selects candidate genes identified in all three methods (GWAS, TWAS and FCT).

**Scripts 11 to 12 are related to GWAS hotspot and candidate genes located in GWAS hotspots.**

11.claim_hotspot.R defines GWAS hotspots using the sliding-window strategy.

12.define_HS_regions.R sets boundaries of GWAS hotspot genomic regions and selects candidate genes allocated in GWAS hotspots.

**Scripts 13 to 14 are related to visualizations of GWAS, TWAS, FCT test  results for WE 49:0 and GWAS hotspots in the maize genome.**

13. Figure 3. Manhattan plots of GWAS, TWAS, and FCT results for WE 49:0.
   
14. Figure 4. Genomic positions of GWAS hotspots in the maize genome.

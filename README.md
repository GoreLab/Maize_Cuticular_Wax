# Maize_Cuticular_Wax

**Abstract**

This is a place holder

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


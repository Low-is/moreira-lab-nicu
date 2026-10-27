# Developmental and Injury-Associated Variation in the Neonatal Pulmonary Immune Compartment. A Donor-Level Spatial Transcriptomic Analysis

<p align="center">
  <img src = "images/Graphical-Abstract.png" alt="Graphical Abstract" width="100%">
</p>

## Overview
This project investigates how developmental age and lung injury shape immune cell transcriptional profiles in the developing human lung using imaging-based spatial transcriptomic data from 16 donors and 122,005 immune cells. Using donor-level pseudobulk analysis, principal component analysis, and spatial characterization, we examined associations with postmenstrual age (PMA), assessed potential transcript spillover from neighboring non-immune cells, and evaluated whether preterm lung injury was associated with transcriptional, compositional, or spatial variation beyond developmental expectations. Developmental age was the dominant measured correlate of immune variation, while no independent association with lung injury was detected in this cohort. These findings highlight the importance of accounting for developmental maturity when investigating immune biology in the preterm lung.

## Study Design and Data Sources
This project uses publicly available human lung imaging-based spatial transcriptomic data from GSE297945, generated using the 10x Genomics Xenium platform with a targeted 343-gene panel. The source study included 17 biological donors from the Human Infant Lung Repository. 


## Cohort Characteristics
| Sample Name | Age at Birth | Life Span | Age at Death | COD | Disease Score (DX) |
|-------------|--------------|-----------|--------------|-----|----|
| PDL001 | 16 weeks | 0 weeks | 16 weeks | Not viable | n/a |
| PDL002 | 18 weeks | 0 weeks | 18 weeks | Not viable | n/a |
| PDL003 | 21 weeks | 0 weeks | 21 weeks | Prematurity | 1 |
| PDL004 | 23 weeks | 7 weeks | 30 weeks | Necrotizing enterocolitis | 1 |
| PDL005 | 24 weeks | 0 weeks | 24 weeks | Prematurity | 3 |
| PDL006 | 24 weeks | 0 weeks | 24 weeks | Prematurity | 0 |
| PDL007 | 26 weeks | 2 weeks | 28 weeks | Pulmonary Hypoplasia | 1 |
| PDL008 | 27 weeks | 16 weeks | 43 weeks | BPD | 2 |
| PDL009 | 28 weeks | 3 weeks | 31 weeks | congenital heart disease | 2 |
| PDL010 | 34 weeks | 4 weeks | 38 weeks | Meningitis | 3 |
| PDL011 | 35 weeks | 0 weeks | 35 weeks | Encephalopathy | 0 |
| PDL012 | 40 weeks | 0 weeks | 40 weeks | Encephalopathy | 0 |
| PDL013 | 36 weeks | 0 weeks | 36 weeks | Encephalitis | 0 |
| PDL014 | 29 weeks | 0 weeks | 29 weeks | CHAOS syndrome | n/a |
| PDL015 | 33 weeks | 0 weeks | 33 weeks | Pulmonary hypoplasia | n/a |
| PDL016 | 35 weeks | 28 weeks | 63 weeks | ARB fetopathy (severe BPD) | n/a |
| PDL017 | 28 weeks | 30 years | n/a | BPD (explant of adult with BPD) | n/a |

*_COD_: **Characteristics of Disease**



Postmenstrual age (PMA) was calculated as gestational age plus postnatal age using source metadata. The source-provided histologic injury score was not used to define study groups because it was incomplete and discordant with the source injury classifications. Ventilation and supplemental oxygen exposures were unavailable.

## Disease Severity Score
The source dataset included a pathologist-assigned histologic disease severity score ranging from 0 to 3 (0, no findings; 1, mild; 2, moderate; 3, severe). The score summarized acute or subacute injury (hyaline membranes, inflammation, and hemorrhage), alveolar architectural changes (hyperinflation and simplification), chronic remodeling (increased lobular smooth muscle and fibrosis), and vascular abnormalities (arterial wall hypertrophy and thrombi). Although available in the source metadata, this score was not used to define study groups or as an analytical variable because it was missing for some donors and discordant with the source investigators’ injury classifications.

## Processed Data Availability
Due to their large file sizes, the processed Seurat object and raw TAR archive are too large to upload directly to this GitHub repository. Both files are available through the supplementary files associated with GSE297945.

- Processed seurat obj: ![GSE297945_final_complete_cniches_tniches_allsamples_2025_05_26.rds](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297945) (RDS; 598.0 Mb)
- Raw data archive: ![GSE974945_RAW.tar](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297945) (TAR; 69.9 Gb)

To access either file, visit the GEO accession page linked above and locate the corresponding file under Supplementary files.


# Citations
Mallapragada S, Lyu R, Williams-Katek AL, Fischer BK, Vannan A, Hadad N, Mee ED, Shirazi SP, Jetter CS, Negretti NM, Hilgendorff A, Eldredge LC, Deutsch GH, McCarthy DJ, Kropski JA, Sucre JMS, Banovich NE. A spatial transcriptomic atlas of acute neonatal lung injury across development and disease severity. bioRxiv [Preprint]. 2025 Jun 5:2025.06.02.656433. doi: 10.1101/2025.06.02.656433. Update in: Am J Physiol Lung Cell Mol Physiol. 2025 Sep 1;329(3):L376-L388. doi: 10.1152/ajplung.00191.2025. PMID: 40502191; PMCID: PMC12157373.

Lal CV, Travers C, Aghai ZH, Eipers P, Jilling T, Halloran B, Carlo WA, Keeley J, Rezonzew G, Kumar R, Morrow C, Bhandari V, Ambalavanan N. The Airway Microbiome at Birth. Sci Rep 6, 2016. doi: 10.1038/SREP31023.

De Biasi S, Neroni A, Nasi M, Lo Tartaro D, Borella R, Gibellini L, Lucaccioni L, Bertucci E, Lugli L, Miselli F, Bedetti L, Neri I, Ferrari F, Facchinetti F, Berardi A, Cossarizza A. Healthy preterm newborns: Altered innate immunity and impaired monocyte function. Eur J Immunol 53, 2023. doi: 10.1002/EJI.202250224.

Melville JM, Moss TJM. The immune consequences of preterm birth. Front Neurosci 7, 2013. doi: 10.3389/FNINS.2013.00079.

Hillman NH, Kallapur SG, Jobe AH. Physiology of transition from intrauterine to extrauterine life. Clin Perinatol 39: 769–783, 2012. doi: 10.1016/j.clp.2012.09.009.

Dirks J, Fortmann I, Marißen J, Pagel J, Reichert L, Kipke H, Dammann MT, Göpel W, Birkner T, Dahm K, Forslund-Startceva SK, Viemann D, Rupp J, Morbach H, Härtel C. Effect of gestational age on lymphocyte phenotypes in hospitalized preterm infants. Journal of Allergy and Clinical Immunology 157, 2025. doi: 10.1016/j.jaci.2025.09.030.

Zaramella P, Munari F, Stocchero M, Molon B, Nardo D, Priante E, Tosato F, Bonadies L, Viola A, Baraldi E. Innate immunity ascertained from blood and tracheal aspirates of preterm newborn provides new clues for assessing bronchopulmonary dysplasia. PLoS One 14, 2019. doi: 10.1371/JOURNAL.PONE.0221206.

Moreira A, Noronha M, Joy J, Bierwirth N, Tarriela A, Naqvi A, Zoretic S, Jones M, Marotta A, Valadie T, Brick J, Winter C, Porter M, Decker I, Bruschettini M, Ahuja SK. Rates of bronchopulmonary dysplasia in very low birth weight neonates: a systematic review and meta-analysis. Respir Res 25: 1–11, 2024. doi: 10.1186/S12931-024-02850-X/TABLES/2.

Thébaud B, Abman SH. Bronchopulmonary Dysplasia. Am J Respir Crit Care Med 175: 978–985, 2007. doi: 10.1164/rccm.200611-1660PP.

Sampah MES, Hackam DJ. Dysregulated Mucosal Immunity and Associated Pathogeneses in Preterm Neonates. Front Immunol 11, 2020. doi: 10.3389/FIMMU.2020.00899.

Humberg A, Fortmann I, Siller B, Kopp MV, Herting E, Göpel W, Härtel C. Preterm birth and sustained inflammation: consequences for the neonate. Semin Immunopathol 42: 451–468, 2020. doi: 10.1007/S00281-020-00803-2.

Honda A, Hoeksema MA, Sakai M, Lund SJ, Lakhdari O, Butcher LD, Rambaldo TC, Sekiya NM, Nasamran CA, Fisch KM, Sajti E, Glass CK, Prince LS. The Lung Microenvironment Instructs Gene Transcription in Neonatal and Adult Alveolar Macrophages. J Immunol 208: 1947–1959, 2022. doi: 10.4049/JIMMUNOL.2101192.

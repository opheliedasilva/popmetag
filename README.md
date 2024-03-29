# Protist population metagenomics

___Genomic differentiation of three pico‐phytoplankton species in the Mediterranean Sea___

Ophélie Da Silva<sup>1,2</sup>, Sakina-Dorothée Ayata<sup>1,2,3</sup>, Enrico Ser-Giacomi<sup>3,4</sup>, Jade Leconte<sup>5</sup>, Eric Pelletier<sup>5,6</sup>, Cécile Fauvelot<sup>1,7</sup>, Mohammed-Amin Madoui<sup>5</sup>, Lionel Guidi<sup>1,6</sup>, Fabien Lombard<sup>1,6,8</sup>, Lucie Bittner<sup>2,8</sup>

1 Sorbonne Université, CNRS, Laboratoire d’Océanographie de Villefranche, LOV, F-06230 Villefranche-sur-Mer, France <br/> 2 Institut de Systématique, Evolution, Biodiversité (ISYEB), Muséum national d’Histoire naturelle, CNRS, Sorbonne Université, EPHE, Université des Antilles CP 50, 57 rue Cuvier, 75005 Paris, France <br/> 3 Sorbonne Université, UMR 7159 CNRS-IRD-MNHN, LOCEAN-IPSL, 75005 Paris, France <br/> 4 Department of Earth, Atmospheric and Planetary Sciences, Massachusetts Institute of Technology, 54-1514 MIT, Cambridge, MA 02139, USA <br/> 5 Génomique Métabolique, Genoscope, Institut François Jacob, CEA, CNRS, Univ Evry, Université Paris-Saclay, Evry, France <br/> 6 Research Federation for the study of Global Ocean Systems Ecology and Evolution, FR2022/Tara Oceans GOSEE, 3 rue Michel-Ange, 75016 Paris, France <br/> 7 Institut de Recherche pour le Développement (IRD), UMR ENTROPIE, Nouméa, New Caledonia <br/> 8 Institut Universitaire de France (IUF), Paris, France

[Article is available here](https://sfamjournals.onlinelibrary.wiley.com/)

## Summary

For more than a decade, high‐throughput sequencing has transformed the study of marine planktonic communities and has highlighted the extent of protist diversity in these ecosystems. Nevertheless, little is known relative to their genomic diversity at the species‐scale as well as their major speciation mechanisms. An increasing number of data obtained from global scale sampling campaigns is becoming publicly available, and we postulate that metagenomic data could contribute to deciphering the processes shaping protist genomic differentiation in the marine realm. As a proof of concept, we developed a FAIR pipeline and focused on the Mediterranean Sea to study three *a priori* abundant protist species: *Bathycoccus prasinos*, *Pelagomonas calceolata* and *Phaeocystis cordata*. We compared the genomic differentiation of each species in light of geographic, environmental and oceanographic distances. We highlighted that isolation‐by‐environment shapes the genomic differentiation of *B. prasinos*, whereas *P. cordata* is impacted by geographic distance (i.e. isolation‐by‐distance). At present time, the use of metagenomics to accurately estimate the genomic differentiation of protists remains challenging since coverages are lower compared to traditional population surveys. However, our approach sheds light on ecological and evolutionary processes occurring within natural marine populations and paves the way for future protist population metagenomic studies. <br/> **Keywords**: Metagenomics; Genetic differentiation; Protists; Marine plankton; Mediterranean Sea.

## Description

The objective of our study was to develop an original, bioinformatic pipeline aiming to exploit the currently available metagenomic data for characterizing genetic differentiation of protists in the ecosystems. We focused on three planktonic protistan species *a priori* abundant in the Mediterranean Sea.  By gathering reference sequences and metagenomic data previously published, we investigated the genetic differentiation at the species scale and then if this one could be related to external drivers (*i.e.,* geography, environmental conditions and oceanographic circulation). The contrasted results for the three species leading us to discuss how current metagenomics could support and provide new resources for population genetics. 
<p align="center">
<img src="https://github.com/opheliedasilva/popmetag/blob/master/additional/global_overview.png" alt="drawing" width="600"/>
</p>

## Data and code access
Codes used to carry out the study and to produce figures are available in this GitHub repository. Genetic data are publicly available [here](https://doi.org/10.5281/zenodo.6434681).

```
                                             + coastline_medit.csv
                                             + coord.txt
  main directory ---+ data  -----+ raw  -----+ mct_protist.txt
                    |            |           + MEDAR-MEDATLAS
                    |            |           + proba_cnx
                    |            |           + snps
                    |            + out
                    |
                    |            + 0_preprocessing_circu_mct.R
                    |            + 0_preprocessing_circu_proba_cnx.R
                    |            + 0_preprocessing_env.R
                    + src  ------+ 0_preprocessing_fst.R
                                 + 0_preprocessing_geo.R
                                 + 1_analysis_environment.R
                                 + 1_analysis_linear_models.R
                                 + Rsources.R
```

For any question, please contact: <oph.dasilva@gmail.com>.

This is a repository for building connections between phenotypes and mass spectrum metabolites. This is a second round of analysis 
attempting to build a simpler system that contains the data and transformation of it and connecting the results to verifiable
observations. The goals are find features that associate with the phenotypic differences in the samples taking into account
the phylogenetic diversity of the species these samples are derived from. The samples are also paired, one from supernatent
and one from whole cell extract and we are interested in the differences between these for a given sample.

Previous findings and results are in ../Rhodotorula_MS2_pheno_explore and we will revisit some of those findings but starting 
with simple data organization and orientation.

The main phenotype we are interested in relates to color and intensity which is captured by L*, a*, and b*. The metadata file 
has these phenotypes encoded. A second set of phenotypes relates to growth rate (Area Under the Curve=AUC) calculated currently
just for Copper media but will be populated with other tested metals. 


There are two sister folder that represent different stages of analysis 
* /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Rodeo/ - functional and gene information for strains - db folder has database of functional domains and the sequences for those individuals sequenced and annotated
* /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_MS2_pheno_explore/ - first attempt at data analysis, has several statistical 


GOALS.md lists analysis goals and things we want to develop and establish statistical investigation

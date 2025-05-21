library("UniprotR")
library(stringr)
# getwd()

all_tissue <- read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/Dataset/Plasma_membrane/tissue_category_rna_Any_Detected_all_tissue.tsv", sep="\t")
many_tissue <- read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/Dataset/Plasma_membrane/tissue_category_rna_Any_Detected_many_tissue.tsv", sep="\t")
some_tissue <- read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/Dataset/Plasma_membrane/tissue_category_rna_Any_Detected_some_tissue.tsv", sep="\t")
single_tissue <- read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/Dataset/Plasma_membrane/tissue_category_rna_Any_Detected_single_tissue.tsv", sep="\t")

#combine datasets
combined_data <- rbind(all_tissue, many_tissue, some_tissue, single_tissue)
combined_data_unique_accession <- subset(combined_data, !grepl(",",combined_data$Uniprot) & !combined_data$Uniprot=="")
GetSubcellular_location(combined_data_unique_accession$Uniprot[], directorypath = "~/GitHub/MVP_code/R_output")
combined_subcellular <- read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Subcellular location Info.csv", sep=",")
combined_subcellular_single_pass <- subset(combined_subcellular,grepl("Single-pass",combined_subcellular$Subcellular.location..CC.))
combined_subcellular_single_pass_extracellular <- subset(combined_subcellular_single_pass, grepl("Extracellular", combined_subcellular_single_pass$Topological.domain))
combined_subcellular_single_pass_extracellular2 <- subset(combined_subcellular_single_pass_extracellular, grepl("Cytoplasmic", sub(".*note=", "", combined_subcellular_single_pass_extracellular$Topological.domain)))
combined_subcellular_single_pass_cytoplasmic <- subset(combined_subcellular_single_pass_extracellular, !grepl("Cytoplasmic", sub(".*note=", "", combined_subcellular_single_pass_extracellular$Topological.domain)))

write.table(combined_subcellular_single_pass_extracellular2, file = "C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/combined_subcellular_single_pass_extracellular.tsv", sep='\t')
write.table(combined_subcellular_single_pass_cytoplasmic, file = "C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/combined_subcellular_single_pass_cytoplasmic.tsv", sep='\t')
write.table(combined_subcellular_single_pass, file = "C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/combined_subcellular_single_pass.tsv", sep='\t')

myFiles = list.files("C:/Users/Zhenh/Documents/GitHub/MVP_code/Dataset", pattern = ".tsv", recursive = T, full.names = T)

combined_subcellular_single_pass_accession <-combined_subcellular$X[grepl("Single-pass",combined_subcellular$Subcellular.location..CC.)]
combined_single_pass <- subset(combined_data_unique_accession, combined_data_unique_accession$Uniprot %in% combined_subcellular_single_pass_accession)

for (file in myFiles) {
  location = str_match(file, "Dataset/\\s*(.*?)\\s*/")[,2]
  filename_extracellular = sub(".tsv", "_single_pass_extracellular.tsv", basename(file))
  filename_cytoplasmic = sub(".tsv", "_single_pass_cytoplasmic.tsv", basename(file))
  filename = sub(".tsv", "_single_pass.tsv", basename(file))
  data1 =read.delim(file)
  data1_unique <- subset(data1, !grepl(",",data1$Uniprot) & !data1$Uniprot=="")
  if(file.exists(paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, sep=""))){
    
  } else {
    dir.create(file.path(paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, sep="")))
  }
  GetSubcellular_location(data1_unique$Uniprot[], directorypath = paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, sep=""))
  data1_subcellular <- read.delim(paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, "/Subcellular location Info.csv", sep=""), sep=",")
  
  data1_subcellular_single_pass <- subset(data1_subcellular,grepl("Single-pass",data1_subcellular$Subcellular.location..CC.))
  data1_subcellular_single_pass_extracellular <- subset(data1_subcellular_single_pass, grepl("Extracellular", data1_subcellular_single_pass$Topological.domain))
  
  data1_subcellular_single_pass_extracellular2 <- subset(data1_subcellular_single_pass_extracellular, grepl("Cytoplasmic", sub(".*note=", "", data1_subcellular_single_pass_extracellular$Topological.domain)))
  data1_subcellular_single_pass_cytoplasmic <- subset(data1_subcellular_single_pass_extracellular, !grepl("Cytoplasmic", sub(".*note=", "", data1_subcellular_single_pass_extracellular$Topological.domain)))
  
  write.table(data1_subcellular_single_pass_extracellular2, file = paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, "/", filename_extracellular, sep=""), sep='\t')
  write.table(data1_subcellular_single_pass_cytoplasmic, file = paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, "/", filename_cytoplasmic, sep=""), sep='\t')
  write.table(data1_subcellular_single_pass, file = paste("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/", location, "/", filename, sep=""), sep='\t')
}

data_cytoplasmic$Protein.primary.name = GetProteinAnnontate(data_cytoplasmic$Accession, "gene_primary")
data_extracellularProtein.primary.name = GetProteinAnnontate(data_extracellular$Accession, "gene_primary")
data_extracellular$Protein.primary.name = data_extracellularProtein.primary.name

data_extracellular =read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Plasma_membrane/Single_pass_plasma_extracellular.csv",sep=",")
data_cytoplasmic =read.delim("C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Plasma_membrane/Single_pass_plasma_cytoplasmic.csv",sep=",")

GetSequences(data_extracellular$Accession, directorypath="C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Plasma_membrane/")
GetSequences(data_cytoplasmic$Accession, directorypath="C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Plasma_membrane/")

write.csv(data_extracellular, file = "C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Plasma_membrane/Single_pass_plasma_extracellular.csv", row.names =F)
write.csv(data_cytoplasmic, file = "C:/Users/Zhenh/Documents/GitHub/MVP_code/R_output/Plasma_membrane/Single_pass_plasma_cytoplasmic.csv", row.names =F)
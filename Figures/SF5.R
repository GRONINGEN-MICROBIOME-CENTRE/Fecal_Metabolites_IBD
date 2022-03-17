

#Heatmap 
library(readxl)
library(reshape2)
library(ggplot2)
library (meta)
library(pheatmap)
sp_metab <- read_excel("~/Desktop/IBD_metabolomics_2022/6.Supplementary_tables/STB_Meta-analysis_microbiome_presence.xlsx")
annot <- read.delim("~/Desktop/IBD_metabolomics_2022//1.Input/Metabolon_annotation_v2.txt", row.names=1)
my_taxa = read.delim("~/Desktop/IBD_metabolomics_2022//1.Input/Taxa_mpa3.txt", row.names=1, check.names = F)

only_sp1=sp_metab[,c("Metabolite","Specie","te_random", "FDR_meta_random")]
only_sp1$te_random[only_sp1$FDR_meta_random>0.05]=0
only_sp2=only_sp1[,c("Metabolite","Specie","te_random")]
only_sp2$te_random=as.numeric(as.character(only_sp2$te_random))

only_sp2=dcast(only_sp2, Metabolite ~ Specie)

row.names(only_sp2)=only_sp2$Metabolite
only_sp2$Metabolite=NULL
only_sp2$sp_Shannon_Index=NULL
colnames(only_sp2)=gsub("sp_","",colnames(only_sp2))


paletteLength <- 9
myColor <- colorRampPalette(c("blue", "white", "firebrick3"))(paletteLength)
mybreaks= c(seq(min(only_sp2), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(only_sp2)/paletteLength, max(only_sp2), length.out=floor(paletteLength/2)))
my_path_annot= list(
  SUPER.PATHWAY = c(Amino_Acid="firebrick3", Carbohydrate = "steelblue2", Cofactors_and_Vitamins = "gold3", Energy = "limegreen", Lipid = "grey77", Nucleotide = "pink1", Peptide = "salmon", SCFA="turquoise", Unknown_compount="black", Xenobiotics="azure"),
  Class = c(Actinobacteria="blue",Bacilli = "tan2", Bacteroidia = "red2", Betaproteobacteria = "gray74", 
            Clostridia = "darkorchid2", Coriobacteriia = "gold", Deltaproteobacteria = "lightblue1", Erysipelotrichia="seagreen2", 
            Firmicutes_unclassified="yellow", Gammaproteobacteria="lightpink1", Methanobacteria="brown", Negativicutes="olivedrab", Verrucomicrobiae="gray89"))


pheatmap(only_sp2, color = myColor, breaks = mybreaks, cutree_rows =2, cutree_cols = 3, fontsize_col = 5, fontsize_row = 3, border_color = NA, scale = "none", annotation_row = annotation_metabolites, 
         annotation_colors = my_path_annot, annotation_col = annotation_taxa)



#Prepare annotation files

annot$names=make.names(row.names(annot))
annotation_metabolites=annot[, c("names","SUPER.PATHWAY")]
row.names(annotation_metabolites)=annotation_metabolites$names
annotation_metabolites$names=NULL
annotation_metabolites$SUPER.PATHWAY=gsub(" ","_",annotation_metabolites$SUPER.PATHWAY)
#my_list=data.frame(Metabolite=unique(only_sp1$metabolite), Bla="X")
#annotation_metabolites=merge(my_list, annotation_metabolites, by.x="Metabolite",by.y="row.names", all.x=T)
#row.names(annotation_metabolites)=annotation_metabolites$Metabolite
annotation_metabolites$Metabolite=NULL
annotation_metabolites$SUPER.PATHWAY[is.na(annotation_metabolites$SUPER.PATHWAY)]="SCFA"
annotation_metabolites$Bla=NULL
annotation_taxa=my_taxa[,c(1:2)]
annotation_taxa=annotation_taxa[grep("s__", rownames(annotation_taxa)), ]
annotation_taxa=annotation_taxa[!grepl("t__", rownames(annotation_taxa)), ]
annotation_taxa$Class=gsub(".*c__","",row.names(annotation_taxa))
annotation_taxa$Class=sapply(strsplit(annotation_taxa$Class, "|", fixed=TRUE), head, 1)
annotation_taxa$Specie=gsub(".*s__","",row.names(annotation_taxa))
annotation_taxa$Specie=sapply(strsplit(annotation_taxa$Specie, "|", fixed=TRUE), head, 1)
row.names(annotation_taxa)=annotation_taxa$Specie
annotation_taxa$Specie=NULL
annotation_taxa$LLDeep_1255=NULL
annotation_taxa$LLDeep_0444=NULL
row.names(annotation_taxa)=gsub("sp_","",rownames(annotation_taxa))



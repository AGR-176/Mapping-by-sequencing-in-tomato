#################################################
#CHANGE THESE PARAMETERS AS NEEDED

# Annovar output files
exonic_variant_function_file<-"mutant_bulk_functional_annotation.exonic_variant_function"
variant_function_file<-"mutant_bulk_functional_annotation.variant_function"

# gff3 file for tomato
gff3_file<-"ITAG2.4_gene_models.gff3"

# vcf file with the region of interest
region_of_interest_vcf<-"all_variants.ch05_0_3000000.vcf"

# Output
output_file<-"candidate_mutations.tsv"

#################################################

# Open Funct Descr from gff3
funct_descr<-system(paste("awk -F\"\t\" '{if($3==\"mRNA\") print $9}' ", gff3_file ,sep=""), intern=T)
gene_name<-gsub("ID=mRNA:","",funct_descr)
gene_name<-gsub(";.*","",gene_name)
gene_name<-gsub(".1$","",gene_name)
annotation<-gsub(".*Note=","",funct_descr)
annotation<-gsub(";.*","",annotation)
funct_descr<-data.frame(cbind(gene_name,annotation))
names(funct_descr)<-c("gene_id","funct_descr")

# Open file with exonic variants
exons<-read.delim(exonic_variant_function_file,header=F)
names(exons)<-c("line","funct_consequ","gene","chr","start","end","base1","base2","chrom","pos","id","ref","alt", "pos_qual", "filter", "info", "format", "genotype")
exons<-exons[,c("funct_consequ","gene","chr","start","end","base1","base2","chrom","pos","id","ref","alt", "pos_qual", "filter", "info", "format", "genotype")]

# Open file with all variants
all_var<-read.delim(variant_function_file,header=F)
names(all_var)<-c("funct_consequ","gene","chr","start","end","base1","base2","chrom","pos","id","ref","alt", "pos_qual", "filter", "info", "format", "genotype")

# Open vcf file (modify path to bcftools)
vcf<-system(paste("bcftools query --print-header --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ",region_of_interest_vcf,sep=""), intern=T)
vcf_header<-unlist(strsplit(gsub("(\\[\\d\\]|# |\\:GT)","", vcf[1]), "\t"))
vcf<-data.frame(matrix(unlist(strsplit(vcf[2:length(vcf)],"\t")), byrow=T, ncol=7))
names(vcf)<-vcf_header

########################
# Merge all variants (exonic and non exonic)
annovar<-rbind(exons,all_var[-which(all_var$funct_consequ=="exonic"),])
annovar<-annovar[, c("chr", "pos", "ref", "alt", "pos_qual","funct_consequ", "gene")]

# Merge variants with human readable genotype from vcf file
annovar<-merge(annovar,vcf,by.x=c("chr", "pos", "ref", "alt"), by.y=c("CHROM", "POS", "REF", "ALT"), all.x=T)

# Merge variants with gene functional description (need to create new gene_id column)
annovar$gene_id<-gsub("gene:","",annovar$gene)
annovar$gene_id<-gsub(":.*","",annovar$gene_id)
annovar$gene_id<-gsub("\\(.*","",annovar$gene_id)
annovar$gene_id[grep("(upstream|downstream|intergenic)",annovar$funct_consequ)] <- ""
annovar<-merge(annovar, funct_descr, all.x=T)
annovar<-annovar[,-which(names(annovar)=="gene_id")]

# Make rank with best variants first 
h1<-paste(annovar$ref,"/",annovar$ref,sep="")
h2<-paste(annovar$alt,"/",annovar$alt,sep="")
het<-paste(annovar$ref,"/",annovar$alt,sep="")

annovar$rank<-4
annovar$rank[which(annovar$mutant_bulk==h2 & annovar$wild_type_bulk==het & annovar$tomato_parent==h1)]<-3
annovar$rank[which(annovar$mutant_bulk==h2 & annovar$wild_type_bulk==het & annovar$tomato_parent==h1 & annovar$funct_consequ%in%c("nonsynonymous SNV"))]<-2
annovar$rank[which(annovar$mutant_bulk==h2 & annovar$wild_type_bulk==het & annovar$tomato_parent==h1 & annovar$funct_consequ%in%c("stopgain", "stoploss", "frameshift_insertion", "frameshift deletion", "splicing"))]<-1
annovar<-annovar[order(annovar$rank, annovar$pos),]

write.table(annovar, output_file, row.names=F, quote=F, sep="\t")

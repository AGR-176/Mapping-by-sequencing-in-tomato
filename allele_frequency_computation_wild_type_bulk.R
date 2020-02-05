
# Open VCF file in R
my_vcf_file<-"variant_calling_filtered.vcf"
vcf<-read.delim("variant_calling_filtered.vcf",comment.char = "#",header=F)
head_vcf <- system(paste("head -n 100 ",my_vcf_file," | grep '#CHROM' ",sep=""),intern=T)
head_vcf<-unlist(strsplit(head_vcf,"\t"))
head_vcf<-gsub("#","",head_vcf)
names(vcf)<-head_vcf
nrow(vcf)


# Obtain names of fields in genotype column
format_names<-as.character(vcf$FORMAT[1])
format_names<-unlist(strsplit(as.character(format_names),":"))

GTcolumn2matrix<-function(dd,format_names){
        temp<-data.frame(matrix(unlist(strsplit(as.character(dd),":")),ncol=5,byrow=T))
        names(temp)<-format_names
        temp$GQ<-as.numeric(as.character(temp$GQ))
        temp$DP<-as.numeric(as.character(temp$DP))
        temp$GT<-as.character(temp$GT)
        temp$PL<-as.character(temp$PL)
        temp$AD<-as.character(temp$AD)
        temp
}

str(vcf)

# Make genotypes character to make next loop possible (otherwise invalid factor level error)
for (i in 10){
        vcf[,i]<-as.character(vcf[,i])
        print(i)
}
str(vcf)
head(vcf)

# Get name columns and select wild_type_bulk genotype for allele frequency calculation
names(vcf)
xx<-as.character(vcf[,12])
wild_type_bulk<-GTcolumn2matrix(xx,format_names)
nrow(wild_type_bulk)
names(wild_type_bulk)
head(wild_type_bulk)
tail(wild_type_bulk)

format_names<-c("AD1","AD2")
ADcolumn2matrix<-function(dd,format_names){
        temp<-data.frame(matrix(unlist(strsplit(as.character(dd),",")),ncol=2,byrow=T))
        names(temp)<-format_names
        temp
}

xx<-wild_type_bulk[,2]
strsplit(as.character(xx),",")
length(strsplit(as.character(xx),","))

wild_type_bulk_AD<-ADcolumn2matrix(xx,format_names)
head(wild_type_bulk_AD)
tail(wild_type_bulk_AD)
nrow(wild_type_bulk_AD)

add_col<-vcf[,c("CHROM","POS")]
head(add_col)
nrow(add_col)

genot<-wild_type_bulk[,c("GT")]
head(genot)

x1<-cbind(add_col,wild_type_bulk_AD)
x<-cbind(x1,genot)
head(x)
nrow(x)
write.csv2(x,file="wild_type_bulk_allele_freq.csv")
x<-read.csv2("wild_type_bulk_allele_freq.csv", header=T)
head(x)
x$ratio<- x$AD2 / (x$AD1 + x$AD2) 
write.csv2(x,file="wild_type_bulk_allele_freq.csv")

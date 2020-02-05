
# Open genotyping file filtered by chromosomal region of interest 
# (in the example shown in Figure 2: chromosome 05, positions from 0 bp to 3000000 bp) 
x<-read.delim("variant_calling_ch09_0_3000000.csv", header=T)
str(x)
head(x)
nrow(x)
names(x)

## Look for candidate mutations 
x2<-x[which(x$tomato_parent=="0/0" & x$wild_type_bulk=="0/1" & x$mutant_bulk=="1/1"),]
head(x2)
nrow(x2)
write.csv2(x2,file="candidate_mutations.csv")

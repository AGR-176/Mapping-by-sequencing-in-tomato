library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape2)

# Open allele frequency computation files
x <- read.delim("wild_type_bulk_allele_freq.csv",sep=";", dec=",")
head(x)
tail(x)
dim(x)

y <- read.delim("mutant_bulk_allele_freq.csv",sep=";", dec=",")
head(y)
tail(y)
dim(y)

# Check that both files contain the same list of variants
which(paste(x$CHROM,x$POS)!=paste(y$CHROM,y$POS))

# Merge both files

names(x)[8]<-"ratio.x"
names(y)[8]<-"ratio.y"

tt<-merge(x[,c("CHROM", "POS", "ratio.x")],y[,c("CHROM", "POS", "ratio.y")])
head(tt)

# Calculate the average allele frequency ratios in 1 Mb windows
summary(tt$POS)
tt$window<-round(tt$POS/1000000)
summary(tt$window)
t2<-ddply(tt, .(CHROM, window), summarize, avg_freq.x=mean(ratio.x), avg_freq.y=mean(ratio.y))
head(t2)

# Plot average average allele frequency ratios along tomato chromosomes
t2<-melt(t2, id.vars=c("CHROM", "window"))
head(t2)
t2<-t2[-which(t2$CHROM=="SL2.50ch00"),] # Delete Chr0
write.csv(t2, "mi_tabla.csv", row.names=FALSE, quote=FALSE)

p1<-ggplot(t2,aes(x=window, y=value, colour=variable)) + geom_line(size = rel(0.3)) + 
        geom_hline(yintercept=0.5, lty=3) + facet_wrap(~CHROM, scales="free_x", ncol=4) + 
        theme(
                strip.background = element_rect(colour="grey50", fill="grey50"),
                strip.text.x = element_text(size=10, colour = "white"),
                axis.text.y = element_text(size=10),
                axis.text.x = element_text(size=10),
                legend.position = "bottom",
                panel.grid.minor.y =  element_blank(),
                panel.grid.major.x = element_line(colour = "grey50", linetype = "dashed", size = rel(0.1)),
                panel.border = element_rect(colour = "grey50", linetype = "solid", fill = NA),
                panel.background = element_rect(fill = NA) 
        ) + scale_colour_manual(values=c("red","darkblue")) + ylim(c(0,1)) + labs(y="Average allele frequency", x="Position (Mb)", title="Mapping-by-sequencing")                        

pdf("allele_frequency_comparison.pdf", height=6, width=8)
p1
dev.off()

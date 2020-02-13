#################################################
#CHANGE THESE PARAMETERS AS NEEDED

# Input file
input_file_name <- "filtered_variants.wild_type_and_mutant_bulks.tsv"

# Output file
output_file_name <- "allele_frequency_comparison.pdf"
pdf_height<-6
pdf_width<-8

# Size of window to calculate average allele frequencies
window_size_in_bp <- 1000000

#################################################

# Load required packages, or install if not present
packages <- c("ggplot2","plyr","reshape2")
package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
    }
})

# Open allele counts 
x <- read.delim(input_file_name)

# Find columns of interest
chr_col<-grep("CHROM$", names(x))
pos_col<-grep("POS$", names(x))
mut_cols<-grep("mutant_bulk", names(x))
wt_cols<-grep("wild_type_bulk", names(x))

# Remove SL2.50ch00
if(length(which(x[,chr_col]=="SL2.50ch00"))>0){
x<-x[-which(x[,chr_col]=="SL2.50ch00"),]
}

# Calculate frequency of alternate allele for each SNP
x$ratio_mut <- x[,mut_cols[2]] / (x[,mut_cols[1]]+x[,mut_cols[2]])
x$ratio_wt <- x[,wt_cols[2]] / (x[,wt_cols[1]]+x[,wt_cols[2]])

# Calculate the average allele frequency ratios in windows: group positions in windows
x$window<-round(x[,pos_col]/window_size_in_bp)
# Calculate the average allele frequency ratios in windows: calculate average per window
names(x)[chr_col]<-"CHROM"
dd<-ddply(x, .(CHROM, window), summarize, avg_freq_mut=mean(ratio_mut), avg_freq_wt=mean(ratio_wt))

# Reformat for plotting
dd<-melt(dd, id.vars=c("CHROM", "window"))

# ggplot
plot_1<-ggplot(dd,aes(x=window, y=value, colour=variable, group=variable)) + geom_line(size = rel(0.3)) + 
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

pdf(output_file_name, height=pdf_height, width=pdf_width)
plot_1
dev.off()


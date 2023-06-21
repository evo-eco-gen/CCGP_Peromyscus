#########################################################
#### TITLE: Plot cumulative coverage for 6 Peromyscus genomes.
#### Script adapted from the plotting code by Meixi Lin 2022 
#### (doi.org/10.1093/jhered/esac031). 
#### In this version, I simplify the code to work with simple vectors of scaffold lengths 
#### parsed from GenBank reports using command line tools.

date()

### Read in input data (parsed from a GenBank Assembly Structure report): 
### each column is a list of scaffold lengths for a genome.
ngx<-read.table("ngx.txt", header=T, sep="\t")
species <- as.vector(colnames(ngx))

### Create a list to store output tables for each species.
datalist <- vector("list", length(species))

### Define the function to format the input data for one species
get_cumsum <- function(input_vector, species_name) {

    # get sum length
    genomesize = sum(input_vector)
    
    # sort by length
    Scaffold_Lengths<-sort(input_vector, decreasing=TRUE)
    outdt <- as.data.frame(Scaffold_Lengths)
    
    # get cumulative sum length
    cumsum = c()
    for (i in 1:nrow(outdt)) {
        tempsum = sum(input_vector[1:i])
        cumsum = c(cumsum,tempsum)
    }
    outdt$Cumulative_Length = cumsum
    
    #Format output table
    outdt$Cumulative_Coverage = outdt$Cumulative_Length/genomesize
    outdt<-cbind(outdt, Genome_Name=species_name)
    outdt = outdt[,c('Scaffold_Lengths','Cumulative_Coverage','Genome_Name')]
    
    # add the starting positions
    startline = outdt[1,]
    startline$Cumulative_Coverage = 0
    outdt = rbind(startline, outdt)
    return(outdt)
}


### Apply the function to all species.
for (i in 1:length(species)) {
    print(species[i])
    datalist[[i]] <- get_cumsum(ngx[,i], species[i])
}

all_species<- do.call(rbind, datalist)
dim(all_species)


######
### Use a palette readable to a colorblind person.
display.brewer.all(colorblindFriendly = TRUE)

### Create the plot
Peromyscus_plot <- ggplot(all_species,aes(x = Cumulative_Coverage, y = Scaffold_Lengths/1e+6, color = Genome_Name)) +
    scale_color_brewer(palette = 'Paired') +
    geom_step(direction = 'vh') +
    geom_vline(xintercept = 0.5, linetype = 'dotted') +
    annotate("text", x=0.54,y=220, label="N50") +
	labs(x = 'Cumulative Coverage',
    y = 'Scaffold Size (Mb)') +
    scale_x_continuous(expand = expansion(mult = c(0, .05))) + 
    theme_bw() + theme(axis.text = element_text(size = 12)) + 
    theme(legend.background = element_rect(fill="lightgrey", linetype="solid")) + 
    guides(color=guide_legend(title="(sub)species"))
    
### Output the plot:

ggsave(filename = 'Peromyscus.Ngx.pdf', plot = Peromyscus_plot, height = 5, width = 5)

pdf(file = 'Peromyscus.Ngx.pdf', height = 5, width = 5)
print(Peromyscus_plot)
dev.off()

tiff(file = 'Peromyscus.Ngx.tiff', height = 5, width = 7, compression="none")
print(Peromyscus_plot)
dev.off()

### Save the workspace
date()
save.image(file = 'Peromyscus.NGx.RData')
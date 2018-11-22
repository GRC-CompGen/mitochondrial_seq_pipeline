---
title: "Mitochondrial genome Circos plot"
author: "Nicholas Harvey"
date: "November 22, 2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```



##Load Required variant frequency data
```{R, echo=TRUE}
require(ggplot2)

SNP_freq <- read.csv('H:/mito_methods_paper/NI_mito_freq.csv')

SNP_freq <- as.data.frame(SNP_freq)
```

Subsection genomic locations and frequency columns from the data (mitodata), create a mitochondrial genome table (lines) for every position. 
We then need to merge the mitodata file with the lines file to get frequency inormation for every position in the mitochondrial genome (this function places NA at unknown values so we replace NA with 0)
```{r, echo=TRUE}
headers <- c("chromStart", "freq")
mitodata <- SNP_freq[headers]

head(mitodata)

lines <- data.frame(x = seq(1,16567,by=1))

df1 <- data.frame(x=lines$x)
#colnames(df1)[1] <- "x"

df2 <- data.frame(x=mitodata$chromStart, y = mitodata$freq)

mitodata <- merge(df1, df2, by.x = "x", by.y = "x", all.x = TRUE)

mitodata[is.na(mitodata)] <- 0

mitodata
```

##Add gene label function to input gene symbols into lines and mitodata tables

This step links mitochondrial genomic position with gene symbols and then gene symbols to colours.
Visible boudries vector used to label circos plot with gene location boundries
```{R, echo = TRUE}
addgenelabel <- function(snp,gene) { gene <- ifelse(snp < 577,gene <- "Control-Region", ifelse(snp < 648,gene <- "tRNA", ifelse(snp < 1602,gene <- "rRNA", ifelse(snp < 1671,gene <- "tRNA", ifelse(snp < 3230,gene <- "rRNA", ifelse(snp < 3305,gene <- "tRNA", ifelse(snp < 3307,gene <- "Non-Coding", ifelse(snp < 4263,gene<- "ND1", ifelse(snp < 4332,gene <- "tRNA", ifelse(snp < 4401,gene <- "tRNA", ifelse(snp < 4402,gene <- "Non-Coding", ifelse(snp < 4470,gene <- "tRNA", ifelse(snp < 5512,gene <- "ND2", ifelse(snp < 5580,gene <- "tRNA", ifelse(snp < 5587,gene <- "Non-Coding", ifelse(snp < 5656,gene <- "tRNA", ifelse(snp < 5657,gene <- "Non-Coding", ifelse(snp < 5730,gene <- "tRNA", ifelse(snp < 5826,gene <- "tRNA", ifelse(snp < 5892,gene <- "tRNA", ifelse(snp < 5904,gene <- "Non-Coding", ifelse(snp < 7446,gene <- "CO1", ifelse(snp < 7515,gene <- "tRNA", ifelse(snp < 7518,gene <- "Non-Coding", ifelse(snp < 7586,gene <- "tRNA", ifelse(snp < 8270,gene <- "CO2", ifelse(snp < 8295,gene <- "Non-Coding", ifelse(snp < 8365,gene <- "tRNA", ifelse(snp < 8366,gene <- "Non-Coding", ifelse(snp < 8573,gene <- "ATP8", ifelse(snp < 9208,gene <- "ATP6", ifelse(snp < 9991,gene <- "CO3", ifelse(snp < 10059,gene <- "tRNA", ifelse(snp < 10405,gene <- "ND3", ifelse(snp < 10470,gene <- "tRNA", ifelse(snp < 10767,gene <- "ND4L", ifelse(snp < 12138,gene <- "ND4", ifelse(snp < 12207,gene <- "tRNA", ifelse(snp < 12266,gene <- "tRNA", ifelse(snp < 12337,gene <- "tRNA", ifelse(snp < 14149,gene <- "ND5", ifelse(snp < 14674,gene <- "ND6", ifelse(snp < 14743,gene <- "tRNA", ifelse(snp < 14747,gene <- "Non-Coding", ifelse(snp < 15888,gene <- "CYB", ifelse(snp < 15954,gene <- "tRNA", ifelse(snp < 15956,gene <- "Non-Coding", ifelse(snp < 16024,gene <- "tRNA", ifelse(snp < 17000,gene <- "Control-Region") ))))))))))))))))))) ))))))))))))))))))) ))))))))) ) }

colours <- c("Control-Region" = "lightblue", "tRNA" = "magenta4", "rRNA" = "mediumaquamarine", "Non-Coding" = "sienna4", "ND1" = "magenta", "ND2" = "mediumblue", "CO1" = "olivedrab", "CO2" = "orange2", "ATP8" = "orchid4", "ATP6" = "red3", "CO3" = "royalblue2", "ND3" = "palegreen4", "ND4L" = "grey0", "ND4" = "pink4", "ND5" = "yellow4", "ND6" = "steelblue4", "CYB" = "tan", "red")

visibleboundaries <- c(576,1601,3229,4262,5511,7445,8269,9207,9990,10404,10766,12137,14148,14673,15887)

bdries <- data.frame(x = visibleboundaries,y= 0.75)

bdries$gene <- addgenelabel(bdries$x,bdries$gene)
```


##Gene labels and frequency data processing
The gene labels above are added to a new "gene" column in the lines and mitodata tables.
as variant frequencies cover a wide range (0-n samples), we need to compress the data. We first subtract 0.5 from the 0 values to avoid logging 0 (mitodata$freq_values). by negative logging the frequency values we are able to get manageable data for the plot (mitodata$actual_values) 

```{R, echo = TRUE}
mitodata$gene <- addgenelabel(mitodata$x,mitodata$gene)

lines$gene <- addgenelabel(lines$x,lines$gene)

mitodata$freq_values <- (mitodata$y - 0.5)

mitodata$actual_values <- -1*log10(mitodata$freq_values)

str(mitodata)
```


vectors for scaling the colour manually are created to avoid long lines of code when we plot
```{R, echo = TRUE}
break_genes <- c("Control-Region","tRNA","rRNA","Non-Coding","ND1","ND2","CO1","CO2","ATP8","ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")
label_genes <- c("Control-Region","tRNA","rRNA","Non-Coding","ND1","ND2","CO1","CO2","ATP8","ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")
```

##Plotting the circos plot with ggplot2
Here we set the theme of the plot and then provide a series or aesthetics and layers to visualise.
Firstly, we call ggplot to use the x and y values from mitodata and colour them according to "gene". geom_point classifies the data as single points and a direction of 1 means the plot will be in a clockwise rotation. the scale_colour_manual function lets us define the colour scheme to use for both labels and gene breaks (classified above). We then insert a layer of second geom_point from the lines data frame. labels on the plot are set using the labs function. A layer using and classifying the gene boudries is then added.

We then create the "solid" mitochondrial genome also using geom_point but making the points thicker (lwd = 5.5), and set the loction using a value instead of the y column in the lines data.

Lastly, we use theme_replace to adjust specific points of the chosen theme. In this instance, we have ensured that the x-axis title and text are blank elements because we dont want any labels at that section.

```{R, echo = TRUE}
theme_set(theme_gray())

p <- ggplot(mitodata, aes(x = mitodata$x, y = mitodata$actual_values, colour = gene)) + 
  geom_point() + coord_polar(direction=1) + 
  scale_colour_manual(values = colours, "Genes", aesthetics = "colour", breaks = break_genes, labels = label_genes) +
  geom_point(aes(x, colour = gene), data=lines, inherit.aes = TRUE, lwd=1.5) + 
  labs(title = "Mitochondrial variant frequency in NI subpopulation", y = "Negative log10 Variant Frequency") +
  layer(mapping = aes(x,y,label = x), data = bdries, geom="text", stat = "identity", position = "identity")

p1 <- p + geom_point(aes(x,0.45, colour = gene), data = lines, inherit.aes = TRUE, lwd=3)

p2 <- p1 + theme_replace(axis.text.x = element_blank(), 
                  axis.title.x = element_blank())
p2
```

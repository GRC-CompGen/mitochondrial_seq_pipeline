---
title: "Untitled"
author: "Nicholas Harvey, Miles Benton"
date: "November 26, 2018"
output: html_document
---

 
# creating circular genome plots for mitocondrial data
## based on the code from Steven Turner: https://github.com/stephenturner/solarplot
### created: 2018-11-25
### modified: 2018-11-26

##Load Required variant frequency data

```{r}
# packages
require('tidyverse')

# load mt freq data
SNP_freq <- read.csv('lines.csv')

# filter out non-variant SNPs
mtSnps <- SNP_freq[SNP_freq$freq != 0,] # no SNP
mtSnps <- mtSnps[mtSnps$freq != 77,] # no SNP# all share these variants

```

Toggle the below to remove singleton variants
mtSnps <- mtSnps[mtSnps$freq != 1,]

## 'simulate' a p value
```{r}
mtSnps$p <- mtSnps$freq / 77
```
## create a snp col
```{r}
mtSnps$snp <- paste0("mt", mtSnps$chromStart)
```
## define new plotting object
```{r}
mitodata <- data.frame(snp = mtSnps$snp, bp = mtSnps$chromStart, freq = mtSnps$freq, p = mtSnps$p)
```
## Sets gene names for bp ranges and colours for each gene name
```{R}
addgenelabel <- function(bp,gene) { gene <- ifelse(bp < 577,gene <- "Control-Region", ifelse(bp < 648,gene <- "tRNA", ifelse(bp < 1602,gene <- "rRNA", ifelse(bp < 1671,gene <- "tRNA", ifelse(bp < 3230,gene <- "rRNA", ifelse(bp < 3305,gene <- "tRNA", ifelse(bp < 3307,gene <- "Non-Coding", ifelse(bp < 4263,gene<- "ND1", ifelse(bp < 4332,gene <- "tRNA", ifelse(bp < 4401,gene <- "tRNA", ifelse(bp < 4402,gene <- "Non-Coding", ifelse(bp < 4470,gene <- "tRNA", ifelse(bp < 5512,gene <- "ND2", ifelse(bp < 5580,gene <- "tRNA", ifelse(bp < 5587,gene <- "Non-Coding", ifelse(bp < 5656,gene <- "tRNA", ifelse(bp < 5657,gene <- "Non-Coding", ifelse(bp < 5730,gene <- "tRNA", ifelse(bp < 5826,gene <- "tRNA", ifelse(bp < 5892,gene <- "tRNA", ifelse(bp < 5904,gene <- "Non-Coding", ifelse(bp < 7446,gene <- "CO1", ifelse(bp < 7515,gene <- "tRNA", ifelse(bp < 7518,gene <- "Non-Coding", ifelse(bp < 7586,gene <- "tRNA", ifelse(bp < 8270,gene <- "CO2", ifelse(bp < 8295,gene <- "Non-Coding", ifelse(bp < 8365,gene <- "tRNA", ifelse(bp < 8366,gene <- "Non-Coding", ifelse(bp < 8573,gene <- "ATP8", ifelse(bp < 9208,gene <- "ATP6", ifelse(bp < 9991,gene <- "CO3", ifelse(bp < 10059,gene <- "tRNA", ifelse(bp < 10405,gene <- "ND3", ifelse(bp < 10470,gene <- "tRNA", ifelse(bp < 10767,gene <- "ND4L", ifelse(bp < 12138,gene <- "ND4", ifelse(bp < 12207,gene <- "tRNA", ifelse(bp < 12266,gene <- "tRNA", ifelse(bp < 12337,gene <- "tRNA", ifelse(bp < 14149,gene <- "ND5", ifelse(bp < 14674,gene <- "ND6", ifelse(bp < 14743,gene <- "tRNA", ifelse(bp < 14747,gene <- "Non-Coding", ifelse(bp < 15888,gene <- "CYB", ifelse(bp < 15954,gene <- "tRNA", ifelse(bp < 15956,gene <- "Non-Coding", ifelse(bp < 16024,gene <- "tRNA", ifelse(bp < 17000,gene <- "Control-Region") ))))))))))))))))))) ))))))))))))))))))) ))))))))) ) }


colours <- c("Control-Region" <- "lightblue4", "tRNA" <- "magenta4", "rRNA" <- "mediumaquamarine", "Non-Coding" <- "sienna4", "ND1" <- "magenta", "ND2" <- "mediumblue", "CO1" <- "olivedrab", "CO2" <- "orange2", "ATP8" <- "orchid4", "ATP6" <- "red3", "CO3" <- "royalblue2", "ND3" <- "palegreen4", "ND4L" <- "grey0", "ND4" <- "pink4", "ND5" <- "yellow4", "ND6" <- "steelblue4", "CYB" <- "tan","red")
```

# Add gene names to each SNP and display internal structure of mitodata
```{r}
mitodata$gene <- addgenelabel(mitodata$bp,mitodata$gene)
str(mitodata)
```
## Creates and stores negative log p as a new variable
```{r}
mitodata$neglogp <- -1*log10(mitodata$p)
```
### (optional) Adds significance threshold lines at negative log of 0.05 and -3 (only add this for actual p-value data)
```{r}
mitodata$neglogpline <- -1*log10(0.05)

mitodata$extraline <- -3
```


## Create gene boundaries and lines
Note: the bdries y value may be changed to place the genomic labels elsewere on the plot
```{r}
visibleboundaries <- c(1,576,1601,3229,4262,5511,7445,8269,9207,9990,10404,10766,12137,14148,14673,15887)

bdries <- data.frame(x = visibleboundaries,y=-.5)

bdries$gene <- addgenelabel(bdries$x,bdries$gene)

lines <- data.frame(x = seq(0,16567,by=1),y = 0)

lines$gene <- addgenelabel(lines$x,lines$gene)
```

vectors for scaling the colour manually are created to avoid long lines of code when we plot

```{R, echo = TRUE}
break_genes <- c("Control-Region","tRNA","rRNA","Non-Coding","ND1","ND2","CO1","CO2","ATP8","ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")
label_genes <- c("Control-Region","tRNA","rRNA","Non-Coding","ND1","ND2","CO1","CO2","ATP8","ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")
```

## Plot everything and GO


Firstly, we call ggplot to use the x and y values from mitodata and colour them according to "gene". geom_point classifies the data as single points and a direction of 1 means the plot will be in a clockwise rotation. the scale_colour_manual function lets us define the colour scheme to use for both labels and gene breaks (classified above). We then insert a layer of second geom_point from the lines data frame. labels on the plot are set using the labs function. A layer using and classifying the gene boudries is then added.the width of the data points may be changed with the lwd= fucntion. 

Option 1:
plot all above and generate solarplot with variant frequencies outside genomic area 

Note: The y value in the third line represents the "significance" threshold. This is not necessary unless handling actual association data

```{r}
ggplot(mitodata, aes(x = bp,y = neglogp,color = gene)) +
  geom_point()+ coord_polar(direction = -1) +
  geom_line(aes(x,1.30,color = "red"),data = lines) +
  #facet_grid(.~pheno) +
  geom_line(aes(y=extraline)) +
  geom_point(aes(x,y,color = gene),data=lines) +
  scale_colour_manual(values = colours,"Genes",breaks = break_genes,labels = label_genes) +
  xlab("Mitochondrial Base-Pair Location") +
  ylab("-log(p-value)") +
  ggtitle("Negative Log P-value of Mitochondrial  Variant Frequencies") +
  geom_text(mapping = aes(x, y, label = x), data = bdries, size = 2.5)

ggsave("solarplot.png", w=6, h=6, dpi=110)
```

Option 2:
Here we set the theme of the plot and then provide a series of aesthetics and layers to visualise.
We then create the "solid" mitochondrial genome also using geom_point but making the points thicker (lwd = 5.5), and set the loction using a value instead of the y column in the lines data.

Lastly, we use theme_replace to adjust specific points of the chosen theme. In this instance, we have ensured that the x-axis title and text are blank elements because we dont want any labels at that section.

```{r}
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
#/END

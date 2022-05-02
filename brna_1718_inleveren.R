# brna eindopdracht 1e kans 1718 
# Dirk van der Torre s1089166
# Robin Zanoni s1080232

######################
## 1: voorbereiding ##
#####################

# 1.1: installeren libraries #
##############################
# tik "sudo R" in in een unix terminal.
# Neem daarna deze commands over:
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("edgeR") (voor cpm())
# biocLite("gplots")
# biocLite("RColorBrewer")
# biocLite("org.Hs.eg.db")
# biocLite("GO.db")

# 1.2 imporeteren libraries #
#############################
library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(org.Hs.eg.db)
library(GO.db)

# 1.3 inlezen data #
####################

annotatie <- read.csv(file="Data/annotables_unique_grch37(1).csv", header = TRUE)
metadata <- read.csv(file="Data/chronische_inflammatie_metadata.csv", header = TRUE)
rawcounts <- read.csv(file="Data/chronische_inflammatie_rawcounts(1).csv", header = TRUE)

# 1.4 bewerken data #
#####################

# cpm wil alleen numerieke colomen, dus halen we eerste kolom weg.
# Namen van genen kunnen nog uit rawcounts worden gehaald.
countdata <- rawcounts[,-1]
rownames(countdata) <- rawcounts[,1]

# Kleine spiekbrief voor alle sample namen.
# groep 1 = SRR1039508 (control) en SRR1039509 (treated) (celtype N61311)
# groep 2 = SRR1039512 (control) en SRR1039513 (treated) (celtype N052611)
# groep 3 = SRR1039516 (control) en SRR1039517 (treated) (celtype N080611)
# groep 4 = SRR1039520 (control) en SRR1039521 (treated) (celtype N061011)

# Kolom namen veranderen.
# c staat voor control en t voor treated.
colnames(countdata) <- c("N61311.c", "N61311.t", "N052611.c", "N052611.t", "N080611.c", "N080611.t", "N061011.c", "N061011.t")


##########################
## 2 inhoudelijke eisen ##
##########################

### 2.1 filteren op lage expressie ###
######################################

# Om verschillen in library grootte tegen te gaan worden log2 counts/10⁶ gebruikt
# cpm (counts per million) = aantal gemapte reads vergeleken met alle reads die zijn gesequenced.
# Nogmaals: je normaliseert dus voor verschil tussen het totaal aantal reads per sample.
CountsPM <- cpm(countdata)

# Er voor zorgen dat er maar 1 figuur in het scherm komt. 
par(mfrow=c(1,1))
# kijken of een threshold van 0.5 cpm overeen komt met ongeveer 10 reads in sample 1.
plot(CountsPM[,1], countdata[,1],ylim=c(0,50),xlim=c(0,3), xlab="cpm", ylab="expressie count")
abline(v=0.5)
# Bij een cpm van minimimCPM (ongeveer 0,5) zijn er komen er ongeveer 10 genen 
# in de eerste sample tot expressie, wat een redelijk minimum is.

# Bepalen welke genen genen in de samples een hogere cpm hebben van 0.5.
thresh <- CountsPM > 0.5 

# Kijken hoeveel TRUE's er zijn per sample, misschien kunnen we samples zonder TRUE's eruit halen.
# Anders gezegd, als genen bijna niet tot uiting zijn gekomen in meerdere samples, dan is het gen niet interesant.
table(rowSums(thresh))
# Alle samples hebben genen die een hogere expressie hebben dan de threshold.
# We kunnen dus alle samples houden.

# Genen die in minder dan 2 samples een 0.5 cpm hebben zijn niet interesant genoeg om te houden.
# Daarom zoeken we genen waarbij minimaal 2 saples boven de 0.5 zitten.
keep = rowSums(thresh) >= 2
# Deze filtering moet nu nog worden toegepast.
counts.keep <- countdata[keep,]
summary(keep)

# Maken DGEList object, dit wordt gebruikt om de count data op te slaan.
countsDGE <- DGEList(counts.keep)
countsDGE$counts
countsDGE$samples

### 2.2 QC op count data ###
############################

# De bedoeling is dat we zoeken wat de meeste variatie veroorzaakt.

# Aantal reads per sample bekijken.
countsDGE$samples$lib.size

# Barplot om de verdeling van samples te bekijken.
barplot(countsDGE$samples$lib.size,names=colnames(countsDGE),las=2)
# Titel geven aan plot.
title("Barplot of library sizes")
# Omdat count data niet normaal verdeeld is moeten de we counts omzetten in logcounts.
logcounts <- cpm(CountsPM,log=TRUE)

# Boxplot maken om de verdeling van de samples te zien.
# Er zijn outliers in een range van 2 tot 13 miljoen, wat de boxplot onleesbaar maakt.
# Outliers kunnen uit de boxplot met "outline=FALSE".
boxplot(logcounts, xlab="", ylab="log2 counts / 10^6",las=2, outline=FALSE)
# Te zien is dat de eerste kwart van de data redelijk wat verschil heeft in expressie niveaus.
# Het tweede kwart heeft ook redelijk wat variatie in niveaus.
# De derde en vierde kwartalen heeft bijna geen variatie. 

# Kijken naar de verdeling van de data.
# Multidimensionale scaling plot (plotMDS)
par(mfrow=c(1,1))
plotMDS(countsDGE)
# kijken naar sample info
levels(metadata$dex)
metadata
# de vergelijkingen:
# id controle id threated
# SRR1039508  SRR1039509
# SRR1039512  SRR1039513
# SRR1039516  SRR1039517
# SRR1039520  SRR1039521
countsDGE$samples

# Plot maken waarbij de samples een kleur krijgen die overeenkomt met het controle en threated.
col.behandeling <- c("purple", "orange")[metadata$dex]
data.frame(metadata$dex,col.behandeling)
# Gelukt, controle wordt straks paars, behandeld oranje
# Nu weer een plot maken van MDS maar dan met de kleuren voor de behandeling er bij
par(mfrow=c(1,1))

plotMDS(countsDGE, col=col.behandeling)
legend("topleft",fill=c("purple", "orange"),legend=levels(metadata$dex))
title("behandeling samples")

# dezelfde plot maken maar dan kleuren voor celltype
col.celtype <- c("red", "green", "black", "blue")[metadata$celltype]
plotMDS(countsDGE, col=col.celtype)
legend("bottomright",fill=c("red", "green", "black", "blue"),legend=levels(metadata$celltype))
title("cel type")
par(mfrow=c(1,1))

# Te zien is dat  de controle en treated steeds hetzelfde verschil geeft tussen de samples.
# Visueel te zien: de paarse en oranje groep (van de behandelde samples plot) lijken erg op elkaar.
# De controle vs treated geeft dus niet zo veel variatie tussen de samples.

# Maar als je gaat kijken naar de cel typen zie je de kleuren totaal niet bij elkaar klusteren.

# Dit geeft aan dat er meer variatie is tussen de verschillende cel typen, dan tussen de behandelingen.
# De cel typen zijn dus de eigenschappen waar we op moeten letten om groot verschil in expressie te vinden.


### 2.3 Hierarchisch clusteren (met heatmaps) ###
#################################################

# Alternatief op plotMDS om clustering van saples te ontdekken.
# moeten kiezen hoeveel genen we willen bekijken voor de heatmap plot.
# laten we nu even van 500 uitgaan.

# countsCPM zijn logcounts (log van cpm van de countdata)
# Variantie bepalen voor elk gen
var_genes <- apply(CountsPM, 1, var)
head(var_genes)
# Maar we willen alleen de gennamen van de eerste 500 meest variabele.
select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
# Even kijken welke bovenaan zitten
head(select_var)
# Logcounts van de 500 meest variabele genen verkrijgen
highly_variable_lcpm <- logcounts[select_var,]
# Als alles goed is hebben we 500 genen en 8 samples
dim(highly_variable_lcpm)
# [1] 500   8
# Dit klopt dus.
# Even kijken naar de data.
head(highly_variable_lcpm)

# Kleuren verkrijgen.
# Volgens de site moet je 1 kleur minder hebben dan dat je samples hebt.
mypalette <- brewer.pal(7, "RdYlBu")
morecols <- colorRampPalette(mypalette)
# We willen kijken naar de behandelingen dus moeten dex selecteren.
col.cell <- c("purple", "orange")[metadata$dex]
# paars = control, oranje = behandeld.

# Heatmap plotten.
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),
          trace="none", main="500 meest variable genen tussen samples",
          ColSideColors=col.cell, scale='row')

# Wanneer een horizontale "strook"(geclusterde genen)
# oranje is in de controle groep, dan is de strook meestal
# blauw in de behandelde groep.
# Dat is een redlijk verschil in variantie.
N61311 <- highly_variable_lcpm[,1]/highly_variable_lcpm[,2]
N052611 <- highly_variable_lcpm[,3]/highly_variable_lcpm[,4]
N080611 <- highly_variable_lcpm[,5]/highly_variable_lcpm[,6]
N061011 <- highly_variable_lcpm[,7]/highly_variable_lcpm[,8]

test = data.frame(N61311, N052611, N080611, N061011)

# Als alles goed is hebben we 500 genen en 4 samples
dim(test)
# [1] 500   4
# Dit klopt dus.
# Even kijken naar de data.
head(test)

# Kleuren verkrijgen
# Volgens de site moet je 1 kleur minder hebben dan dat je samples hebt(??)
mypalette <- brewer.pal(3, "RdYlBu")
morecols <- colorRampPalette(mypalette)
# Ee willen kijken naar de behandelingen dus moeten dex selecteren
col.cell <- c("purple", "orange", "green", "blue")[metadata$celltype]
col.cell <- c("blue", "purple", "green", "orange")
# paars = control, oranje = behandeld.

# heatmap plotten
heatmap.2(as.matrix(test),col=rev(morecols(50)),
          trace="none", main="500 meest variable genen tussen celtypen",
          ColSideColors=col.cell, scale='row')
# De samples hebben een duidelijk verschil in logcounts per million.

### 2.4 Normalisatie ###
########################

# De grootte van de library's wordt genormaliseerd
# (We waren al genormaliseerd op (library grootte), en nu worden ze genormaliseerd op (composition bias)).

# TMM normalisatie
countsDGEnorm <- calcNormFactors(countsDGE) 
# De normalisatie factoren zitten redelijk rond de 1.
countsDGEnorm$samples
# Dat komt omdat de library grootten al redelijk en de buurt van elkaar kwamen.
countsDGE$samples

# Scatterplots maken van library's met hoogste en laagste normalisatie factor
# Nog een keer kijken naar de normalisatie factoren
countsDGEnorm$samples[3]
# Grootste en kleinste normalisatiefactor eruit halen.
# Code kan nog wat mooier zodat automatisch de ID van min en max eruit wordt gehaald
max(countsDGEnorm$samples[3]) # N61311.c
min(countsDGEnorm$samples[3]) # N052611.t
maxNum = which(colnames(logcounts)=="N61311.c")
minNum = which(colnames(logcounts)=="N052611.t")
# Kijken hoe de data van de grootste en kleinste norm.factor eruit ziet,
# logcounts gebruikt omdat deze nog niet genormaliseerd zijn voor compositie bias, wel voor library grootte.
par(mfrow=c(1,2))
plotMD(logcounts, column=maxNum)
abline(h=0,col="grey")
plotMD(logcounts, column=minNum)
abline(h=0,col="grey")
# zo goed als geen verschil

par(mfrow=c(1,2))
plotMD(countsDGEnorm, column=maxNum)
abline(h=0,col="grey")
plotMD(countsDGEnorm, column=minNum)
abline(h=0,col="grey")
# Vergeleken met logcounts heeft countsDGEnorm een miniscule verbetering gegeven.
# Dit geeft aan dat de samples nu iets meer bij elkaar in de buurt komen.

### 2.5 Differential Expressed Gene testing ###
###############################################

# Clue van QC uit 2.2: Er is meer variatie tussen cel typen, dan tussen behandeling
# Daarom kijken naar verschil tussen cellen.

### 2.5.1 Design matrix maken ###
##################################

# Design matrix maken
subject <- factor(metadata$celltype)
treat <- factor(metadata$dex)
design <- model.matrix(~0+subject+treat)

# Namen zijn raar, maakt niet uit, het gaat er om dat de samplesbij elkaar
# worden gezet en wordt aangegven welke treated is (laaste kolom).
# Maar we gaan de namen toch lekker veranderen, want dat maakt alles duidelijker.
colnames(design) <- c("N052611", "N061011", "N080611", "N61311", "treated")


### 2.5.2 Testen voor diff. gen expr. ###
#########################################

# Voom transformatie van de data:
# Voom past de libraray grootes aan aan de normalisatie
# factoren die al zijn berekend in countsDGEnorm$samples.
# Voom transformatie geeft een Elist object met hulp
# van de design matrix. De gemiddelde variantie trend
# kan worden gevisualiseerd met een lijn.
# De plot geeft aan of er genen zijn die variabel eruit zien 
# in de data en of lage counts er goed uit zijn gefilterd.
par(mfrow=c(1,1))
v <- voom(countsDGEnorm, design, plot=TRUE)
v$design


fit <- lmFit(v)
names(fit)

# empirical Bayes shrinkage on the variances, and estimates moderated t-statistics and the associated p-values.
fit <- eBayes(fit)
# generate a quick summary of DE genes for the contrasts.
results <- decideTests(fit)
summary(results)

### 2.5.3 Annotatie toevoegen ###
#################################

head(annotatie)
# Genen uit annotatie halen die we nodig hebben.
useme = annotatie[annotatie$ensgene %in% rownames(fit),]
# Kijken of alle waarden TRUE zijn (zijn alle ensgene ID's gevonden?).
all(useme$ensgene==rownames(fit))
# Ja, allemaal gevonden, opslaan in fit.
fit$genes <- useme

# Even kijken naar enzymen die het significants aanwezig zijn.
topTable(fit, coef="treated",sort.by="p")

## Alle resultaten opslaan.
# limma.res.treated <- topTable(fit, coef="treated",sort.by="p",n="Inf")
## resultaten naar een bestand schrijven
#write.csv(limma.res.treated,file="treatedResults.csv",row.names=FALSE)

# Plots maken na DGE testing
# MA en vulcano plots maken.
par(mfrow=c(1,2))
# MA plot.
plotMD(fit,coef="treated",status=results[,"treated"], values = c(-1, 1))
# Genen die up gereguleerd zijn (1) zijn groen, down regulatie (-1) rood, geen verschil (0) zwart.

# Volcanoplot om significantie van fold change te bekijken.
# top 100 significante genen krijgen een hightlight (naam staat bij punt)
volcanoplot(fit,coef="treated",highlight=100,names=fit$genes$SYMBOL)
# Een hogere y coordinaat geeft een hogere significantie aan
# Deze 100 genen (geen idee hoe je de id's ophaald uit de vulcanoplot) zijn het meest significant.

#### 2.5.4 Testing relative to a threshold ####
###############################################

# When there is a lot of differential expression, sometimes we may want to cut-off 
# on a fold change threshold as well as a p-value threshold so that we follow up on 
# the most biologically significant genes. However, it is not recommended to simply 
# rank by p-value and then discard genes with small logFC’s, as this has been shown 
# to increase the false discovery rate.
# In other words, you are not controlling the false discovery rate at 5% any more.
# There is a function called treat in the limma package that performs this style of
# analysis correctly (McCarthy and Smyth 2009).
# treat will simply take our fit.cont object, as well as a user-specified log fold 
# change cut-off, and recalculate the moderated t-statistics and p-values with
# the new information about logFC.

# Let's decide that we are only interested in genes that have a absolute logFC of 1.
# This corresponds to a fold change of 2, or 0.5 (i.e. double or half).
# We can perform a treat analysis which ranks our genes according to p-value AND logFC.
# This is easy to do after our analysis, we just give the treat function the fit.cont object and specify our cut-off.
fit.treat <- treat(fit,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)

# Top genes from linear model fit
topTreat(fit.treat,coef=3)

par(mfrow=c(1,2))
# Kijken naar oude MA plot zonder cut-off.
plotMD(fit,coef="treated",status=results[,"treated"], values = c(-1, 1))
# Kijken naar oude MA plot met cut-off.
plotMD(fit.treat,coef="treated",status=res.treat[,"treated"])
abline(h=0,col="grey")
# Dit heeft veel ruis verwijderd.

### 2.6 GSEA/Pathway level analyse ###
######################################

# GSEA = gen/eiwit klassen detecteren die overrepresentatief zijn in de samples.
# pathway level analyse = de overrepresentatieve genen/eiwitten koppelen aan pathways om te kijken welke pathways anders zijn.

### 2.6.1 GSEA: Gene ontology testing with goana ###
####################################################
# Self-contained tests, which include the ROAST procedure, ask the question “Are the genes
# in the set/pathway differentially expressed as a whole?” 

# Competitive gene set tests, like goana and  camera ask the question whether 
# the differentially expressed genes tend to be over-represented in the gene set, 
# compared to all the other genes in the experiment. 
# These different questions use different statistical methodology.

# 10 meest significante up gereguleerde eiwit onthologien die gerelateerd zijn aan een biologisch proces.
gen_lengten = fit$genes$end - fit$genes$start + 1
go_length <- goana(fit, coef="treated",species = "Hs", covariate=gen_lengten, geneid="entrez")
topGO(go_length, ontology="BP", n=20)
# resultaat:
# Term Ont     N   Up Down         P.Up       P.Down
# GO:0042221                                             response to chemical  BP  2706  628  475 8.184667e-16 2.918390e-02
# GO:0007167                 enzyme linked receptor protein signaling pathway  BP   729  219  137 9.331050e-16 7.043901e-02
# GO:0050896                                             response to stimulus  BP  5551 1172  987 1.101079e-15 2.297995e-04
# GO:0007169 transmembrane receptor protein tyrosine kinase signaling pathway  BP   520  167   94 8.242853e-15 2.366843e-01
# GO:0007165                                              signal transduction  BP  3726  829  694 9.704034e-15 4.191474e-05
# GO:0023052                                                        signaling  BP  4051  891  759 1.437951e-14 5.662438e-06
# GO:0007154                                               cell communication  BP  4064  893  760 1.458908e-14 6.798134e-06
# GO:0070887                           cellular response to chemical stimulus  BP  2065  491  354 4.063544e-14 1.430072e-01
# GO:0051716                                    cellular response to stimulus  BP  4702 1003  842 1.398039e-13 6.222644e-04
# GO:0006629                                          lipid metabolic process  BP   996  267  169 2.232218e-13 3.289951e-01
# GO:0009888                                               tissue development  BP  1204  312  247 3.955324e-13 7.882289e-05
# GO:0048513                                         animal organ development  BP  2150  507  429 6.049820e-13 2.268003e-06
# GO:0044255                                 cellular lipid metabolic process  BP   771  215  127 9.425349e-13 4.996529e-01
# GO:0008150                                               biological_process  BP 11472 2161 1899 9.841516e-13 3.442495e-02
# GO:0009653                               anatomical structure morphogenesis  BP  1747  427  347 1.903690e-12 1.444220e-04
# GO:0010033                                    response to organic substance  BP  2125  490  368 5.365737e-12 8.365117e-02
# GO:0051179                                                     localization  BP  4417  937  731 6.479997e-12 4.362609e-01
# GO:0048856                                 anatomical structure development  BP  3712  807  718 1.134552e-11 2.373758e-08
# GO:0006810                                                        transport  BP  3492  757  546 1.440822e-11 9.286914e-01
# GO:0003012                                            muscle system process  BP   267   94   55 1.531733e-11 5.296110e-02

# Het hele onderzoek ging om chronische inflammatie in verschillende weefsels

# Van https://www.erasmusmc.nl/dermatologie/Inflammatoire-dermatologie/4387596/
# Hoe ontstaan inflammatoire dermatosen?
# Het afweersysteem of immuunsysteem bestaat uit een groep cellen en weefsels die het  
# lichaam beschermt tegen vreemde indringers zoals virussen, bacteriën of schimmels. 
# Het immuunsysteem voorkomt dat de indringers schade aanbrengen. Ook lichaamsvreemde stoffen zoals tumoren
# worden opgeruimd om het lichaam gaaf te houden. Als het afweersysteem niet alleen beschermt tegen
# vreemde indringers, maar ook lichaamseigen cellen aanvalt, dan gaat er iets mis. 
# Psoriasis is hier een goed voorbeeld van.

# Reacties die te maken hebben met het reactie van het lichaam/weefsel/cel zijn dus interesant.
# Genen die differentieel tot expressie komen in de volgende pathways kunnen interesant zijn: 
# GO:0042221 (response to chemical) 
# GO:0007167 (enzyme linked receptor protein signaling pathway)
# GO:0050896 (response to stimulus)

sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
# 
# Matrix products: default
# BLAS: /usr/lib/libblas/libblas.so.3.6.0
# LAPACK: /usr/lib/lapack/liblapack.so.3.6.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=nl_NL.UTF-8        LC_COLLATE=nl_NL.UTF-8    
#   [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
#   [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] GO.db_3.5.0          org.Hs.eg.db_3.5.0   AnnotationDbi_1.40.0 IRanges_2.12.0       S4Vectors_0.16.0    
#   [6] Biobase_2.38.0       BiocGenerics_0.24.0  RColorBrewer_1.1-2   gplots_3.0.1         edgeR_3.20.9        
#   [11] limma_3.34.9        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.16       bit_1.1-12         lattice_0.20-35    BiasedUrn_1.07     blob_1.1.1         caTools_1.17.1    
#   [7] tools_3.4.3        grid_3.4.3         KernSmooth_2.23-15 DBI_0.8            gtools_3.5.0       digest_0.6.15     
#   [13] bit64_0.9-7        bitops_1.0-6       memoise_1.1.0      RSQLite_2.1.0      gdata_2.18.0       compiler_3.4.3    
#   [19] locfit_1.5-9.1     pkgconfig_2.0.1   
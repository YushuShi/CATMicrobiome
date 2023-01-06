# CATMicrobiome

## General info
CATMicrobiome package performs conditional association test for the list of the taxa provided. The main function of the package is CAT.
	
## Usage

```
CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,
outcomeVar,tree=NULL,method="PERMANOVA",
nperm=9999,parallel=TRUE,nCore=2)
```

* **testList** A list of taxa to be tested.
* **otutable** An OTU table, rows correspond to OTUs, while columns correspond to observations. Row names must match the row names in taxonomy, while the column names must match the row names in metaData.
* **taxonomy** A taxonomy table, rows correspond to OTUs, while columns correspond to taxonomy assignment at different levels. Row names must match the row names in the OTU table.
* **metric** Can be one single metric or a vector of metrics(for MiRKAT only). It accepts "Weighted UniFrac", "Unweighted UniFrac" and all metrics that can be calculated from the "vegdist"" function in the R package "vegan". The default is "Weighted UniFrac".
* **metaData** The dataset containing the information of the study objects. Row names must match the column names in the OTU table.
* **outcomeVar** The outcome variable in metaData. If it is a vector of two elements, the function will consider the outcome is time-to-event (for now only compatible with MiRKAT) and use the first element in the list as the time variable and the second element in the list as the event indicator. For now it only accepts right censored data.
* **tree** An object of the class "tree", which is needed for calculating weighted and unweighted UniFrac distances.
* **method** Can be "PERMANOVA" or "MiRKAT". The default choice is "PERMANOVA".
* **nperm** Number of permutations used in the test. The default is 9999.
* **parallel** Whether or not use parallel computing.
* **nCore** How many cores to use for parallel computing. This is only relevant when parallel=TRUE.

## Output
A vector of p-values for the taxa list given.

## Details

CAT implements a novel permutation-based conditional association test, which can account for other features and phylogenetic relatedness when testing the association between a feature and an outcome. CAT adopts a leave-out method, measuring the importance of a feature in predicting the outcome by removing that feature from the data and quantifying how the association with the outcome is weakened through a two-proportion z-test. By pairing with PERMANOVA and MiRKAT-based methods, the package allows association testing for continuous, binary, and survival outcomes.

## Reference
Shi Y, Zhang L, Do KA, Jenq RR, Peterson CB (2022) _CAT: a conditional association test for microbiome data using a permutation-based approach_

## Examples

```

library(devtools)
install_github("YushuShi/CATMicrobiome")
library(CATMicrobiome) 
 
# Example 1 with PERMANOVA
otuPath<-system.file("extdata","GopalakrishnanOTUtable.csv", 
package = "CATMicrobiome")
otutable<-read.csv(otuPath,header=TRUE,row.names = 1)
taxonomyPath<-system.file("extdata","GopalakrishnanTaxonomy.csv", 
package = "CATMicrobiome")
taxonomy<-read.csv(taxonomyPath,header=TRUE,row.names = 1)
metaPath<-system.file("extdata","GopalakrishnanMeta.csv", 
package = "CATMicrobiome")
metaData<-read.csv(metaPath,header=TRUE,row.names = 1)
treePath<-system.file("extdata","GopalakrishnanTree.tree", 
package = "CATMicrobiome")
tree<-read.tree(treePath)

testList<-c("p__Firmicutes",
            "c__Clostridia",
            "o__Clostridiales",
            "f__Ruminococcaceae",
            "s__prausnitzii",
            "g__Faecalibacterium",
            "s__bromii",
            "g__Ruminococcus")

testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="ContOutcomes",
tree,nperm=9999,parallel=TRUE,nCore=2)
testResult
testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="ContOutcomes",
tree,nperm=9999,parallel=FALSE)
testResult

testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="BinOutcomes",
tree,nperm=9999,parallel=TRUE,nCore=2)
testResult
testResult<-CAT(testList,otutable,taxonomy,
metric="Weighted UniFrac",metaData,outcomeVar="BinOutcomes",
tree,nperm=9999,parallel=FALSE)
testResult

#Example 2 with survival outcomes
otuPath<-system.file("extdata","RiquelmeOTUtable.csv", 
package = "CATMicrobiome")
otutable<-read.csv(otuPath,header=TRUE,row.names = 1)
taxonomyPath<-system.file("extdata","RiquelmeTaxonomy.csv", 
package = "CATMicrobiome")
taxonomy<-read.csv(taxonomyPath,header=TRUE,row.names = 1)
metaPath<-system.file("extdata","RiquelmeMeta.csv", 
package = "CATMicrobiome")
metaData<-read.csv(metaPath,header=TRUE,row.names = 1)
treePath<-system.file("extdata","RiquelmeTree.tree", 
package = "CATMicrobiome")
tree<-read.tree(treePath)

testList<-c("Clostridia",
            "Clostridiales",
            "Lachnospiraceae",
            "Corynebacteriales",
            "Bacteroidia",
            "Bacteroidales",
            "Corynebacteriaceae",
            "Corynebacterium",
            "Pseudomonadales",
            "Pseudomonas")


testResult<-CAT(testList,otutable,taxonomy,
metric=c("bray","Unweighted UniFrac"),metaData,
outcomeVar=c("OS.Months","event"),tree=tree,
method="MiRKAT",parallel=TRUE,nCore=2)
testResult

testResult<-CAT(testList,otutable,taxonomy,
metric=c("bray","Weighted UniFrac","Unweighted UniFrac"),
metaData,outcomeVar=c("OS.Months","event"),
tree=tree,method="MiRKAT",parallel=FALSE)
testResult
```

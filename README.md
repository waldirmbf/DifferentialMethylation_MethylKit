# **Differential methylation analysis using the R package methylKit**
____
 The following documentation describes the R script used to estimate differentially-methylated cytosines between genetically-identical *K. marmoratus* individuals reared in enriched and poor environments with the R package methylKit. Methods are descibed in detail in **Berbel-Filho et al. (2020)**

```
> library('methylKit')
> packageVersion('methylKit')

### Create a list of all BAM files to get raw methylation data. ###

> file.list1 = (list(system.file("extdata", "PE-R03_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "PE-R04_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "PE-R08_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "PE-R09_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "PE-R12_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "PE-R14_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R05_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R06_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R07_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R08_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R10_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R11_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R12_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R13_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R14_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "RE-R15_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
 )


### Extract raw methylation call files (proportions of Cs and Ts per base) from .bam files aligned with Bismark. ###

> objs = processBismarkAln(location=file.list1,
  sample.id=list("PE-R03","PE-R04","PE-R08","PE-R09","PE-R12","PE-R14","RE-R05","RE-R06","RE-R07","RE-R08","RE-R10","RE-R11","RE-R12","RE-R13","RE-R14","RE-R15"),
  save.folder=NULL,save.context=NULL,read.context="CpG", nolap=FALSE,mincov=10,minqual=20,phred64=FALSE)

### Input methylation call files into methylKit. ###

> file.list2 = list(system.file("extdata", "PE-R03_CpG.txt", package = "methylKit"),
    system.file("extdata", "PE-R04_CpG.txt", package = "methylKit"),
    system.file("extdata", "PE-R08_CpG.txt", package = "methylKit"),
    system.file("extdata", "PE-R09_CpG.txt", package = "methylKit"),
	system.file("extdata", "PE-R12_CpG.txt", package = "methylKit"),
	system.file("extdata", "PE-R14_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R05_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R06_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R07_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R08_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R10_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R11_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R12_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R13_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R14_CpG.txt", package = "methylKit"),
	system.file("extdata", "RE-R15_CpG.txt", package = "methylKit"),
 )

### Read the files as a methylRawList object called "myobj" with minimum coverage requriment (at least 10 reads per base).
Treatment '0' for individuals in poor environments; Treatment '1' for individuals in enriched environments:

> myobj = methRead(file.list2,
    sample.id=list("PE-R03","PE-R04","PE-R08","PE-R09","PE-R12","PE-R14","RE-R05","RE-R06","RE-R07","RE-R08","RE-R10","RE-R11","RE-R12","RE-R13","RE-R14","RE-R15")),
    treatment=c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
    context="CpG",
    mincov = 10)

### Merging all samples to have coverage information across all individuals:

> united_all = unite(myobj, destrand=FALSE)

### Finding differentially methylated bases within 20% difference:

> myDiff20p = getMethylDiff(myDiff,difference=20,qvalue=0.01)

### Writing results

> write(myDiff20p,file="results_DMCs_parents.txt",append=TRUE)
```

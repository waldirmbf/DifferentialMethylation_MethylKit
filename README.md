# **Differential methylation analysis using the R package methylKit**

 The following documentation describes the R script used to estimate differentially-methylated cytosines between genetically-identical *K. marmoratus* individuals reared in enriched and poor environments with the R package methylKit. Methods are descibed in detail in [Berbel-Filho et al. (2020)](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15481)


Setting the local directory with all BAM files:
```
> setwd("/data/home/waldir/Desktop/2nd_generation/Analysis/Permutation/BAM")
> library('methylKit')
> packageVersion('methylKit')
```
Creating a list of all BAM files to get raw methylation data:
```
> file.list1 <- (list(system.file("extdata", "PE-R03_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
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
    system.file("extdata", "2nd_PE_R02_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_PE_R03_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_PE_R04_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_PE_R05_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_PE_R07_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_PE_R15_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_RE_R02_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_RE_R03_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_RE_R04_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd_RE_R11_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd-RE_R12_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd-RE_R13_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd-RE_R14_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
    system.file("extdata", "2nd-RE_R15_trimmed_bismarkbt2.bam.sorted", package = "methylKit"),
 )
```

Extracting raw methylation call files (proportions of Cs and Ts per base) from .bam files aligned with Bismark. For quantitation, a minimum coverage of 10 reads across all individuals (parental and offspring) was estabilished.
```
> objs <- processBismarkAln(location=file.list1,
                            sample.id=list("PE-R03","PE-R04","PE-R08","PE-R09","PE-R12","PE-R14","RE-R05",
                            "RE-R06","RE-R07","RE-R08","RE-R10","RE-R11","RE-R12","RE-R13","RE-R14","RE-R15",
                            "2nd_PE_R02","2nd_PE_R03","2nd_PE_R04","2nd_PE_R05","2nd_PE_R07","2nd_PE_R15",
                            "2nd_RE_R02","2nd_RE_R03","2nd_RE_R04","2nd_RE_R11","2nd_RE_R12","2nd_RE_R13","2nd_RE_R14","2nd_RE_R15"),
                            save.folder=NULL,
                            save.context=NULL,
                            read.context="CpG",
                            nolap=FALSE,
                            mincov=10,
                            minqual=20,
                            phred64=FALSE)
```

Inputting methylation call files into methylKit:
```
> file.list2 <- list(system.file("extdata", "PE-R03_CpG.txt", package = "methylKit"),
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
  system.file("extdata", "2nd_PE_R02_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_PE_R03_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_PE_R04_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_PE_R05_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_PE_R07_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_PE_R15_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R02_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R03_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R04_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R11_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R12_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R13_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R14_CpG.txt", package = "methylKit"),
  system.file("extdata", "2nd_RE_R15_CpG.txt", package = "methylKit"),
 )
```
Reading the files as a methylRawList object called "myobj" with minimum coverage requirement (at least 10 reads per base). We were only interested  to check DMCs  between parents,so there was no need here to dinstiguish the offspring treatments. Treatment '0' for individuals in poor environments; Treatment '1' for individuals in enriched environments.  We set Treatment '2' for all offspring individuals.
```
> myobj <- methRead(file.list2,
                    sample.id=list("PE-R03","PE-R04","PE-R08","PE-R09","PE-R12","PE-R14",
                    "RE-R05","RE-R06","RE-R07","RE-R08","RE-R10","RE-R11","RE-R12","RE-R13",
                    "RE-R14","RE-R15","2nd_PE_R02","2nd_PE_R03","2nd_PE_R04","2nd_PE_R05",
                    "2nd_PE_R07","2nd_PE_R15","2nd_RE_R02","2nd_RE_R03","2nd_RE_R04",
                    "2nd_RE_R11","2nd_RE_R12","2nd_RE_R13","2nd_RE_R14","2nd_RE_R15"),
                    assembly="ASM164957v1",
                    dbtype = NA,
                    pipeline = "bismark",  
                    header = TRUE,
                    skip = 0,
                    sep = "\t",
                    context = "CpG",  
                    resolution = "base",
                    treatment=c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                    mincov = 10)
```
Merging all samples to have coverage information across all individuals:
```
> united_all <- unite(myobj, destrand=FALSE)
```
Now, the object ``` united_all``` contains coverage information for all cytosines with > 10 reads common to all individuals. As we only wanted to check for epigenetic differences in the parents, we reorganised the data to create a dataset only including parental individuals.Reorganising dataset for the parental generation only:
```
> reorganized_all <- reorganize(united_all,
                            sample.id=c("PE-R03","PE-R04","PE-R08","PE-R09","PE-R12","PE-R14","RE-R05",
                            "RE-R06","RE-R07","RE-R08","RE-R10","RE-R11","RE-R12","RE-R13","RE-R14","RE-R15"),
                            treatment=c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1))
```
Getting differentially cytosines between treatments:
```
> myDiff <- calculateDiffMeth(reorganized_all)
```
Finding differentially methylated bases within 20% difference:
```
> myDiff20p <- getMethylDiff(myDiff,difference=20,qvalue=0.01)
```
Writing results:
```
> write(myDiff20p,file="results_DMCs_parents.txt",append=TRUE)
```
We found __1854 DMCs__ between parents living in poor or enriched environments. We need to test whether the number of DMCs is higher or lower than expected by chance.

Creating list wtih 4000 permutations of the original treatment:
```
> treatment <- c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)
  mytreatments <- list()
  for(i in 1:4000){
	   mytreatments[[i]]<-sample(treatment, 16, replace = FALSE)
  }
```
Defining the function ```calculateDMCs```, which receives each treatment permutation and recalculates the number of DMCs:
```
> calculateDMCs <- function(x) {
	reorganized_all <- reorganize(united_all,
                              sample.ids=c("PE-R03","PE-R04","PE-R08","PE-R09",
                              "PE-R12","PE-R14","RE-R05","RE-R06","RE-R07","RE-R08",
                              "RE-R10","RE-R11","RE-R12","RE-R13","RE-R14","RE-R15"),
                              treatment=x)
	myDiffParent <- calculateDiffMeth(reorganized_all)
	myDiffParent_20 <- getMethylDiff(myDiffParent,
                                    difference=20,
                                    qvalue=0.01)
	write(nrow(myDiffParent_20),file="results_PERMUTATIONS.txt",append=TRUE)
}
```
Calling the function ```calculateDMCs``` for each permutation.  The resulting number of DMCs for each treatment permutation is written into the file __results_PERMUTATIONS.txt.__
```
> lapply(mytreatments,calculateDMCs)
```


### Citation ###

Berbel‐Filho, Waldir M., et al. "Environmental enrichment induces intergenerational behavioural and epigenetic effects on fish." Molecular Ecology. https://doi.org/10.1111/mec.15481

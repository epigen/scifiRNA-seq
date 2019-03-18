sciRNA-seq/scifi-RNA-seq
===================

Pipeline and analysis scripts for combinatorial indexing single-cell RNA-seq (sciRNA-seq) and scifi-RNA-seq.

Check the [Makefile](Makefile) and [source files in src](src/) for more.

Due to the volatility of development of the project in the lab, each experiment can differ wildely from the previous.
I decided not to keep all the code specific to each version/run, but simply to git tag each version and modify the code to the latest version. A brief description of the various versions is below:

|experiment|description|outcome|
|-|-|-|
|**SCI004**|First test of sciRNA-seq, but run in bulk mode|Success|
|**BSF_0477_HJ7J7BGX2**|First single-cell sciRNA-seq. A species-mixing experiment|Complete failure|
|**SCI007/010/011**|sciRNA-seq optimizations: nuclei extraction/preparation/sorting|Mildly improved|
|**SCI012-013**|sciRNA-seq optimizations: various enzymes, cycles, cells vs nuclei, anchored vs unanchored primers|Greatly improved|
|**SCI016**|first scifiRNA-seq|Failure, but promising|
|**SCI017**|scifiRNA-seq|Success, but somewhat innefficient|
|**SCI017**-reanalysis|scifiRNA-seq|Success, found out major bottlenecks, understand system better|

the experiment IDs correspond to the [metadata annotations](metadata/annotation.csv)

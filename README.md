sciRNA-seq/scifi-RNA-seq
===================

Pipeline and analysis scripts for combinatorial indexing single-cell RNA-seq (sciRNA-seq) and scifi-RNA-seq.

Check the [Makefile](Makefile) and [source files in src](src/) for more.

Due to the volatility of development of the project in the lab, each experiment can differ wildely from the previous.
I decided not to keep all the code specific to each version/run, but simply to git tag each version and modify the code to the latest version. A brief description of the various versions is below:

|experiment|description|outcome|
|-|-|-|
|**PD187/188**|scifiRNA-seq - new version|Amaaaazing!|
|**SCI023/24**|scifiRNA-seq - T cells from 4 donors|Some problems remain|
|**SCI022**|scifiRNA-seq - Mixture cells overtook primary cells; first 384 well experiment; learned a few things about going into primary cells; barcodes need further inspection|Failure|
|**SCI021**|scifiRNA-seq| Large scale 125K mixture experiment, hitting native 10X limitations, scifi still with a lot of space to grow, barcodes need inspection though|Success|
|**SCI020-resequencing**|scifiRNA-seq, publication 4K mixture experiment"|Success|
|**SCI020**|scifiRNA-seq| found out optimal tagmentation conditions"|Success|
|**SCI019-resequencing**|scifiRNA-seq|Success|
|**SCI019**|scifiRNA-seq| found out 10M reads are enough for rough estimate of experiment pass/fail"|Success|
|**SCI017-reanalysis**|scifiRNA-seq| found out major bottlenecks| understand system better"|Success|
|**SCI017**|scifiRNA-seq| but somewhat innefficient"|Success|
|**SCI016**|first scifiRNA-seq|Failure, but promising|
|**SCI012-013**|sciRNA-seq optimizations: various enzymes, cycles, cells vs nuclei, anchored vs unanchored primers|Greatly improved|
|**SCI007/010/011**|sciRNA-seq optimizations: nuclei extraction/preparation/sorting|Mildly improved|
|**BSF_0477_HJ7J7BGX2**|First single-cell sciRNA-seq. A species-mixing experiment|Complete failure|
|**SCI004**|First test of sciRNA-seq, but run in bulk mode|Success|

the experiment IDs correspond to the [metadata annotations](metadata/annotation.csv)

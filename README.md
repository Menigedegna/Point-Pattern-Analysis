# Point pattern analysis

The jupyter scripts in these repository are necessary to process and run data analysis on files produced by Imaris plugins that can be downloaded from the Image processing repository. Although the pattern anlysis was designed for RNA PolII foci distrubution, in theory same script could be adapted to any foci distribution in nuclei, as long as plugins in image processing repository were used to process the 3D nuclei images.
![RNA PolII Data Processing goal scheme](DistanceScheme.png?=250x)



## Setting up the  environment

Install <b>jupyter</b> and set up <b>R kernel</b> for RNA PolII Data processing and <b>python kernel</b> for Machine_Learning_And_Pattern_Analysis_of_RNA_PolII_Distribution. 
### Prerequisites

R 3.4 and higher versions and the following packages:
```
matrixStats
spam
data.table
ggplot2
reshape2
```

Python 2.7 and higher versions and the following packages:
```
pandas
numpy
matplotlib
sklearn
seaborn
scipy 
```

## Built With

* Python 2.7
* R 3.4

## Acknowledgments

* Prof. Reinhard Furrer and Dr. Peter Majer for their recommendations
* Users of Stack Overflow for their questions and usefull solutions.

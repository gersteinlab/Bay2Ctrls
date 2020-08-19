## Bay2Ctrls: bayesian model to integrate two types of ChIP-seq controls for binding site detection

This R package is based on spp-1.14. And the wrapper script (run_Bay2Ctrls.R) for the R package is adapted from the ENCODE ChIP-seq pipeline. 

### INSTALLATION
`R CMD INSTALL Bay2Ctrls-0.1.tar.gz`

### DATA INPUT
This software require one IP (i) and two set of controls:
The first control is the DNA input (di) generated for the IP experiment. This control (di) and the IP (i) are the results from conventional ChIP-seq.
The second control contains a mock IP (m) and its corresponding DNA input (dm)

### DATA OUTPUT
Binding peaks are in the narrow peak file format

### USAGE
`R run_Bay2Ctrls.R` 

### Paramters specific for Bay2Ctrls
<pre>
-ip            a *.tagAlign.gz file containing mapped reads from the IP experiment (i)  
-mock          a *.tagAlign.gz file containing mapped reads from the mock IP experiment (m)  
-input4ip      a *.tagAlign.gz file containing mapped reads from the DNA input control for the IP experiment (di)  
-input4mock    a *.tagAlign.gz file containing mapped reads from the DNA input control for the mock IP experiment (dm)  
-totReads      scaling the numbers of reads of the four experiments to the same level, e.g. -totReads=10000000
-mcstep        simulation steps for the Bayesian model, e.g. -mcstep=1000000
</pre>

### Intrinsic parameters of SPP
<pre>
-npeak         the number of total narrow peaks to be output, e.g. -npeak=30000  
-x             avoid phantom peaks, e.g. -x=-500:85  
-s             genomic distance for signal correlation, e.g. -s=0:5:1200  
-out           output log file  
-rf            remove previous results  
-odir          output directory  
-p             number of processors to be used
</pre>
More options for output format please refer to the spp R package


**An example of the command:**  
`Rscript ./run_Bay2Ctrls.R -ip=test_data/a.rep0.tagAlign.gz -mock=test_data/EMb1.rep0.tagAlign.gz -input4ip=test_data/ap.rep0.tagAlign.gz -input4mock=test_data/EMb1p.rep0.tagAlign.gz -npeak=30000 -x=-500:85 -s=0:5:1200 -odir=./ -p=1 -filtchr='.*_[CD].*' -savr -savp -rf -out=test.cc -npeak=30000 -totReads=10000000 -mcstep=1000000`

### Data download
Download the test data [here](http://archive2.gersteinlab.org/proj/MockOrNot/Bay2Ctrls/test_data/)  
All raw and processed data are available [here](http://archive2.gersteinlab.org/proj/MockOrNot/Data/)

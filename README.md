## Bay2Ctrls: Bayesian model to integrate two types of ChIP-seq controls for binding site detection
Two types of controls were widely used for ChIP-seq: DNA input controls and mock IP controls. This software was developed to combine both for better binding site detection.

This R package is based on SPP, and the wrapper script (run_Bay2Ctrls.R) for the R package is adapted from the ENCODE ChIP-seq pipeline. 

### DATA INPUT
The software accecpt four files as input:  
The first set is the DNA input controlled IP experiment. This DNA control (di) and the IP (i) are the results from conventional ChIP-seq.  
The second set is the mock IP experiment. This includes the mock IP (m) and its corresponding DNA input (dm).

### DATA OUTPUT
Binding peaks are in the narrow peak file format

### DATA DOWNLOAD
Download the test data [here](http://archive2.gersteinlab.org/proj/MockOrNot/Bay2Ctrls/test_data/)  
All raw and processed data are available [here](http://archive2.gersteinlab.org/proj/MockOrNot/Data/)

### DEPENDENCE
- [SPP](https://cran.r-project.org/web/packages/spp/index.html)
- [caTools](https://cran.r-project.org/web/packages/caTools/index.html)

Tested with
- R 3.6.1
- SPP 1.16.0
- caTools 1.18.0

Please note that you might encounter a lot of issues if you want to install SPP from source. We recommend you install it via R `install.package("spp")`.

### INSTALLATION
`devtools::install_github("gersteinlab/Bay2Ctrls")`

### USAGE
Run the wrapper script **run_Bay2Ctrls.R** as the following:  
`Rscript run_Bay2Ctrls.R -ip=<IP_file> -mock=<mockIP_control_file> -input4ip=<DNA_input_control_file> -input4mock=<DNA_input_control_for_mockIP_file> [Bay2Ctrls_parameters, ...] [SPP_parameters, ...]`

##### Example  
`Rscript run_Bay2Ctrls.R -ip=test_data/a.rep0.tagAlign.gz -mock=test_data/EMb1.rep0.tagAlign.gz -input4ip=test_data/ap.rep0.tagAlign.gz -input4mock=test_data/EMb1p.rep0.tagAlign.gz -npeak=30000 -x=-500:85 -s=0:5:1200 -odir=./ -filtchr='.*_[CD].*' -savr -savp -rf -out=test.cc -npeak=30000 -totReads=10000000 -mcstep=1000000`

You might also use the R package Bay2Ctrls in your scripts. See man pages for individual functions(TODO).  

##### Paramters specific for Bay2Ctrls
<pre>
-ip            a *.tagAlign.gz file containing mapped reads from the IP experiment (i)  
-mock          a *.tagAlign.gz file containing mapped reads from the mock IP experiment (m)  
-input4ip      a *.tagAlign.gz file containing mapped reads from the DNA input control for the IP experiment (di)  
-input4mock    a *.tagAlign.gz file containing mapped reads from the DNA input control for the mock IP experiment (dm)  
-totReads      scaling the numbers of reads of the four experiments to the same level, e.g. -totReads=10000000
-mcstep        simulation steps for the Bayesian model, e.g. -mcstep=1000000
</pre>

##### Intrinsic parameters of SPP
<pre>
-npeak         the number of total narrow peaks to be output, e.g. -npeak=30000  
-x             avoid phantom peaks, e.g. -x=-500:85  
-s             genomic distance for signal correlation, e.g. -s=0:5:1200  
-out           output log file  
-rf            remove previous results  
-odir          output directory  
-p             number of processors to be used
</pre>
More options for output format please refer to the spp R package [github](https://github.com/hms-dbmi/spp) or [CRAN](https://cran.r-project.org/web/packages/spp/index.html)

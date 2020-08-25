library(spp)
library(Bay2Ctrls)

args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

get.file.parts <- function(file.fullpath) {
    if (! is.character(file.fullpath)) {
        stop('File name must be a string')
    }
    
    file.parts <- strsplit(as.character(file.fullpath), .Platform$file.sep, fixed=TRUE)[[1]] # split on file separator
    
    if (length(file.parts) == 0) { # if empty file name
        return(list(path='',
            fullname='',
            name='',
            ext='')
            )
    } else {
        if (length(file.parts) == 1) { # if no path then just the file name itself
            file.path <- '.'
            file.fullname <- file.parts
        } else {
            file.path <- paste(file.parts[1:(length(file.parts)-1)], collapse=.Platform$file.sep) # 1:last-1 token is path
            file.fullname <- file.parts[length(file.parts)] # last token is filename
        }
        file.fullname.parts <- strsplit(file.fullname,'.',fixed=TRUE)[[1]] # split on .
        if (length(file.fullname.parts) == 1) { # if no extension
            file.ext <- ''
            file.name <- file.fullname.parts
        } else {
            file.ext <- paste('.', file.fullname.parts[length(file.fullname.parts)], sep="") # add the . to the last token
            file.name <- paste(file.fullname.parts[1:(length(file.fullname.parts)-1)], collapse=".")
        }
        return(list(path=file.path,
            fullname=file.fullname,
            name=file.name,
            ext=file.ext))
    }
} # end: get.file.parts()


parse.arguments <- function(args) {
    # Set arguments to default values
    ip.file <- NA  # main ChIP tagAlign/BAM file name
    isurl.ip.file <- FALSE # flag indicating whether ChIP file is a URL
    input4ip.file <- NA # control tagAlign/BAM file name
    isurl.input4ip.file <- FALSE # flag indicating whether control file is a URL   

    mock.file <- NA
    isurl.mock.file <- FALSE 
    input4mock.file <- NA 
    isurl.input4mock.file <- FALSE 

    totReads <- 1e+7
    wdfold <- 0
    mcstep <- 1e+6
    flagShift <- 0
    rgn <- "mock"
    weightspp <- 1
    weightmc <- 1
    rmAbnormal <- 1 # 0 no remove, 1 w/ remove
    mthd <- 'mcbin'
    rankby <- 'spp'
    # mcbin parameters: spp, pvaap, pvab, pvmc, 
    # poisson parameters: signal, pvalue

    sep.min <- -100  # min strand shift
    sep.max <- 600  # max strand shift
    sep.bin <- 5    # increment for strand shift
    sep.peak <- NA # user-defined peak shift
    exclude.min <- 10 # lowerbound of strand shift exclusion region
    exclude.max <- NaN # upperbound of strand shift exclusion region
    n.nodes <- NA # number of parallel processing nodes
    fdr <- 0.01 # false discovery rate threshold for peak calling
    npeak <- NA # threshold on number of peaks to call
    temp.dir <- tempdir() # temporary directory
    chrname.rm.pattern <- NA # chromosome name pattern used to remove tags
    output.odir <- NA # Output directory name
    output.npeak.file <- NA # Output narrowPeak file name
    output.rpeak.file <- NA # Output regionPeak file name
    output.rdata.file <- NA # Rdata file
    output.plot.file <- NA  # cross correlation plot file
    output.result.file <- NA # result file
    replace.flag <- FALSE # replace file flag
    clean.files.flag <- FALSE # file deletion flag
    
    # Parse arguments   
    for (each.arg in args) {
        if (grepl('^-ip=',each.arg)) { #-c=<chip.file>  IP data, a
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                ip.file <- arg.split[2] # second part is chip.file
            } else {
                stop('No tagAlign/BAM file name provided for parameter -ip=')
            }
            
        } else if (grepl('^-input4ip=',each.arg)) { #-i=<control.file>, control of IP, ap
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                input4ip.file <- arg.split[2] # second part is control.file
        }

        } else if (grepl('^-mock=',each.arg)) { # mock IP data, b
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                mock.file <- arg.split[2] # second part is chip.file
            }

        } else if (grepl('^-input4mock=',each.arg)) { # control of mock IP
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                input4mock.file <- arg.split[2] # second part is control.file
            }

        } else if (grepl('^-mcstep=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                mcstep <- as.numeric(arg.split[2]) # second part is chip.file
            }
        } else if (grepl('^-mthd=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                mthd <- as.character(arg.split[2]) # second part is chip.file
            }

        } else if (grepl('^-totReads=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                totReads <- as.numeric(arg.split[2]) # second part is chip.file
            } else {
                   totReads <- 1e+7
            }

        } else if (grepl('^-rankby=',each.arg)) { #-rankby= for poisson "pvalue or ratio"; for 
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                   rankby <- as.character(arg.split[2]) # second part is chip.file
            } else {
                   rankby <- "spp"
            }

        } else if (grepl('^-wdfold=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                wdfold <- as.numeric(arg.split[2]) # second part is chip.file
            } else {
              wdfold <- 0
            }
        } else if (grepl('^-flagShift=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                flagShift <- as.numeric(arg.split[2]) # second part is chip.file
            } else {
              flagShift <- 0
            }
        } else if (grepl('^-rgn=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                   rgn <- as.character(arg.split[2]) # second part is chip.file
            } else {
                   rgn <- "mock"
            }

        } else if (grepl('^-weightmc=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                   weightmc <- as.numeric(arg.split[2]) # second part is chip.file
            } else {
                   weightmc <- 1
            }
        } else if (grepl('^-weightspp=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                   weightspp <- as.numeric(arg.split[2]) # second part is chip.file
            } else {
                   weightspp <- 1
            }
        } else if (grepl('^-rmAbnormal=',each.arg)) { #-c=<chip.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                   rmAbnormal <- as.numeric(arg.split[2]) # second part is chip.file
            } else {
                   rmAbnormal <- 1 # 0 do not remove; 1 remove
            }
        } else if (grepl('^-s=',each.arg)) { #-s=<sep.min>:<sep.bin>:<sep.max>          
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                sep.vals <- arg.split[2] # second part is sepmin:sepbin:sepmax
                sep.vals.split <- strsplit(sep.vals,':',fixed=TRUE)[[1]] # split on :                
                if (length(sep.vals.split) != 3) { # must have 3 parts
                    stop('Strand shift limits must be specified as -s=sepmin:sepbin:sepmax')                    
                } else {
                    if (any(is.na(as.numeric(sep.vals.split)))) { # check that sep vals are numeric
                        stop('Strand shift limits must be numeric values')
                    }
                    sep.min <- round(as.numeric(sep.vals.split[1]))
                    sep.bin <- round(as.numeric(sep.vals.split[2]))
                    sep.max <- round(as.numeric(sep.vals.split[3]))
                    if ((sep.min > sep.max) || (sep.bin > (sep.max - sep.min)) || (sep.bin < 0)) {
                        stop('Illegal separation values -s=sepmin:sepbin:sepmax')
                    }
                }                                    
            } else {
                stop('Strand shift limits must be specified as -s=sepmin:sepbin:sepmax')
            }
            
        } else if (grepl('^-speak=',each.arg)) { #-speak=<sep.peak> , user-defined cross-correlation peak strandshift
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                sep.peak <- arg.split[2] # second part is <sep.peak>
                if (is.na(as.numeric(sep.peak))) { # check that sep.peak is numeric
                    stop('-speak=<sep.peak>: User defined peak shift must be numeric')
                }
                sep.peak <- as.numeric(sep.peak)
            } else {
                stop('User defined peak shift must be provided as -speak=<sep.peak>')
            }
            
        } else if (grepl('^-x=',each.arg)) { #-x=<exclude.min>:<exclude.max>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                exclude.vals <- arg.split[2] # second part is excludemin:excludemax
                exclude.vals.split <- strsplit(exclude.vals,':',fixed=TRUE)[[1]] # split on :                
                if (length(exclude.vals.split) != 2) { # must have 2 parts
                    stop('Exclusion limits must be specified as -x=excludemin:excludemax')
                } else {
                    if (any(is.na(as.numeric(exclude.vals.split)))) { # check that exclude vals are numeric
                        stop('Exclusion limits must be numeric values')
                    }
                    exclude.min <- round(as.numeric(exclude.vals.split[1]))                    
                    exclude.max <- round(as.numeric(exclude.vals.split[2]))
                    if (exclude.min > exclude.max) {
                        stop('Illegal exclusion limits -x=excludemin:excludemax')
                    }                    
                }                                    
            } else {
                stop('Exclusion limits must be specified as -x=excludemin:excludemax')
            }
            
        } else if (grepl('^-p=',each.arg)) { #-p=<n.nodes> , number of parallel processing nodes, default=NULL            
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                n.nodes <- arg.split[2] # second part is numnodes
                if (is.na(as.numeric(n.nodes))) { # check that n.nodes is numeric
                    stop('-p=<numnodes>: numnodes must be numeric')
                }
                n.nodes <- round(as.numeric(n.nodes))
            } else {
                stop('Number of parallel nodes must be provided as -p=<numnodes>')
            }
            
        } else if (grepl('^-fdr=',each.arg)) { #-fdr=<fdr> , false discovery rate, default=0.01            
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                fdr <- arg.split[2] # second part is fdr
                if (is.na(as.numeric(fdr))) { # check that fdr is numeric
                    stop('-fdr=<falseDiscoveryRate>: false discovery rate must be numeric')
                }
                fdr <- as.numeric(fdr)
            } else {
                stop('False discovery rate must be provided as -fdr=<fdr>')
            }
            
        } else if (grepl('^-npeak=',each.arg)) { #-npeak=<numPeaks> , number of peaks threshold, default=NA            
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                npeak <- arg.split[2] # second part is npeak
                if (is.na(as.numeric(npeak))) { # check that npeak is numeric
                    stop('-npeak=<numPeaks>: threshold on number of peaks must be numeric')
                }
                npeak <- round(as.numeric(npeak))
            } else {
                stop('Threshold on number of peaks must be provided as -npeak=<numPeaks>')
            }
            
        } else if (grepl('^-tmpdir=',each.arg)) { #-tmpdir=<temp.dir>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                temp.dir <- arg.split[2] # second part is temp.dir
            } else {
                stop('No temporary directory provided for parameter -tmpdir=')
            }

        } else if (grepl('^-filtchr=',each.arg)) { #-filtchr=<chrname.rm.pattern>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                chrname.rm.pattern <- arg.split[2] # second part is chrname.rm.pattern
            } else {
                stop('No pattern provided for parameter -filtchr=')
            }
            
        } else if (grepl('^-odir=',each.arg)) { #-odir=<output.odir>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                output.odir <- arg.split[2] # second part is output.odir
            } else {
                stop('No output directory provided for parameter -odir=')
            }
            
        } else if (grepl('^-savn',each.arg)) { # -savn=<output.npeak.file> OR -savn , save narrowpeak
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2])) {                
                output.npeak.file <- arg.split[2] #-savn=
            } else if (each.arg=='-savn') {
                output.npeak.file <- NULL # NULL indicates get the name from the main file name
            } else {
                stop('Argument for saving narrowPeak file must be -savn or -savn=<filename>')
            }
            
        } else if (grepl('^-savr',each.arg)) { # -savr=<output.rpeak.file> OR -savr , save regionpeak
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2])) {                
                output.rpeak.file <- arg.split[2] #-savr=
            } else if (each.arg=='-savr') {
                output.rpeak.file <- NULL # NULL indicates get the name from the main file name
            } else {
                stop('Argument for saving regionPeak file must be -savr or -savr=<filename>')
            }
            
        } else if (grepl('^-savd',each.arg)) { # -savd=<output.rdata.file> OR -savd , save Rdata file
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2])) {                
                output.rdata.file <- arg.split[2] #-savd=
            } else if (each.arg=='-savd') {
                output.rdata.file <- NULL # NULL indicates get the name from the main file name
            } else {
                stop('Argument for saving Rdata file must be -savd or -savd=<filename>')
            }
            
        } else if (grepl('^-savp',each.arg)) { # -savp=<output.plot.file> OR -savp , save cross-correlation plot                       
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2])) {                
                output.plot.file <- arg.split[2] #-savp=
            } else if (each.arg=='-savp') {
                output.plot.file <- NULL # NULL indicates get the name from the main file name
            } else {
                stop('Argument for saving Rdata file must be -savp or -savp=<filename>')
            }
            
        } else if (grepl('^-out=',each.arg)) { #-out=<output.result.file>
            
            arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
            if (! is.na(arg.split[2]) ) {
                output.result.file <- arg.split[2] # second part is output.result.file
            } else {
                stop('No result file provided for parameter -out=')
            }
        
        } else if (each.arg == '-rf') {
            
            replace.flag <- TRUE
        
        } else if (each.arg == '-clean') {
            
            clean.files.flag <- TRUE
            
        } else {
            
            stop('Illegal argument ',each.arg)
        }        
    }
    # End: for loop
    
    # Check arg combination
    if (wdfold != 0 & flagShift == 0) {
           stop('wdfold anc flagShift problem')
    }

    # Check mandatory arguments
    if (is.na(ip.file)) {
        stop('-ip=<tagAlign/BAMFileName> is a mandatory argument')
    }
    
    # Correct other arguments
    if (is.na(output.odir)) { # Reconstruct output.odir if not provided
        output.odir <- get.file.parts(ip.file)$path
    }

    chip.file <- ip.file
    if (! is.na(mock.file)) {
        control.file <- mock.file
    }else if(! is.na(input4ip.file)) {
        control.file <- input4ip.file
    }
    
    if (is.null(output.npeak.file)) { # Reconstruct output.npeak.file if NULL
        output.npeak.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '_VS_', get.file.parts(control.file)$name,'.narrowPeak', sep=""))
    }

    if (is.null(output.rpeak.file)) { # Reconstruct output.rpeak.file if NULL
        output.rpeak.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '_VS_', get.file.parts(control.file)$name,'.regionPeak', sep=""))
    }
    
    if (is.null(output.rdata.file)) { # Reconstruct output.rdata.file if NULL
        output.rdata.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '.Rdata', sep=""))
    }
    
    if (is.null(output.plot.file)) { # Reconstruct output.plot.file if NULL
        output.plot.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '.pdf', sep=""))
    }
    
    return(list(ip.file=ip.file,
                    input4ip.file=input4ip.file,
                    mock.file=mock.file,
                    input4mock.file=input4mock.file,

                    totReads=totReads,
                    wdfold=wdfold,
                    mcstep=mcstep,
                    flagShift=flagShift,
                    rgn=rgn,
                    rankby=rankby,
                    weightspp=weightspp,
                    weightmc=weightmc,
                    rmAbnormal=rmAbnormal,
                    mthd=mthd,

                    sep.range=c(sep.min,sep.bin,sep.max),
                    sep.peak=sep.peak,
                    ex.range=c(exclude.min,exclude.max),
                    n.nodes=n.nodes,
                    fdr=fdr,
                    npeak=npeak,
                    temp.dir=temp.dir,
                    chrname.rm.pattern=chrname.rm.pattern,
                    output.odir=output.odir,
                    output.npeak.file=output.npeak.file,
                    output.rpeak.file=output.rpeak.file,
                    output.rdata.file=output.rdata.file,
                    output.plot.file=output.plot.file,
                    output.result.file=output.result.file,
                    replace.flag=replace.flag,
                    clean.files.flag=clean.files.flag))          
} # end: parse.arguments()

read.align <- function(align.filename) {
# ===================================
# Function will read a tagAlign or BAM file
# ===================================   
    if (grepl('(\\.bam)?.*(\\.tagAlign)',align.filename)) { # if tagalign file
        chip.data <- read.tagalign.tags(align.filename)
        # get readlength info
        tmpDataRows <- read.table(align.filename,nrows=500)
        chip.data$read.length <- round(median(tmpDataRows$V3 - tmpDataRows$V2))
    } else if (grepl('(\\.tagAlign)?.*(\\.bam)',align.filename)) { # if bam file
        # create BAM file name
        bam2align.filename <- sub('\\.bam','.tagAlign',align.filename)
        # generate command to convert bam to tagalign
        command <- vector(length=2)
        command[1] <- sprintf("samtools view -F 0x0204 -o - %s",align.filename)
        command[2] <- paste("awk 'BEGIN{FS=" , '"\t"' , ";OFS=", '"\t"} {if ($2==16) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"}}', "' 1> ", bam2align.filename, sep="")
        # command[2] <- paste("awk 'BEGIN{OFS=", '"\t"} {if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"}}', "' 1> ", bam2align.filename, sep="")   
        command <- paste(command,collapse=" | ")
        cat(command,"\n")
        # Run command
        status <- system(command,intern=FALSE,ignore.stderr=FALSE)
        if ((status != 0) || !file.exists(bam2align.filename)) {
            cat(sprintf("Error converting BAM to tagalign file: %s\n",align.filename),file=stderr())
            q(save="no",status=1)
        }
        # read converted BAM file
        chip.data <- read.tagalign.tags(bam2align.filename)
        # get readlength info
        tmpDataRows <- read.table(bam2align.filename,nrows=500)
        chip.data$read.length <- round(median(tmpDataRows$V3 - tmpDataRows$V2))
        # delete temporary tagalign file
        file.remove(bam2align.filename) 
    } else {
        cat(sprintf("Error:Unknown file format for file:%s\n",align.fname),file=stderr())
        q(save="no",status=1)   
    }
    return(chip.data)
} # end: read.align()

print.run.params <- function(iparams){
# ===================================
# Output run parameters
# ===================================       
    cat('################\n',file=stdout())
    cat(iparams$ip.file,
            iparams$input4ip.file,
            iparams$sep.range,
            iparams$sep.peak,
            iparams$ex.range,
            iparams$n.nodes,
            iparams$fdr,
            iparams$npeak,
            iparams$output.odir,
            iparams$output.npeak.file,
            iparams$output.rpeak.file,
            iparams$output.rdata.file,
            iparams$output.plot.file,
            iparams$output.result.file,
            iparams$replace.flag,  
            labels=c(
                'ChIP data:',
                'Control data:', 
                'strandshift(min):',
                'strandshift(step):',
                'strandshift(max)',
                'user-defined peak shift',
                'exclusion(min):',
                'exclusion(max):',
                'num parallel nodes:',
                'FDR threshold:',
                'NumPeaks Threshold:',
                'Output Directory:',
                'narrowPeak output file name:',
                'regionPeak output file name:',
                'Rdata filename:', 
                'plot pdf filename:',
                'result filename:',
                'Overwrite files?:'),
            fill=18,
            file=stdout())  
    cat('\n',file=stdout()) 
} # end: print.run.parameters()


check.replace.flag <- function(iparams){
# ===================================
# Check if files exist
# ===================================
# If replace.flag is NOT set, check if output files exist and abort if necessary
    if (! iparams$replace.flag) {
        if (! is.na(iparams$output.npeak.file)) {
            if (file.exists(iparams$output.npeak.file)) {
                cat('narrowPeak file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
                q(save="no",status=1)
            }
        }
        if (! is.na(iparams$output.rpeak.file)) {
            if (file.exists(iparams$output.rpeak.file)) {
                cat('regionPeak file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
                q(save="no",status=1)
            }
        }    
        if (! is.na(iparams$output.plot.file)) {
            if (file.exists(iparams$output.plot.file)) {
                cat('Plot file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
                q(save="no",status=1)
            }
        }
        if (! is.na(iparams$output.rdata.file)) {
            if (file.exists(iparams$output.rdata.file)) {
                cat('Rdata file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
                q(save="no",status=1)
            }
        }
    }   
}

# #############################################################################
# MAIN FUNCTION
# #############################################################################

#print("here2\n");
# Check number of arguments
minargs = 1;
maxargs = 40;

iparams <- parse.arguments(args)

# Print run parameters
print.run.params(iparams)

# Check if output files exist 
#check.replace.flag(iparams)

# curr.chip.file and curr.control.file always point to the original ChIP and control files on disk
# ta.chip.filename & ta.control.filename always point to the final but temporary versions of the ChIP and control files that will be passed to read.align

# Download ChIP and control files if necessary to temp.dir
curr.ip.file <- iparams$ip.file
curr.input4ip.file <- iparams$input4ip.file
curr.mock.file <- iparams$mock.file
curr.input4mock.file <- iparams$input4mock.file

# unzip ChIP and input files if required AND copy to temp directory
if (get.file.parts(curr.ip.file)$ext == '.gz') {
    ta.ip.filename <- tempfile(get.file.parts(curr.ip.file)$name, tmpdir=iparams$temp.dir) # unzip file to temp.dir/[filename with .gz removed][randsuffix]
    cat('Decompressing ChIP file\n',file=stdout())
    if (system(paste("gunzip -c",curr.ip.file,">",ta.ip.filename)) != 0) {
        stop('Unable to decompress file:', iparams$ip.file)
    }
    if (iparams$clean.files.flag) { # Remove original file if clean.files.flag is set
        file.remove(curr.ip.file)
    }
} else {
    ta.ip.filename <- tempfile(get.file.parts(curr.ip.file)$fullname, tmpdir=iparams$temp.dir)
    if (iparams$clean.files.flag) {
        file.rename(curr.ip.file,ta.ip.filename) # move file to temp.dir/[filename][randsuffix]
    } else {
        file.copy(curr.ip.file,ta.ip.filename) # copy file to temp.dir/[filename][randsuffix]       
    }   
}

if (! is.na(iparams$mock.file)) {
    if (get.file.parts(curr.mock.file)$ext == '.gz') {
        ta.mock.filename <- tempfile(get.file.parts(curr.mock.file)$name, tmpdir=iparams$temp.dir) # unzip file to temp.dir/[filename with .gz removed][randsuffix]
        cat('Decompressing mockIP file\n',file=stdout())
        if (system(paste("gunzip -c",curr.mock.file,">",ta.mock.filename)) != 0) {
            stop('Unable to decompress file:', iparams$mock.file)
        }
        if (iparams$clean.files.flag) { # Remove original file if clean.files.flag is set
            file.remove(curr.mock.file)
        }
    } else {
        ta.mock.filename <- tempfile(get.file.parts(curr.mock.file)$fullname, tmpdir=iparams$temp.dir)
        if (iparams$clean.files.flag) {
            file.rename(curr.mock.file,ta.mock.filename) # move file to temp.dir/[filename][randsuffix]
        } else {
            file.copy(curr.mock.file,ta.mock.filename) # copy file to temp.dir/[filename][randsuffix]
        }
    }
}

if (! is.na(iparams$input4ip.file)) {
    if (get.file.parts(curr.input4ip.file)$ext == '.gz') {
        ta.input4ip.filename <- tempfile(get.file.parts(curr.input4ip.file)$name, tmpdir=iparams$temp.dir) # unzip file to temp.dir/[filename with .gz removed][randsuffix]
        cat('Decompressing control file\n',file=stdout())        
        if (system(paste("gunzip -c",curr.input4ip.file,">",ta.input4ip.filename)) != 0) {
            stop('Unable to decompress file:', iparams$input4ip.file)
        }
        if (iparams$clean.files.flag) { # Remove original file if clean.files.flag is set
            file.remove(curr.input4ip.file)
        }               
    } else {
        ta.input4ip.filename <- tempfile(get.file.parts(curr.input4ip.file)$fullname, tmpdir=iparams$temp.dir) # copy file to temp.dir/[filename][randsuffix]
        
        if (iparams$clean.files.flag) {
            file.rename(curr.input4ip.file,ta.input4ip.filename) # move file to temp.dir/[filename][randsuffix]
        } else {
            file.copy(curr.input4ip.file,ta.input4ip.filename) # copy file to temp.dir/[filename][randsuffix]       
        }           
    }
}

if (! is.na(iparams$input4mock.file)) {
    if (get.file.parts(curr.input4mock.file)$ext == '.gz') {
        ta.input4mock.filename <- tempfile(get.file.parts(curr.input4mock.file)$name, tmpdir=iparams$temp.dir) # unzip file to temp.dir/[filename with .gz removed][randsuffix]
        cat('Decompressing controlMockIP file\n',file=stdout())        
        if (system(paste("gunzip -c",curr.input4mock.file,">",ta.input4mock.filename)) != 0) {
            stop('Unable to decompress file:', iparams$input4mock.file)
        }
        if (iparams$clean.files.flag) { # Remove original file if clean.files.flag is set
            file.remove(curr.input4mock.file)
        }               
    } else {
        ta.input4mock.filename <- tempfile(get.file.parts(curr.input4mock.file)$fullname, tmpdir=iparams$temp.dir) # copy file to temp.dir/[filename][randsuffix]
        
        if (iparams$clean.files.flag) {
            file.rename(curr.input4mock.file,ta.input4mock.filename) # move file to temp.dir/[filename][randsuffix]
        } else {
            file.copy(curr.input4mock.file,ta.input4mock.filename) # copy file to temp.dir/[filename][randsuffix]       
        }           
    }
}

# Read ChIP tagAlign/BAM files
#cat("Reading ChIP tagAlign/BAM file",iparams$ip.file,"\n",file=stdout())
ip.data <- read.align(ta.ip.filename)
#cat("ChIP data read length",chip.data$read.length,"\n",file=stdout())
#cat("ChIP data read length",ip.data$read.length,"\n",file=stdout())
file.remove(ta.ip.filename) # Delete temporary file
if (length(ip.data$tags)==0) {
    stop('Error in ChIP file format:', iparams$ip.file)
}

# Remove illegal chromosome names
if (! is.na(iparams$chrname.rm.pattern)) {
    selectidx <- which(grepl(iparams$chrname.rm.pattern,names(ip.data$tags))==FALSE)
    ip.data$tags <- ip.data$tags[selectidx]
    ip.data$quality <- ip.data$quality[selectidx]
}
ip.data$num.tags <- sum(unlist(lapply(ip.data$tags,function(d) length(d))))

# Read Control tagAlign/BAM files
if (! is.na(iparams$input4ip.file)) {
    cat("Reading Control tagAlign/BAM file",iparams$input4ip.file,"\n",file=stdout())
    input4ip.data <- read.align(ta.input4ip.filename)
    file.remove(ta.input4ip.filename) # Delete temporary file    
    if (length(input4ip.data$tags)==0) {
        stop('Error in control file format:', iparams$input4ip.file)
    }    
    cat("Control data read length",input4ip.data$read.length,"\n",file=stdout())
    # Remove illegal chromosome names
    if (! is.na(iparams$chrname.rm.pattern)) {
        selectidx <- which(grepl(iparams$chrname.rm.pattern,names(input4ip.data$tags))==FALSE)
        input4ip.data$tags <- input4ip.data$tags[selectidx]
        input4ip.data$quality <- input4ip.data$quality[selectidx]
    }
    input4ip.data$num.tags <- sum(unlist(lapply(input4ip.data$tags,function(d) length(d))))
}


if (! is.na(iparams$mock.file)) {
    cat("Reading Control tagAlign/BAM file",iparams$mock.file,"\n",file=stdout())
    mock.data <- read.align(ta.mock.filename)
    file.remove(ta.mock.filename) # Delete temporary file    
    if (length(mock.data$tags)==0) {
        stop('Error in control file format:', iparams$mock.file)
    }    
    cat("mockIP data read length",mock.data$read.length,"\n",file=stdout())
    # Remove illegal chromosome names
    if (! is.na(iparams$chrname.rm.pattern)) {
        selectidx <- which(grepl(iparams$chrname.rm.pattern,names(mock.data$tags))==FALSE)
        mock.data$tags <- mock.data$tags[selectidx]
        mock.data$quality <- mock.data$quality[selectidx]
    }
    mock.data$num.tags <- sum(unlist(lapply(mock.data$tags,function(d) length(d))))
}

if (! is.na(iparams$input4mock.file)) {
    cat("Reading Control tagAlign/BAM file",iparams$input4mock.file,"\n",file=stdout())
    input4mock.data <- read.align(ta.input4mock.filename)
    file.remove(ta.input4mock.filename) # Delete temporary file    
    if (length(input4mock.data$tags)==0) {
        stop('Error in control file format:', iparams$input4mock.file)
    }    
    cat("Control of mockIP data read length",input4mock.data$read.length,"\n",file=stdout())
    # Remove illegal chromosome names
    if (! is.na(iparams$chrname.rm.pattern)) {
        selectidx <- which(grepl(iparams$chrname.rm.pattern,names(input4mock.data$tags))==FALSE)
        input4mock.data$tags <- input4mock.data$tags[selectidx]
        input4mock.data$quality <- input4mock.data$quality[selectidx]
    }
    input4mock.data$num.tags <- sum(unlist(lapply(input4mock.data$tags,function(d) length(d))))
}

# Open multiple processes if required
if (is.na(iparams$n.nodes)) {
    cluster.nodes <- NULL
} else {
    cat("loading snow\n")
    library(snow)
    cat("loaded snow\n")
    cluster.nodes <- makeCluster(iparams$n.nodes,type="SOCK")
    cat("multiple nodes ",cluster.nodes$rank,"\n")
}

# #################################    
# Calculate cross-correlation for various strand shifts
# #################################    
cat("Calculating peak characteristics\n",file=stdout())
# crosscorr
# $cross.correlation : Cross-correlation profile as an $x/$y data.frame
# $peak : Position ($x) and height ($y) of automatically detected cross-correlation peak.
# $whs: Optimized window half-size for binding detection (based on the width of the cross-correlation peak) 
crosscorr <- get.binding.characteristics(ip.data,
        srange=iparams$sep.range[c(1,3)],
        bin=iparams$sep.range[2],
        accept.all.tags=T,
        cluster=cluster.nodes)

if (!is.na(iparams$n.nodes)) {
    stopCluster(cluster.nodes)
}

# Smooth the cross-correlation curve if required
cc <- crosscorr$cross.correlation
crosscorr$min.cc <- crosscorr$cross.correlation[ which.min(crosscorr$cross.correlation$y) , ] # minimum value and shift of cross-correlation
cat("Minimum cross-correlation value", crosscorr$min.cc$y,"\n",file=stdout())
cat("Minimum cross-correlation shift", crosscorr$min.cc$x,"\n",file=stdout())
sbw <- 2*floor(ceiling(5/iparams$sep.range[2]) / 2) + 1 # smoothing bandwidth
cc$y <- caTools::runmean(cc$y,sbw,alg="fast")

# Compute cross-correlation peak
bw <- ceiling(2/iparams$sep.range[2]) # crosscorr[i] is compared to crosscorr[i+/-bw] to find peaks
peakidx <- (diff(cc$y,bw)>=0) # cc[i] > cc[i-bw]
peakidx <- diff(peakidx,bw)
peakidx <- which(peakidx==-1) + bw        

# exclude peaks from the excluded region
if ( is.nan(iparams$ex.range[2]) ) {
    iparams$ex.range[2] <- chip.data$read.length+10
}
peakidx <- peakidx[(cc$x[peakidx] < iparams$ex.range[1]) | (cc$x[peakidx] > iparams$ex.range[2])]    
cc <- cc[peakidx,]

# Find max peak position and other peaks within 0.9*max_peakvalue that are further away from maxpeakposition   
maxpeakidx <- which.max(cc$y)
maxpeakshift <- cc$x[maxpeakidx]
maxpeakval <- cc$y[maxpeakidx]
peakidx <-which((cc$y >= 0.9*maxpeakval) & (cc$x >= maxpeakshift)) 
cc <- cc[peakidx,]

# sort the peaks and get the top 3
sortidx <- order(cc$y,decreasing=TRUE)
sortidx <- sortidx[c(1:min(3,length(sortidx)))]
cc.peak <- cc[sortidx,]

# Override peak shift if user supplies peak shift
if (! is.na(iparams$sep.peak)) {
    cc.peak <- approx(crosscorr$cross.correlation$x,crosscorr$cross.correlation$y,iparams$sep.peak,rule=2)
}
cat("Peak cross-correlation value", paste(cc.peak$y,collapse=","),"\n",file=stdout())
cat("Peak strand shift",paste(cc.peak$x,collapse=","),"\n",file=stdout())

# Reset values in crosscorr
crosscorr$peak$x <- cc.peak$x[1]
crosscorr$peak$y <- cc.peak$y[1]

# Compute window half size
whs.thresh <- crosscorr$min.cc$y + (crosscorr$peak$y - crosscorr$min.cc$y)/3
crosscorr$whs <- max(crosscorr$cross.correlation$x[crosscorr$cross.correlation$y >= whs.thresh])
cat("Window half size",crosscorr$whs,"\n",file=stdout())

# Compute phantom peak coefficient
ph.peakidx <- which( ( crosscorr$cross.correlation$x >= ( ip.data$read.length - round(2*iparams$sep.range[2]) ) ) & 
                     ( crosscorr$cross.correlation$x <= ( ip.data$read.length + round(1.5*iparams$sep.range[2]) ) ) )
ph.peakidx <- ph.peakidx[ which.max(crosscorr$cross.correlation$y[ph.peakidx]) ]
crosscorr$phantom.cc <- crosscorr$cross.correlation[ph.peakidx,]
cat("Phantom peak location",crosscorr$phantom.cc$x,"\n",file=stdout())
cat("Phantom peak Correlation",crosscorr$phantom.cc$y,"\n",file=stdout())
crosscorr$phantom.coeff <- crosscorr$peak$y / crosscorr$phantom.cc$y
crosscorr$phantom.coeff <- crosscorr$peak$y / crosscorr$min.cc$y
cat("Normalized cross-correlation coefficient (NCCC)",crosscorr$phantom.coeff,"\n",file=stdout())
crosscorr$rel.phantom.coeff <- (crosscorr$peak$y - crosscorr$min.cc$y) / (crosscorr$phantom.cc$y - crosscorr$min.cc$y)
cat("Relative Cross correlation Coefficient (RCCC)",crosscorr$rel.phantom.coeff,"\n",file=stdout())
crosscorr$phantom.quality.tag <- NA
if ( (crosscorr$rel.phantom.coeff >= 0) & (crosscorr$rel.phantom.coeff < 0.25) ) {
    crosscorr$phantom.quality.tag <- -2
} else if ( (crosscorr$rel.phantom.coeff >= 0.25) & (crosscorr$rel.phantom.coeff < 0.5) ) {
    crosscorr$phantom.quality.tag <- -1
} else if ( (crosscorr$rel.phantom.coeff >= 0.5) & (crosscorr$rel.phantom.coeff < 1) ) {
    crosscorr$phantom.quality.tag <- 0
} else if ( (crosscorr$rel.phantom.coeff >= 1) & (crosscorr$rel.phantom.coeff < 1.5) ) {
    crosscorr$phantom.quality.tag <- 1
} else if ( (crosscorr$rel.phantom.coeff >= 1.5) ) {
    crosscorr$phantom.quality.tag <- 2
}
cat("Phantom Peak Quality Tag",crosscorr$phantom.quality.tag,"\n",file=stdout())

# Output result to result file if required
#Filename\tnumReads\tPeak_shift\tPeak_Correlation\tRead_length\tPhantomPeak_Correlation\tMin_Correlation_Shift\tMin_Correlation\tNormalized_CrossCorrelation_Coefficient\tRelative_CrossCorrelation_Coefficient\tQualityTag)
if (! is.na(iparams$output.result.file)) {
    cat(get.file.parts(iparams$ip.file)$fullname,
            ip.data$num.tags,
            paste(cc.peak$x,collapse=","),
            paste(cc.peak$y,collapse=","),
            crosscorr$phantom.cc$x,
            crosscorr$phantom.cc$y,
            crosscorr$min.cc$x,
            crosscorr$min.cc$y,
            crosscorr$phantom.coeff,
            crosscorr$rel.phantom.coeff,
            crosscorr$phantom.quality.tag,
            sep="\t",
            file=iparams$output.result.file,
            append=TRUE)
    cat("\n",
            file=iparams$output.result.file,
            append=TRUE)    
}

# Save figure if required
if (! is.na(iparams$output.plot.file)) {
    pdf(file=iparams$output.plot.file,width=5,height=5)
    par(mar = c(4,3.5,2,0.5), mgp = c(1.5,0.5,0), cex = 0.8);
    plot(crosscorr$cross.correlation,
            type='l',
            xlab=sprintf("strand-shift (%s)",paste(cc.peak$x,collapse=",")),
            ylab="cross-correlation")
    abline(v=cc.peak$x,lty=2,col=2)
    abline(v=crosscorr$phantom.cc$x,lty=2,col=4)
    title(main=get.file.parts(iparams$ip.file)$fullname,
          sub=sprintf("NSC=%g,RSC=%g,Qtag=%d",crosscorr$phantom.coeff,crosscorr$rel.phantom.coeff,crosscorr$phantom.quality.tag))
    dev.off();    
}

# Save RData file if required
if (! is.na(iparams$output.rdata.file)) {
    save(iparams,
            crosscorr,
            cc.peak,
            file=iparams$output.rdata.file);    
}

# #################################    
# Call peaks
# #################################
if ( !is.na(iparams$output.npeak.file) || !is.na(iparams$output.rpeak.file) ) {
    
    # Remove local tag anomalies
    if(iparams$rmAbnormal == 1){
        cat('Removing read stacks\n',file=stdout())
        ip.data <- remove.local.tag.anomalies(ip.data$tags)
        if (! is.na(iparams$input4ip.file)) {
            input4ip.data <- remove.local.tag.anomalies(input4ip.data$tags)
        }
        if (! is.na(iparams$mock.file)) {
            mock.data <- remove.local.tag.anomalies(mock.data$tags)
        }
        if (! is.na(iparams$input4mock.file)) {
            input4mock.data <- remove.local.tag.anomalies(input4mock.data$tags)
        }
    }

    if(iparams$flagShift == 0){
        tagShift <- 0
    }else{
        tagShift <- as.numeric(round(crosscorr$peak$x/2))
        tagShift <- mean(tagShift)
    }

    # Open multiple processes if required
    if (is.na(iparams$n.nodes)) {
        cluster.nodes <- NULL
    } else {
        cluster.nodes <- makeCluster(iparams$n.nodes,type="SOCK")
    }
    
    # Find peaks
    cat('Finding peaks\n',file=stdout())
    if (!is.na(iparams$npeak)) {
        iparams$fdr <- 0.99999
    }


    if(iparams$mthd == 'mcbin'){
        if(is.na(iparams$input4mock.file) | is.na(iparams$input4ip.file) | is.na(iparams$mock.file) | is.na(iparams$ip.file) ){
            stop ("missing chip-seq data files for mcbin\n")
        }
        data.a <- ip.data
        data.b <- mock.data
        data.ap <- input4ip.data
        data.bp <- input4mock.data
        if(iparams$rgn == "mock"){
            narrow.peaks <- find.binding.positions(signal.data=data.a,control.data=data.b,fdr=iparams$fdr,method=tag.lwcc,whs=crosscorr$whs,cluster=cluster.nodes,tec.filter=T,enrichment.z=0,min.thr=0,background.density.scaling = T,min.mle.threshold=0,e.value= 1000000000,enrichment.background.scales=c(1))
        }else if( iparams$rgn == "input"){
            narrow.peaks <- find.binding.positions(signal.data=data.a,control.data=data.ap,fdr=iparams$fdr,method=tag.lwcc,whs=crosscorr$whs,cluster=cluster.nodes,tec.filter=T,enrichment.z=0,min.thr=0,background.density.scaling = T,min.mle.threshold=0,e.value  =1000000000,enrichment.background.scales=c(1))
        }else{
        stop("-rgn= parameter is missing mock or input?")
    }

        if (!is.na(iparams$n.nodes)) {
            stopCluster(cluster.nodes)
        }
        cat(paste("Detected",sum(unlist(lapply(narrow.peaks$npl,function(d) length(d$x)))),"peaks"),"\n",file=stdout())

        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, pv = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.diff = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.a = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.ap = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.bp = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.b = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.ratio = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.diffRatio = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, pvab = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, pvaap = NA);

        # Compute and write regionPeak file
        if (!is.na(iparams$output.rpeak.file)) {
            if(iparams$rgn == "mock"){
                region.peaks <- add.broad.peak.regions(data.a,data.b,narrow.peaks,window.size=max(50,round(crosscorr$whs/4)),z.thr=10)
            }else{
                region.peaks <- add.broad.peak.regions(data.a,data.ap,narrow.peaks,window.size=max(50,round(crosscorr$whs/4)),z.thr=10)
            }
            region.peaks$npl <- add.rsre(region.peaks,margin=round(crosscorr$whs/2))
            region.peaks$npl <- add.count4(region.peaks,data.a,data.b,data.ap,data.bp,iparams$totReads,tagShift,iparams$wdfold)
            region.peaks$npl <- add.pois(region.peaks,"both")
            write.narrowpeak.mcbin(region.peaks,iparams$output.rpeak.file,rby=iparams$rankby,mcstep = iparams$mcstep,coe.mc = iparams$weightmc, coe.spp = iparams$weightspp,npeaks=iparams$npeak)
            system(paste('gzip -f ',iparams$output.rpeak.file))
        }
    }else if(iparams$mthd == 'poisson'){
        data.a <- ip.data
        if(! is.na(iparams$input4ip.file) & is.na(iparams$mock.file)){
            data.ap <- input4ip.data
        }else if(is.na(iparams$input4ip.file) & ! is.na(iparams$mock.file)){
            data.ap <- mock.data
        }else{
            stop('chip-seq data file prolem')
        }
        narrow.peaks <- find.binding.positions(signal.data=data.a,control.data=data.ap,fdr=iparams$fdr,method=tag.lwcc,whs=crosscorr$whs,cluster=cluster.nodes,tec.filter=T,enrichment.z=0,min.thr=0,background.density.scaling = T,min.mle.threshold=0,e.value=1000000000,enrichment.background.scales=c(1))
        if (!is.na(iparams$n.nodes)) {
            stopCluster(cluster.nodes)
        }
        cat(paste("Detected",sum(unlist(lapply(narrow.peaks$npl,function(d) length(d$x)))),"peaks"),"\n",file=stdout())

        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, pv = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.diff = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.a = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.ap = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.bp = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.b = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.ratio = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, nt.diffRatio = NA)
        narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, pvab = NA); narrow.peaks$npl <- lapply(narrow.peaks$npl, cbind, pvaap = NA);

        # Compute and write regionPeak file
        if (!is.na(iparams$output.rpeak.file)) {
            region.peaks <- add.broad.peak.regions(data.a,data.ap,narrow.peaks,window.size=max(50,round(crosscorr$whs/4)),z.thr=10)
            region.peaks$npl <- add.rsre(region.peaks,margin=round(crosscorr$whs/2))
            region.peaks$npl <- add.count2(region.peaks,data.a,data.ap,iparams$totReads,tagShift,iparams$wdfold)
            region.peaks$npl <- add.pois(region.peaks,"aap")
            write.narrowpeak.pois(region.peaks,iparams$output.rpeak.file,rby=iparams$rankby,npeaks=iparams$npeak)
            system(paste('gzip -f ',iparams$output.rpeak.file))
        }
    }else{
        stop('-mthd= poisson or mcbin')
    }
    
    # Write to narrowPeak file
    if (!is.na(iparams$output.npeak.file)) {
        write.narrowpeak.binding(narrow.peaks,iparams$output.npeak.file,margin=round(crosscorr$whs/2),npeaks=iparams$npeak)
        system(paste('gzip -f ',iparams$output.npeak.file))
    }
        
    # Save Rdata file    
    if (! is.na(iparams$output.rdata.file)) {
        save(iparams,
                crosscorr,
                cc.peak,
                narrow.peaks,
                region.peaks,
                file=iparams$output.rdata.file);    
    }
    
}

tag.shift <- round(crosscorr$peak$x/2)
cat("tag.shift: ")
cat(tag.shift)
cat("\n")

whs <- round(crosscorr$whs)
cat("whs: ")
cat(whs)
cat("\n")

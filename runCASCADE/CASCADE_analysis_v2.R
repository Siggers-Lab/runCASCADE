#!/usr/bin/env Rscript

### CASCADE ANALYSIS command-line program

### DESCRIPTION ------------------------------------------------------------------------

# analyzes and generates plots for CASCADE-formatted nextPBM experiments
#
# modules that are run here:
#     sorting data matrices and saving them as intermediates
#     integration of sorted data matrices into a unified CASCADE table
#     output of ALL possible scatterplot comparisons

### INSTRUCTIONS -----------------------------------------------------------------------

# DEPENDENCIES
#     R/3.5.2 (on the SCC cluster)
#     RColorBrewer
#     ggplot2
#     cowplot
#     ggseqlogo
#     MotifDb
#     gridExtra
#
# USAGE:
# Rscript CASCADE_analysis.R args[1] args[2] args[3] args[4] args[5] args[6] args[7]
#
# ARGUMENTS:
#     args[1]     output directory name
#     args[2]     experiment name prefix
#     args[3]     path to annotation table
#     args[4]     path to directory with data matrices
#     args[5]     path to text file with the array conditions (to be used as column names)
#     args[6]     path to text file with previously annotated TF sites
#     args[7]     TRUE or FALSE for whether PWM similarity suite should be run

### SORT THE DATA MATRICES BY PROBE NAME AND SAVE AN INTERMEDIATE ----------------------

# set up environment
rm(list=ls())
options(scipen=999)
options(digits=22)

# check installation of packages
library.path <- "C:/Users/ykoga07/Documents/R/win-library/3.6"
list.of.packages <- c("ggplot2", "cowplot","ggseqlogo","MotifDb","gridExtra", "crayon")
 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc = library.path)[,"Package"])]
if(length(new.packages)) install.packages(new.packages, 
                                           lib = library.path,
                                           repos = "http://cran.us.r-project.org")


# parse the command line arguments
args <- commandArgs(trailingOnly=TRUE)

# check if files exist
if(!file.exists(args[3])){
stop(paste0("The annotation table ", args[3], " does not exist. Please check if this is the correct file/directory name."))}

if(!file.exists(args[4])){
stop(paste0("The directory ", args[4], " does not exist. Please check if this is the correct directory name."))}

if(!file.exists(args[5])){
stop(paste0("The text file ", args[5], " does not exist. Please check if this is the correct file/directory name."))}

if(!file.exists(args[6])){
stop(paste0("The text file ", args[6], " does not exist. Please check if this is the correct file/directory name."))}

# output status message
message("sorting nextPBM data matrices...")

# create a directory for the analysis using the user-specified prefix
dir.create(paste(getwd(), args[1], sep="/"))

# create a directory for the sorted data files
dir.create(paste(getwd(), args[1], "sorted_data", sep="/"))

# declare orientations
or <- c("o1", "o2", "or", "br")

# for each orientation
for (o in or) {
  
  # open the data matrix corresponding to the current matrix
  x <- read.table(paste(args[4], "/", args[2], "_", o, ".dat", sep=""))
  
  # sort by the probe name in the first column
  x <- x[order(x$V1),]
  
  # write the sorted file
  write.table(x, file=paste(getwd(), "/", args[1], "/", "sorted_data/", args[2], "_", o, "_sorted.dat", sep=""), sep='\t', quote=F, row.names=F, col.names=F)
  
}



### INTEGRATE PBM RESULTS INTO CASCADE ANNOTATION TABLE --------------------------------------

# output status message
message("integrating nextPBM data into CASCADE annotation table...")

# read in the full annotation table (via args[3])
df_name <- args[3]
df <- read.table(df_name, header=T, sep='\t')

# generate a probeID column using the probe_name (to be able to merge data into this table)
df$probeID <- gsub("_o[1-5]_r[1-5]", "", df$probe_name)

# collapse annotation by the number of replicates (one entry per unique target sequence)
df <- df[which(df$probe_repl=="r1" & df$probe_or=="o1"),]

# declare the different orientations for the PBM probes (to be read/integrated into the larger table)
or <- c("o1", "o2", "or", "br")

# declare the different PBM results conditions (in order, via args[5])
PBM_cond <- read.table(args[5], sep="\t")
PBM_cond <- as.character(PBM_cond$V1)

# for each orientation and each condition read and integrate the nextPBM results
for (o in or) {
  
  # read the results table for the given orientation
  x <- read.table(paste(getwd(), "/", args[1], "/", "sorted_data/", args[2], "_", o, "_sorted.dat", sep=""))
  
  # change names of the columns
  names(x) <- c("probeID", "probe_seq", "dummy", PBM_cond)
  
  # modify the probeID to not contain the orientation tag
  x$probeID <- gsub("_o[1-5]", "", x$probeID)
  
  # bind every PBM condition to the growing table
  for (i in PBM_cond) {
    
    # keep only the columns needed to bind to the larger annotation
    y <- x[, c("probeID", i)]
    
    # change name of the data column to contain the run and the orientation tag
    names(y)[2] <- paste(args[1], o, names(y)[2], sep="_")
    
    # merge results into the larger annotation
    df <- merge(df, y, by="probeID", sort=F, all=T)
    
  }
}

# save the table as an intermediate
write.table(df, file=paste(args[1], "/", args[2], "_FULL_ANNOT_", args[1], ".bed", sep=""), col.names=T, row.names=F, sep='\t', quote=F)



### ORIENTATION (o1 vs. o2) SCATTERS -----------------------------------------------------------

# output a status message
message("generating orientation (o1 vs. o2) scatterplots...")

# import libraries
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# create a directory to house the orientation scatterplots
dir.create(paste(getwd(), args[1], "o1_vs_o2_scatters", sep="/"))

# read in the full results table for the current experiment
df <- read.table(paste(args[1], "/", args[2], "_FULL_ANNOT_", args[1], ".bed", sep=""), header=T, sep='\t')

# collect the unique IDs of all the background probes (needed for z-score calculations)
bg_IDs <- droplevels(df[which(df$seed_names=="BG"), "probeID"])
bg_rows <- which(df$probeID %in% bg_IDs)

# create a data frame with just the data columns
PBM_dat <- df[,18:ncol(df)]

# helper function to compute column-wise z-scores
fluo_to_zscore <- function(fluo) {
  # log transform the fluorescence value
  log_fluo <- suppressWarnings(log(fluo))
  
  # save the mean and standard deviation for the background probe fluorescence values
  bg_mean <- mean(log_fluo[bg_rows], na.rm=T)
  bg_sd <- sd(log_fluo[bg_rows], na.rm=T)
  
  # compute zscore for the vector of log(fluo) values
  zscores <- (log_fluo-bg_mean)/bg_sd
  
  # return the vector of zscores to the parent function
  return(zscores)
}

# apply the function to each column in the PBM
PBM_dat <- as.data.frame(lapply(PBM_dat, fluo_to_zscore))

# bind the unique IDs and probe sequences back to the numerical data prior to merging
df <- cbind(df[,1:17], PBM_dat)
rm(PBM_dat)

# declare names of seed regions of interest
to_keep <- c("BACKGROUND", "STATIC", "TILE")
x <- droplevels(df[which(df$probe_type %in% to_keep),])
seeds <- unique(as.character(x$seed_names))
rm(x)

# read in the PBM conditions again
PBM_cond <- read.table(args[5], sep="\t")
PBM_cond <- as.character(PBM_cond$V1)

# use the number of seeds to determine which colors to use for the plot
if (length(seeds)<=8) {
  
  # there are enough colors in the default palette to plot
  PBM_colors <- brewer.pal(8, "Set2")
  
} else {
  
  # construct a palette with the correct number of colors (equal to the number of seeds to be plotted) using the Set2 as a reference
  PBM_colors <- colorRampPalette(brewer.pal(8, "Set2"))
  PBM_colors <- PBM_colors(length(seeds))
}

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# generate different plots for the seeds and seeds+SNVs plots
for (cond in c("SEED", "SNV") ) {
  
  # keep only probes that need to be plotted
  x <- droplevels(df[which(df$seed_names %in% seeds),])
  
  # drop the SNV probes if plotting the seeds only
  if (cond == "SEED") {
    x <- droplevels(x[which(x$SNV_pos_offset==0),])
  }
  
  # for each of the PBM experimental conditions tested in this set of experiments
  for (i in 1:length(PBM_cond)) {
    
    # declare a plot variable for the current scatter. Needs to be done this way because of cowplot
    ggp <- ggplot(x, aes_string(y=paste(args[1], "o1", PBM_cond[i], sep="_"), x=paste(args[1], "o2", PBM_cond[i], sep="_"), color="seed_names")) + 
      geom_point(size=2.5, alpha=0.4) + scale_color_manual(values=PBM_colors) +
      scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
      ggtitle(PBM_cond[i]) +
      ylab("o1 z-score") + xlab("o2 z-score")
    
    # assign as a separate variable
    assign(paste("ggp_", i, sep=""), ggp)
    
  }
  
  # use cowplot to save the scatterplots
  pdf_name <- paste("o1_vs_o2", cond, "scatters.pdf", sep="_")
  ggp_grid <- plot_grid(plotlist=mget(paste0("ggp_", 1:length(PBM_cond))), ncol = 2, nrow = length(PBM_cond)/2, label_size = 8)
  save_plot(paste(getwd(), args[1], "o1_vs_o2_scatters", pdf_name, sep="/"), ggp_grid, ncol = 2, nrow = length(PBM_cond)/2, base_height=4, base_aspect_ratio=1.4, limitsize=FALSE)
  
}



### PAIRWISE COMPARISION SCATTERS -----------------------------------------------------------

# output a status message
message("generating pairwise comparison scatterplots...")

# create a directory to house the orientation scatterplots
dir.create(paste(getwd(), args[1], "pairwise_scatters", sep="/"))

# declare names of seed regions of interest
to_keep <- c("BACKGROUND", "STATIC", "TILE")
x <- droplevels(df[which(df$probe_type %in% to_keep),])
seeds <- unique(as.character(x$seed_names))
rm(x)

# read in the PBM conditions again
PBM_cond <- read.table(args[5], sep="\t")
PBM_cond <- as.character(PBM_cond$V1)

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# generate different plots for the seeds and seeds+SNVs plots
for (cond in c("SEED", "SNV") ) {
  
  # keep only probes that need to be plotted
  x <- droplevels(df[which(df$seed_names %in% seeds),])
  
  # drop the SNV probes if plotting the seeds only
  if (cond == "SEED") {
    x <- droplevels(x[which(x$SNV_pos_offset==0),])
  }
  
  # plot the pairwise scatters across all PBM conditions (compared against all PBM conditions in a grid)
  for (i in 1:length(PBM_cond)) {
    for (j in 1:length(PBM_cond)) {
      
      # declare a plot variable for the current scatter. Needs to be done this way because of cowplot
      ggp <- ggplot(x, aes_string(y=paste(args[1], "br", PBM_cond[i], sep="_"), x=paste(args[1], "br", PBM_cond[j], sep="_"), color="seed_names")) +
        geom_point(size=2.5, alpha=0.4) + scale_color_manual(values=PBM_colors) +
        scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
        ggtitle(paste(PBM_cond[i], " vs. ", PBM_cond[j], sep="")) +
        ylab(paste(PBM_cond[i], " z-score", sep="")) + xlab(paste(PBM_cond[j], " z-score", sep="")) +
        theme(legend.position="right")
      
      # assign as a separate variable
      assign(paste("ggp_", j, sep=""), ggp)
      
    }
    
    # use cowplot to save the scatterplots
    pdf_name <- paste("pairwise", cond, PBM_cond[i], "scatters.pdf", sep="_")
    ggp_grid <- plot_grid(plotlist=mget(paste0("ggp_", 1:length(PBM_cond))), ncol = length(PBM_cond), nrow = 1, label_size = 8)
    save_plot(paste(getwd(), args[1], "pairwise_scatters", pdf_name, sep="/"), ggp_grid, ncol = length(PBM_cond), nrow = 1, base_height=4, base_aspect_ratio=1.4, limitsize=FALSE)
    
  }
}



### DECLARE FUNCTIONS THAT NEED TO BE USED BY PWM SIMILARITY SUITE ----------------------

# checks pwm arguments for PWMSimilarity()
normargPwm <- function(pwm, argname="pwm") {
  if (!is.matrix(pwm) || !is.numeric(pwm))
    stop("'", argname, "' must be a numeric matrix")
  if (!identical(rownames(pwm), DNA_BASES))
    stop("'rownames(", argname, ")' must be the 4 DNA bases ('DNA_BASES')")
  if (!is.double(pwm))
    storage.mode(pwm) <- "double"
  if (any(is.na(pwm)))
    stop("'", argname, "' contains NAs")
  pwm
}

# source: https://github.com/ge11232002/TFBSTools/blob/master/R/PWM-methods.r (TFBSTools bioconductor package)

# implemented from (Harbison et al. 2004)
PWMEuclidean = function(pwm1, pwm2) {
  # now the pwm1 and pwm2 must have same widths
  stopifnot(isConstant(c(ncol(pwm1), ncol(pwm2))))
  pwm1 = normargPwm(pwm1)
  pwm2 = normargPwm(pwm2)
  width = ncol(pwm1)
  diffMatrix = (pwm1 - pwm2)^2
  PWMDistance = sum(sqrt(colSums(diffMatrix))) / sqrt(2) / width
  return(PWMDistance) 
}

PWMModPearson = function(pwm1, pwm2) {
  # now the pwm1 and pwm2 must have the same widths
  stopifnot(isConstant(c(ncol(pwm1), ncol(pwm2))))
  pwm1 = normargPwm(pwm1)
  pwm2 = normargPwm(pwm2)
  # obtain diagonal of the correlation matrix between pwm1 and pwm2
  r <- diag(cor(pwm1, pwm2, method="pearson"), names=F)
  # take 1 - the pearson
  r <- 1 - r
  # obtain product across all columns
  r <- prod(r)
  return(r)
}

PWMPearson = function(pwm1, pwm2) {
  # now the pwm1 and pwm2 must have the same widths
  stopifnot(isConstant(c(ncol(pwm1), ncol(pwm2))))
  pwm1 = normargPwm(pwm1)
  pwm2 = normargPwm(pwm2)
  top = colSums((pwm1 - 0.25) * (pwm2 - 0.25))
  bottom = sqrt(colSums((pwm1 - 0.25)^2) * colSums((pwm2 - 0.25)^2))
  r = 1 / ncol(pwm1) * sum((top / bottom))
  return(r)
}

PWMKL = function(pwm1, pwm2) {
  # now the pwm1 and pwm2 must have the same widths
  stopifnot(isConstant(c(ncol(pwm1), ncol(pwm2))))
  pwm1 = normargPwm(pwm1)
  pwm2 = normargPwm(pwm2)
  KL = 0.5 / ncol(pwm1) * sum(colSums(pwm1 * log(pwm1 / pwm2) +
                                        pwm2 * log(pwm2 / pwm1)))
  return(KL)
}

PWMSimilarity <- function(pwmSubject, pwmQuery, method=c("Euclidean", "ModPearson", "Pearson", "KL")) {
  pwm1 = pwmSubject
  pwm2 = pwmQuery
  method = match.arg(method)
  widthMin = min(ncol(pwm1), ncol(pwm2))
  ans = c()
  for(i in 1:(1+ncol(pwm1)-widthMin)){
    for(j in 1:(1+ncol(pwm2)-widthMin)){
      pwm1Temp = pwm1[ ,i:(i+widthMin-1)]
      pwm2Temp = pwm2[ ,j:(j+widthMin-1)]
      ansTemp = switch(method,
                       "Euclidean"=PWMEuclidean(pwm1Temp, pwm2Temp),
                       "ModPearson"=PWMModPearson(pwm1Temp, pwm2Temp),
                       "Pearson"=PWMPearson(pwm1Temp, pwm2Temp),
                       "KL"=PWMKL(pwm1Temp, pwm2Temp))
      ans = c(ans, ansTemp)
    }
  }
  ans = switch(method,
               "Euclidean"=min(ans),
               "ModPearson"=min(ans),
               "Pearson"=max(ans),
               "KL"=min(ans)
  )
  return(ans)
}

### PLOT ALL THE LOCUS LOGOS AND Z-SCORE CONFETTI PLOTS -----------------------------------------------------------------

# print a status message
message("plotting locus logos and z-score confetti plots...")

# set up environment
library(ggseqlogo)
library(MotifDb)

# function to obtain reverse complement PWM
rev_compl_pwm <- function(m) {
  # reverse the input pwm
  rev_pwm <- m[,ncol(m):1]
  
  # reverse complement the pwm (A to T, C to G)
  rc_pwm <- rev_pwm
  rc_pwm[1,] <- rev_pwm[4,]
  rc_pwm[2,] <- rev_pwm[3,]
  rc_pwm[3,] <- rev_pwm[2,]
  rc_pwm[4,] <- rev_pwm[1,]
  rownames(rc_pwm) <- rownames(m)
  
  # return reversed pwm
  return(rc_pwm)
}

# function to obtain reverse complement of a DNA sequence
rev_compl <- function(DNA_seq){
  # returns the input DNA sequence in reverse
  temp_seq <- unlist(strsplit(DNA_seq, split=''))
  temp_seq <- rev(temp_seq)
  temp_seq <- paste(temp_seq, collapse='')
  temp_seq <- chartr("ACGTN", "TGCAN", temp_seq)
  return(temp_seq)
}

# function that takes a SNV matrix, converts to PWM using the supplied beta parameter, and performs motif similarity
SNV_to_PWM_sim <- function(SNV, SNV_trans, beta) {
  
  # assign the y axis label based on the unit being plotted
  y_axis_label <- expression(paste(Delta, "z-score", sep=""))
  
  # scale transformation function for rounding axes
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  # declare a spacer length to use for within a motif
  spacer <- 3
  
  # declare a minimum motif length to enforce
  min_motif_length <- 6
  
  # declare minimum stack heigh for high scoring positions
  stack_thresh <- 0.5
  
  # initialize high-scoring position vector
  high_score_pos <- vector(mode="integer")
  
  # determine the positive stack heights of each position
  SNV_trans_pos <- SNV_trans
  SNV_trans_pos[SNV_trans_pos<0] <- 0
  
  # get the colSums to take a cumulative stack height (for just the positive values)
  colsums_SNV <- colSums(SNV_trans_pos)
  
  # for each position in the SNV matrix
  for (k in 1:ncol(SNV_trans)) {
    
    # if the current position is high-scoring
    if (colsums_SNV[k] > stack_thresh) {
      
      # add the current position as a high-scoring position
      high_score_pos <- c(high_score_pos, k)
    }
  }
  
  # if there are enough high-scoring positions to start a motif
  if (length(high_score_pos)>1) {
    
    # initialize motif start and end vectors to keep track
    motif_starts <- vector(mode="integer")
    motif_ends <- vector(mode="integer")
    
    # initialize a current motif start and end using the first high-scoring position
    curr_start <- high_score_pos[1]
    curr_end <- high_score_pos[1]
    
    # for the remaining high-scoring positions from the SNV matrix
    for (k in 2:length(high_score_pos)) {
      
      # if the current position is beyond the declared spacer length
      if (high_score_pos[k]-high_score_pos[k-1] > spacer) {
        
        # if the current motif start and end represent a motif
        if (length(curr_start:curr_end) >= min_motif_length) {
          
          # push the current motif
          motif_starts <- c(motif_starts, curr_start)
          motif_ends <- c(motif_ends, curr_end)
          
          # use the current position as the new motif start and end
          curr_start <- high_score_pos[k]
          curr_end <- high_score_pos[k]
          
        }
        
        # if the current motif start and end DO NOT represent a motif
        if (length(curr_start:curr_end) < min_motif_length) {
          
          # DO NOT update the current motifs vectors
          
          # use the current position as the new motif start and end
          curr_start <- high_score_pos[k]
          curr_end <- high_score_pos[k]
          
        }
        
      }
      
      # if the current position is within the spacer length allowance
      if (high_score_pos[k]-high_score_pos[k-1] <= spacer) {
        
        # update the current motif end to extend the current motif forward
        curr_end <- high_score_pos[k]
        
        # if this current position happens to be the final high-scoring position, motif needs to be pushed as-is
        if (k == length(high_score_pos)) {
          
          # if the current motif start and end represent a motif
          if (length(curr_start:curr_end) >= min_motif_length) {
            
            # push the current motif
            motif_starts <- c(motif_starts, curr_start)
            motif_ends <- c(motif_ends, curr_end)
          }
        }
        
      }
    }
    
    # if there are motifs within the locus SNV matrix to plot
    if (length(motif_starts)>0) {
      
      # for each of the motifs uncovered
      for (x in 1:length(motif_starts)) {
        
        # isolate the SNV matrix positions relevant to the current motif
        SNV_motif <- SNV_trans[,(motif_starts[x]:motif_ends[x])]
        rownames(SNV_motif) <- c("A", "C", "G", "T")
        
        # reverse the motif if direction passed to the function is REVERSE
        if (dir == "REVERSE") {
          SNV_motif <- rev_compl_pwm(SNV_motif)
        }
        
        # save the isolated SNV matrix to file
        SNV_name <- paste(seed, cond, dir, "MOTIF", as.character(x), "SNV_matrix.txt", sep="_")
        write.table(SNV_motif, paste(getwd(), args[1], "motif_similarity", SNV_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        # assign the y axis label based on the unit being plotted
        y_axis_label <- expression(paste(Delta, "z-score", sep=""))
        
        # save the isolated SNV motif logo as a variable
        ggp_SNV <- ggseqlogo(SNV_motif, method="custom")
        my.ggp.yrange <- ggplot_build(ggp_SNV)$layout$panel_params[[1]]$y.range
        ggp_SNV <- ggp_SNV +
          scale_x_continuous(expand = c(0.01,0.01)) +
          scale_y_continuous(expand = c(0.01,0.01), labels=scaleFUN) +
          annotate('rect', xmin = 0.5, xmax = ncol(SNV_motif)+0.5, ymin = my.ggp.yrange[1], ymax = 0.0, alpha = 0.75, col='white', fill='white') +
          annotate('rect', xmin = 0.5, xmax = ncol(SNV_motif)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=1.1) +
          annotate('rect', xmin = 0.5, xmax = ncol(SNV_motif)+0.5, ymin = my.ggp.yrange[1], ymax = 0, col="lightgray", fill=NA, size=1.1) +
          ylab(y_axis_label) +
          theme(axis.text.x=element_blank()) +
          theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain"))
        
        # convert the SNV matrix to a PWM using the beta parameter (USING THE RELEVANT POSITIONS ONLY FROM MOTIF FINDING)
        beta <- 2.0
        SNV_pwm_pos <- SNV[,(motif_starts[x]:motif_ends[x])]
        SNV_pwm <- exp(beta * SNV_pwm_pos)
        
        # compute probabilities for PWM using the column-wise totals
        PWM <- matrix(nrow=4, ncol=ncol(SNV_pwm))
        for (i in 1:nrow(PWM)) {
          for (j in 1:ncol(PWM)) {
            PWM[i, j] <- SNV_pwm[i, j]/sum(SNV_pwm[,j])
          }
        }
        rownames(PWM) <- c("A", "C", "G", "T")
        
        # reverse the PWM if necessary
        if (dir == "REVERSE") {
          PWM <- rev_compl_pwm(PWM)
        }
        
        # save the isolated PWM matrix to file
        PWM_name <- paste(seed, cond, dir, "MOTIF", as.character(x), "PWM_matrix.txt", sep="_")
        write.table(PWM, paste(getwd(), args[1], "motif_similarity", PWM_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
        
        # save PWM logo as a variable to plot later
        ggp_PWM <- ggseqlogo(PWM)
        ggp_PWM <- ggp_PWM + 
          ylim(0, 2) +
          scale_x_continuous(expand = c(0.01,0.01)) +
          scale_y_continuous(expand = c(0.01,0.01)) +
          annotate('rect', xmin = 0.5, xmax = ncol(PWM)+0.5, ymin = 0, ymax = 1, col="lightgray", fill=NA, size=1.1) +
          annotate('rect', xmin = 0.5, xmax = ncol(PWM)+0.5, ymin = 1, ymax = 2, col="lightgray", fill=NA, size=1.1) +
          theme(axis.text.x=element_blank()) +
          theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain"))
        
        # plot the matching SNV motif and PWM motif logos to the same file using cowplot
        pdf_name <- paste(seed, cond, dir, "MOTIF", as.character(x), "logos.pdf", sep="_")
        ggp_logo <- plot_grid(ggp_SNV, ggp_PWM,
                              ncol = 1, nrow = 2, rel_heights = c(1, 1), align = "v")
        save_plot(paste(getwd(), args[1], "motif_similarity", pdf_name, sep="/"), ggp_logo, ncol = 1, nrow = 2, base_height = 2.5, base_width = ncol(PWM)/2, base_aspect_ratio = 1, limitsize=FALSE)
        
        # only perform the motif similarity during the FORWARD run. It will do both FORWARD and REVERSE in the same step
        if (dir == "FORWARD") {
          
          # check for and apply pseudocounts if a uniform column exists in the PWM
          for (col_idx in 1:ncol(PWM)) {
            if (isTRUE(all.equal(as.vector(PWM[,col_idx]), c(0.25, 0.25, 0.25, 0.25)))) {
              PWM[,col_idx] <- c(0.250000000000001, 0.249999999999999, 0.250000000000001, 0.249999999999999)
            }
          }
          
          # compare SNV against motif database
          motif_sims_forward <- lapply(motif_db, function(m) PWMSimilarity(PWM, m, method="ModPearson"))
          motif_sims_reverse <- lapply(motif_db, function(m) PWMSimilarity(rev_compl_pwm(PWM), m, method="ModPearson"))
          motif_sims <- c(motif_sims_forward, motif_sims_reverse)
          
          # collect the top motifs
          top_motifs <- motif_sims[order(unlist(motif_sims))]
          top_motifs <- top_motifs[which(!duplicated(top_motifs))]
          top_motifs <- top_motifs[1:25]
          top_scores <- top_motifs
          top_motifs <- as.list(motif_db[names(top_motifs)])
          
          
          # create plot with the top motifs
          p_list <- lapply(1:length(top_motifs), function(m) {
            ggseqlogo(top_motifs[m]) +
              ggtitle(paste("", as.character(formatC(top_scores[m][[1]], format="e", digits=5)), sep='\n')) +
              scale_x_continuous(expand = c(0.01,0.01)) +
              scale_y_continuous(expand = c(0.01,0.01)) +
              ylim(0, 2) +
              theme(plot.title=element_text(size=6), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
          })
          
          # save database motif matches to file
          pdf_name <- paste(seed, cond, beta, "MOTIF", as.character(x), "MotifDb_matches.pdf", sep="_")
          pdf(file=paste(getwd(), args[1], "motif_similarity", pdf_name, sep="/"), width=3.5, height=1.2*length(top_motifs))
          do.call(gridExtra::grid.arrange, c(p_list, nrow=length(top_motifs)))
          dev.off()
        }
        
      }
      
    }
    
  }
  
}

# function to obtain SNV value matrix (using seed-weighted mean - gain-of-function variants are also flagged and reverted to seed z-score value)
locus_SNV_matrix <- function(PBM, PBM_col, genomic_window) {
  
  # initialize a matrix to hold the SNV probe z-score values
  nucs <- c("A", "C", "G", "T")
  SNV <- matrix(nrow=4, ncol=length(genomic_window))
  
  # collect fluo (treatment/UT) values for each position in the matrix
  for (i in 1:nrow(SNV)) {
    for (j in 1:ncol(SNV)) {
      
      # determine the current genomic position
      genomic_pos <- genomic_window[j]
      
      # collect data relevant to the SNV probes at the current genomic position
      SNV_probes <- PBM[which(PBM$SNV_pos==genomic_pos),]
      
      # determine the current genomic nucleotide
      curr_seed_nuc <- unique(as.character(SNV_probes$seed_nuc))
      
      # if the current SNV nucleotide is the genomic nucleotide
      if (curr_seed_nuc == nucs[i]) {
        
        # collect tile numbers that contain the given SNV position
        tiles <- unique(SNV_probes$tile_order)
        
        # collect info about the tile probes that contain the current nucleotide at the SNV position
        curr_nuc_TILE_probes <- PBM[which(PBM$tile_order %in% tiles & PBM$SNV_pos_offset==0),]
        
        # determine weights of genomic tile probes to use for the weighted average
        weights <- curr_nuc_TILE_probes[, PBM_col]
        
        # save the z-scores of the current TILE probes for use in the weighted mean formula (in a way that preserves tile order used in prev step)
        probe_z_scores <- curr_nuc_TILE_probes[, PBM_col]
        
        # replace non-finite values with 0
        weights[!is.finite(weights)] <- 0
        probe_z_scores[!is.finite(probe_z_scores)] <- 0
        
        # replace non-positive values with a pseudo-count
        weights[weights<=0] <- 0.000001
        probe_z_scores[probe_z_scores<=0] <- 0.000001
        
        # assign the SNV value as a function of the PBM values for the genomic tiles containing the given SNV position
        SNV[i,j] <- weighted.mean(probe_z_scores, weights)
        
      } else { # if the current nucleotide corresponds to a SNV of the genomic seed nucleotide
        
        # collect tile numbers that contain the given SNV position
        tiles <- unique(SNV_probes$tile_order)
        
        # collect info about the tile probes that contain the current nucleotide at the SNV position
        curr_nuc_TILE_probes <- PBM[which(PBM$tile_order %in% tiles & PBM$SNV_pos_offset==0),]
        
        # collect info about the SNV probes that contain the current nucleotide at the SNV position
        curr_nuc_SNV_probes <- SNV_probes[which(SNV_probes$SNV_nuc==nucs[i]),]
        
        # determine weights of genomic tile probes to use
        weights <- vector(mode="numeric")
        for (k in 1:nrow(curr_nuc_SNV_probes)) {
          weight_k <- curr_nuc_TILE_probes[which(curr_nuc_TILE_probes$tile_order==curr_nuc_SNV_probes[k, "tile_order"]), PBM_col]
          weights <- c(weights, weight_k)
        }
        
        # save the z-scores of the current SNV probes for use in the weighted mean formula (in a way that preserves tile order used in prev step)
        probe_z_scores <- curr_nuc_SNV_probes[, PBM_col]
        
        # replace non-finite values with 0
        weights[!is.finite(weights)] <- 0
        probe_z_scores[!is.finite(probe_z_scores)] <- 0
        
        # replace non-positive values with a pseudo-count
        weights[weights<=0] <- 0.000001
        probe_z_scores[probe_z_scores<=0] <- 0.000001
        
        # modify the z-scores if the seed value is below threshold
        for (k in 1:length(probe_z_scores)) {
          
          # if the variant percentile is greater than 0.95 and the seed percentile is less than 0.95
          if (probe_z_scores[k] > 1.645 & weights[k] <= 1.645) {
            
            # nullify effects of the variant by setting it as the seed value
            probe_z_scores[k] <- weights[k]
            
          }
        }
        
        # assign the SNV value as a function of the PBM values for the genomic tiles containing the given SNV position
        SNV[i,j] <- weighted.mean(probe_z_scores, weights)
        
      }
    }
  }
  
  # transform fluo value matrix prior to returning to parent function to plot
  SNV_trans <- matrix(nrow=4, ncol=ncol(SNV))
  rownames(SNV_trans) <- c("A", "C", "G", "T")
  for (i in 1:nrow(SNV_trans)) {
    for (j in 1:ncol(SNV_trans)) {
      # subtract the positional column-wise median if looking at condition-specific z-scores relative to the background set of probes
      SNV_trans[i, j] <- SNV[i, j] - median(SNV[,j])
    }
  }
  
  # convert the SNV matrix to PWM format and perform PWM similarity for different beta values
  if (args[7]=="TRUE") {
    SNV_to_PWM_sim(SNV, SNV_trans, 2.0)
  }
  
  # reverse complement the pwm if necessary
  if (dir == "REVERSE") {
    SNV_trans <- rev_compl_pwm(SNV_trans)
  }
  
  # return the matrix to be plotted as a sequence logo
  return(SNV_trans)
}

# function to obtain z-score data frame for each position and SNV in the CASCADE locus
locus_zscore_df <- function(PBM, PBM_col, genomic_window) {
  
  # initialize nucleotide vector needed to iterate through
  nucs <- c("A", "C", "G", "T")
  
  # initialize vectors needed to be included in the z-score data frame
  z_genomic_pos <- vector(mode="integer")
  z_nuc <- vector(mode="character")
  z_SNV_pos_offset <- vector(mode="integer")
  z_zscore <- vector(mode="numeric")
  z_probe_type <- vector(mode="character")
  z_order <- vector(mode="integer")
  
  # collect z-score fluo values for each position in the genomic sequence
  for (i in 1:length(genomic_window)) {
    for (j in 1:length(nucs)) {
      
      # determine the current genomic position
      genomic_pos <- genomic_window[i]
      
      # collect data relevant to the SNV probes at the current genomic position
      SNV_probes <- PBM[which(PBM$SNV_pos==genomic_pos),]
      
      # determine the current genomic nucleotide
      curr_seed_nuc <- unique(as.character(SNV_probes$seed_nuc))
      
      # if the current SNV nucleotide is the genomic nucleotide
      if (curr_seed_nuc == nucs[j]) {
        
        # collect tiles that contain the given SNV position
        tiles <- unique(SNV_probes$tile_order)
        
        # collect z-scores associated with the tiles that contain the current genomic position
        curr_zscores <- PBM[which(PBM$tile_order %in% tiles & PBM$SNV_pos_offset==0), PBM_col]
        
        # update the vectors with the current values
        n <- length(curr_zscores)
        z_genomic_pos <- c(z_genomic_pos, rep(genomic_pos, n))
        z_nuc <- c(z_nuc, rep(nucs[j], n))
        z_SNV_pos_offset <- c(z_SNV_pos_offset, rep(0, n))
        z_zscore <- c(z_zscore, curr_zscores)
        z_probe_type <- c(z_probe_type, rep("genomic", n))
        z_order <- c(z_order, rep(i, n))
        
      } else { # if the current nucleotide corresponds to a SNV of the genomic seed nucleotide
        
        # collect z-scores from the SNV probes that have the given nucleotide as a SNV
        curr_zscores <- SNV_probes[which(SNV_probes$SNV_nuc==nucs[j]), PBM_col]
        
        # collect the SNV pos offsets from the SNV probes that have the given nucleotide as a SNV
        curr_SNV_pos_offsets <- SNV_probes[which(SNV_probes$SNV_nuc==nucs[j]), "SNV_pos_offset"]
        
        # update the vectors with the current values
        n <- length(curr_zscores)
        z_genomic_pos <- c(z_genomic_pos, rep(genomic_pos, n))
        z_nuc <- c(z_nuc, rep(nucs[j], n))
        z_SNV_pos_offset <- c(z_SNV_pos_offset, curr_SNV_pos_offsets)
        z_zscore <- c(z_zscore, curr_zscores)
        z_probe_type <- c(z_probe_type, rep("SNV", n))
        z_order <- c(z_order, rep(i, n))
        
      }
    }
  }
  
  # construct the zscore data frame using the vectors produced
  zscore_df <- data.frame(z_genomic_pos, z_nuc, z_SNV_pos_offset, z_zscore, z_probe_type, z_order)
  
  # return the matrix to be plotted as a sequence logo
  return(zscore_df)
}

# declare the different PBM results conditions (in order, via args[5])
PBM_cond <- read.table(args[5], sep="\t")
PBM_cond <- as.character(PBM_cond$V1)

# declare names of seed regions of interest
to_keep <- c("TILE")
x <- droplevels(df[which(df$probe_type %in% to_keep),])
seeds <- unique(as.character(x$seed_names))
rm(x)

# read in the known TF binding site annotation
TF_sites <- read.table(args[6], header=T, sep='\t')
TF_sites$hex_color <- paste("#", TF_sites$hex_color, sep="")

# intialize a list of human TF motifs (MotifDb)
motif_db <- MotifDb[grep("Hsapiens", values(MotifDb)$organism)]

# filter the database for motifs greater than or equal to 6 bases in length
motif_db <- motif_db[lapply(motif_db, ncol) >= 6]

# create a directory to house the locus plots
dir.create(paste(getwd(), args[1], "locus_plots", sep="/"))

# create a directory to house the locus SNV matrices
dir.create(paste(getwd(), args[1], "locus_SNV_matrices", sep="/"))

# if user wants to run the motif similarity module
if (args[7]=="TRUE") {
  
  # create a directory to house the motif similarity results
  dir.create(paste(getwd(), args[1], "motif_similarity", sep="/"))
  message("and performing PWM similarity suite...")
}

# for each PBM condition in the current series of experiments
for (cond in PBM_cond) {
  
  # modify the cond to contain the run prefix and the "br" orientation tag
  cond <- paste(args[1], "br", cond, sep="_")
  
  # start next iteration if greater than 10% of entries are repeated values (low quality experiment)
  if (sum(duplicated(df[,cond]))>(0.10*nrow(df))) {
    
    # print a warning message to console
    message(paste(cond, " has been excluded from subsequent analysis due to low quality", sep=""))
    
    # proceed to next cond experiment column
    next
  }
  
  # for each seed sequence
  for (seed in seeds) {
    
    # keep only data that matches the seed
    x <- droplevels(df[which(df$seed_names==seed),])
    
    # grab the actual genomic window (nucleotide positions within a given chromosome)
    genomic_window <- min(x$tile_starts):max(x$tile_ends)
    
    # construct the genomic sequence of the given locus
    genomic_seq <- unique(x[,c("SNV_pos", "seed_nuc")]) # collect the genomic nucleotides for each genomic position
    genomic_seq <- genomic_seq[order(genomic_seq$SNV_pos),] # ensure SNV positions are sorted
    genomic_seq <- as.character(genomic_seq[which(!is.na(genomic_seq$seed_nuc)),"seed_nuc"]) # get rid of the NA from the tile probe
    genomic_seq <- paste(genomic_seq, collapse="") # collapse the character vector into a single character string
    
    # restrict TF binding site annotation to given seed
    seed_TF_sites <- droplevels(TF_sites[which(TF_sites$seed_names==seed),])
    
    # declare the possible strand directions (FORWARD is 5' to 3' on the + strand, REVERSE is 5' to 3' on the - strand)
    directions <- c("FORWARD",
                    "REVERSE")
    
    # for both possible directions (FORWARD is 5' to 3' on the + strand, REVERSE is 5' to 3' on the - strand)
    for (dir in directions) {
      
      # reverse complement the current genomic sequence and reverse the window if REVERSE
      curr_genomic_seq <- genomic_seq
      curr_genomic_window <- genomic_window
      if (dir == "REVERSE") {
        curr_genomic_seq <- rev_compl(curr_genomic_seq)
        curr_genomic_window <- rev(genomic_window)
      }
      
      # scale transformation function
      scaleFUN <- function(x) sprintf("%.1f", x)
      
      
      
      ### GENOMIC LADDER
      
      # plot genomic ladder
      seq_info <- data.frame(x = 1:length(genomic_window),
                             y = 1,
                             nuc = strsplit(curr_genomic_seq, "")[[1]])
      seq_info$nuc <- as.character(seq_info$nuc)
      ggp_1 <- ggplot(seq_info, aes(x, y, label=nuc)) +
        geom_text() +
        labs(title = paste("hg38: ", as.character(unique(x$tile_chroms, sep="")))) +
        theme(plot.title = element_text(lineheight=0.8, face="bold", hjust=0)) +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
        theme(axis.line=element_blank()) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
        # extra space in limits ensures perfect alignment with motif plotted in ggp_2 below
        scale_x_continuous(limits=c(0.5, length(genomic_window)+0.5), expand = c(0.01,0.01), position = "top", breaks = which(curr_genomic_window %% 10 == 0), labels = curr_genomic_window[which(curr_genomic_window %% 10 == 0)]) +
        theme(axis.text.x=element_text(hjust=0.15)) +
        scale_y_continuous(limits=c(-2.5, 2.5), expand = c(0.01,0.01))
      
      # plot known TF sites onto the genomic ladder
      for (i in 1:nrow(seed_TF_sites)) {
        min_x <- which(curr_genomic_window==(seed_TF_sites[i,"start"]-1))-0.5 # determines the start of the current TF binding site window
        max_x <- which(curr_genomic_window==(seed_TF_sites[i,"end"]+1))+0.5 # determines the end of the current TF binding site window
        
        # change the min_x and max_x if looking at the FORWARD logo
        if (dir == "FORWARD") {
          min_x <- which(curr_genomic_window==(seed_TF_sites[i,"start"]))-0.5 # determines the start of the current TF binding site window
          max_x <- which(curr_genomic_window==(seed_TF_sites[i,"end"]))+0.5 # determines the end of the current TF binding site window
        }
        
        # annotate the known TF site by highlighting it with a semi-transparent box
        ggp_1 <- ggp_1 + annotate('rect', xmin = min_x, xmax = max_x, ymin=0, ymax=2, fill=seed_TF_sites[i,"hex_color"], alpha=0.35)
        # label the current TF site using the name listed in the file (centered beneath the locus)
        ggp_1 <- ggp_1 + annotate('text', x = (min_x + max_x)/2, y = -1.5, label = seed_TF_sites[i,"TF_names"])
      }
      
      
      
      ### SNV ENERGY LOGO
      
      # build the given SNV energy matrix
      pwm <-  suppressMessages(locus_SNV_matrix(x, cond, genomic_window))
      
      # replace flat 0s (very rare but possible) with a pseudocount to prevent ggseqlogo from failing while trying to plot
      pwm[pwm==0] <- 0.000001
      
      # save the SNV energy matrix to file
      mat_name <- paste(seed, cond, dir, "SNV_matrix.txt", sep="_")
      write.table(pwm, paste(getwd(), args[1], "locus_SNV_matrices", mat_name, sep="/"), quote=F, row.names=T, col.names=F, sep='\t')
      
      # assign the y axis label based on the unit being plotted
      y_axis_label <- suppressMessages(expression(paste(Delta, "z-score", sep="")))
      
      # plot the SNV energy logo
      ggp_2 <- ggseqlogo(pwm, method="custom")
      my.ggp.yrange <- ggplot_build(ggp_2)$layout$panel_params[[1]]$y.range
      ggp_2 <- ggp_2 +
        scale_x_continuous(expand = c(0.01,0.01)) +
        scale_y_continuous(expand = c(0.01,0.01), labels=scaleFUN) +
        annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = my.ggp.yrange[1], ymax = 0.0, alpha = 0.75, col='white', fill='white') +
        annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=1.1) +
        annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = my.ggp.yrange[1], ymax = 0, col="lightgray", fill=NA, size=1.1) +
        ylab(y_axis_label) +
        theme(axis.text.x=element_blank()) +
        theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain"))
      
      
      
      ### RAW Z-SCORE TRACK
      
      # build the raw z-score data frame
      zscore_df <- locus_zscore_df(x, cond, genomic_window)
      
      # reverse complement the data frame if REVERSE
      if (dir == "REVERSE") {
        # determine the max z_order value
        max_z <- max(zscore_df$z_order)
        
        # re-map the z_order (equivalent to reversing it)
        zscore_df$z_order <- -zscore_df$z_order + max_z + 1
        
        # reverse complement the nucleotide
        zscore_df$z_nuc <- chartr("ACGTN", "TGCAN", zscore_df$z_nuc)
      }
      
      # assign the y axis label based on the unit being plotted
      y_axis_label <- "z-score"
      
      # assign different colors for different nucleotides
      nuc_colors <- c("#109648", # A
                      "#255c99", # C
                      "#f7b32b", # G
                      "#d62839") # T
      
      # use ggplot to plot a scatter of the z-score vs. the position
      ggp_3 <- ggplot(zscore_df, aes(x=z_order, y=z_zscore, shape=z_probe_type, color=z_nuc)) + geom_point()
      my.ggp.yrange <- ggplot_build(ggp_3)$layout$panel_params[[1]]$y.range
      ggp_3 <- ggp_3 +
        scale_color_manual(values=nuc_colors) +
        scale_x_continuous(limits=c(0.5, length(genomic_window)+0.5), expand = c(0.01,0.01), name="") +
        scale_y_continuous(expand = c(0.01, 0.01), labels=scaleFUN) +
        annotate('rect', xmin = 0.5, xmax = max(zscore_df$z_order)+0.5, ymin = my.ggp.yrange[1], ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=1.1) +
        theme(axis.line=element_blank()) +
        ylab(y_axis_label) +
        theme(legend.position="none") +
        theme(axis.text.x=element_blank(), axis.ticks=element_blank()) +
        theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain"))
      
      
      
      ### USE COWPLOT TO SAVE THE PLOTS
      pdf_name <- paste(seed, cond, dir, "locus", "plots.pdf", sep="_")
      ggp <- plot_grid(ggp_1, ggp_2, ggp_3,
                       ncol = 1, nrow = 3, rel_heights = c(1, 1.75, 1.75), align = "v")
      save_plot(paste(getwd(), args[1], "locus_plots", pdf_name, sep="/"), ggp, ncol = 1, nrow = 3, base_height = 2, base_width = length(genomic_window)/7, base_aspect_ratio = 1, limitsize=FALSE)
      
    }
  }
}



### PLOT ALL OF THE TILE SEED LOGOS ----------------------------------

# print a status message
message("plotting tile seed logos...")

# create a directory to house the orientation scatterplots
dir.create(paste(getwd(), args[1], "tile_logos", sep="/"))

# function to obtain SNV value matrix
SNV_to_energy <- function(PBM, PBM_col, target_seq_col, target_seq) {
  
  # initialize a matrix to hold the SNV probe z-score values
  nucs <- c("A", "C", "G", "T")
  SNV <- matrix(nrow=4, ncol=nchar(as.character(target_seq)))
  
  # collect z-score values for each position in the matrix
  for (i in 1:nrow(SNV)) {
    for (j in 1:ncol(SNV)) {
      # construct the current sequence of interest using the current SNV
      SNV_seq <- paste(substr(target_seq, 1, j-1), nucs[i], substr(target_seq, j+1, ncol(SNV)), sep="")
      
      # collect the z-score for that given SNV probe sequence and condition
      SNV[i, j] <- PBM[which(as.character(PBM[, target_seq_col])==SNV_seq), PBM_col]
    }
  }
  
  # transform z-score values to reflect deviation from column-wise median
  SNV_trans <- matrix(nrow=4, ncol=nchar(as.character(target_seq)))
  rownames(SNV_trans) <- c("A", "C", "G", "T")
  for (i in 1:nrow(SNV_trans)) {
    for (j in 1:ncol(SNV_trans)) {
      SNV_trans[i, j] <- SNV[i, j] - median(SNV[,j])
    }
  }
  
  # return the PWM matrix to be plotted as a sequence logo
  return(SNV_trans)
}

# declare names of seed regions of interest
to_keep <- c("STATIC", "TILE")
x <- droplevels(df[which(df$probe_type %in% to_keep),])
seeds <- unique(as.character(x$seed_names))
rm(x)

# for each seed locus of interest
for (seed in seeds) {
  
  # keep only data that matches the current seed
  x <- df[which(df$seed_names==seed),]
  
  # collect tile order for the given seed
  tiles <- unique(x$tile_order)
  
  # for each PBM test condition
  for (cond in PBM_cond) {
    
    # modify cond to contain run name and "br" tag
    cond <- paste(args[1], "br", cond, sep="_")
    
    # start next iteration if greater than 10% of entries are repeated values (low quality experiment)
    if (sum(duplicated(df[,cond]))>(0.10*nrow(df))) {
      
      # print a warning message to console
      message(paste(cond, " has been excluded from subsequent analysis due to low quality", sep=""))
      
      # proceed to next cond experiment column
      next
    }
    
    # for both possible directions (FORWARD is 5' to 3' on the + strand, REVERSE is 5' to 3' on the - strand)
    for (dir in c("FORWARD", "REVERSE")) {
      
      # open pdf file
      pdf_name <- paste(seed, cond, dir, "tile_logos.pdf", sep="_")
      pdf(file=paste(getwd(), args[1], "tile_logos", pdf_name, sep="/"), width=0.35*26, height=2.5*length(tiles))
      
      # initialize a general pwm to have dimensions on hand
      pwm <- matrix(nrow=4, ncol=26)
      rownames(pwm) <- c("A", "C", "G", "T")
      
      # plot each SNV model energy logo
      p_list <- lapply(tiles, function(f) {
        
        # determine the reference seq for the given seed and tile
        ref_seed <- as.character(x[which(x$tile_order==f & x$SNV_pos_offset==0),"target_seq"])
        
        # scale transformation function
        scaleFUN <- function(x) sprintf("%.1f", x)
        
        # build the SNV matrix using the current data
        pwm <- SNV_to_energy(x, cond, "target_seq", ref_seed)
        
        # reverse complement the pwm and reference seed sequence if necessary
        if (dir == "REVERSE") {
          pwm <- rev_compl_pwm(pwm)
          ref_seed <- rev_compl(ref_seed)
        }
        
        # assign the y axis label based on the unit being plotted
        y_axis_label <- expression(paste(Delta, "z-score", sep=""))
        
        # plot the given SNV energy logo
        ggp <- ggseqlogo(pwm, method="custom") + ggtitle(paste(seed, f, dir, ref_seed, as.character(round(x[which(x$tile_order==f & x$SNV_pos_offset==0), cond], 3)), sep="    "))
        my.ggp.yrange <- ggplot_build(ggp)$layout$panel_params[[1]]$y.range
        ggp <- ggp +
          annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = my.ggp.yrange[1], ymax = 0.0, alpha = 0.75, col='white', fill='white') +
          scale_y_continuous(labels=scaleFUN)
        ggp <- ggp + 
          if (x[which(x$tile_order==f & x$SNV_pos_offset==0), cond] <= 1.5) {
            annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = my.ggp.yrange[1], ymax = my.ggp.yrange[2], alpha = 0.75, col='white', fill='white')
          }
        ggp <- ggp +
          annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=1.1) + 
          annotate('rect', xmin = 0.5, xmax = ncol(pwm)+0.5, ymin = my.ggp.yrange[1], ymax = 0, col="lightgray", fill=NA, size=1.1) +
          ylab(y_axis_label) +
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.title.y=element_text(size=17)) +
          theme(axis.text.y=element_text(size=14))
      })
      
      # arrange a grid of motifs
      do.call(gridExtra::grid.arrange, c(p_list, ncol=1))
      
      # close the pdf file
      dev.off()
    }
  }
}

# save the sessionInfo for reproducibility purposes
writeLines(capture.output(sessionInfo()), paste(getwd(), "/", args[1], "/RsessionInfo_", args[1],".txt", sep=""))

# print final message
message("CASCADE ANALYSIS COMPLETE")



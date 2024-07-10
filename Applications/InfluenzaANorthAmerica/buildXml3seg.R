library(stringr)
library(seqinr)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Make a directory to store the xml files
if (dir.exists("xmls")) {
  unlink("xmls", recursive = TRUE)
}
dir.create("xmls")

# build a wgs xml file for each virus and year between 2015 and 2018
virus = c('H1N1', 'H3N2')
year = seq(2011,2018,1)

# get all fasta files in the data directory that contain one of the year in the 
# name
fasta = list.files(path = "./data", pattern = "fasta", full.names = TRUE)
fasta = fasta[grepl(paste(year, collapse="|"), fasta)]
file.copy(fasta, "xmls")


# define the root height for Ne and reassortment variant rates
timediff = 5
rateshiftvals = c(seq(0, 1, length.out=10), seq(2, timediff, length.out=4))
# loop over all combinations of files
for (a in 1:length(virus)) {
  for (b in 1:length(year)){
    base = paste(virus[[a]], "_", year[[b]], sep="")
    
    fasta1 = read.fasta(file = paste("./data/", base, "_HA.fasta", sep=""), seqtype = "DNA")
    isolates = names(fasta1)
    # require at least 100 isolates
    if (length(isolates) < 50) {
      next
    }
    
    # read in the log file to get the network height
    filename = paste(virus[[a]], "_", year[[b]] , ".constant", sep="")
    # build an inference xml files
    f <- file(sprintf('xmls/%s.xml',filename), 'w')
    # Open the template file
    template <- file('../H5N1NorthAmerica/inference_template_3seg.xml', 'r')
    
    
  
    while (length(line <- readLines(template, n = 1)) > 0) {
      if (grepl('insert_name', line)) {
        writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
      } else if (grepl('insert_tips', line)) {
        # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
        for (i in seq(1, length(isolates))) {
          writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
        }
      }else if (grepl('insert_clock_rate', line)){
        if (virus[[a]] == 'H1N1') {
          writeLines(gsub('insert_clock_rate', '0.0025', line), f)
        } else {
          writeLines(gsub('insert_clock_rate', '0.0021', line), f)
        }

      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', paste(base, "HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(base, "NS", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(base, "PB1", sep="_"), line), f)
      }else if (grepl('insert_seg4', line)){
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
  
      # mask all the non isConstant entries, i.e. isinfectedSkyline
      } else if (grepl('isinfectedSkyline-->', line)) {
        writeLines(gsub('isinfectedSkyline-->', '', line), f)
      } else if (grepl('<!--isinfectedSkyline', line)){
        writeLines(gsub('<!--isinfectedSkyline', '', line), f)      
      } else if (grepl('isNeSkyline-->', line)) {
        writeLines(gsub('isNeSkyline-->', '', line), f)
      } else if (grepl('<!--isNeSkyline', line)){
        writeLines(gsub('<!--isNeSkyline', '', line), f)  
      } else if (grepl('isISkyline-->', line)) {
        writeLines(gsub('isISkyline-->', '', line), f)
      } else if (grepl('<!--isISkyline', line)){
        writeLines(gsub('<!--isISkyline', '', line), f)  
        
      } else if (grepl('isVariable-->', line)) {
        writeLines(gsub('isVariable-->', '', line), f)
      } else if (grepl('<!--isVariable', line)){
        writeLines(gsub('<!--isVariable', '', line), f)  
        
      } else if (grepl('insert_heights', line)) {
        # write the date of each isolate
        value=''
        for (i in seq(1, length(isolates))) {
          date = strsplit(isolates[i], "\\|")[[1]][[3]]
          
          if (i == length(isolates)) {
            value = paste(value, isolates[i], '=', date, sep = '')
          } else {
            value = paste(value, isolates[i], '=', date, sep = '')
            value = paste(value, ',', sep = '')
          }
        }  
        writeLines(gsub('insert_heights', value, line), f)   
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        seq1 = read.fasta(file = paste("./data/", base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = read.fasta(file = paste("./data/", base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq3 = read.fasta(file = paste("./data/", base, "_PB1.fasta", sep=""), seqtype = "DNA")

        #get the number of characters in the first segment
        sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
        sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
        sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
        writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
                                                sequence_length_seg3), line), f)
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)

    # make a second xml where the .trees is replaced by .infected
    filename = paste(virus[[a]], "_", year[[b]] , ".variable", sep="")
    f <- file(sprintf('xmls/%s.xml', filename), 'w')
    # Open the template file
    template <- file('../H5N1NorthAmerica/inference_template_3seg.xml', 'r')
    while (length(line <- readLines(template, n = 1)) > 0) {
      if (grepl('insert_name', line)) {
        writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
      } else if (grepl('insert_tips', line)) {
        # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
        for (i in seq(1, length(isolates))) {
          writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
        }
      }else if (grepl('insert_clock_rate', line)){
        if (virus[[a]] == 'H1N1') {
          writeLines(gsub('insert_clock_rate', '0.0025', line), f)
        } else {
          writeLines(gsub('insert_clock_rate', '0.0021', line), f)
        }
        
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', paste(base, "HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(base, "NS", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(base, "PB1", sep="_"), line), f)
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
        
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
        
      # mask all the non isConstant entries, i.e. isinfectedSkyline
      } else if (grepl('isconstant-->', line)) {
        writeLines(gsub('isconstant-->', '', line), f)
      } else if (grepl('<!--isconstant', line)){
        writeLines(gsub('<!--isconstant', '', line), f)
      } else if (grepl('isNeSkyline-->', line)) {
        writeLines(gsub('isNeSkyline-->', '', line), f)
      } else if (grepl('<!--isNeSkyline', line)){
        writeLines(gsub('<!--isNeSkyline', '', line), f)
        
      } else if (grepl('isISkyline-->', line)) {
        writeLines(gsub('isISkyline-->', '', line), f)
      } else if (grepl('<!--isISkyline', line)){
        writeLines(gsub('<!--isISkyline', '', line), f)
        
      } else if (grepl('isinfectedSkyline-->', line)) {
        writeLines(gsub('isinfectedSkyline-->', '', line), f)
      } else if (grepl('<!--isinfectedSkyline', line)){
        writeLines(gsub('<!--isinfectedSkyline', '', line), f)
        
        
      } else if (grepl('insert_heights', line)) {
        # write the height of each tip using heights seperated by a = and , between
        value=''
        for (i in seq(1, length(isolates))) {
          date = strsplit(isolates[i], "\\|")[[1]][[3]]
          
          if (i == length(isolates)) {
            value = paste(value, isolates[i], '=', date, sep = '')
          } else {
            value = paste(value, isolates[i], '=', date, sep = '')
            value = paste(value, ',', sep = '')
          }
        }
        writeLines(gsub('insert_heights', value, line), f)
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        # get the length of both alignments
        seq1 = read.fasta(file = paste("./data/", base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = read.fasta(file = paste("./data/", base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq3 = read.fasta(file = paste("./data/", base, "_PB1.fasta", sep=""), seqtype = "DNA")

        #get the number of characters in the first segment
        sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
        sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
        sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
        writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
                                                sequence_length_seg3), line), f)      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)

    # make a second xml where the .trees is replaced by .infected
    filename = paste(virus[[a]], "_", year[[b]] , ".ne", sep="")
    f <- file(sprintf('xmls/%s.xml', filename), 'w')
    # Open the template file
    template <- file('../H5N1NorthAmerica/inference_template_3seg.xml', 'r')
    while (length(line <- readLines(template, n = 1)) > 0) {
      if (grepl('insert_name', line)) {
        writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
      } else if (grepl('insert_tips', line)) {
        # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
        for (i in seq(1, length(isolates))) {
          writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
        }
      }else if (grepl('insert_clock_rate', line)){
        if (virus[[a]] == 'H1N1') {
          writeLines(gsub('insert_clock_rate', '0.0025', line), f)
        } else {
          writeLines(gsub('insert_clock_rate', '0.0021', line), f)
        }
        
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', paste(base, "HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(base, "NS", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(base, "PB1", sep="_"), line), f)
        
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
        
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
        # mask all the non isConstant entries, i.e. isinfectedSkyline
      } else if (grepl('isconstant-->', line)) {
        writeLines(gsub('isconstant-->', '', line), f)
      } else if (grepl('<!--isconstant', line)){
        writeLines(gsub('<!--isconstant', '', line), f)      
      } else if (grepl('isISkyline-->', line)) {
        writeLines(gsub('isISkyline-->', '', line), f)
      } else if (grepl('<!--isISkyline', line)){
        writeLines(gsub('<!--isISkyline', '', line), f)  
        
      } else if (grepl('isVariable-->', line)) {
        writeLines(gsub('isVariable-->', '', line), f)
      } else if (grepl('<!--isVariable', line)){
        writeLines(gsub('<!--isVariable', '', line), f)  
        
        
      } else if (grepl('insert_heights', line)) {
        # write the height of each tip using heights seperated by a = and , between
        value=''
        for (i in seq(1, length(isolates))) {
          date = strsplit(isolates[i], "\\|")[[1]][[3]]
          
          if (i == length(isolates)) {
            value = paste(value, isolates[i], '=', date, sep = '')
          } else {
            value = paste(value, isolates[i], '=', date, sep = '')
            value = paste(value, ',', sep = '')
          }
        }  
        writeLines(gsub('insert_heights', value, line), f)  
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        # get the length of both alignments
        seq1 = read.fasta(file = paste("./data/", base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = read.fasta(file = paste("./data/", base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq3 = read.fasta(file = paste("./data/", base, "_PB1.fasta", sep=""), seqtype = "DNA")

        #get the number of characters in the first segment
        sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
        sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
        sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
        writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
                                                sequence_length_seg3), line), f)      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)
  }
}



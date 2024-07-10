library(stringr)
library(seqinr)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get all the aligned fasta files (the fasta files of the form HA.fasta, PB1.fasta)
# from the data directory
fasta <- list.files("./data/", pattern="*.fasta", full.names=TRUE)
# remove all files that contain the word aligned
fasta <- fasta[!grepl("aligned", fasta)]
# remove the H5N1.... fasta file as well
fasta <- fasta[!grepl("H5N1", fasta)]

# read in the first fasta file and keep all the names
fasta1 = read.fasta(file = fasta[1], seqtype = "DNA")
isolates = names(fasta1)
# define the root height for Ne and reassortment variant rates
timediff = 5
# copy all the files in fasta to the xmls directory
file.copy(fasta, "xmls")
# loop over all combinations of files
for (a in 1:(length(fasta)-1)) {
  for (b in (a+1):length(fasta)){
    seg1 = gsub(".fasta", "", basename(fasta[a]))
    seg2 = gsub(".fasta", "", basename(fasta[b]))
    
    # read in the log file to get the network height
    filename = paste("H5N1_", seg1, "_", seg2, ".constant", sep="")
    # build an inference xml files
    f <- file(sprintf('xmls/%s.xml',filename), 'w')
    # Open the template file
    template <- file('inference_template.xml', 'r')
  
    while (length(line <- readLines(template, n = 1)) > 0) {
      if (grepl('insert_name', line)) {
        writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
      } else if (grepl('insert_tips', line)) {
        # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
        for (i in seq(1, length(isolates))) {
          writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
        }
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', seg1, line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', seg2, line), f)
        
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(seq(0,timediff, length.out=15), collapse=' '), line), f)
  
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
        dates = c()
        for (i in seq(1, length(isolates))) {
          date = strsplit(isolates[i], "\\|")[[1]][[3]]
          dates = c(dates, as.Date(date))
          if (i == length(isolates)) {
            value = paste(value, isolates[i], '=', date, sep = '')
          } else {
            value = paste(value, isolates[i], '=', date, sep = '')
            value = paste(value, ',', sep = '')
          }
        }  
        writeLines(gsub('insert_heights', value, line), f)   
      }else if (grepl('insert_EOS', line)){
        # get the most recent sampling time
        EOS = max(as.numeric(strsplit(isolates[length(isolates)], "\\|")[[1]][[3]]))
        # compute EOS minus the august first of the previous year
        EOS = max(dates) - min(dates)
        writeLines(gsub('insert_EOS', EOS/365, line), f)
        
        
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        seq1 = read.fasta(file = fasta[a], seqtype = "DNA")
        seq2 = read.fasta(file = fasta[b], seqtype = "DNA")
        #get the number of characters in the first segment
        sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
        sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
        writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2), line), f)
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)
    
    
    
    # # make a second xml where the .trees is replaced by .infected
    # filename = paste("H5N1_", seg1, "_", seg2, ".variable", sep="")
    # f <- file(sprintf('xmls/%s.xml', filename), 'w')
    # # Open the template file
    # template <- file('inference_template.xml', 'r')
    # while (length(line <- readLines(template, n = 1)) > 0) {
    #   if (grepl('insert_name', line)) {
    #     writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
    #   } else if (grepl('insert_tips', line)) {
    #     # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
    #     for (i in seq(1, length(isolates))) {
    #       writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
    #     }
    #   }else if (grepl('insert_seg1', line)){
    #     writeLines(gsub('insert_seg1', seg1, line), f)
    #   }else if (grepl('insert_seg2', line)){
    #     writeLines(gsub('insert_seg2', seg2, line), f)
    #     
    #   } else if (grepl('insert_times', line)) {
    #     writeLines(gsub('insert_times', paste(seq(0,timediff, length.out=15), collapse=' '), line), f)
    #     
    #   # mask all the non isConstant entries, i.e. isinfectedSkyline
    #   } else if (grepl('isconstant-->', line)) {
    #     writeLines(gsub('isconstant-->', '', line), f)
    #   } else if (grepl('<!--isconstant', line)){
    #     writeLines(gsub('<!--isconstant', '', line), f)
    #   } else if (grepl('isNeSkyline-->', line)) {
    #     writeLines(gsub('isNeSkyline-->', '', line), f)
    #   } else if (grepl('<!--isNeSkyline', line)){
    #     writeLines(gsub('<!--isNeSkyline', '', line), f)
    #     
    #   } else if (grepl('isISkyline-->', line)) {
    #     writeLines(gsub('isISkyline-->', '', line), f)
    #   } else if (grepl('<!--isISkyline', line)){
    #     writeLines(gsub('<!--isISkyline', '', line), f)
    #     
    #   } else if (grepl('isinfectedSkyline-->', line)) {
    #     writeLines(gsub('isinfectedSkyline-->', '', line), f)
    #   } else if (grepl('<!--isinfectedSkyline', line)){
    #     writeLines(gsub('<!--isinfectedSkyline', '', line), f)
    #     
    #     
    #   } else if (grepl('insert_heights', line)) {
    #     # write the height of each tip using heights seperated by a = and , between
    #     value=''
    #     dates = c()
    #     for (i in seq(1, length(isolates))) {
    #       date = strsplit(isolates[i], "\\|")[[1]][[3]]
    #       dates = c(dates, as.numeric(date))
    #       if (i == length(isolates)) {
    #         value = paste(value, isolates[i], '=', date, sep = '')
    #       } else {
    #         value = paste(value, isolates[i], '=', date, sep = '')
    #         value = paste(value, ',', sep = '')
    #       }
    #     }
    #     writeLines(gsub('insert_heights', value, line), f)
    #   }else if (grepl('insert_EOS', line)){
    #     # get the most recent sampling time
    #     EOS = max(as.numeric(strsplit(isolates[length(isolates)], "\\|")[[1]][[3]]))
    #     # compute EOS minus the august first of the previous year
    #     EOS = max(dates) - min(dates)
    #     writeLines(gsub('insert_EOS', EOS/365, line), f)
    #     
    #   }else if (grepl('insert_weights', line)) {
    #     # get the length of both alignments
    #     seq1 = read.fasta(file = fasta[a], seqtype = "DNA")
    #     seq2 = read.fasta(file = fasta[b], seqtype = "DNA")
    #     #get the number of characters in the first segment
    #     sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
    #     sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
    #     writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2), line), f)
    #   } else {
    #     writeLines(line, f)
    #   }
    # }
    # close(f)
    # close(template)

    # make a second xml where the .trees is replaced by .infected
    filename = paste("H5N1_", seg1, "_", seg2, ".ne", sep="")
    f <- file(sprintf('xmls/%s.xml', filename), 'w')
    # Open the template file
    template <- file('inference_template.xml', 'r')
    while (length(line <- readLines(template, n = 1)) > 0) {
      if (grepl('insert_name', line)) {
        writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
      } else if (grepl('insert_tips', line)) {
        # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
        for (i in seq(1, length(isolates))) {
          writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
        }
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', seg1, line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', seg2, line), f)
        
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(seq(0,timediff, length.out=15), collapse=' '), line), f)
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
        dates = c()
        for (i in seq(1, length(isolates))) {
          date = strsplit(isolates[i], "\\|")[[1]][[3]]
          dates = c(dates, as.Date(date))
          if (i == length(isolates)) {
            value = paste(value, isolates[i], '=', date, sep = '')
          } else {
            value = paste(value, isolates[i], '=', date, sep = '')
            value = paste(value, ',', sep = '')
          }
        }  
        writeLines(gsub('insert_heights', value, line), f)  
      }else if (grepl('insert_EOS', line)){
        # get the most recent sampling time
        EOS = max(as.numeric(strsplit(isolates[length(isolates)], "\\|")[[1]][[3]]))
        # compute EOS minus the august first of the previous year
        EOS = max(dates) - min(dates)
        writeLines(gsub('insert_EOS', EOS/365, line), f)
        
        
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        seq1 = read.fasta(file = fasta[a], seqtype = "DNA")
        seq2 = read.fasta(file = fasta[b], seqtype = "DNA")
        #get the number of characters in the first segment
        sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
        sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
        writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2), line), f)
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)
  }
}



library(stringr)
library(seqinr)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# delete all *wgs*.xml files in xmls
file.remove(list.files("./xmls/", pattern="*.xml", full.names=TRUE))

# define the root height for Ne and reassortment variant rates
first_intro = 5

rateshiftvals = c(seq(0, first_intro, length.out=15), seq(first_intro, 10, length.out=6))
rateshiftvals = unique(rateshiftvals)
rateshiftvals2 = rateshiftvals


for (timeaveraged in c('averaged', 'proportional')){
  # read in the first fasta file and keep all the names
  fasta1 = read.fasta(file = paste("./xmls/",timeaveraged, "_HA.fasta", sep=""), seqtype = "DNA")
  isolates = names(fasta1)
  
  # # read in the log file to get the network height
  # filename = paste("H5N1_wgs_", timeaveraged , ".constant", sep="")
  # # build an inference xml files
  # f <- file(sprintf('xmls/%s.rep0.xml',filename), 'w')
  # # Open the template file
  # template <- file('../H5N1NorthAmerica/inference_template_wgs_cr.xml', 'r')
  # 
  # base=""
  # 
  # while (length(line <- readLines(template, n = 1)) > 0) {
  #   if (grepl('insert_name', line)) {
  #     writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
  #   } else if (grepl('insert_tips', line)) {
  #     # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
  #     for (i in seq(1, length(isolates))) {
  #       writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
  #     }
  #   }else if (grepl('insert_clock_rate', line)){
  #     writeLines(gsub('insert_clock_rate', '0.0035', line), f)
  # 
  #   }else if (grepl('insert_seg1', line)){
  #     writeLines(gsub('insert_seg1', paste(timeaveraged,"HA", sep="_"), line), f)
  #   }else if (grepl('insert_seg2', line)){
  #     writeLines(gsub('insert_seg2', paste(timeaveraged,"NA", sep="_"), line), f)
  #   }else if (grepl('insert_seg3', line)){
  #     writeLines(gsub('insert_seg3', paste(timeaveraged,"MP", sep="_"), line), f)
  #   }else if (grepl('insert_seg4', line)){
  #     writeLines(gsub('insert_seg4', paste(timeaveraged,"NS", sep="_"), line), f)
  #   }else if (grepl('insert_seg5', line)){
  #     writeLines(gsub('insert_seg5', paste(timeaveraged,"NP", sep="_"), line), f)
  #   }else if (grepl('insert_seg6', line)){
  #     writeLines(gsub('insert_seg6', paste(timeaveraged,"PB1", sep="_"), line), f)
  #   }else if (grepl('insert_seg7', line)){
  #     writeLines(gsub('insert_seg7', paste(timeaveraged,"PB2", sep="_"), line), f)
  #   }else if (grepl('insert_seg8', line)){
  #     writeLines(gsub('insert_seg8', paste(timeaveraged,"PA", sep="_"), line), f)
  #   } else if (grepl('insert_times', line)) {
  #     writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
  #   } else if (grepl('insert_ratetimes', line)) {
  #     writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
  #     
  #   # mask all the non isConstant entries, i.e. isinfectedSkyline
  #   } else if (grepl('isinfectedSkyline-->', line)) {
  #     writeLines(gsub('isinfectedSkyline-->', '', line), f)
  #   } else if (grepl('<!--isinfectedSkyline', line)){
  #     writeLines(gsub('<!--isinfectedSkyline', '', line), f)      
  #   } else if (grepl('isNeSkyline-->', line)) {
  #     writeLines(gsub('isNeSkyline-->', '', line), f)
  #   } else if (grepl('<!--isNeSkyline', line)){
  #     writeLines(gsub('<!--isNeSkyline', '', line), f)  
  #   } else if (grepl('isISkyline-->', line)) {
  #     writeLines(gsub('isISkyline-->', '', line), f)
  #   } else if (grepl('<!--isISkyline', line)){
  #     writeLines(gsub('<!--isISkyline', '', line), f)  
  #     
  #   } else if (grepl('isVariable-->', line)) {
  #     writeLines(gsub('isVariable-->', '', line), f)
  #   } else if (grepl('<!--isVariable', line)){
  #     writeLines(gsub('<!--isVariable', '', line), f)  
  #     
  #   } else if (grepl('insert_heights', line)) {
  #     # write the date of each isolate
  #     value=''
  #     dates=c()
  #     for (i in seq(1, length(isolates))) {
  #       date = strsplit(isolates[i], "\\|")[[1]][[3]]
  #       dates=c(dates, as.Date(date))
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
  #     # make 3 points between 0 and EOS/365
  #     writeLines(gsub('insert_EOS', EOS/365, line), f)
  #     
  #   }else if (grepl('insert_weights', line)) {
  #     # get the length of both alignments
  #     seq1 = read.fasta(file = paste("./xmls/",timeaveraged, "_HA.fasta", sep=""), seqtype = "DNA")
  #     seq2 = read.fasta(file = paste("./xmls/",timeaveraged, "_NA.fasta", sep=""), seqtype = "DNA")
  #     seq3 = read.fasta(file = paste("./xmls/",timeaveraged, "_MP.fasta", sep=""), seqtype = "DNA")
  #     seq4 = read.fasta(file = paste("./xmls/",timeaveraged, "_NS.fasta", sep=""), seqtype = "DNA")
  #     seq5 = read.fasta(file = paste("./xmls/",timeaveraged, "_NP.fasta", sep=""), seqtype = "DNA")
  #     seq6 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB1.fasta", sep=""), seqtype = "DNA")
  #     seq7 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB2.fasta", sep=""), seqtype = "DNA")
  #     seq8 = read.fasta(file = paste("./xmls/",timeaveraged, "_PA.fasta", sep=""), seqtype = "DNA")
  #     
  #     #get the number of characters in the first segment
  #     sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
  #     sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
  #     sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
  #     sequence_length_seg4 <- nchar(getSequence(seq4[[1]], as.string = TRUE))
  #     sequence_length_seg5 <- nchar(getSequence(seq5[[1]], as.string = TRUE))
  #     sequence_length_seg6 <- nchar(getSequence(seq6[[1]], as.string = TRUE))
  #     sequence_length_seg7 <- nchar(getSequence(seq7[[1]], as.string = TRUE))
  #     sequence_length_seg8 <- nchar(getSequence(seq8[[1]], as.string = TRUE))
  #     writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
  #                                             sequence_length_seg3, sequence_length_seg4,
  #                                             sequence_length_seg5, sequence_length_seg6,
  #                                             sequence_length_seg7, sequence_length_seg8), line), f)
  #   } else {
  #     writeLines(line, f)
  #   }
  # }
  # close(f)
  # close(template)
  
  
  
  
  # read in the log file to get the network height
  filename = paste("H5N1_wgs_", timeaveraged , ".constant", sep="")
  # build an inference xml files
  f <- file(sprintf('xmls/%s.rep0.xml',filename), 'w')
  # Open the template file
  template <- file('../H5N1NorthAmerica/inference_template_wgs_cr.xml', 'r')
  
  base=""
  
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_name', line)) {
      writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
    } else if (grepl('insert_tips', line)) {
      # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
      for (i in seq(1, length(isolates))) {
        writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
      }
    }else if (grepl('insert_clock_rate', line)){
      writeLines(gsub('insert_clock_rate', '0.0035', line), f)
    }else if (grepl('<f idref="reassortmentRate"/>', line)){
      writeLines(line, f)
    }else if (grepl('<log idref="reassortmentRate"/>', line)){
      writeLines(line, f)
    }else if (grepl('insert_seg1', line)){
      writeLines(gsub('insert_seg1', paste(timeaveraged,"HA", sep="_"), line), f)
    }else if (grepl('insert_seg2', line)){
      writeLines(gsub('insert_seg2', paste(timeaveraged,"NA", sep="_"), line), f)
    }else if (grepl('insert_seg3', line)){
      writeLines(gsub('insert_seg3', paste(timeaveraged,"MP", sep="_"), line), f)
    }else if (grepl('insert_seg4', line)){
      writeLines(gsub('insert_seg4', paste(timeaveraged,"NS", sep="_"), line), f)
    }else if (grepl('insert_seg5', line)){
      writeLines(gsub('insert_seg5', paste(timeaveraged,"NP", sep="_"), line), f)
    }else if (grepl('insert_seg6', line)){
      writeLines(gsub('insert_seg6', paste(timeaveraged,"PB1", sep="_"), line), f)
    }else if (grepl('insert_seg7', line)){
      writeLines(gsub('insert_seg7', paste(timeaveraged,"PB2", sep="_"), line), f)
    }else if (grepl('insert_seg8', line)){
      writeLines(gsub('insert_seg8', paste(timeaveraged,"PA", sep="_"), line), f)
      
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
    } else if (grepl('insert_ratetimes', line)) {
      writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
      
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
      dates=c()
      for (i in seq(1, length(isolates))) {
        date = strsplit(isolates[i], "\\|")[[1]][[3]]
        dates=c(dates, as.Date(date))
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
      #make 3 points between 0 and EOS/365
      points = seq(0, EOS/365, length.out=10)
      writeLines(gsub('insert_EOS', paste(rateshiftvals2, collapse = " "), line), f)
    }else if (grepl('insert_weights', line)) {
      # get the length of both alignments
      seq1 = read.fasta(file = paste("./xmls/",timeaveraged, "_HA.fasta", sep=""), seqtype = "DNA")
      seq2 = read.fasta(file = paste("./xmls/",timeaveraged, "_NA.fasta", sep=""), seqtype = "DNA")
      seq3 = read.fasta(file = paste("./xmls/",timeaveraged, "_MP.fasta", sep=""), seqtype = "DNA")
      seq4 = read.fasta(file = paste("./xmls/",timeaveraged, "_NS.fasta", sep=""), seqtype = "DNA")
      seq5 = read.fasta(file = paste("./xmls/",timeaveraged, "_NP.fasta", sep=""), seqtype = "DNA")
      seq6 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB1.fasta", sep=""), seqtype = "DNA")
      seq7 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB2.fasta", sep=""), seqtype = "DNA")
      seq8 = read.fasta(file = paste("./xmls/",timeaveraged, "_PA.fasta", sep=""), seqtype = "DNA")

      #get the number of characters in the first segment
      sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
      sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
      sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
      sequence_length_seg4 <- nchar(getSequence(seq4[[1]], as.string = TRUE))
      sequence_length_seg5 <- nchar(getSequence(seq5[[1]], as.string = TRUE))
      sequence_length_seg6 <- nchar(getSequence(seq6[[1]], as.string = TRUE))
      sequence_length_seg7 <- nchar(getSequence(seq7[[1]], as.string = TRUE))
      sequence_length_seg8 <- nchar(getSequence(seq8[[1]], as.string = TRUE))
      writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
                                              sequence_length_seg3, sequence_length_seg4,
                                              sequence_length_seg5, sequence_length_seg6,
                                              sequence_length_seg7, sequence_length_seg8), line), f)
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)
  
  # make a second xml where the .trees is replaced by .infected
  filename = paste("H5N1_wgs_", timeaveraged , ".variable", sep="")
  f <- file(sprintf('xmls/%s.rep0.xml', filename), 'w')
  # Open the template file
  template <- file('../H5N1NorthAmerica/inference_template_wgs_cr.xml', 'r')
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_name', line)) {
      writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
    } else if (grepl('insert_tips', line)) {
      # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
      for (i in seq(1, length(isolates))) {
        writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
      }
    }else if (grepl('insert_clock_rate', line)){
      writeLines(gsub('insert_clock_rate', '0.0035', line), f)

    }else if (grepl('insert_seg1', line)){
      writeLines(gsub('insert_seg1', paste(timeaveraged,"HA", sep="_"), line), f)
    }else if (grepl('insert_seg2', line)){
      writeLines(gsub('insert_seg2', paste(timeaveraged,"NA", sep="_"), line), f)
    }else if (grepl('insert_seg3', line)){
      writeLines(gsub('insert_seg3', paste(timeaveraged,"MP", sep="_"), line), f)
    }else if (grepl('insert_seg4', line)){
      writeLines(gsub('insert_seg4', paste(timeaveraged,"NS", sep="_"), line), f)
    }else if (grepl('insert_seg5', line)){
      writeLines(gsub('insert_seg5', paste(timeaveraged,"NP", sep="_"), line), f)
    }else if (grepl('insert_seg6', line)){
      writeLines(gsub('insert_seg6', paste(timeaveraged,"PB1", sep="_"), line), f)
    }else if (grepl('insert_seg7', line)){
      writeLines(gsub('insert_seg7', paste(timeaveraged,"PB2", sep="_"), line), f)
    }else if (grepl('insert_seg8', line)){
      writeLines(gsub('insert_seg8', paste(timeaveraged,"PA", sep="_"), line), f)
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
    }else if (grepl('insert_EOS', line)){
      # get the most recent sampling time
      writeLines(gsub('insert_EOS', paste(rateshiftvals2[1:length(rateshiftvals2)], collapse = " "), line), f)
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
    } else if (grepl('insert_ratetimes', line)) {
      writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
    } else if (grepl('x="@reassortmentRate"',line)){
      # skip 1 line
      line <- readLines(template, n = 1)
      
      # write the following lines
      #<distribution spec="Prior">
      #  <x spec="coalre.dynamics.LogDifference" arg="@InfectedToRho" rateShift="@rateShifts"/>
      #    <distr spec="beast.base.inference.distribution.Normal" mean="0" sigma="1.0"/>
      writeLines('\t\t\t\t\t<distribution spec="Prior">', f)
      writeLines('\t\t\t\t\t\t<x spec="coalre.dynamics.LogDifference" arg="@reassortmentRate"/>', f)
      writeLines('\t\t\t\t\t\t<distr spec="beast.base.inference.distribution.Normal" mean="0" sigma="1.0"/>', f)
    
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
      seq1 = read.fasta(file = paste("./xmls/",timeaveraged, "_HA.fasta", sep=""), seqtype = "DNA")
      seq2 = read.fasta(file = paste("./xmls/",timeaveraged, "_NA.fasta", sep=""), seqtype = "DNA")
      seq3 = read.fasta(file = paste("./xmls/",timeaveraged, "_MP.fasta", sep=""), seqtype = "DNA")
      seq4 = read.fasta(file = paste("./xmls/",timeaveraged, "_NS.fasta", sep=""), seqtype = "DNA")
      seq5 = read.fasta(file = paste("./xmls/",timeaveraged, "_NP.fasta", sep=""), seqtype = "DNA")
      seq6 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB1.fasta", sep=""), seqtype = "DNA")
      seq7 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB2.fasta", sep=""), seqtype = "DNA")
      seq8 = read.fasta(file = paste("./xmls/",timeaveraged, "_PA.fasta", sep=""), seqtype = "DNA")
      
      #get the number of characters in the first segment
      sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
      sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
      sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
      sequence_length_seg4 <- nchar(getSequence(seq4[[1]], as.string = TRUE))
      sequence_length_seg5 <- nchar(getSequence(seq5[[1]], as.string = TRUE))
      sequence_length_seg6 <- nchar(getSequence(seq6[[1]], as.string = TRUE))
      sequence_length_seg7 <- nchar(getSequence(seq7[[1]], as.string = TRUE))
      sequence_length_seg8 <- nchar(getSequence(seq8[[1]], as.string = TRUE))
      writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
                                              sequence_length_seg3, sequence_length_seg4,
                                              sequence_length_seg5, sequence_length_seg6,
                                              sequence_length_seg7, sequence_length_seg8), line), f)      } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)
  
  # make a second xml where the .trees is replaced by .infected
  # filename = paste("H5N1_wgs_", timeaveraged , ".ne", sep="")
  # f <- file(sprintf('xmls/%s.rep0.xml', filename), 'w')
  # # Open the template file
  # template <- file('../H5N1NorthAmerica/inference_template_wgs.xml', 'r')
  # while (length(line <- readLines(template, n = 1)) > 0) {
  #   if (grepl('insert_name', line)) {
  #     writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
  #   } else if (grepl('insert_tips', line)) {
  #     # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
  #     for (i in seq(1, length(isolates))) {
  #       writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
  #     }
  #   }else if (grepl('insert_clock_rate', line)){
  #     writeLines(gsub('insert_clock_rate', '0.0035', line), f)
  # 
  #     writeLines(gsub('insert_seg1', paste(timeaveraged,"HA", sep="_"), line), f)
  #   }else if (grepl('insert_seg2', line)){
  #     writeLines(gsub('insert_seg2', paste(timeaveraged,"NA", sep="_"), line), f)
  #   }else if (grepl('insert_seg3', line)){
  #     writeLines(gsub('insert_seg3', paste(timeaveraged,"MP", sep="_"), line), f)
  #   }else if (grepl('insert_seg4', line)){
  #     writeLines(gsub('insert_seg4', paste(timeaveraged,"NS", sep="_"), line), f)
  #   }else if (grepl('insert_seg5', line)){
  #     writeLines(gsub('insert_seg5', paste(timeaveraged,"NP", sep="_"), line), f)
  #   }else if (grepl('insert_seg6', line)){
  #     writeLines(gsub('insert_seg6', paste(timeaveraged,"PB1", sep="_"), line), f)
  #   }else if (grepl('insert_seg7', line)){
  #     writeLines(gsub('insert_seg7', paste(timeaveraged,"PB2", sep="_"), line), f)
  #   }else if (grepl('insert_seg8', line)){
  #     writeLines(gsub('insert_seg8', paste(timeaveraged,"PA", sep="_"), line), f)
  #   } else if (grepl('insert_times', line)) {
  #     writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
  #     
  #   } else if (grepl('insert_times', line)) {
  #     writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
  #     # mask all the non isConstant entries, i.e. isinfectedSkyline
  #   } else if (grepl('isconstant-->', line)) {
  #     writeLines(gsub('isconstant-->', '', line), f)
  #   } else if (grepl('<!--isconstant', line)){
  #     writeLines(gsub('<!--isconstant', '', line), f)      
  #   } else if (grepl('isISkyline-->', line)) {
  #     writeLines(gsub('isISkyline-->', '', line), f)
  #   } else if (grepl('<!--isISkyline', line)){
  #     writeLines(gsub('<!--isISkyline', '', line), f)  
  #     
  #   } else if (grepl('isVariable-->', line)) {
  #     writeLines(gsub('isVariable-->', '', line), f)
  #   } else if (grepl('<!--isVariable', line)){
  #     writeLines(gsub('<!--isVariable', '', line), f)  
  #     
  #     
  #   } else if (grepl('insert_heights', line)) {
  #     # write the height of each tip using heights seperated by a = and , between
  #     value=''
  #     dates=c()
  #     for (i in seq(1, length(isolates))) {
  #       date = strsplit(isolates[i], "\\|")[[1]][[3]]
  #       dates=c(dates, as.Date(date))
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
  #     # get the length of both alignments
  #     seq1 = read.fasta(file = paste("./xmls/",timeaveraged, "_HA.fasta", sep=""), seqtype = "DNA")
  #     seq2 = read.fasta(file = paste("./xmls/",timeaveraged, "_NA.fasta", sep=""), seqtype = "DNA")
  #     seq3 = read.fasta(file = paste("./xmls/",timeaveraged, "_MP.fasta", sep=""), seqtype = "DNA")
  #     seq4 = read.fasta(file = paste("./xmls/",timeaveraged, "_NS.fasta", sep=""), seqtype = "DNA")
  #     seq5 = read.fasta(file = paste("./xmls/",timeaveraged, "_NP.fasta", sep=""), seqtype = "DNA")
  #     seq6 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB1.fasta", sep=""), seqtype = "DNA")
  #     seq7 = read.fasta(file = paste("./xmls/",timeaveraged, "_PB2.fasta", sep=""), seqtype = "DNA")
  #     seq8 = read.fasta(file = paste("./xmls/",timeaveraged, "_PA.fasta", sep=""), seqtype = "DNA")
  #     
  #     
  #     #get the number of characters in the first segment
  #     sequence_length_seg1 <- nchar(getSequence(seq1[[1]], as.string = TRUE))
  #     sequence_length_seg2 <- nchar(getSequence(seq2[[1]], as.string = TRUE))
  #     sequence_length_seg3 <- nchar(getSequence(seq3[[1]], as.string = TRUE))
  #     sequence_length_seg4 <- nchar(getSequence(seq4[[1]], as.string = TRUE))
  #     sequence_length_seg5 <- nchar(getSequence(seq5[[1]], as.string = TRUE))
  #     sequence_length_seg6 <- nchar(getSequence(seq6[[1]], as.string = TRUE))
  #     sequence_length_seg7 <- nchar(getSequence(seq7[[1]], as.string = TRUE))
  #     sequence_length_seg8 <- nchar(getSequence(seq8[[1]], as.string = TRUE))
  #     writeLines(gsub('insert_weights', paste(sequence_length_seg1, sequence_length_seg2,
  #                                             sequence_length_seg3, sequence_length_seg4,
  #                                             sequence_length_seg5, sequence_length_seg6,
  #                                             sequence_length_seg7, sequence_length_seg8), line), f)      } else {
  #     writeLines(line, f)
  #   }
  # }
  # close(f)
  # close(template)
  
  
  
}
# read in the H5N1 varying xml file and add/change some lines

# make two copies of the rep0 files into rep1 and rep2
files = list.files("xmls", pattern="*.rep0.xml", full.names=TRUE)
for (file in files) {
  file.copy(file, gsub("rep0", "rep1", file))
  file.copy(file, gsub("rep0", "rep2", file))
}



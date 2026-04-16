library(stringr)
library(seqinr)
# Clear workspace
rm(list=ls())

# This script creates XML files for BEAST analysis of H5N1 North America data
# The first type of XML created uses all H5 HPAI and LPAI data
# Option to subsample LPAI data while keeping all HPAI data is available

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- getwd()
setwd(this.dir)

# Configuration options
SUBSAMPLE_LPAI <- TRUE  # Set to TRUE to subsample LPAI data, FALSE to use all data
LPAI_SUBSAMPLE_FRACTION <- 0.25  # Fraction of LPAI sequences to keep (0.5 = 50%)

# Load HPAI/LPAI classification data
clade_file = read.csv("./tables/HPAI_LPAI.csv", stringsAsFactors = FALSE, header = TRUE)
# remove all ' ( and ) from names
clade_file$taxa = gsub("'", "", clade_file$taxa)
clade_file$taxa = gsub("\\(", "", clade_file$taxa)
clade_file$taxa = gsub("\\)", "", clade_file$taxa)

# delete all *wgs*.xml files in xmls
file.remove(list.files("./xmls/", pattern="*.xml", full.names=TRUE))
file.remove(list.files("./xmls/", pattern="*.fasta", full.names=TRUE))

# Get all fasta files
fasta_files <- list.files("./data/", pattern = "\\.fasta$", full.names = TRUE)
# remove files containing H5
fasta_files <- fasta_files[!grepl("H5", fasta_files)]

# Quote-like characters to remove (apostrophes, curly quotes, etc.)
quote_chars <- "['‘’`´ʼʹʽˈˊ′‵‶‷]"

# Loop through each file
for (file in fasta_files) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  
  # Clean headers
  lines_cleaned <- ifelse(
    grepl("^>", lines),
    gsub(paste0(quote_chars, "|[\\(\\)]"), "", lines),  # remove quotes + parens
    lines
  )
  
  out_file <- file.path("xmls", basename(file))
  writeLines(lines_cleaned, out_file, useBytes = TRUE)
}

# Apply LPAI subsampling to fasta files if enabled
if (SUBSAMPLE_LPAI) {
  # Get HPAI and LPAI sequences from the classification file
  hpai_sequences = clade_file$taxa[clade_file$status == "HPAI"]
  lpai_sequences = clade_file$taxa[clade_file$status == "LPAI"]
  
  # Subsample LPAI sequences
  if (length(lpai_sequences) > 0) {
    n_lpai_to_keep = max(1, round(length(lpai_sequences) * LPAI_SUBSAMPLE_FRACTION))
    lpai_to_keep = sample(lpai_sequences, n_lpai_to_keep)
  } else {
    lpai_to_keep = character(0)
  }
  
  # Combine HPAI and subsampled LPAI sequences
  sequences_to_keep = c(hpai_sequences, lpai_to_keep)
  
  cat(sprintf("Subsampling LPAI data: keeping %d/%d LPAI sequences (%.1f%%), %d HPAI sequences\n", 
              length(lpai_to_keep), length(lpai_sequences), 
              length(lpai_to_keep)/max(1,length(lpai_sequences))*100, 
              length(hpai_sequences)))
  
  # Apply subsampling to all fasta files in xmls directory
  xmls_fasta_files <- list.files("./xmls/", pattern = "\\.fasta$", full.names = TRUE)
  for (fasta_file in xmls_fasta_files) {
    # Read the fasta file
    fasta_data = seqinr::read.fasta(file = fasta_file, seqtype = "DNA")
    
    # Filter to keep only the selected sequences
    fasta_data_filtered = fasta_data[names(fasta_data) %in% sequences_to_keep]
    
    # Write the filtered fasta file back
    seqinr::write.fasta(sequences = fasta_data_filtered, 
                        names = names(fasta_data_filtered), 
                        file.out = fasta_file)
  }
}



fasta_files <- list.files("./xmls/", pattern = "\\.fasta$", full.names = TRUE)

# remove all fasta files that do not contain HA
fasta_files = fasta_files[grep("HA", fasta_files)]

for (fastafile in fasta_files){
  # read in the first fasta file and keep all the names
  fasta1 = seqinr::read.fasta(file = fastafile, seqtype = "DNA")
  isolates = names(fasta1)
  
  # collect the dates as the last group after splitting on | 
  dates = sapply(isolates, function(x) strsplit(x, "\\|")[[1]][[3]])
  min = min(as.Date(dates))
  max = max(as.Date(dates))
  start=as.Date("2021-01-01")
  first_intro = as.numeric(max-start)/365
  

  # define the root height for Ne and reassortment variant rates 
  # this is the time of the first introduction
  rateshiftvals = c(seq(0,  ceiling(first_intro), length.out=5), seq(ceiling(first_intro+1),  50, length.out=2), 200)
  # rateshiftvals = seq(0,  first_intro, length.out=5)
  rateshiftvals = unique(rateshiftvals)
  rateshiftvals2 = rateshiftvals

  # get the filename as everything before the first _
  filename_base = strsplit(basename(fastafile), "_")[[1]][1]
  
  # read in the log file to get the network height
  filename = paste(filename_base, ".skygrowth", sep="")
  # build an inference xml files
  f <- file(sprintf('xmls/%s.rep0.xml',filename), 'w')
  # Open the template file
  template <- file('../H5N1NorthAmerica/inference_template_wgs_cr_fixed_evol.xml', 'r')
  
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
      writeLines(gsub('insert_seg1', paste(filename_base,"HA", sep="_"), line), f)
    }else if (grepl('insert_seg2', line)){
      writeLines(gsub('insert_seg2', paste(filename_base,"NA", sep="_"), line), f)
    }else if (grepl('insert_seg3', line)){
      writeLines(gsub('insert_seg3', paste(filename_base,"MP", sep="_"), line), f)
    }else if (grepl('insert_seg4', line)){
      writeLines(gsub('insert_seg4', paste(filename_base,"NS", sep="_"), line), f)
    }else if (grepl('insert_seg5', line)){
      writeLines(gsub('insert_seg5', paste(filename_base,"NP", sep="_"), line), f)
    }else if (grepl('insert_seg6', line)){
      writeLines(gsub('insert_seg6', paste(filename_base,"PB1", sep="_"), line), f)
    }else if (grepl('insert_seg7', line)){
      writeLines(gsub('insert_seg7', paste(filename_base,"PB2", sep="_"), line), f)
    }else if (grepl('insert_seg8', line)){
      writeLines(gsub('insert_seg8', paste(filename_base,"PA", sep="_"), line), f)
      
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
    } else if (grepl('insert_ratetimes', line)) {
      writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
      
      # mask all the non isConstant entries, i.e. isinfectedSkyline
    } else if (grepl('<f idref="sigma2"/>', line)) {
      # skip this line for constant models
    } else if (grepl('<operator id="Sigma2Scaler.s:"', line)) {
      # skip the sigma2 operator block for constant models
      while (length(line <- readLines(template, n = 1)) > 0 && !grepl('</operator>', line)) {
        # skip lines until we find the closing operator tag
      }
      # the </operator> line is also skipped
    } else if (grepl('isconstant-->', line)) {
      writeLines(gsub('isconstant-->', '', line), f)
    } else if (grepl('<!--isconstant', line)){
      writeLines(gsub('<!--isconstant', '', line), f)
    # } else if (grepl('isinfectedSkyline-->', line)) {
    #   writeLines(gsub('isinfectedSkyline-->', '', line), f)
    # } else if (grepl('<!--isinfectedSkyline', line)){
    #   writeLines(gsub('<!--isinfectedSkyline', '', line), f)
    } else if (grepl('isNeSkyline-->', line)) {
      writeLines(gsub('isNeSkyline-->', '', line), f)
    } else if (grepl('<!--isNeSkyline', line)){
      writeLines(gsub('<!--isNeSkyline', '', line), f)
    # } else if (grepl('isISkyline-->', line)) {
    #   writeLines(gsub('isISkyline-->', '', line), f)
    # } else if (grepl('<!--isISkyline', line)){
    #   writeLines(gsub('<!--isISkyline', '', line), f)
      
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
      seq1 = seqinr::read.fasta(file = paste("./data/",filename_base, "_HA.fasta", sep=""), seqtype = "DNA")
      seq2 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NA.fasta", sep=""), seqtype = "DNA")
      seq3 = seqinr::read.fasta(file = paste("./data/",filename_base, "_MP.fasta", sep=""), seqtype = "DNA")
      seq4 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NS.fasta", sep=""), seqtype = "DNA")
      seq5 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NP.fasta", sep=""), seqtype = "DNA")
      seq6 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB1.fasta", sep=""), seqtype = "DNA")
      seq7 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB2.fasta", sep=""), seqtype = "DNA")
      seq8 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PA.fasta", sep=""), seqtype = "DNA")

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
    }else if (grepl('insert_independent_after', line)) {
      writeLines(gsub('independentAfter="insert_independent_after"', '', line), f)
      # writeLines(gsub('insert_independent_after', length(rateshiftvals)-1, line), f)
    }else if (grepl('insert_double_discount', line)){
      writeLines(gsub('insert_double_discount', 'true', line), f)
      
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)
  next
  
  # make a second xml where the .trees is replaced by .infected
  filename = paste(filename_base, ".independent", sep="")
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
      writeLines(gsub('insert_seg1', paste(filename_base,"HA", sep="_"), line), f)
    }else if (grepl('insert_seg2', line)){
      writeLines(gsub('insert_seg2', paste(filename_base,"NA", sep="_"), line), f)
    }else if (grepl('insert_seg3', line)){
      writeLines(gsub('insert_seg3', paste(filename_base,"MP", sep="_"), line), f)
    }else if (grepl('insert_seg4', line)){
      writeLines(gsub('insert_seg4', paste(filename_base,"NS", sep="_"), line), f)
    }else if (grepl('insert_seg5', line)){
      writeLines(gsub('insert_seg5', paste(filename_base,"NP", sep="_"), line), f)
    }else if (grepl('insert_seg6', line)){
      writeLines(gsub('insert_seg6', paste(filename_base,"PB1", sep="_"), line), f)
    }else if (grepl('insert_seg7', line)){
      writeLines(gsub('insert_seg7', paste(filename_base,"PB2", sep="_"), line), f)
    }else if (grepl('insert_seg8', line)){
      writeLines(gsub('insert_seg8', paste(filename_base,"PA", sep="_"), line), f)
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

    # } else if (grepl('isISkyline-->', line)) {
    #   writeLines(gsub('isISkyline-->', '', line), f)
    # } else if (grepl('<!--isISkyline', line)){
    #   writeLines(gsub('<!--isISkyline', '', line), f)

    # } else if (grepl('isinfectedSkyline-->', line)) {
    #   writeLines(gsub('isinfectedSkyline-->', '', line), f)
    # } else if (grepl('<!--isinfectedSkyline', line)){
    #   writeLines(gsub('<!--isinfectedSkyline', '', line), f)


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
    }else if (grepl('insert_independent_after', line)) {
      writeLines(gsub('independentAfter="insert_independent_after"', '', line), f)
      # writeLines(gsub('insert_independent_after', length(rateshiftvals)-1, line), f)
    }else if (grepl('insert_double_discount', line)){
      writeLines(gsub('insert_double_discount', 'false', line), f)
      
      
    }else if (grepl('insert_weights', line)) {
      # get the length of both alignments
      # get the length of both alignments
      seq1 = seqinr::read.fasta(file = paste("./data/",filename_base, "_HA.fasta", sep=""), seqtype = "DNA")
      seq2 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NA.fasta", sep=""), seqtype = "DNA")
      seq3 = seqinr::read.fasta(file = paste("./data/",filename_base, "_MP.fasta", sep=""), seqtype = "DNA")
      seq4 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NS.fasta", sep=""), seqtype = "DNA")
      seq5 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NP.fasta", sep=""), seqtype = "DNA")
      seq6 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB1.fasta", sep=""), seqtype = "DNA")
      seq7 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB2.fasta", sep=""), seqtype = "DNA")
      seq8 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PA.fasta", sep=""), seqtype = "DNA")
      
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
}


rateshiftvals = c(seq(0,  first_intro, length.out=40), seq(ceiling(first_intro+1),  50, length.out=5), 200)
rateshiftvals = unique(rateshiftvals)
rateshiftvals2 = rateshiftvals

clade_file = read.csv("./tables/HPAI_LPAI.csv", stringsAsFactors = FALSE, header = TRUE)
# remove all ' ( and ) from names
clade_file$taxa = gsub("'", "", clade_file$taxa)
clade_file$taxa = gsub("\\(", "", clade_file$taxa)
clade_file$taxa = gsub("\\)", "", clade_file$taxa)

for (clade in unique(clade_file$status)){
  next
  
  isolates = clade_file[clade_file$status == "HPAI", "taxa"]
  
  
  # read in all fasta files in data
  fasta = list.files("./xmls/", pattern="*.fasta", full.names=TRUE)
  # only keep the ones starting with HLHxNx
  fasta = fasta[startsWith(basename(fasta), "HL")]
  for (fastafile in fasta){
    # read in the fasta file
    fasta1 = seqinr::read.fasta(file = fastafile, seqtype = "DNA")
    # remove all sequences not in isolates
    if (clade=="HPAI"){
      fasta1 = fasta1[names(fasta1) %in% isolates]
      
    }else{
      fasta1 = fasta1[!names(fasta1) %in% isolates]
    }
    # write the fasta file to the xmls3seg directory
    seqinr::write.fasta(sequences = fasta1, 
                        names = names(fasta1), 
                        file.out = paste0("./xmls/", clade, "_", basename(fastafile)))
  }
  
  for (fastafile in fasta_files){
    # read in the first fasta file and keep all the names
    isolates = names(fasta1)
    # collect the dates as the last group after splitting on | 
    dates = sapply(isolates, function(x) strsplit(x, "\\|")[[1]][[3]])
    min = min(as.Date(dates))
    max = max(as.Date(dates))
    first_intro = as.numeric(max-min)/365+0.5
    

    # get the filename as everything before the first _
    filename_base = strsplit(basename(fastafile), "_")[[1]][1]

    # read in the log file to get the network height
    filename = paste(clade, "_", filename_base, ".dependent", sep="")
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
        writeLines(gsub('insert_seg1', paste(clade, filename_base,"HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(clade, filename_base,"NA", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(clade, filename_base,"MP", sep="_"), line), f)
      }else if (grepl('insert_seg4', line)){
        writeLines(gsub('insert_seg4', paste(clade, filename_base,"NS", sep="_"), line), f)
      }else if (grepl('insert_seg5', line)){
        writeLines(gsub('insert_seg5', paste(clade, filename_base,"NP", sep="_"), line), f)
      }else if (grepl('insert_seg6', line)){
        writeLines(gsub('insert_seg6', paste(clade, filename_base,"PB1", sep="_"), line), f)
      }else if (grepl('insert_seg7', line)){
        writeLines(gsub('insert_seg7', paste(clade, filename_base,"PB2", sep="_"), line), f)
      }else if (grepl('insert_seg8', line)){
        writeLines(gsub('insert_seg8', paste(clade, filename_base,"PA", sep="_"), line), f)
        
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
      } else if (grepl('insert_ratetimes', line)) {
        writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
        
        # mask all the non isConstant entries, i.e. isinfectedSkyline
      } else if (grepl('isconstant-->', line)) {
        writeLines(gsub('isconstant-->', '', line), f)
      } else if (grepl('<!--isconstant', line)){
        writeLines(gsub('<!--isconstant', '', line), f)
        # } else if (grepl('isinfectedSkyline-->', line)) {
        #   writeLines(gsub('isinfectedSkyline-->', '', line), f)
        # } else if (grepl('<!--isinfectedSkyline', line)){
        #   writeLines(gsub('<!--isinfectedSkyline', '', line), f)      
        # } else if (grepl('isNeSkyline-->', line)) {
        #   writeLines(gsub('isNeSkyline-->', '', line), f)
        # } else if (grepl('<!--isNeSkyline', line)){
        #   writeLines(gsub('<!--isNeSkyline', '', line), f)  
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
        seq1 = seqinr::read.fasta(file = paste("./data/",filename_base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NA.fasta", sep=""), seqtype = "DNA")
        seq3 = seqinr::read.fasta(file = paste("./data/",filename_base, "_MP.fasta", sep=""), seqtype = "DNA")
        seq4 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq5 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NP.fasta", sep=""), seqtype = "DNA")
        seq6 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB1.fasta", sep=""), seqtype = "DNA")
        seq7 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB2.fasta", sep=""), seqtype = "DNA")
        seq8 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PA.fasta", sep=""), seqtype = "DNA")
        
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
      }else if (grepl('insert_independent_after', line)) {
        writeLines(gsub('independentAfter="insert_independent_after"', '', line), f)
        # writeLines(gsub('insert_independent_after', length(rateshiftvals)-1, line), f)
      }else if (grepl('insert_double_discount', line)){
        writeLines(gsub('insert_double_discount', 'true', line), f)
        
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)
    
    # make a second xml where the .trees is replaced by .infected
    filename = paste(clade, "_", filename_base, ".independent", sep="")
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
        writeLines(gsub('insert_seg1', paste(clade, filename_base,"HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(clade, filename_base,"NA", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(clade, filename_base,"MP", sep="_"), line), f)
      }else if (grepl('insert_seg4', line)){
        writeLines(gsub('insert_seg4', paste(clade, filename_base,"NS", sep="_"), line), f)
      }else if (grepl('insert_seg5', line)){
        writeLines(gsub('insert_seg5', paste(clade, filename_base,"NP", sep="_"), line), f)
      }else if (grepl('insert_seg6', line)){
        writeLines(gsub('insert_seg6', paste(clade, filename_base,"PB1", sep="_"), line), f)
      }else if (grepl('insert_seg7', line)){
        writeLines(gsub('insert_seg7', paste(clade, filename_base,"PB2", sep="_"), line), f)
      }else if (grepl('insert_seg8', line)){
        writeLines(gsub('insert_seg8', paste(clade, filename_base,"PA", sep="_"), line), f)
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
        
        # } else if (grepl('isISkyline-->', line)) {
        #   writeLines(gsub('isISkyline-->', '', line), f)
        # } else if (grepl('<!--isISkyline', line)){
        #   writeLines(gsub('<!--isISkyline', '', line), f)
        
        # } else if (grepl('isinfectedSkyline-->', line)) {
        #   writeLines(gsub('isinfectedSkyline-->', '', line), f)
        # } else if (grepl('<!--isinfectedSkyline', line)){
        #   writeLines(gsub('<!--isinfectedSkyline', '', line), f)
        
        
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
      }else if (grepl('insert_independent_after', line)) {
        writeLines(gsub('independentAfter="insert_independent_after"', '', line), f)
        # writeLines(gsub('insert_independent_after', length(rateshiftvals)-1, line), f)
      }else if (grepl('insert_double_discount', line)){
        writeLines(gsub('insert_double_discount', 'false', line), f)
        
        
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        # get the length of both alignments
        seq1 = seqinr::read.fasta(file = paste("./data/",filename_base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NA.fasta", sep=""), seqtype = "DNA")
        seq3 = seqinr::read.fasta(file = paste("./data/",filename_base, "_MP.fasta", sep=""), seqtype = "DNA")
        seq4 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq5 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NP.fasta", sep=""), seqtype = "DNA")
        seq6 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB1.fasta", sep=""), seqtype = "DNA")
        seq7 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB2.fasta", sep=""), seqtype = "DNA")
        seq8 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PA.fasta", sep=""), seqtype = "DNA")
        
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
  }

}


############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
#
# Estimate for case positivity
#
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################



# make two copies of the rep0 files into rep1 and rep2
clade_file = read.csv("./tables/HPAI_LPAI.csv", stringsAsFactors = FALSE, header = TRUE)
# replace all ' in isolates
clade_file$taxa = gsub("'", "", clade_file$taxa)
isolates = clade_file[clade_file$status == "HPAI", "taxa"]
clade ="HPAI"
# read in all fasta files in data
fasta = list.files("./data/", pattern="*.fasta", full.names=TRUE)
for (fastafile in fasta){
  # read in the fasta file
  fasta1 = seqinr::read.fasta(file = fastafile, seqtype = "DNA")
  # remove all sequences not in isolates
  if (clade=="HPAI"){
    fasta1 = fasta1[names(fasta1) %in% isolates]
  }else{
    fasta1 = fasta1[!names(fasta1) %in% isolates]
  }
  # write the fasta file to the xmls3seg directory
  seqinr::write.fasta(sequences = fasta1, 
                      names = names(fasta1), 
                      file.out = paste0("./xmls/", clade, "_", basename(fastafile)))
}

library(lubridate)



for (fastafile in fasta_files){
  # read in the first fasta file and keep all the names
  isolates = names(fasta1)
  # collect the dates as the last group after splitting on | 
  dates = sapply(isolates, function(x) strsplit(x, "\\|")[[1]][[3]])
  min = min(as.Date(dates))
  max = max(as.Date(dates))
  first_intro = as.numeric(max-min)/365

  rateshiftvals = c(seq(0,  first_intro, length.out=40), seq(ceiling(first_intro),  10, length.out=2))
  rateshiftvals = unique(rateshiftvals)
  rateshiftvals2 = rateshiftvals


  
  cases = read.csv("./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard_with_flyways.csv")
  # convert dates to decimal
  cases$date = as.Date(cases$Date_Collected, format="%Y-%m-%d")
  # now, convert to decimal
  cases$decimal_date <- decimal_date(cases$date)
  
  # get the first rate shift value above first_intro
  first_intro_index = which(rateshiftvals >= first_intro)[1]
  
  
  
  # compute the smoothed average for the cases as at the time points of the rate
  # shifts
  # loop over all days
  max_date = decimal_date(max)
  smoothed_case_data = data.frame()
  # smoothing area (expanded window)
  diff = 2 * min(diff(rateshiftvals))
  for (d in max_date-rateshiftvals){
    # get all instances within that time window
    window = cases[cases$decimal_date >= d-diff & cases$decimal_date <= d+diff, ]
    window_alt = cases[cases$decimal_date >= d-diff-7/365 & cases$decimal_date <= d+diff-7/365, ]
    
    
    # get how many instances of Final_IAV are Detected
    total_AIV = sum(window$Final_IAV == "Detected", na.rm=TRUE)
    total_H5 = sum(window$Final_H5 == "Detected", na.rm=TRUE)
    
    # get the number of high path cases
    total_HPAI = sum(window$Final_H5 == "Detected" &  window$Final_Pathogenicity == "High Path AI", na.rm=TRUE)
    
    lpai = window[window$Final_H5 == "Detected" &  window$Final_Pathogenicity != "High Path AI", ]
    hpai = window[window$Final_H5 == "Detected" &  window$Final_Pathogenicity == "High Path AI", ]
    
    lpai_alt = window_alt[window_alt$Final_H5 == "Detected" &  window_alt$Final_Pathogenicity != "High Path AI", ]
    hpai_alt = window_alt[window_alt$Final_H5 == "Detected" &  window_alt$Final_Pathogenicity == "High Path AI", ]
    
    print(lpai$flyway)
    # compute overlap as function of LPAI prevalence * HPAI prevalence
    if (nrow(lpai) > 0 && nrow(hpai) > 0){
      # Calculate total samples per flyway in this window
      total_samples_by_flyway = aggregate(rep(1, nrow(window)), by=list(flyway=window$flyway), FUN=sum)
      
      # Calculate LPAI prevalence by flyway
      lpai_counts = aggregate(rep(1, nrow(lpai)), by=list(flyway=lpai$flyway), FUN=sum)
      names(lpai_counts)[2] = "lpai_count"
      lpai_prevalence = merge(total_samples_by_flyway, lpai_counts, by="flyway", all.x=TRUE)
      names(lpai_prevalence)[2] = "total_count"
      lpai_prevalence$lpai_count[is.na(lpai_prevalence$lpai_count)] = 0
      lpai_prevalence$lpai_prevalence = lpai_prevalence$lpai_count / lpai_prevalence$total_count
      
      # Calculate HPAI prevalence by flyway
      hpai_counts = aggregate(rep(1, nrow(hpai)), by=list(flyway=hpai$flyway), FUN=sum)
      names(hpai_counts)[2] = "hpai_count"
      hpai_prevalence = merge(total_samples_by_flyway, hpai_counts, by="flyway", all.x=TRUE)
      names(hpai_prevalence)[2] = "total_count"
      hpai_prevalence$hpai_count[is.na(hpai_prevalence$hpai_count)] = 0
      hpai_prevalence$hpai_prevalence = hpai_prevalence$hpai_count / hpai_prevalence$total_count
      
      # Calculate overlap as sum of flyways that have both LPAI and HPAI cases
      # Get flyways with LPAI cases
      lpai_flyways = unique(lpai_prevalence$flyway[lpai_prevalence$lpai_count > 0])
      # Get flyways with HPAI cases  
      hpai_flyways = unique(hpai_prevalence$flyway[hpai_prevalence$hpai_count > 0])
      # calculate the product of the prevalence of the two flyways
      product = lpai_prevalence$lpai_prevalence * hpai_prevalence$hpai_prevalence
      overlap = sum(product)
    }else{
      overlap = 0
    }
    
    
    # smoothed_case_data = rbind(smoothed_case_data, data.frame(
    #   date = d,
    #   positivity = (total_H5-total_HPAI),
    #   type = "lpai_nosummer"
    # ))
    
    # smoothed_case_data = rbind(smoothed_case_data, data.frame(
    #   date = d,
    #   positivity = (total_H5-total_HPAI),
    #   type = "h5_lpai"
    # ))
    smoothed_case_data = rbind(smoothed_case_data, data.frame(
      date = d,
      positivity = (total_AIV-total_HPAI),
      type = "lpai"
    ))
    smoothed_case_data = rbind(smoothed_case_data, data.frame(
      date = d,
      positivity = total_HPAI,
      type = "hpai"
    ))
    smoothed_case_data = rbind(smoothed_case_data, data.frame(
      date = d,
      positivity = total_AIV,
      type = "total"
    ))

    # smoothed_case_data = rbind(smoothed_case_data, data.frame(
    #   date = d,
    #   positivity = overlap,
    #   type = "overlap"
    # ))
    
  }
  
  # set the values for summer 2022 to min value for the no_summer predictor
  smoothed_case_data$positivity[smoothed_case_data$type == "lpai_nosummer" & smoothed_case_data$date >= decimal_date(as.Date("2022-07-01")) & smoothed_case_data$date <= decimal_date(as.Date("2023-03-01"))] = 0 
  
  # plot the smoothed case data
  library(ggplot2)
  l=length(unique(smoothed_case_data$type))
  # set the first values of each type that where NaN to the second values
  smoothed_case_data$positivity[1:l] = smoothed_case_data$positivity[(l+1):(2*l)]  # set all values that are Na to 0
  smoothed_case_data$positivity[is.na(smoothed_case_data$positivity)] <- 0
  # add the mimum case positivity above 0 to each value
  
  # set independetafter as the rate shift corresponding to the first_intro
  independentafter = first_intro_index
  


  # get the filename as everything before the first _
  filename_base = strsplit(basename(fastafile), "_")[[1]][1]
  
  # read in the log file to get the network height
  filename = paste(clade, "_", filename_base, ".glm" , sep="")
  # build an inference xml files
  f <- file(sprintf('xmls/%s.rep0.xml',filename), 'w')
  # Open the template file
  template <- file('../H5N1NorthAmerica/inference_template_wgs_cr.xml', 'r')
  
  base=""
  plot_data = data.frame()
  
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_name', line)) {
      writeLines(gsub('insert_name', sprintf('%s', filename), line), f)
    }else if (grepl('spec="coalre.dynamics.SkygrowthReassortmentRatesFromSkygrowthNe" logNe="@logNe"', line)){
      writeLines(sprintf('\t\t\t\t\t<timeVaryingReassortmentRates id="splinePop2" spec="coalre.dynamics.GLMReassortmentRates" independentAfter="%d" neToReassortment="@InfectedToRho"  rateShifts="@rateShifts" predictorIsActive="@predictorActive">\n', independentafter), f)  
      
      for (case_type in  unique(smoothed_case_data$type)){
        # make a vector out of the HPAI cases
        writeLines(sprintf('\t\t\t\t\t\t<predictor idref="%s"/>\n', case_type), f)
      }
      writeLines(sprintf('\t\t\t\t\t\t<predictor idref="logNe"/>\n'), f)
      writeLines(sprintf('\t\t\t\t\t\t<effectSize idref="effectSize"/>\n'), f)
      writeLines('\t\t\t\t\t</timeVaryingReassortmentRates>', f)
    
    } else if (grepl('insert_tips', line)) {
      # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
      for (i in seq(1, length(isolates))) {
        writeLines(sprintf('\t\t<taxon spec="Taxon" id="%s"/>', isolates[i]), f)
      }
    }else if (grepl('insert_clock_rate', line)){
      writeLines(gsub('insert_clock_rate', '0.0035', line), f)
      for (case_type in  unique(smoothed_case_data$type)){
        # make a vector out of the HPAI cases
        cases_vector_all = smoothed_case_data$positivity[smoothed_case_data$type == case_type]
        cases_vector = cases_vector_all[1:independentafter]
        time_all=smoothed_case_data$date[smoothed_case_data$type == case_type]
        
        # log standardize the cases vector
        log_cases_vector = log(cases_vector+min(cases_vector[cases_vector>0]))
        log_cases_vector = log_cases_vector - mean(log_cases_vector)
        log_cases_vector = log_cases_vector / sd(log_cases_vector)
        
        plot_data = rbind(plot_data, data.frame(
          time = time_all[1:independentafter],
          positivity = log_cases_vector,
          type = case_type
        ))
        
        # fill up with 0's
        log_cases_vector = c(log_cases_vector, rep(0, length(rateshiftvals)-length(log_cases_vector)))
        
        cases_string = paste(log_cases_vector, collapse=' ')
        writeLines(sprintf('\t\t\t<stateNode id="%s" spec="parameter.RealParameter" value="%s"/>\n', case_type, cases_string), f)
      }
      writeLines('\t\t\t<parameter id="predictorActive" spec="parameter.IntegerParameter" name="stateNode" upper="4" lower="0">4</parameter>\n', f)
      writeLines('\t\t\t<parameter id="effectSize" spec="parameter.RealParameter" name="stateNode" lower="0">1</parameter>\n', f)
    }else if (grepl(' <operator id="FixMeanMutationRatesOperator"', line)){
      writeLines(sprintf('\t\t\t<operator id="PredictorOperator" spec="ChangePredictorOperator" weight="5"  independentAfter="%d" predictorIsActive="@predictorActive" neToReassortment="@InfectedToRho">',independentafter), f)
      for (case_type in  unique(smoothed_case_data$type)){
        # make a vector out of the HPAI cases
        writeLines(sprintf('\t\t\t\t<predictor idref="%s"/>', case_type), f)
      }
      writeLines(sprintf('\t\t\t\t<predictor idref="logNe"/>'), f)
      writeLines(sprintf('\t\t\t\t<effectSize idref="effectSize"/>'), f)
      writeLines('\t\t\t</operator>\n', f);
      
      writeLines(sprintf('\t\t\t<operator id="EffectSizeOperator" spec="EffectSizePredictorOperator" weight="1"  independentAfter="%d" predictorIsActive="@predictorActive" neToReassortment="@InfectedToRho">',independentafter), f)
      for (case_type in  unique(smoothed_case_data$type)){
        # make a vector out of the HPAI cases
        writeLines(sprintf('\t\t\t\t<predictor idref="%s"/>', case_type), f)
      }
      writeLines(sprintf('\t\t\t\t<predictor idref="logNe"/>'), f)
      writeLines(sprintf('\t\t\t\t<effectSize idref="effectSize"/>'), f)
      writeLines('\t\t\t</operator>\n', f);
      
      writeLines('\t\t\t<operator id="EffectiveSize.ScalerX.t" spec="beast.base.inference.operator.kernel.BactrianRandomWalkOperator" parameter="@effectSize" scaleFactor="0.5" weight="1.0"/>\n',f);
      writeLines('\t\t\t<operator id="Uni.ScalerX.t" spec="beast.base.inference.operator.UniformOperator" parameter="@predictorActive" weight="1.0"/>\n',f);        

      writeLines(line, f)
    }else if (grepl('<f idref="reassortmentRate"/>', line)){
      writeLines(line, f)
    }else if (grepl('<log idref="treeLikelihood.seg8"/>', line)){
      writeLines(line, f)
      writeLines('\t\t\t<log idref="predictorActive"/>\n', f)
      writeLines('\t\t\t<log idref="effectSize"/>\n', f)
      writeLines('\t\t\t<log idref="splinePop2"/>\n', f)
    }else if (grepl('insert_seg1', line)){
      writeLines(gsub('insert_seg1', paste(clade, filename_base,"HA", sep="_"), line), f)
    }else if (grepl('insert_seg2', line)){
      writeLines(gsub('insert_seg2', paste(clade, filename_base,"NA", sep="_"), line), f)
    }else if (grepl('insert_seg3', line)){
      writeLines(gsub('insert_seg3', paste(clade, filename_base,"MP", sep="_"), line), f)
    }else if (grepl('insert_seg4', line)){
      writeLines(gsub('insert_seg4', paste(clade, filename_base,"NS", sep="_"), line), f)
    }else if (grepl('insert_seg5', line)){
      writeLines(gsub('insert_seg5', paste(clade, filename_base,"NP", sep="_"), line), f)
    }else if (grepl('insert_seg6', line)){
      writeLines(gsub('insert_seg6', paste(clade, filename_base,"PB1", sep="_"), line), f)
    }else if (grepl('insert_seg7', line)){
      writeLines(gsub('insert_seg7', paste(clade, filename_base,"PB2", sep="_"), line), f)
    }else if (grepl('insert_seg8', line)){
      writeLines(gsub('insert_seg8', paste(clade, filename_base,"PA", sep="_"), line), f)
      
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
    } else if (grepl('insert_ratetimes', line)) {
      writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
      
      # mask all the non isConstant entries, i.e. isinfectedSkyline
    } else if (grepl('isconstant-->', line)) {
      writeLines(gsub('isconstant-->', '', line), f)
    } else if (grepl('<!--isconstant', line)){
      writeLines(gsub('<!--isconstant', '', line), f)
      # } else if (grepl('isinfectedSkyline-->', line)) {
      #   writeLines(gsub('isinfectedSkyline-->', '', line), f)
      # } else if (grepl('<!--isinfectedSkyline', line)){
      #   writeLines(gsub('<!--isinfectedSkyline', '', line), f)      
      # } else if (grepl('isNeSkyline-->', line)) {
      #   writeLines(gsub('isNeSkyline-->', '', line), f)
      # } else if (grepl('<!--isNeSkyline', line)){
      #   writeLines(gsub('<!--isNeSkyline', '', line), f)  
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
      seq1 = seqinr::read.fasta(file = paste("./data/",filename_base, "_HA.fasta", sep=""), seqtype = "DNA")
      seq2 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NA.fasta", sep=""), seqtype = "DNA")
      seq3 = seqinr::read.fasta(file = paste("./data/",filename_base, "_MP.fasta", sep=""), seqtype = "DNA")
      seq4 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NS.fasta", sep=""), seqtype = "DNA")
      seq5 = seqinr::read.fasta(file = paste("./data/",filename_base, "_NP.fasta", sep=""), seqtype = "DNA")
      seq6 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB1.fasta", sep=""), seqtype = "DNA")
      seq7 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PB2.fasta", sep=""), seqtype = "DNA")
      seq8 = seqinr::read.fasta(file = paste("./data/",filename_base, "_PA.fasta", sep=""), seqtype = "DNA")
      
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
    }else if (grepl('insert_independent_after', line)) {
      writeLines(gsub('insert_independent_after"', paste0(independentafter, '"'), line), f)
    }else if (grepl('insert_from_to_prior', line)) {
      writeLines('\t\t\t\t\t<distribution spec="Prior" x="@effectSize">\n', f)
      writeLines('\t\t\t\t\t\t<distr spec="beast.base.inference.distribution.Normal" mean="1" sigma="1"/>\n', f)
      writeLines('\t\t\t\t\t</distribution>\n', f)
    }else if (grepl('insert_double_discount', line)){
      writeLines(gsub('insert_double_discount', 'false', line), f)
      
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)
  
  p=ggplot(plot_data, aes(x=time, y=positivity, color=type)) +
    geom_line() +
    labs(title = "Smoothed Case Positivity", x = "Date", y = "Positivity") +
    theme_minimal() +
    facet_wrap(~type, scales = "free_y") +
    geom_vline(xintercept = max_date-rateshiftvals[independentafter], linetype="dashed", color = "red") +
    xlim(c(2020,2025.5))
  plot(p)
  
  # Create flyway-level plot
  # Aggregate cases by flyway and type
  flyway_data <- data.frame()
  
  for (d in max_date-rateshiftvals){
    # get all instances within that time window
    window = cases[cases$decimal_date >= d-diff & cases$decimal_date <= d+diff, ]
    
    # get LPAI and HPAI cases by flyway
    lpai_window = window[window$Final_H5 == "Detected" & window$Final_Pathogenicity != "High Path AI", ]
    hpai_window = window[window$Final_H5 == "Detected" & window$Final_Pathogenicity == "High Path AI", ]
    
    lpai_by_flyway = data.frame()
    hpai_by_flyway = data.frame()
    
    if(nrow(lpai_window) > 0) {
      lpai_by_flyway = aggregate(rep(1, nrow(lpai_window)), 
                                 by=list(flyway=lpai_window$flyway), 
                                 FUN=sum, na.rm=TRUE)
    }
    
    if(nrow(hpai_window) > 0) {
      hpai_by_flyway = aggregate(rep(1, nrow(hpai_window)), 
                                 by=list(flyway=hpai_window$flyway), 
                                 FUN=sum, na.rm=TRUE)
    }
    
    # Add type labels and combine
    if(nrow(lpai_by_flyway) > 0) {
      lpai_by_flyway$type <- "LPAI"
      names(lpai_by_flyway)[2] <- "count"
      flyway_data <- rbind(flyway_data, data.frame(date=d, lpai_by_flyway))
    }
    
    if(nrow(hpai_by_flyway) > 0) {
      hpai_by_flyway$type <- "HPAI"
      names(hpai_by_flyway)[2] <- "count"
      flyway_data <- rbind(flyway_data, data.frame(date=d, hpai_by_flyway))
    }
  }
  
  # Create the flyway plot
  if(nrow(flyway_data) > 0) {
    p_flyway <- ggplot(flyway_data, aes(x=date, y=count, color=type)) +
      geom_line() +
      labs(title = "LPAI and HPAI Cases by Flyway", x = "Date", y = "Case Count") +
      theme_minimal() +
      facet_wrap(~flyway, scales = "free_y") +
      geom_vline(xintercept = max_date-rateshiftvals[independentafter], linetype="dashed", color = "red") +
      xlim(c(2020,2025.5)) +
      scale_color_manual(values = c("LPAI" = "blue", "HPAI" = "red"))
    
    plot(p_flyway)
  }
  
  # Create overlap prevalence plot by flyway
  overlap_data <- data.frame()
  
  for (d in max_date-rateshiftvals){
    # get all instances within that time window
    window = cases[cases$decimal_date >= d-diff & cases$decimal_date <= d+diff, ]
    
    # get LPAI and HPAI cases by flyway
    lpai_window = window[window$Final_H5 == "Detected" & window$Final_Pathogenicity != "High Path AI", ]
    hpai_window = window[window$Final_H5 == "Detected" & window$Final_Pathogenicity == "High Path AI", ]
    
    if(nrow(lpai_window) > 0 && nrow(hpai_window) > 0) {
      # Calculate total samples per flyway in this window
      total_samples_by_flyway = aggregate(rep(1, nrow(window)), by=list(flyway=window$flyway), FUN=sum)
      
      # Calculate LPAI prevalence by flyway
      lpai_counts = aggregate(rep(1, nrow(lpai_window)), by=list(flyway=lpai_window$flyway), FUN=sum)
      names(lpai_counts)[2] = "lpai_count"
      lpai_prevalence = merge(total_samples_by_flyway, lpai_counts, by="flyway", all.x=TRUE)
      names(lpai_prevalence)[2] = "total_count"
      lpai_prevalence$lpai_count[is.na(lpai_prevalence$lpai_count)] = 0
      lpai_prevalence$lpai_prevalence = lpai_prevalence$lpai_count / lpai_prevalence$total_count
      
      # Calculate HPAI prevalence by flyway
      hpai_counts = aggregate(rep(1, nrow(hpai_window)), by=list(flyway=hpai_window$flyway), FUN=sum)
      names(hpai_counts)[2] = "hpai_count"
      hpai_prevalence = merge(total_samples_by_flyway, hpai_counts, by="flyway", all.x=TRUE)
      names(hpai_prevalence)[2] = "total_count"
      hpai_prevalence$hpai_count[is.na(hpai_prevalence$hpai_count)] = 0
      hpai_prevalence$hpai_prevalence = hpai_prevalence$hpai_count / hpai_prevalence$total_count
      
      # Calculate overlap prevalence (LPAI prevalence * HPAI prevalence) for each flyway
      merged_prevalence = merge(lpai_prevalence[,c("flyway", "lpai_prevalence")], 
                               hpai_prevalence[,c("flyway", "hpai_prevalence")], 
                               by="flyway", all.x=TRUE)
      merged_prevalence$hpai_prevalence[is.na(merged_prevalence$hpai_prevalence)] = 0
      merged_prevalence$overlap_prevalence = merged_prevalence$lpai_prevalence * merged_prevalence$hpai_prevalence
      
      # Add to plot data
      for(i in 1:nrow(merged_prevalence)) {
        overlap_data <- rbind(overlap_data, data.frame(
          date = d,
          flyway = merged_prevalence$flyway[i],
          overlap_prevalence = merged_prevalence$overlap_prevalence[i]
        ))
      }
    }
  }
  
  # Create the overlap prevalence plot
  if(nrow(overlap_data) > 0) {
    p_overlap <- ggplot(overlap_data, aes(x=date, y=overlap_prevalence, color=flyway)) +
      geom_line() +
      labs(title = "Overlap Prevalence by Flyway (LPAI × HPAI)", x = "Date", y = "Overlap Prevalence") +
      theme_minimal() +
      facet_wrap(~flyway, scales = "free_y") +
      geom_vline(xintercept = max_date-rateshiftvals[independentafter], linetype="dashed", color = "red") +
      xlim(c(2020,2025.5)) +
      scale_color_viridis_d()
    
    # Save the overlap prevalence plot
    ggsave("overlap_prevalence_by_flyway.png", p_overlap, width = 12, height = 8, dpi = 300)
    plot(p_overlap)
  }



}
  

files = list.files("xmls", pattern="*.rep0.xml", full.names=TRUE)


for (file in files) {
  # Check if filename contains 'constant'
  if (grepl("skygrowth", basename(file))) {
    # For files with 'constant', create all 10 copies (rep0-rep9)
    file.copy(file, gsub("rep0", "rep1", file))
    file.copy(file, gsub("rep0", "rep2", file))
    file.copy(file, gsub("rep0", "rep3", file))
    file.copy(file, gsub("rep0", "rep4", file))
    file.copy(file, gsub("rep0", "rep5", file))
    file.copy(file, gsub("rep0", "rep6", file))
    file.copy(file, gsub("rep0", "rep7", file))
    file.copy(file, gsub("rep0", "rep8", file))
    file.copy(file, gsub("rep0", "rep9", file))
  } else {
    # For files without 'constant', only create 2 copies (rep1 and rep2)
    file.copy(file, gsub("rep0", "rep1", file))
    file.copy(file, gsub("rep0", "rep2", file))
  }
}

# remove anything that doens't contain HLHxNx_
files = list.files("xmls", pattern="*", full.names=TRUE)
files = files[!grepl("HLHxNx_", files)]
# remove the files that are not in the list of files

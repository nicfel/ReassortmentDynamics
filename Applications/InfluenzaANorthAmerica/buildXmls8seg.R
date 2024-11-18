library(stringr)
library(seqinr)
library(ggplot2)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# build a wgs xml file for each virus and year between 2015 and 2018
virus = c('H1N1', 'H3N2')
year = seq(2015,2023,1)

# delete all xml files in xmls  
file.remove(list.files("./xmls", pattern="*.xml", full.names=TRUE))

# define the root height for Ne and reassortment variant rates
timediff = 5

# read in the who data
who_NA = read.csv("./whoData/NA_who_sentinel.csv", header = TRUE)
# remove all rows for which COUNTRY_CODE is not USA
who_NA = who_NA[who_NA$COUNTRY_CODE == "USA",]
who_NA$date = as.Date(who_NA$ISO_SDATE)

# loop over all combinations of files
for (a in 1:length(virus)) {
  for (b in 1:length(year)){
    base = paste(virus[[a]], "_", year[[b]], sep="")
    # check if the fasta file exists otherwise next
    if (!file.exists(paste("./xmls/", base, "_HA.fasta", sep=""))) {
      next
    }
    
    fasta1 <- seqinr::read.fasta(file = paste("./xmls/", base, "_HA.fasta", sep=""), seqtype = "DNA")
    isolates = names(fasta1)
    # require at least 100 isolates
    if (length(isolates) < 50) {
      next
    }
    
    peak = who_NA[who_NA$date > as.Date(paste(year[[b]], "06", "01", sep="-")) &
                    who_NA$date < as.Date(paste(year[[b]]+1, "05", "31", sep="-")),]
    if (virus[[a]] == 'H1N1') {
      peak$cases = peak$AH1N12009
    } else {
      peak$cases = peak$AH3
    }
    peak_date = peak[which.max(peak$cases),"date"]
    peak_prior = peak[peak$date <= peak_date,]
    # get the last time, befor cases where 1/20 of the peak
    start = peak_prior[peak_prior$cases < peak[which.max(peak$cases),"cases"]/8,]
    start_time = start[length(start$date),"date"]
    if (length(start_time)==0){
      start_time = peak_date-3*365
    }
    
    dates = c()
    for (i in seq(1, length(isolates))) {
      date = strsplit(isolates[i], "\\|")[[1]][[3]]
      dates = c(dates, as.Date(date))
    }  
    
    date_vals = data.frame(dates=as.Date(dates))

    independentAfterTime = max(min(dates), start_time)
    timediff = (max(dates)-independentAfterTime)/365

    print(timediff)
    
    rateshiftvals = unique(c(seq(0, 1, length.out=20), 1.5, 3, 5))
    rateshiftvals2 = rateshiftvals

    # find the index of the first rateshiftvals2 > timediff
    independentIndex = which(rateshiftvals2 > timediff)[1]
    
    
    p1=ggplot(data=peak, aes(x=date, y=cases)) + geom_line() + geom_point() + theme_minimal() + 
      geom_vline(xintercept = peak_date, linetype="dashed")+
      geom_vline(xintercept = start_time, linetype="dashed")+
      geom_vline(xintercept =  min(dates), linetype="dashed")+
      geom_vline(xintercept = max(dates)-rateshiftvals2[independentIndex+1]*365, color="red")+
      ggtitle(paste(virus[[a]], " ", year[[b]], " peak date: ", peak_date, " start date: ", start_time))+
      coord_cartesian(xlim = c(min(peak$date), max(peak$date))) +
      geom_text(label=timediff, x=max(dates)-rateshiftvals2[independentIndex+1]*365, y=10, color="red")
    p2=ggplot(date_vals, aes(x=dates)) + geom_histogram(binwidth=7) + theme_minimal() + 
      geom_vline(xintercept = peak_date, linetype="dashed")+
      geom_vline(xintercept = start_time, linetype="dashed")+
      geom_vline(xintercept =  min(dates), linetype="dashed")+
      ggtitle(paste(virus[[a]], " ", year[[b]], " peak date: ", peak_date, " start date: ", start_time))+
      coord_cartesian(xlim = c(min(peak$date), max(peak$date)))
    
    library(gridExtra)
    p = grid.arrange(p1, p2, ncol=1)
    plot(p)
    

    # save plot p to samplingFigs/
    ggsave(plot=p, filename=paste("./samplingFigs/", base, "_sampling.pdf", sep=""), width=6, height=5)
    
    # read in the log file to get the network height
    filename = paste(virus[[a]], "_", year[[b]] , ".constant", sep="")
    # build an inference xml files
    f <- file(sprintf('xmls/%s.rep0.xml',filename), 'w')
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
        if (virus[[a]] == 'H1N1') {
          writeLines(gsub('insert_clock_rate', '0.0025', line), f)
        } else {
          writeLines(gsub('insert_clock_rate', '0.0025', line), f)
        }
        
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', paste(base, "HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(base, "NA", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(base, "MP", sep="_"), line), f)
      }else if (grepl('insert_seg4', line)){
        writeLines(gsub('insert_seg4', paste(base, "NS", sep="_"), line), f)
      }else if (grepl('insert_seg5', line)){
        writeLines(gsub('insert_seg5', paste(base, "NP", sep="_"), line), f)
      }else if (grepl('insert_seg6', line)){
        writeLines(gsub('insert_seg6', paste(base, "PB1", sep="_"), line), f)
      }else if (grepl('insert_seg7', line)){
        writeLines(gsub('insert_seg7', paste(base, "PB2", sep="_"), line), f)
      }else if (grepl('insert_seg8', line)){
        writeLines(gsub('insert_seg8', paste(base, "PA", sep="_"), line), f)
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
        dates = c()
        for (i in seq(1, length(isolates))) {
          date = strsplit(isolates[i], "\\|")[[1]][[3]]
          
          if (i == length(isolates)) {
            value = paste(value, isolates[i], '=', date, sep = '')
          } else {
            value = paste(value, isolates[i], '=', date, sep = '')
            value = paste(value, ',', sep = '')
          }
          dates = c(dates, as.Date(date))
        }  
        writeLines(gsub('insert_heights', value, line), f)   
      }else if (grepl('insert_EOS', line)){
        # get the most recent sampling time
        # EOS = max(as.numeric(strsplit(isolates[length(isolates)], "\\|")[[1]][[3]]))
        # # get the peak of that season
        # EOS = max(dates) - as.numeric(start_time)
        pretime = rateshiftvals2
        writeLines(gsub('insert_EOS', paste(pretime, collapse = " "), line), f)
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        seq1 = seqinr::read.fasta(file = paste("./xmls/", base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = seqinr::read.fasta(file = paste("./xmls/", base, "_NA.fasta", sep=""), seqtype = "DNA")
        seq3 = seqinr::read.fasta(file = paste("./xmls/", base, "_MP.fasta", sep=""), seqtype = "DNA")
        seq4 = seqinr::read.fasta(file = paste("./xmls/", base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq5 = seqinr::read.fasta(file = paste("./xmls/", base, "_NP.fasta", sep=""), seqtype = "DNA")
        seq6 = seqinr::read.fasta(file = paste("./xmls/", base, "_PB1.fasta", sep=""), seqtype = "DNA")
        seq7 = seqinr::read.fasta(file = paste("./xmls/", base, "_PB2.fasta", sep=""), seqtype = "DNA")
        seq8 = seqinr::read.fasta(file = paste("./xmls/", base, "_PA.fasta", sep=""), seqtype = "DNA")
        
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
        writeLines(gsub('insert_independent_after', as.character(independentIndex), line), f)
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)
    
    # make a second xml where the .trees is replaced by .infected
    filename = paste(virus[[a]], "_", year[[b]] , ".variable", sep="")
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
        if (virus[[a]] == 'H1N1') {
          writeLines(gsub('insert_clock_rate', '0.003', line), f)
        } else {
          writeLines(gsub('insert_clock_rate', '0.003', line), f)
        }
        
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', paste(base, "HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(base, "NA", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(base, "MP", sep="_"), line), f)
      }else if (grepl('insert_seg4', line)){
        writeLines(gsub('insert_seg4', paste(base, "NS", sep="_"), line), f)
      }else if (grepl('insert_seg5', line)){
        writeLines(gsub('insert_seg5', paste(base, "NP", sep="_"), line), f)
      }else if (grepl('insert_seg6', line)){
        writeLines(gsub('insert_seg6', paste(base, "PB1", sep="_"), line), f)
      }else if (grepl('insert_seg7', line)){
        writeLines(gsub('insert_seg7', paste(base, "PB2", sep="_"), line), f)
      }else if (grepl('insert_seg8', line)){
        writeLines(gsub('insert_seg8', paste(base, "PA", sep="_"), line), f)
      } else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
        
      }else if (grepl('insert_EOS', line)){
        # get the most recent sampling time
        writeLines(gsub('insert_EOS', paste(rateshiftvals2, collapse = " "), line), f)
        
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
        seq1 = seqinr::read.fasta(file = paste("./xmls/", base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = seqinr::read.fasta(file = paste("./xmls/", base, "_NA.fasta", sep=""), seqtype = "DNA")
        seq3 = seqinr::read.fasta(file = paste("./xmls/", base, "_MP.fasta", sep=""), seqtype = "DNA")
        seq4 = seqinr::read.fasta(file = paste("./xmls/", base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq5 = seqinr::read.fasta(file = paste("./xmls/", base, "_NP.fasta", sep=""), seqtype = "DNA")
        seq6 = seqinr::read.fasta(file = paste("./xmls/", base, "_PB1.fasta", sep=""), seqtype = "DNA")
        seq7 = seqinr::read.fasta(file = paste("./xmls/", base, "_PB2.fasta", sep=""), seqtype = "DNA")
        seq8 = seqinr::read.fasta(file = paste("./xmls/", base, "_PA.fasta", sep=""), seqtype = "DNA")
        
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
        timediff = (max(dates)-independentAfterTime)/365
        # find the index of the first rateshiftvals2 > timediff
        index = which(rateshiftvals2 > timediff)[1]
        writeLines(gsub('insert_independent_after', as.character(independentIndex), line), f)
      }else if (grepl('insert_cases_cases', line)) {
        weekly_cases = peak$cases
        week_start_offset = (as.Date(max(dates)) - peak$date)/365
        
        # get the first index where week_start_offset is smaller than 0
        index = which(week_start_offset <= 0.0)[1]
        # log standardize the cases
        weekly_cases = log(c(weekly_cases[index:1]+1, 1))
        weekly_cases = (weekly_cases - mean(weekly_cases))/sd(weekly_cases)
        week_start_offset = week_start_offset
        week_start_offset[index]=0
        week_start_offset = c(week_start_offset[index:1], rateshiftvals2[length(rateshiftvals2)])
        print(week_start_offset)
        if (length(week_start_offset)!=length(weekly_cases)){
          print("error")
        }
        
        writeLines(gsub('insert_cases_cases', paste(weekly_cases, collapse=" "), line), f)
      } else if (grepl('insert_cases_time', line)){
        writeLines(gsub('insert_cases_time', paste(week_start_offset, collapse=" "), line), f)
      } else if (grepl('<stateNode id="InfectedToRho" spec="RealParameter"', line)){
        writeLines('\t\t\t\t<stateNode id="InfectedToRho" spec="RealParameter" lower="-6" value="-1 -1.5"/>', f)
      
      } else {
        writeLines(line, f)
      }
    }
    close(f)
    close(template)

    # skip this file
    next
    
    # make a second xml where the .trees is replaced by .infected
    filename = paste(virus[[a]], "_", year[[b]] , ".ne", sep="")
    f <- file(sprintf('xmls/%s.xml', filename), 'w')
    # Open the template file
    template <- file('../H5N1NorthAmerica/inference_template.xml', 'r')
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
          writeLines(gsub('insert_clock_rate', '0.003', line), f)
        } else {
          writeLines(gsub('insert_clock_rate', '0.003', line), f)
        }
        
      }else if (grepl('insert_seg1', line)){
        writeLines(gsub('insert_seg1', paste(base, "HA", sep="_"), line), f)
      }else if (grepl('insert_seg2', line)){
        writeLines(gsub('insert_seg2', paste(base, "NA", sep="_"), line), f)
      }else if (grepl('insert_seg3', line)){
        writeLines(gsub('insert_seg3', paste(base, "MP", sep="_"), line), f)
      }else if (grepl('insert_seg4', line)){
        writeLines(gsub('insert_seg4', paste(base, "NS", sep="_"), line), f)
      }else if (grepl('insert_seg5', line)){
        writeLines(gsub('insert_seg5', paste(base, "NP", sep="_"), line), f)
      }else if (grepl('insert_seg6', line)){
        writeLines(gsub('insert_seg6', paste(base, "PB1", sep="_"), line), f)
      }else if (grepl('insert_seg7', line)){
        writeLines(gsub('insert_seg7', paste(base, "PB2", sep="_"), line), f)
      }else if (grepl('insert_seg8', line)){
        writeLines(gsub('insert_seg8', paste(base, "PA", sep="_"), line), f)
      }else if (grepl('insert_times', line)) {
        writeLines(gsub('insert_times', paste(rateshiftvals, collapse=' '), line), f)
        
      } else if (grepl('insert_ratetimes', line)) {
        writeLines(gsub('insert_ratetimes', paste(rateshiftvals2, collapse=' '), line), f)
        
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
        dates =c()
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
        # compute EOS as one month prior to the most distant sample
        EOS = max(dates) - as.numeric(start_time)
        writeLines(gsub('insert_EOS', EOS/365, line), f)
        
      }else if (grepl('insert_weights', line)) {
        # get the length of both alignments
        # get the length of both alignments
        seq1 = seqinr::read.fasta(file = paste("./xmls/", base, "_HA.fasta", sep=""), seqtype = "DNA")
        seq2 = seqinr::read.fasta(file = paste("./xmls/", base, "_NA.fasta", sep=""), seqtype = "DNA")
        seq3 = seqinr::read.fasta(file = paste("./xmls/", base, "_MP.fasta", sep=""), seqtype = "DNA")
        seq4 = seqinr::read.fasta(file = paste("./xmls/", base, "_NS.fasta", sep=""), seqtype = "DNA")
        seq5 = seqinr::read.fasta(file = paste("./xmls/", base, "_NP.fasta", sep=""), seqtype = "DNA")
        seq6 = seqinr::read.fasta(file = paste("./xmls/", base, "_PB1.fasta", sep=""), seqtype = "DNA")
        seq7 = seqinr::read.fasta(file = paste("./xmls/", base, "_PB2.fasta", sep=""), seqtype = "DNA")
        seq8 = seqinr::read.fasta(file = paste("./xmls/", base, "_PA.fasta", sep=""), seqtype = "DNA")
        
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
  }
}


# for all *.rep0.xml files in xmls, make a copy with the same name but .rep1.xml
# and .rep2.xml
for (file in list.files("./xmls", pattern="*.rep0.xml", full.names=TRUE)) {
  file.copy(file, gsub("rep0", "rep1", file))
  file.copy(file, gsub("rep0", "rep2", file))
}
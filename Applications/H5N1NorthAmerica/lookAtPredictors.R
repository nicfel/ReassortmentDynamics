library(stringr)
library(seqinr)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)




# read in all fasta files in data
fasta = list.files("./xmls/", pattern="*.fasta", full.names=TRUE)
# only keep the ones starting with HLHxNx
fasta = fasta[startsWith(basename(fasta), "HL")]
for (fastafile in fasta){
  # read in the fasta file
  fasta1 = seqinr::read.fasta(file = fastafile, seqtype = "DNA")
  # remove all sequences not in isolates
  # write the fasta file to the xmls3seg directory
  seqinr::write.fasta(sequences = fasta1, 
                      names = names(fasta1), 
                      file.out = paste0("./xmls/", clade, "_", basename(fastafile)))
}


# read in the first fasta file and keep all the names
isolates = names(fasta1)
# collect the dates as the last group after splitting on | 
dates = sapply(isolates, function(x) strsplit(x, "\\|")[[1]][[3]])
min = min(as.Date(dates))
max = max(as.Date(dates))
first_intro = as.numeric(max-min)/365

rateshiftvals = c(seq(0,  first_intro, length.out=40))
rateshiftvals = unique(rateshiftvals)
rateshiftvals2 = rateshiftvals


cases = read.csv("./tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv")
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
# smoothing area
diff = min(diff(rateshiftvals))/2
for (d in max_date-rateshiftvals){
  # get all instances within that time window
  window = cases[cases$decimal_date >= d-diff & cases$decimal_date <= d+diff, ]
  window_alt = cases[cases$decimal_date >= d-diff-14/365 & cases$decimal_date <= d+diff-14/365, ]
  
  
  # get how many instances of Final_IAV are Detected
  total_AIV = sum(window$Final_IAV == "Detected", na.rm=TRUE)
  total_H5 = sum(window$Final_H5 == "Detected", na.rm=TRUE)
  
  # get the number of high path cases
  total_HPAI = sum(window$Final_H5 == "Detected" &  window$Final_Pathogenicity == "High Path AI", na.rm=TRUE)
  
  lpai = window[window$Final_H5 == "Detected" &  window$Final_Pathogenicity != "High Path AI", ]
  hpai = window[window$Final_H5 == "Detected" &  window$Final_Pathogenicity == "High Path AI", ]
  
  lpai_alt = window_alt[window_alt$Final_H5 == "Detected" &  window_alt$Final_Pathogenicity != "High Path AI", ]
  hpai_alt = window_alt[window_alt$Final_H5 == "Detected" &  window_alt$Final_Pathogenicity == "High Path AI", ]
  
  # check for states that are in both
  if (nrow(lpai) > 0 && nrow(hpai) > 0){
    overlap = length(intersect(lpai$County, hpai$County))
  }else{
    overlap= 0
  }
  
  smoothed_case_data = rbind(smoothed_case_data, data.frame(
    date = d,
    positivity = (total_H5-total_HPAI),
    type = "lpai_nosummer"
  ))
  
  smoothed_case_data = rbind(smoothed_case_data, data.frame(
    date = d,
    positivity = (total_H5-total_HPAI),
    type = "h5_lpai"
  ))
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

  smoothed_case_data = rbind(smoothed_case_data, data.frame(
    date = d,
    positivity = overlap,
    type = "overlap"
  ))
  
}


# set the values for summer 2022 to min value for the no_summer predictor
smoothed_case_data$positivity[smoothed_case_data$type == "lpai_nosummer" & smoothed_case_data$date >= decimal_date(as.Date("2022-07-01")) & smoothed_case_data$date <= decimal_date(as.Date("2023-03-01"))] = 0
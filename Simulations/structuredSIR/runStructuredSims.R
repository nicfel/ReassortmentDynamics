library(stringr)
# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# 
# # Make a directory to store the xml files
# if (dir.exists("master")) {
#   unlink("master", recursive = TRUE)
# }
# dir.create("master")
# 
# # Initialize the parameters for the SIR simulations
# param_file <- file('SIR_simulations.txt', 'w')
# 
# # define the number of states
# nstates = 20
# 
# # Initialize the parameters for the SIR simulations
# param_file <- file('structuredSIR_simulations.txt', 'w')
# 
# # write one header with a different rate for each state from 1...nstates
# 
# writeLines("transmission\trecovery\tpopulation_size\td", param_file)
# recovery_rate <- 1
# 
# # Loop over 100 runs
# for (i in 1:100) {
#   print(i)
#   # Make a new xml file for each run
#   f <- file(sprintf('master/structuredSIR_simulations_%d.xml', i), 'w')
#   # Open the template file
#   template <- file('structuredSimulation_template.xml', 'r')
# 
#   # Randomly sample one transmission rate from a lognormal distribution with mean in real space of 1.5 and S 0.25
#   transmission <- rlnorm(1, meanlog = 0.4, sdlog = 0.25)
#   while (transmission < 1.01) {
#     transmission <- rlnorm(1, meanlog = 0.4, sdlog = 0.25)
#   }
# 
#   # Randomly sample 10 sampling rates from a lognormal distribution with mean in real space of 0.01 and S 0.5
#   sampling <- rlnorm(nstates, meanlog = -3.605170185988091, sdlog = 0.5)
#   # Randomly sample 10 population sizes from a lognormal distribution with mean 10000 and S 0.5
#   population_size <- round(rlnorm(nstates, meanlog = 6.510340371976182, sdlog = 0.5))
#   # Randomly sample the k of the negative binomial distribution from a lognormal distribution with mean 1 and S 0.5
#   k <- rlnorm(1, meanlog = 0, sdlog = 1)
#   # Randomly sample nstates*(nstates-1) migration rates from a log normal with mean in real space of 0.1 and S of 0.25
#   if (i>50){
#     migration <- rlnorm(nstates*(nstates-1), meanlog = -2.5, sdlog = 1)
#   }else{
#     migration <- rlnorm(nstates*(nstates-1), meanlog = -4.5, sdlog = 1)
#   }
# 
#   # replace transmission, recover and k and waning (just 0's) nstates times
#   transmission = rep(transmission, nstates)
#   recovery = rep(recovery_rate, nstates)
#   waning = rep(0,nstates)
# 
#   # Write the parameters to the file
#   cat(sprintf('%f\t%f\t%f\t%f\n', transmission[[1]], recovery_rate[[1]], sum(population_size), k), file=param_file)
# 
#   # Write the parameters to the xml file
#   while (length(line <- readLines(template, n = 1)) > 0) {
#     if(grepl('insert_transmission', line)) {
#       writeLines(gsub('insert_transmission', paste(transmission, collapse=" "), line), f)
#     } else if (grepl('insert_recovery', line)) {
#       writeLines(gsub('insert_recovery', paste(recovery, collapse=" "), line), f)
#     } else if (grepl('insert_sampling', line)) {
#       writeLines(gsub('insert_sampling', paste(sampling, collapse=" "), line), f)
#     } else if (grepl('insert_population_size', line)) {
#       writeLines(gsub('insert_population_size', paste(population_size, collapse=" "), line), f)
#     } else if (grepl('insert_migration', line)) {
#       writeLines(gsub('insert_migration', paste(migration, collapse=" "), line), f)
#     } else if (grepl('insert_waning', line)) {
#       writeLines(gsub('insert_waning', paste(waning, collapse=" "), line), f)
#     }else if (grepl('insert_k', line)){
#       writeLines(gsub('insert_k', as.character(k), line), f)
#     } else {
#       writeLines(line, f)
#     }
#   }
#   close(f)
#   close(template)
# 
# 
#   # Run the xml using BEAST and the system command, while preventing any logging to screen
#   system(sprintf('/Applications/BEAST\\ 2.7.6/bin/beast -seed %d -overwrite master/structuredSIR_simulations_%d.xml', i, i))
# }
# close(param_file)

# Read tree files and simulate sequence alignment
library(ape)
library(phytools)


# Make a directory to store the xml files
if (dir.exists("xmls")) {
  unlink("xmls", recursive = TRUE)
}
dir.create("xmls")

trees <- list.files(path = "master", pattern = "\\.trees$", full.names = TRUE)
for (tree_file in trees) {
  lines <- readLines(tree_file)
  network <- lines[1]
  for (j in 1:2) {
    tree <- read.tree(text = lines[j + 1])

    # replace all the name to deal with the seq gen same length concat
    for (k in 1:length(tree$tip.label)) {
      tree$tip.label[k] <- gsub("sample", "s",tree$tip.label[k])
    }
    
    # Write the tree to a file in NEWICK format for Seq-Gen
    tree_filename <- sprintf('master/temp_tree_%d.newick', j)
    write.tree(tree, file = tree_filename)
    
    # Call Seq-Gen to simulate sequences
    output_filename <- sprintf('xmls/%s_%d.nexus', basename(tree_file), j)
    
    seq_gen_command <- sprintf("./seq-gen -mHKY -on -s0.001 < %s > %s", tree_filename, output_filename)
    system(seq_gen_command)

    # Optionally, you can remove the temporary tree file
    file.remove(tree_filename)

    # read in the nexus files line by line and remove all z's and replace them by a \t
    nex_lines = readLines(output_filename)
    nex_lines = gsub("s_", "sample_", nex_lines)
    # replace the fasta file
    write(nex_lines, file = output_filename)
  }

  # get the heights of all the tips in the tree
  tree <- read.tree(text = lines[2])
  # calculate the distances between all tips from the root node
  dists <- nodeHeights(tree)
  # get the distance of each tip from the root from the dists matrix
  # using the index of each tip and the root
  heights = c()
  for (i in seq(1, length(tree$tip.label))) {
    heights <- c(heights, dists[tree$edge[,2]==i, 2])
  } 
  
  heights = abs(heights-max(heights))

  # close all connections
  closeAllConnections()
  
  # read in the log file to get the network height
  filename = gsub(".trees", ".log", basename(tree_file))  
  lines = readLines(sprintf('master/%s',filename))
  data = strsplit(lines[[2]], split="\t")[[1]]
  lins <- str_split(str_replace_all(data[[3]], "\\[|\\]", ""), ', ')[[1]]
  
  timediff = as.numeric(strsplit(lins[[1]], split=":")[[1]][[2]]) - 
    as.numeric(strsplit(lins[[length(lins)]], split=":")[[1]][[2]])

  # get the basename of the xml by replacing the .trees. with .constant.
  filename = gsub(".trees", ".constant", basename(tree_file))  
  # build an inference xml files
  f <- file(sprintf('xmls/%s.xml',filename), 'w')
  # Open the template file
  template <- file('inference_template.xml', 'r')

  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_name', line)) {
      writeLines(gsub('insert_name', sprintf('%s', basename(tree_file)), line), f)
    } else if (grepl('insert_tips', line)) {
      # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
      for (i in seq(1, length(tree$tip.label))) {
        writeLines(sprintf('<taxon spec="Taxon" id="%s"/>', tree$tip.label[i]), f)
      }  
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(seq(0,timediff*1.1, length.out=15), collapse=' '), line), f)

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
      # write the height of each tip using heights seperated by a = and , between
      value=''
      for (i in seq(1, length(tree$tip.label))) {
        if (i == length(tree$tip.label)) {
          value = paste(value, tree$tip.label[i], '=', heights[i], sep = '')
        } else {
          value = paste(value, tree$tip.label[i], '=', heights[i], sep = '')
          value = paste(value, ',', sep = '')
        }
      }  
      writeLines(gsub('insert_heights', value, line), f)      
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)

  # make a second xml where the .trees is replaced by .infected
  filename = gsub(".trees", ".variable", basename(tree_file))
  f <- file(sprintf('xmls/%s.xml', filename), 'w')
  # Open the template file
  template <- file('inference_template.xml', 'r')
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_name', line)) {
      writeLines(gsub('insert_name', sprintf('%s', basename(tree_file)), line), f)
    } else if (grepl('insert_tips', line)) {
      # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
      for (i in seq(1, length(tree$tip.label))) {
        writeLines(sprintf('<taxon spec="Taxon" id="%s"/>', tree$tip.label[i]), f)
      }
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(seq(0,timediff*1.1, length.out=15), collapse=' '), line), f)

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
      for (i in seq(1, length(tree$tip.label))) {
        if (i == length(tree$tip.label)) {
          value = paste(value, tree$tip.label[i], '=', heights[i], sep = '')
        } else {
          value = paste(value, tree$tip.label[i], '=', heights[i], sep = '')
          value = paste(value, ',', sep = '')
        }
      }
      writeLines(gsub('insert_heights', value, line), f)
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)
  
  
  # make a second xml where the .trees is replaced by .infected
  filename = gsub(".trees", ".ne", basename(tree_file))
  f <- file(sprintf('xmls/%s.xml', filename), 'w')
  # Open the template file
  template <- file('inference_template.xml', 'r')
  while (length(line <- readLines(template, n = 1)) > 0) {
    if (grepl('insert_name', line)) {
      writeLines(gsub('insert_name', sprintf('%s', basename(tree_file)), line), f)
    } else if (grepl('insert_tips', line)) {
      # write the name of each tip, adding a , for all but the last one to the end <taxon spec="Taxon" id="t1"/>
      for (i in seq(1, length(tree$tip.label))) {
        writeLines(sprintf('<taxon spec="Taxon" id="%s"/>', tree$tip.label[i]), f)
      }   
    } else if (grepl('insert_times', line)) {
      writeLines(gsub('insert_times', paste(seq(0,timediff*1.1, length.out=15), collapse=' '), line), f)
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
      for (i in seq(1, length(tree$tip.label))) {
        if (i == length(tree$tip.label)) {
          value = paste(value, tree$tip.label[i], '=', heights[i], sep = '')
        } else {
          value = paste(value, tree$tip.label[i], '=', heights[i], sep = '')
          value = paste(value, ',', sep = '')
        }
      }  
      writeLines(gsub('insert_heights', value, line), f)      
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(template)

}

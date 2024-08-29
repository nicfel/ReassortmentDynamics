# Clear workspace
rm(list=ls())

# Set random seed for reproducibility
set.seed(6465546)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# # Make a directory to store the xml files
# if (dir.exists("master")) {
#   unlink("master", recursive = TRUE)
# }
# dir.create("master")
# 
# # Initialize the parameters for the SIR simulations
# param_file <- file('SIR_simulations.txt', 'w')
# writeLines('transmission\trecovery\tsampling\tpopulation_size\tk', param_file)
# 
# recovery_rate <- 1
# 
# # Loop over 100 runs
# for (i in 1:100) {
#   # Make a new xml file for each run
#   f <- file(sprintf('master/SIR_simulations_%d.xml', i), 'w')
#   # Open the template file
#   template <- file('simulation_template.xml', 'r')
#   
#   # Randomly sample the transmission rate from a lognormal distribution with mean 3 and S 0.25
#   transmission <- rlnorm(1, meanlog = 1.0986122886681098, sdlog = 0.25)
#   while (transmission < 1.2) {
#     transmission <- rlnorm(1, meanlog = 1.0986122886681098, sdlog = 0.25)
#   }
#   
#   # Randomly sample the sampling rate from a lognormal distribution with mean 0.01 and S 0.5
#   sampling <- rlnorm(1, meanlog = -4.605170185988091, sdlog = 0.25)
#   # Randomly sample the population size from a lognormal distribution with mean 10000 and S 0.5
#   population_size <- round(rlnorm(1, meanlog = 9.5, sdlog = 0.1))
#   # Randomly sample the k of the negative binomial distribution from a lognormal distribution with mean 1 and S 0.5
#   k <- rlnorm(1, meanlog = 0, sdlog = 1)
# 
#   
#   # Write the parameters to the file
#   cat(sprintf('%f\t%f\t%f\t%f\t%f\n', transmission, recovery_rate, sampling, population_size, k), file=param_file)
#   
#   # Write the parameters to the xml file
#   if (i>50){
#     while (length(line <- readLines(template, n = 1)) > 0) {
#       if (grepl('spec="SIRwithReassortment"', line)) {
#         writeLines(gsub('spec="SIRwithReassortment"', 'spec="SuperspreadingSIRwithReassortment"', line), f)
#       } else if(grepl('insert_transmission', line)) {
#         writeLines(gsub('insert_transmission', as.character(transmission), line), f)
#       } else if (grepl('insert_recovery', line)) {
#         writeLines(gsub('insert_recovery', as.character(recovery_rate), line), f)
#       } else if (grepl('insert_sampling', line)) {
#         writeLines(gsub('insert_sampling', as.character(sampling), line), f)
#       } else if (grepl('insert_population_size', line)) {
#         writeLines(gsub('insert_population_size', as.character(population_size), line), f)
#         writeLines(gsub('populationSize="insert_population_size"', paste('k="', as.character(k), '"', sep=""), line), f)
#       } else {
#         writeLines(line, f)
#       }
#     }
#     close(f)
#     close(template)
#   }else{
#     while (length(line <- readLines(template, n = 1)) > 0) {
#       if (grepl('insert_transmission', line)) {
#         writeLines(gsub('insert_transmission', as.character(transmission), line), f)
#       } else if (grepl('insert_recovery', line)) {
#         writeLines(gsub('insert_recovery', as.character(recovery_rate), line), f)
#       } else if (grepl('insert_sampling', line)) {
#         writeLines(gsub('insert_sampling', as.character(sampling), line), f)
#       } else if (grepl('insert_population_size', line)) {
#         writeLines(gsub('insert_population_size', as.character(population_size), line), f)
#       } else {
#         writeLines(line, f)
#       }
#     }
#     close(f)
#     close(template)
#   }
#   
#   # Run the xml using BEAST and the system command
#   system(sprintf('/Applications/BEAST\\ 2.7.6/bin/beast -seed %d -overwrite master/SIR_simulations_%d.xml', i, i))
# }
# 
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
  for (j in 1:4) {
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
  
  rateshifts = seq(0,max(dists),length.out = 20)

  # close all connections
  closeAllConnections()
    
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
    }else if (grepl('insert_rateShifts', line)){
      writeLines(gsub('insert_rateShifts', paste(rateshifts, collapse = ' '), line), f)
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
  filename = gsub(".trees", ".infected", basename(tree_file))
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
    }else if (grepl('insert_rateShifts', line)){
      writeLines(gsub('insert_rateShifts', paste(rateshifts, collapse = ' '), line), f)
      
    # mask all the non isConstant entries, i.e. isinfectedSkyline
    } else if (grepl('isconstant-->', line)) {
      writeLines(gsub('isconstant-->', '', line), f)
    } else if (grepl('<!--isconstant', line)){
      writeLines(gsub('<!--isconstant', '', line), f)      
    } else if (grepl('isNeSkyline-->', line)) {
      writeLines(gsub('isNeSkyline-->', '', line), f)
    } else if (grepl('<!--isNeSkyline', line)){
      writeLines(gsub('<!--isNeSkyline', '', line), f)  
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
    }else if (grepl('insert_rateShifts', line)){
      writeLines(gsub('insert_rateShifts', paste(rateshifts, collapse = ' '), line), f)
      # mask all the non isConstant entries, i.e. isinfectedSkyline
    } else if (grepl('isconstant-->', line)) {
      writeLines(gsub('isconstant-->', '', line), f)
    } else if (grepl('<!--isconstant', line)){
      writeLines(gsub('<!--isconstant', '', line), f)      
    } else if (grepl('isISkyline-->', line)) {
      writeLines(gsub('isISkyline-->', '', line), f)
    } else if (grepl('<!--isISkyline', line)){
      writeLines(gsub('<!--isISkyline', '', line), f)  
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

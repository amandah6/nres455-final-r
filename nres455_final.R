setwd("C:/Users/amand/Downloads/IL_DATA")

install.packages("dplyr")
library("dplyr")

il.tree <- read.table("IL_TREE.csv", sep=",", header=TRUE)
il.cond <- read.table("IL_COND.csv", sep=",", header=TRUE)
il.plot <- read.table("IL_PLOT.csv", sep=",", header=TRUE)

tree.cn <- il.tree$CN
tree.plot.cn <- il.tree$PLT_CN
tree.invyr <- il.tree$INVYR
tree.species <-il.tree$SPCD
tree.cause <- il.tree$AGENTCD

# replace all NA or 0 entries with -999 (before 1999, 0 meant alive trees)
tree.cause.final <- ifelse(is.na(tree.cause) | tree.cause == 0, -999, tree.cause)

il.tree.df <- data.frame(tree.cn, tree.plot.cn, tree.invyr, tree.species, tree.cause.final)
il.tree.df.final <- na.exclude(il.tree.df)
print(nrow(il.tree.df) - nrow(il.tree.df.final))

head(il.tree.df.final)

# create condition data frame
cond.cn <- il.cond$CN
cond.plot.cn <- il.cond$PLT_CN
cond.invyr <- il.cond$INVYR
cond.canopy.cover <- il.cond$LIVE_CANOPY_CVR_PCT

il.cond.df <- data.frame(cond.cn, cond.plot.cn, cond.invyr, cond.canopy.cover)
# remove all rows without canopy cover information and arrange by sequence number
il.cond.df.final <- arrange(filter(il.cond.df, !is.na(cond.canopy.cover)), cond.cn)
print(nrow(il.cond.df) - nrow(il.cond.df.final))
head(il.cond.df.final)


# create plot data frame
plot.cn <- il.plot$CN
plot.invyr <- il.plot$INVYR
plot.manual <- il.plot$MANUAL
plot.prev <- il.plot$PREV_PLT_CN
plot.county <- il.plot$COUNTYCD
plot.unit <- il.plot$UNITCD
plot.plot <- il.plot$PLOT
plot.lat <- il.plot$LAT
plot.lon <- il.plot$LON
plot.status <- il.plot$PLOT_STATUS_CD
plot.rem <- il.plot$REMPER

il.plot.df.final <- data.frame(plot.cn, plot.invyr, plot.manual, plot.prev, plot.county, plot.unit, plot.plot, plot.lat, plot.lon, plot.status, plot.rem)
#print(nrow(il.plot.df) - nrow(il.plot.df.final))
head(il.plot.df.final)
# sum(is.na(il.plot.df.final$plot.prev))
# sum(duplicated(il.tree.df.final$tree.plot.cn))


tree.plot.cn.unique <- unique(il.tree.df.final$tree.plot.cn)
length(tree.plot.cn.unique)

# average number of trees per survey plot
treeperplot <- numeric(length(tree.plot.cn.unique))
for (i in 1:length(tree.plot.cn.unique)) {
  treeperplot[i] <- sum(il.tree.df.final$tree.plot.cn == tree.plot.cn.unique[i])
}
mean(treeperplot)

# using the tree.df, calculate the shannon diversity index for each plot and 
# add it as a column to il.plot.df.final

# join tree and plot table
# for each plot sequence number and unique year: (aggregate?)
# create a data frame of just trees in plot in 1985, 1998, 2010, 2023
years.keep <- c(2001, 2006, 2011, 2016, 2023)
#years.keep <- unique(il.plot.df.final$plot.invyr)
il.plot.df.years<- filter(il.plot.df.final, plot.invyr %in% years.keep)
il.tree.df.years <- filter(il.tree.df.final, tree.invyr %in% years.keep)

# https://www.geeksforgeeks.org/apply-function-to-each-row-in-r-dataframe/
# given a data frame, calculate the shannon index for each plot

# add a tree count column to il.plot.df.years
tree.counter <- function(x) {
  count <- sum(il.tree.df.years$tree.plot.cn == x[1])
  return(count)
}

plot.tree.counts = apply(il.plot.df.years, 1, tree.counter)
plot.tree.counts[1:20]

# plots in our desired years with plot.tree.counts column
il.plot.df.years.counts <- cbind(il.plot.df.years, plot.tree.counts)

# filter out all plots with no trees, but keep ones that were not resampled for a reason
il.plot.df.years.counts.filtered <- filter(il.plot.df.years.counts, plot.tree.counts > 0 | plot.status == 3)

# rename tree.plot.cn to plot.cn so that we can merge it correctly
colnames(il.tree.df.years)[2] <- "plot.cn"
# join valid plots with trees table
plottree <- merge(il.plot.df.years.counts.filtered, il.tree.df.years)
# remove tree.invyr since it's a duplicate of plot.invyr
plottree$tree.invyr <- NULL
# column index of tree.species
tree.species.idx <- which(colnames(plottree[1,]) == "tree.species")
tree.count.idx <- which(colnames(plottree[1,]) == "plot.tree.counts")

shannon.index <- function(x) {
  plottree.subset <- subset(plottree, plot.cn == x[1])
  if (nrow(plottree.subset) > 0) {
    # create a matrix, species.count, where the index is the species code and the value is the count
    species.list <- unique(plottree.subset$tree.species)
    species.count <- matrix(numeric(length(species.list)))
    rownames(species.count) <- species.list
    
    
    # count the number of trees in each species
    for (i in 1:nrow(plottree.subset)) {
      species.count[as.character(plottree.subset[i, tree.species.idx]),] <- species.count[as.character(plottree.subset[i, tree.species.idx]),] + 1
    }
    
    # divide all elements by t.count
    species.count.pct <- species.count/x[tree.count.idx]
    
    # sum up all p_i ln(p_i) for species.count.pct
    h <- sum(species.count.pct * log(species.count.pct)) * -1
    
    return(h)
  } else {
    return(NA)
  }
}


il.plot.shannon.index <- apply(il.plot.df.years.counts.filtered, 1, shannon.index)
il.plot.df.final.shannon <- cbind(il.plot.df.years.counts.filtered, il.plot.shannon.index)
il.plot.df.final.shannon
il.plot.df.final.shannon.export <- select(il.plot.df.final.shannon, plot.cn, plot.invyr, plot.lat, plot.lon, plot.plot, il.plot.shannon.index, )
write.csv(il.plot.df.final.shannon.export, "IL_PLOT_SHANNON_2001_2023.csv", row.names=FALSE)


# find the most common tree species by year
count <- function(g) {
  species.list <- unique(g)
  species.count <- matrix(numeric(length(species.list)))
  rownames(species.count) <- species.list
  print(species.list)
  for (i in 1:length(g)) {
    species.count[as.character(g[i])] <- species.count[as.character(g[i])] + 1
  }
  return(species.count)
}

# group the species column by year
year.groups <- split(il.tree.df.final$tree.species, il.tree.df.final$tree.invyr)
typeof(year.groups)
s.counts <- il.tree.df.final %>% group_by(tree.invyr, tree.species) %>% summarise(Count = n())
s.counts
typeof(s.counts)
arrange(filter(s.counts, tree.invyr==1985), desc(Count))

# create a table of number of trees per year
p.counts <- il.tree.df.final %>% group_by(tree.invyr) %>% summarise(t_count = n())
head(p.counts)

# convert all numbers into common name
master.tree.species.list <- read.table("v9-5_2024-10_Natl_MasterTreeSpeciesList.csv", sep=",", header=TRUE)
head(master.tree.species.list)
tree.species <- master.tree.species.list$FIA.Code
common.name <- master.tree.species.list$Common.Name
scientific.name <- master.tree.species.list$SCIENTIFIC_NAME
tree.species.list <- data.frame(tree.species, common.name, scientific.name)
sp.count <- merge(s.counts, tree.species.list)
spp.counts <- merge(sp.count, p.counts, by="tree.invyr")

head(spp.counts)
tree.pct.calc <- function(x) {
  return((as.numeric(x[3])/as.numeric(x[6]))*100)
}
tree.pct <- apply(spp.counts, 1, tree.pct.calc)
head(tree.pct)

tree.pct.final <- cbind(spp.counts, tree.pct)
tree.pct.final.sorted <- arrange(tree.pct.final, desc(tree.pct)) %>% group_by(tree.invyr)
write.csv(tree.pct.final, "tree_species.csv", row.names=FALSE)




# forest loss: rate of change of canopy cover

# remove all plots that don't have a prev plot and are not a prev plot themselves
# or inversely, keep all plots that have a prev plot or are in the prev plot list
prev.plots <- unique(il.plot.df.final$plot.prev)
length(prev.plots)
il.plot.df.remer <- filter(il.plot.df.final, !is.na(plot.prev) | plot.cn %in% prev.plots)


plot.county.idx <- which(colnames(il.plot.df.remer[1,]) == "plot.county")
plot.plot.idx <- which(colnames(il.plot.df.remer[1,]) == "plot.plot")

# function that pastes county, unit, and plot together to uniquely identify a plot
create.plot.id <- function(x) {
  return(paste(x[plot.county.idx:plot.plot.idx], collapse=""))
}

# apply above function and add ids to il.plot.df.remer
plot.id <- apply(il.plot.df.remer, 1, create.plot.id)
il.plot.df.remer.id <- cbind(il.plot.df.remer, plot.id)

# group by newly created id, only keep records that have exactly 5 resamples over the years
# https://dplyr.tidyverse.org/articles/grouping.html
plot.group <- il.plot.df.remer.id %>% group_by(plot.id)
plot.group.freq <- filter(plot.group %>% tally(sort=TRUE), n==5)
# https://stackoverflow.com/questions/21618423/extract-a-dplyr-tbl-column-as-a-vector
plot.five <- pull(plot.group.freq, plot.id)

# a data frame of all plots that were measured sequentially in 2001, 2006, 2011, 2016, 2023
il.plot.df.remer.final <- filter(il.plot.df.remer.id, plot.id %in% plot.five)

# calculate change in canopy cover
# merge cond with this
# rename tree.plot.cn to plot.cn so that we can merge it correctly
# https://stackoverflow.com/questions/1299871/how-to-join-merge-data-frames-inner-outer-left-right
colnames(il.cond.df.final)[2] <- "plot.cn"
il.plotcond.df <- filter(merge(il.plot.df.remer.final, il.cond.df.final, by="plot.cn"), cond.canopy.cover > 0)
# il.plotcond.df.more <- merge(il.plot.df.remer.final, il.cond.df.final, by="plot.cn", all.x=TRUE)

plot.prev.idx <- which(colnames(il.plotcond.df) == "plot.prev")
plot.canopy.idx <- which(colnames(il.plotcond.df) == "cond.canopy.cover")
# for each unique plot id, we can now calculate percent change from prev 
canopy.cover.pct.change <- function(x) {
  prev <- x[plot.prev.idx]
  prev.record <- filter(il.plotcond.df, plot.cn == prev)
  # don't divide by 0
  if (nrow(prev.record) > 0 & as.numeric(prev.record[1, plot.canopy.idx]) > 0) {
    # percent change is (new - old)/old * 100
    pct.change <- (as.numeric(x[plot.canopy.idx]) - as.numeric(prev.record[1, plot.canopy.idx]))/as.numeric(prev.record[1, plot.canopy.idx])
    if (is.finite(pct.change)) {
      return(pct.change*100)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
  
}
# cast sequence numbers as strings, this solves automatic casting from big int to e string
# https://stackoverflow.com/questions/32339636/long-numbers-as-a-character-string
il.plotcond.df$plot.cn <- as.character(il.plotcond.df$plot.cn)
il.plotcond.df$plot.prev <- as.character(il.plotcond.df$plot.prev)

# apply above function to add pct.change column to df
pct.change <- apply(il.plotcond.df, 1, canopy.cover.pct.change)
il.plotcond.df.final <- cbind(il.plotcond.df, pct.change)

write.csv(il.plotcond.df.final, "IL_CANOPY_COVER_CHANGE.csv", row.names=FALSE)

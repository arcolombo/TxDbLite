source("repPrep.R")
source("repBaseToRmsk.R")

humanRepeats <- c("humrep.fa", "humsub.fa")
human <- suppressWarnings(repPrep(humanRepeats))

repClass <- readRDS("repClass.rds")
humanRepClass <- repBaseToRmsk(human, repClass, method='lv')

## same for mouse
#
# mouseRepeats <- c("mousub.fa", "rodrep.fa")
# mouse <- suppressWarnings(repPrep(mouseRepeats))


# export mice data for stata
library(mice)
# don't really need data.table but makes adding the indices easier
library(data.table)

Data <-
  "B:/Kim/02-Inflammation-mediation/02-analysis/data/"

setwd(Data)

## Load data and convert to long
load("HPC/06b-imp-all-sub.RData")

# Function to export mice imputed datasets to Stata ice datasets
# https://stackoverflow.com/questions/49965155/importing-mice-object-to-stata-for-analysis

mice2stata <- function(imp, path = "stata", type = "dta") {
  completed <-
    lapply(seq_len(imp_all_sub$m), function(i)
      complete(imp_all_sub, i))
  data_out <- rbindlist(completed, idcol = "_mj")
  data_out <-
    bind_rows(imp_all_sub$data, data_out, .id = imp_all_sub$aln)
  setDT(data_out)
  data_out[, `_mj` := replace(`_mj`, is.na(`_mj`), 0L)]
  data_out[, `_mi` := rowid(`_mj`)]
  if (type == "dta") {
    foreign::write.dta(data_out, file = paste(path, type, sep = "."))

  } else {
    write.csv(
      data_out,
      file = paste(path, type, sep = "."),
      na = "",
      row.names = FALSE
    )
  }
}

mice2stata(imp_all_sub, path = "imp_for_stata", type = "dta")

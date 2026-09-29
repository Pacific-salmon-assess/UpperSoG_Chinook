source("01-Salmon-camp-CM.R")
source("01-Woss-camp-CM.R")
source("01-Adam-camp-CM.R")
source("01-QuinsamCampbell-camp-CM.R")


run_QC_CM(prod = "base")
run_QC_CM(prod = "low")
run_QC_CM(prod = "high")
run_Woss_CM(prod = "base")
run_Woss_CM(prod = "low")
run_Woss_CM(prod = "high")
run_Salmon_CM(prod = "base")
run_Salmon_CM(prod = "low")
run_Salmon_CM(prod = "high")
run_Adam_CM(prod = "base")
run_Adam_CM(prod = "low")
run_Adam_CM(prod = "high")



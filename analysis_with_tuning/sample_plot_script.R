
if (.Platform$OS.type == "windows"){
  directory_where_code_is = "C:\\Users\\kapelner\\workspace\\CGMBART_GPL"
}
setwd(directory_where_code_is)
source("r_scripts/bart_package_inits.R")
source("r_scripts/bart_package_data_preprocessing.R")
source("r_scripts/bart_package_builders.R")
source("r_scripts/bart_package_plots.R") ##altered to trash label at the top
source("r_scripts/bart_package_variable_selection.R")
source("r_scripts/bart_package_f_tests.R")

setwd("C:/Users/kapelner/workspace/bart_gene/analysis_with_tuning/")
      

source("simulation_params.R")

setwd(directory_where_code_is)

y = gene_train[,1]
x = data.frame(tf_train)


bart_machine = build_bart_machine(X=x, y=y)

investigate_var_importance(bart_machine, num_replicates_for_avg = 1, num_var_plot = 15)
counts = get_var_counts_over_chain(bart_machine)

sel = var_selection_by_permute_response_three_methods(bart_machine, num_reps_for_avg = 5, num_var_plot = 20)
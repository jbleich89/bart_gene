
get_sparsity_by_method_and_plot_num_tfs = function(method){
	gene_by_tf = matrix(0, nrow = MAX_GENE_NUM, ncol = ncol(tf_train))
	colnames(gene_by_tf) = colnames(tf_train)
	for (g in 1 : MAX_GENE_NUM){
		for (t in 1 : ncol(tf_train)){
			names_of_all_tfs = colnames(tf_train)
			tfs = names(all_results[[g]][['0']][['20']][["0.05"]][[method]])
			cols = which(names_of_all_tfs %in% tfs)
			gene_by_tf[g, cols] = 1
		}
	}	
	
	sort_gene_by_tf = t(t(sort(colSums(gene_by_tf), decr = TRUE)))
	xtable(sort_gene_by_tf)
	sum(colSums(gene_by_tf))
	sum(rowSums(gene_by_tf) == 0) / nrow(gene_by_tf)
	
	par(mar = c(5,5,5,5))
	poissonity_pval = poisson.mtest(rowSums(gene_by_tf), R = 999)$p.value
	barplot(table(rowSums(gene_by_tf)), main = paste("# of TF's Affecting", MAX_GENE_NUM, "Genes for method", method), xlab = "# TF's", ylab = "# Genes")
	
	sum(gene_by_tf)/(nrow(gene_by_tf) * ncol(gene_by_tf))
}


get_sparsity_by_method_and_plot_num_tfs("important_tfs_at_alpha_simul_max")
get_sparsity_by_method_and_plot_num_tfs("important_tfs_at_alpha_simul_se")
get_sparsity_by_method_and_plot_num_tfs("important_tfs_at_alpha_pointwise")


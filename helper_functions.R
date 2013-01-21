##Setup function
##This function creates:
#PriorWeights: rows are genes, cols are TFs
##gene.train and gene.test rows are obs and cols are genes
##TF train and TF test: rows are obs and cols are TFs
setup=function(){
	##remove first two columns
	prior_nums=priors[,3:ncol(priors)]
	geneNames=as.character(gene.exp[,2])
	
#  ##get weights- rownames are genes. colnames are TFs
#  priorWeights=apply(prior_nums,2,ChangePriorVec)
#  rownames(priorWeights)=geneNames
#  
#  ##prior col 1 goes with gene col 2 and vice versa
#  priorWeights[,1]==gene.exp[,2]
	
	
	##need transpose of tf mat
	tf.names=tf.exp[,1]
	temp.tf=tf.exp[,5:ncol(tf.exp)]
	tf.full=t(temp.tf) ##This is training matrix
	colnames(tf.full)=tf.names
	rownames(tf.full)=colnames(tf.exp)[5:ncol(tf.exp)]
	dim(tf.full) #check 314 conditions and 39 TFs
	
	##randomize later
	cv_hold=ceiling(.1*nrow(tf.full))
	test_hold=ceiling(.1*nrow(tf.full))
	
	#get training and test
	train.end=nrow(tf.full)-(cv_hold + test_hold)
	tf.train=tf.full[1:train.end,]
	tf.cv=tf.full[(train.end+1):(train.end+cv_hold),]
	tf.test=tf.full[(train.end+cv_hold+1):nrow(tf.full),]
	
	##Clean gene matrix
	geneMatClean=t(gene.exp[,5:ncol(gene.exp)])
	rownames(geneMatClean)=rownames(tf.full)
	colnames(geneMatClean)=geneNames
	
	##Get gene test and train
	gene.train=geneMatClean[1:train.end,]
	gene.cv=geneMatClean[(train.end+1):(train.end+cv_hold),]
	gene.test=geneMatClean[(train.end+cv_hold+1):nrow(geneMatClean),]
	list(
			tf.train=tf.train,
			tf.cv=tf.cv, 
			tf.test=tf.test,
			gene.train=gene.train, 
			gene.cv=gene.cv, 
			gene.test=gene.test
	)
}

prior_adj = function(prior_mat, lower_bound = .05, upper_bound = .95){
	adj_mat = prior_mat
	for (i in 1 : nrow(adj_mat)){
		row = prior_mat[i, ]
		row = ifelse(row < lower_bound, lower_bound, row)
		row = ifelse(row > upper_bound, upper_bound, row)
		adj_mat[i,] = row
	}
	adj_mat
}

get_gene_response_as_dataframe = function(gene, data){
	as.numeric(data[,gene])
}
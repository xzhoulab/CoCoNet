#' @title The main function that leverages gene co-expression pattern to infer trait-relevant tissues in genome-wide association studies
#' @description CoCoNet applies a composite likelihood-based covariance regression network model to identify trait-relevant tissues or cell types 
#' @param gene_effect_size: A length m vector of gene level effect sizes for m genes
#' @param max_path: maximum number of paths linking between two genes considered in the model. The default value is 1. 
#' @param adjacency_mat: A m by m adjacency matrix
#' @return A list of estimated parameters
#' \item{loglikelihood}{A value of log likelihood of selected disease-tissue pair}
#' \item{sigma0}{The estimated sigma0 square in our model, where sigma0 square multiply adjacency_mat measures gene-level effect size correlation due to direct connections among genes}
#' \item{sigma1}{The estimated sigma1 square in our model, where sigma1 square multiply adjacency_mat measures gene-level effect size correlation due to indirect connections among genes through two-paths}
#' \item{sigma2}{If A2 model is selected, this is the estimated sigma2 square in our model, where sigma2 square multiply adjacency_mat measures gene-level effect size correlation due to indirect connections among genes through two-paths}
#' \item{sigma3}{If A3 model is selected, this is the estimated sigma3 square in our model, where sigma3 square multiply adjacency_mat measures gene-level effect size correlation due to indirect connections among genes through three-paths}
#' \item{sigma4}{If A4 model is selected, this is the estimated sigma4 square in our model, where sigma4 square multiply adjacency_mat measures gene-level effect size correlation due to indirect connections among genes through four-paths}
#' \item{rho}{A value measuring relative signal strength of one-path gene co-expression pattern on gene-level effect sizes}
CoCoNet <- function(gene_effect_size, max_path = 1, adjacency_mat){ 
	
if(max_path == 1){
	y = scale(gene_effect_size) 
	A = adjacency_mat
	res = optim(c(1, 0, 0), coconetA1, A=A, y=y, method = "Nelder-Mead",hessian = FALSE)
	result = list()
	result$loglikelihood = -res$value
	result$sigma2_0 = res$par[1]
	result$sigma2_1 = res$par[2]
	result$rho = res$par[2]/(res$par[1]+res$par[2])
	return(result)
}else if(max_path == 2){
	y = scale(gene_effect_size) 
	A = adjacency_mat
	A2 = t(A)%*%A 
	res = optim(c(1,0,0,0), coconetA2, A=A,A2=A2, y=y, method = "Nelder-Mead",hessian = FALSE)
	result = list()
	result$loglikelihood = -res$value
	result$sigma2_0 = res$par[1]
	result$sigma2_1 = res$par[2]
	result$sigma2_2 = res$par[3]
	result$rho = res$par[2]/(res$par[1]+res$par[2])
	return(result)
}else if(max_path == 3){
	y = scale(gene_effect_size) 
	A = adjacency_mat
	A2 = t(A)%*%A 
	diag(A2) = 0
	A3 = t(A2)%*%A 
	res = optim(c(1,0,0,0,0), coconetA3, A=A,A2=A2,A3=A3, y=y, method = "Nelder-Mead",hessian = FALSE)
	result = list()
	result$loglikelihood = -res$value
	result$sigma2_0 = res$par[1]
	result$sigma2_1 = res$par[2]
	result$sigma2_2 = res$par[3]
	result$sigma2_3 = res$par[4]
	result$rho = res$par[2]/(res$par[1]+res$par[2])
	return(result)
}else if(max_path == 4){
	y = scale(gene_effect_size) 
	A = adjacency_mat
	A2 = t(A)%*%A 
	diag(A2) = 0
	A3 = t(A2)%*%A 
	diag(A3) = 0
	A4 = t(A4)%*%A 
	res = optim(c(1,0,0,0,0,0), coconetA4, A=A,A2=A2,A3=A3,A4=A4, y=y, method = "Nelder-Mead",hessian = FALSE)
	result = list()
	result$loglikelihood = -res$value
	result$sigma2_0 = res$par[1]
	result$sigma2_1 = res$par[2]
	result$sigma2_2 = res$par[3]
	result$sigma2_3 = res$par[4]
	result$sigma2_4 = res$par[5]
	result$rho = res$par[2]/(res$par[1]+res$par[2])
	return(result)
}
}

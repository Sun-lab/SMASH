# Minor Functions
poss_mult = function(x,alloc){
	x2 = unique(c(1,x))
	output = x2[x2 > 0]
	if(alloc != 1){
		output = 1
	}
	
	return(output)
	
}
calc_maf = function(purity,vec_q,mult,tCN,vec_alloc){
	purity * mult * sum(vec_q * vec_alloc) /
		(tCN * purity + 2*(1-purity))
}

#' A collection of pre-defined subclone configurations.
#' 
#' A R list containing subclone configurations in matrix 
#'	form for 1 to 5 subclones. For each matrix, each column 
#'	corresponds to a subclone and each row corresponds to a 
#'	variant's allocation across all subclones. For example, 
#'	the first row of each matrix is a vector of 1's to 
#'	represent clonal variants, variants present in all subclones.
#'
"eS"

#' @title gen_subj_truth
#' @description Simulates copy number states, multiplicities, allocations
#' @param mat_eS A subclone configuration matrix pre-defined in R list \code{eS}
#' @param maxLOCI A positive integer number of simulated somatic variant calls
#' @param nCN A positive integer for the number of allelic copy number pairings 
#'	to sample from. If \code{NULL}, it will be randomly sampled between 1 and 5.
#' @return A list containing the following components:
#' \describe{
#' \item{\code{subj_truth}}{dataframe of each variant's simulated minor 
#'	(\code{CN_1}) and major (\code{CN_2}) copy number states, total copy 
#'	number (\code{tCN}), subclone allocation (\code{true_A}), multiplicity 
#'	(\code{true_M}), mutant allele frequency (\code{true_MAF}), and cellular 
#'	prevalence (\code{true_CP})}
#' \item{\code{purity}}{tumor purity}
#' \item{\code{eta}}{the product of tumor purity and subclone proportions}
#' \item{\code{q}}{vector of subclone proportions}
#' }
#' 
#' @export
gen_subj_truth = function(mat_eS,maxLOCI,nCN = NULL){

	# Enumerate possible CN states
	all_CN_states = c()
	for(tCN in seq(5)){
		CN_1 = seq(0,floor(tCN/2))
		CN_2 = tCN - CN_1
		all_CN_states = rbind(all_CN_states,
			smart_df(tCN = tCN,CN_1 = CN_1,CN_2 = CN_2))
		rm(CN_1,CN_2)
	}
	all_CN_states$prob = dgeom(abs(all_CN_states[,2] - 1) 
		+ abs(all_CN_states[,3] - 1),prob = 0.45)
	all_CN_states$prob = all_CN_states$prob / sum(all_CN_states$prob)
	tmp_all_CN_states = all_CN_states[which(all_CN_states$CN_1 == all_CN_states$CN_2),]
	tmp_all_CN_states$prob = tmp_all_CN_states$prob / sum(tmp_all_CN_states$prob)
  
	true_purity = runif(1,0.2,1)
	while( TRUE ){
		# Calculate eta and purity and q
		true_q = stats::runif(ncol(mat_eS),0,2)
		true_q = true_q / sum(true_q)
		true_eta = true_q * true_purity
		if( ncol(mat_eS) == 1 ) break
		min_CP_diff = min(diff(sort(mat_eS %*% true_q)))
		if( min(true_q) >= 0.025 && min_CP_diff >= 0.025 ) break
	}
  
	# Generate CN states
	if( !is.null(nCN) ){
		num_CN_states = nCN
	} else {
		num_CN_states = sample(seq(1,5),1)
	}
	if( num_CN_states == 1 ){
		# sample either tCN = 2 or 4 with allelic balance
		poss_CN = tmp_all_CN_states[sample(x = seq(nrow(tmp_all_CN_states)),
			size = 1,prob = tmp_all_CN_states$prob),c("CN_1","CN_2","prob")]
	} else if( num_CN_states > 1 ){
		poss_CN = all_CN_states[sample(x = seq(nrow(all_CN_states)),
			size = num_CN_states,replace = FALSE,
			prob = all_CN_states$prob),c("CN_1","CN_2","prob")]
	}
  poss_CN$prob = poss_CN$prob / sum(poss_CN$prob)
  
  # Generate multiplicity, allocation
  subj_truth = poss_CN[sample(x = seq(nrow(poss_CN)),
		size = maxLOCI,replace = TRUE,prob = poss_CN$prob),
		c("CN_1","CN_2")]
	rownames(subj_truth) = NULL
	subj_truth$tCN = rowSums(subj_truth)
  
	while(TRUE){
		subj_truth$true_A = sample(seq(nrow(mat_eS)),maxLOCI,replace = TRUE)
		tmp_tab = table(subj_truth$true_A)
		min_alloc = min(tmp_tab)
		if( min_alloc >= 2 && length(tmp_tab) == ncol(mat_eS) ) break
	}
  
	subj_truth$true_M = apply(subj_truth[,c("CN_1","CN_2","true_A")],1,
		function(x) sample(poss_mult(x = x[1:2],alloc = x[3]),1))
	subj_truth$true_MAF = apply(subj_truth[,c("true_M","tCN","true_A")],1,
		function(x) calc_maf(purity = true_purity,vec_q = true_q,
			mult = x[1],tCN = x[2],vec_alloc = mat_eS[x[3],]))
	subj_truth$true_CP = sapply(subj_truth$true_A,function(x) sum(mat_eS[x,] * true_q))
	uniq_data2 = unique(subj_truth[,c("CN_1","CN_2")])
  
  # Output
  return(list(subj_truth = subj_truth,
		purity = true_purity,eta = true_eta,
		q = true_q))
	
}

#' @title gen_ITH_RD
#' @description Simulates observed alternate and reference read counts
#' @param DATA The output data.frame from \code{gen_subj_truth} 
#' @param RD A positive integer for the mean read depth generated from 
#'	the negative binomial distribution
#' @return A matrix of simulated alternate and reference read counts.
#' 
#' @export
gen_ITH_RD = function(DATA,RD){
	
	BB = nrow(DATA)
	while(TRUE){
		DATA$tDP = rnbinom(n = BB,mu = RD,size = 2) + 30
		DATA$tAD = rbinom(n = BB,size = DATA$tDP,prob = DATA$true_MAF)
		if( mean(DATA$tAD == 0) < 0.1 ) break
	}
	DATA$tRD = DATA$tDP - DATA$tAD
	return(as.matrix(DATA[,c("tAD","tRD")]))
}

#' @title ITH_optim
#' @description Performs EM algorithm for a given configuration matrix
#' @param my_data A R dataframe containing the following columns:
#' \describe{
#' \item{\code{tAD}}{tumor alternate read counts}
#' \item{\code{tRD}}{tumor reference read counts}
#' \item{\code{CN_1}}{minor allele count}
#' \item{\code{CN_2}}{major allele count, where \code{CN_1 <= CN_2}}
#' \item{\code{tCN}}{\code{CN_1 + CN_2}}
#' }
#' @param my_purity A single numeric value of known/estimated purity
#' @param init_eS A subclone configuration matrix pre-defined in R 
#'	list \code{eS}
#' @param pi_eps0 A user-specified parameter denoting the proportion 
#'	of loci not explained by the combinations of purity, copy number, 
#'	multiplicity, and allocation. If \code{NULL}, it is initialized at 
#'	1e-3. If set to 0.0, the parameter is not estimated.
#' @param my_unc_q An optimal initial vector for the unconstrained 
#'	\code{q} vector, useful after running \code{grid_ITH_optim}. If 
#'	this variable is \code{NULL}, then the subclone proportions, 
#'	\code{q}, are randomly initialized. For instance, if 
#'	\code{my_unc_q = ( x1 , x2 )}, then \code{q = ( exp(x1) / (1 + exp(x1) + exp(x2)) , exp(x2) / (1 + exp(x1) + exp(x2)) , 1 / (1 + exp(x1) + exp(x2))}.
#' @param max_iter Positive integer, preferably 1000 or more, setting 
#'	the maximum number of iterations
#' @param my_epsilon Convergence criterion threshold for changes in 
#'	the log likelihood, preferably 1e-6 or smaller
#' @return If the EM algorithm converges, the output will be a list containing
#' \describe{
#' \item{\code{iter}}{number of iterations}
#' \item{\code{converge}}{convergence status}
#' \item{\code{unc_q0}}{initial unconstrained subclone proportions parameter}
#' \item{\code{unc_q}}{unconstrained estimate of \code{q}}
#' \item{\code{q}}{estimated subclone proportions among cancer cells}
#' \item{\code{CN_MA_pi}}{estimated mixture probabilities of multiplicities 
#'	and allocations given copy number states}
#' \item{\code{eta}}{estimated subclone proportion among tumor cells}
#' \item{\code{purity}}{user-inputted tumor purity}
#' \item{\code{entropy}}{estimated entropy}
#' \item{\code{infer}}{A R dataframe containing inferred variant allocations 
#'	(\code{infer_A}), multiplicities (\code{infer_M}), cellular prevalences 
#'	(\code{infer_CP}).}
#' \item{\code{ms}}{model size, number of parameters within parameter space}
#' \item{\code{LL}}{The observed log likelihood evaluated at maximum likelihood 
#'	estimates.}
#' \item{\code{AIC = 2 * LL - 2 * ms}}{Negative AIC, used for model selection}
#' \item{\code{BIC = 2 * LL - ms * log(LOCI)}}{Negative BIC, used for model selection}
#' \item{\code{LOCI}}{The number of inputted somatic variants.}
#' }
#'
#' @export
ITH_optim = function(my_data,my_purity,init_eS,
	pi_eps0 = NULL,my_unc_q = NULL,max_iter = 4e3,
	my_epsilon = 1e-6){

  if( !is.data.frame(my_data) ){
    stop("my_data is not a data.frame")
  }
  req_vars = c("tCN","CN_1","CN_2","tAD","tRD")
  miss_vars = req_vars[!req_vars %in% names(my_data)]
  if( length(miss_vars) > 0 ){
		stop(sprintf("The following columns are missing from my_data: %s",
			paste(miss_vars,collapse=", ")))
  }
  if( min(my_data$tAD) <= 0 ){
		stop("Input variants with positive integer alternate read counts")
	}
	if( !all(my_data$tAD == round(my_data$tAD)) || !all(my_data$tRD == round(my_data$tRD)) 
			|| !all(my_data$tCN == round(my_data$tCN)) || !all(my_data$CN_1 == round(my_data$CN_1))
			|| !all(my_data$CN_2 == round(my_data$CN_2)) ){
		stop("Alternate and reference read counts, total, minor, and major copy numbers should all be integers.")
	}
	if( my_purity <= 0 || my_purity > 1 ){
		stop("Tumor purity should be strictly > 0 and <= 1.")
	}
	
  # Calculate some fixed values
  my_data$tDP = rowSums(my_data[,c("tAD","tRD")])
  my_data$log_DP = log(my_data$tDP)
  my_data$LBC = apply(my_data[,c("tDP","tAD")],1,
		function(xx) lchoose(xx[1],xx[2]))
  
  # Initialize Params
  rownames(my_data) = NULL
  uniq_CN = unique(my_data[,c("tCN","CN_2","CN_1")])
  uniq_CN = uniq_CN[order(uniq_CN$tCN,uniq_CN$CN_2),]
  uniq_CN_MA = c()
  for(ii in seq(nrow(uniq_CN))){
    one_row = uniq_CN[ii,c("tCN","CN_2","CN_1")]
    for(kk in seq(nrow(init_eS))){
      tmp_mult = poss_mult(as.numeric(one_row[,c("CN_1","CN_2")]),
				alloc = kk)
      for(jj in seq(length(tmp_mult))){
        uniq_CN_MA = rbind(uniq_CN_MA,smart_df(one_row,
					M = tmp_mult[jj],A = kk))
      }
    }
  }
  rownames(uniq_CN_MA) = NULL
  init_pi_CN_MA = rep(NA,nrow(uniq_CN_MA))
  num_pi = 0
  for(cn in unique(uniq_CN$tCN)){
    for(CN_2 in uniq_CN$CN_2[which(uniq_CN$tCN == cn)]){
      tmp_index = which(uniq_CN_MA$tCN == cn & uniq_CN_MA$CN_2 == CN_2)
      tmp_len = length(tmp_index)
      init_pi_CN_MA[tmp_index] = rep(1/tmp_len,tmp_len)
      num_pi = num_pi + tmp_len - 1
    }
  }
  
  # Initialize vector q, cancer subclone proportions
  if( is.null(my_unc_q) ){
    rr = sample(c(0,1),1)
    if(rr == 0){
      init_unc_q = runif(ncol(init_eS) - 1,-3,1)
    } else if(rr == 1){
      init_unc_q = runif(ncol(init_eS) - 1,1,2)
    }
  } else {
    init_unc_q = my_unc_q
  }
  unc_q0 = init_unc_q
  
	if( is.null(pi_eps0) ){
		pi_eps0 = 1e-3
	}
  
	# Create B matrix: Stores 0s and 1s for what's feasible given copy number state
	num_LOCI = nrow(my_data)
	my_mat_B = matrix(0,num_LOCI,nrow(uniq_CN_MA))
	for(bb in seq(num_LOCI)){
		one_row = my_data[bb,]
		one_index = which(uniq_CN_MA$tCN == one_row$tCN 
			& uniq_CN_MA$CN_2 == one_row$CN_2)
		my_mat_B[bb,one_index] = 1
	}
  
  # Optimize
  if( ncol(init_eS) == 1 ){
		ith_out = Rcpp_ITH_opt(RD = as.matrix(my_data[,c("tAD","tRD")]),
			log_DP = my_data$log_DP,LBC = my_data$LBC,BB = my_mat_B,
			uniq_BB = unique(my_mat_B),uniq_CN_MA = as.matrix(uniq_CN_MA),
			eS = init_eS,purity = my_purity,tCN = my_data$tCN,
			pi_eps0 = pi_eps0,pi0 = init_pi_CN_MA,unc_qq0 = unc_q0,
			mstep = FALSE,max_iter = max_iter,eps = my_epsilon,show = FALSE)
  } else {
		ith_out = Rcpp_ITH_opt(RD = as.matrix(my_data[,c("tAD","tRD")]),
			log_DP = my_data$log_DP,LBC = my_data$LBC,BB = my_mat_B,
			uniq_BB = unique(my_mat_B),uniq_CN_MA = as.matrix(uniq_CN_MA),
			eS = init_eS,purity = my_purity,tCN = my_data$tCN,
			pi_eps0 = pi_eps0,pi0 = init_pi_CN_MA,unc_qq0 = unc_q0,
			mstep = TRUE,max_iter = max_iter,eps = my_epsilon,show = FALSE)
  }
  # ith_out[c("converge","iter","LL","ms","BIC","pi_eps","CN_MA_pi","qq")]
  
  # Output
  if( ith_out$converge == 1 ){
    uniq_CN_MA$infer_MAF = apply(uniq_CN_MA,1,function(x) 
      calc_maf(my_purity,ith_out$qq,x[4],x[1],init_eS[x[5],]))
    uniq_CN_MA$infer_CP = sapply(uniq_CN_MA$A,function(x) 
      sum(init_eS[x,] * ith_out$qq))
    noise_df = smart_df(t(rep(NA,ncol(uniq_CN_MA))))
    names(noise_df) = names(uniq_CN_MA)
    uniq_CN_MA_noise = rbind(noise_df,uniq_CN_MA)
    infer_out = uniq_CN_MA_noise[apply(ith_out$Z,1,which.max),
			c("M","A","infer_MAF","infer_CP")]
		rownames(infer_out) = NULL
		names(infer_out)[1:2] = paste0("infer_",names(infer_out)[1:2])
		list(converge = TRUE,unc_q0 = unc_q0,iter = ith_out$iter,
			q = ith_out$qq,CN_MA_pi = smart_df(uniq_CN_MA,pi = ith_out$CN_MA_pi),
			pi_eps = ith_out$pi_eps,eta = ith_out$qq * my_purity,
			purity = my_purity,entropy = round(ith_out$entropy,2),
			infer = infer_out,ms = ith_out$ms,LL = ith_out$LL,
			AIC = ith_out$AIC,BIC = ith_out$BIC,LOCI = num_LOCI)
  } else {
    list(converge = FALSE,BIC = -Inf)
  }
  
}

#' @title grid_ITH_optim
#' @description This function performs a grid search over enumerated 
#' 	configurations within the pre-defined list \code{eS}
#' @inheritParams ITH_optim
#' @param list_eS A nested list of subclone configuration matrices
#' @param trials Positive integer, number of random initializations 
#'	of subclone proportions
#' @return A R list containing two objects. \code{GRID} is a 
#'	dataframe where each row denotes a feasible subclone configuration 
#'	with corresponding subclone proportion estimates \code{q} and 
#'	somatic variant allocations \code{alloc}. \code{INFER} is a list 
#'	where \code{INFER[[i]]} corresponds to the \code{i}-th row or 
#'	model of \code{GRID}.
#' 
#' @export
grid_ITH_optim = function(my_data,my_purity,list_eS,pi_eps0 = NULL,
	trials = 20,max_iter = 4e3,my_epsilon = 1e-6){
	
	# Run EM on all provided subclone configurations
	all_infer = list() # new code
	all_ck = c()
	
	for(cc in seq(length(list_eS))){
		cat(sprintf("cc = %s; ",cc))
		tmp_ck0 = c()
		tmp_infer0 = list() # new code
	
	for(kk in seq(length(list_eS[[cc]]))){
		cat(sprintf("kk = %s ",kk))
		tmp_ck = c()
		tmp_infer = list() # new code
		
		for(count in seq(trials)){
			if(count %% 5 == 0) cat("*")
			ith_out = list(converge = FALSE,BIC = -Inf)
			
			ith_out1 = ITH_optim(my_data = my_data,my_purity = my_purity,
				init_eS = list_eS[[cc]][[kk]],pi_eps0 = 0,my_unc_q = NULL, 
				max_iter = max_iter,my_epsilon = my_epsilon)
			if( ith_out1$converge ){
				ith_out = ith_out1
				rm(ith_out1)
			}
			
			ith_out2 = ITH_optim(my_data = my_data,my_purity = my_purity,
				init_eS = list_eS[[cc]][[kk]],pi_eps0 = pi_eps0,my_unc_q = NULL, 
				max_iter = max_iter,my_epsilon = my_epsilon)
			if( ith_out2$converge && ith_out2$BIC >= ith_out$BIC ){
				ith_out = ith_out2
			}
			
			if( ith_out$converge ){
				# Even though we have convergence, check all subclones have an allocation
				tab_A = table(ith_out$infer$infer_A)
				len_A = length(tab_A)
				if( len_A == cc && min(round(ith_out$q,2)) > 0 ){
					tmp_ck = rbind(tmp_ck,smart_df(cc = cc,kk = kk,
						ms = ith_out$ms,entropy = ith_out$entropy,
						LL = ith_out$LL,AIC = ith_out$AIC,BIC = ith_out$BIC,
						q = paste(round(ith_out$q,2),collapse = ","),
						unc_q0 = paste(round(ith_out$unc_q0,6),collapse = ","),
						alloc = paste(paste0(names(tab_A),"|",tab_A),collapse = ";"),
						pi_eps = ith_out$pi_eps))
					tmp_infer[[length(tmp_infer) + 1]] = ith_out$infer # new code
				}
			}
			
		}
		cat(" ")
		
		if( !is.null(tmp_ck) && nrow(tmp_ck) > 0 ){
			keep_idx = !duplicated(tmp_ck[,c("cc","kk","ms","entropy","LL")]) # new code
			tmp_ck = tmp_ck[keep_idx,]
			tmp_ck0 = rbind(tmp_ck0,tmp_ck)
			
			# new code
			tmp_infer = tmp_infer[keep_idx]
			for(jk in seq(length(tmp_infer))){
				tmp_infer0[[length(tmp_infer0) + 1]] = tmp_infer[[jk]]
			}
			
		}
		
	}
    
    if( !is.null(tmp_ck0) && nrow(tmp_ck0) > 0 ){
      all_ck = rbind(all_ck,tmp_ck0)
			
			# new code
			for(jk in seq(length(tmp_infer0))){
				all_infer[[length(all_infer) + 1]] = tmp_infer0[[jk]]
			}
			
    }
    cat("\n")
  }
	# all_ck = all_ck[order(all_ck$cc,all_ck$kk,all_ck$ms),]
	rownames(all_ck) = NULL
	# all_ck
  
	if(FALSE){ # Obtain details of inferred results
		cat("Obtain inferred results: ")
		thres = 0; num_res = nrow(all_ck)
		infer_list = list()
		for(ii in seq(num_res)){
			
			if( ii / num_res * 100 >= thres ){
				cat(paste0(thres,"% "))
				thres = thres + 5
			}
			
			if( all_ck$cc[ii] == 1 ){
				ith_out = ITH_optim(my_data = my_data,
					my_purity = my_purity,
					init_eS = list_eS[[all_ck$cc[ii]]][[all_ck$kk[ii]]],
					pi_eps0 = pi_eps0,my_unc_q = NULL,
					max_iter = max_iter,
					my_epsilon = my_epsilon)
			} else if( all_ck$cc[ii] > 1 ){
				ith_out = ITH_optim(my_data = my_data,
					my_purity = my_purity,
					init_eS = list_eS[[all_ck$cc[ii]]][[all_ck$kk[ii]]],
					pi_eps0 = pi_eps0,my_unc_q = as.numeric(strsplit(all_ck$unc_q0[ii],",")[[1]]),
					max_iter = max_iter,
					my_epsilon = my_epsilon)
				if( !ith_out$converge ){
					while(TRUE){
						ith_out = ITH_optim(my_data = my_data,
							my_purity = my_purity,
							init_eS = list_eS[[all_ck$cc[ii]]][[all_ck$kk[ii]]],
							pi_eps0 = pi_eps0,my_unc_q = NULL,
							max_iter = max_iter,
							my_epsilon = my_epsilon)
						if( ith_out$converge 
								&& round(ith_out$LL,2) == all_ck$LL[ii] 
								&& ith_out$ms == all_ck$ms[ii] ){
							break
						}
					}
				}
			}
			
			# After recovering inferred result, save to list
			infer_list[[ii]] = ith_out$infer
			rm(ith_out)
		}
		cat("\n")
  }
	
  return(list(GRID = all_ck,INFER = all_infer))
	
}

#' @title vis_GRID
#' @description A simple visualization of SMASH's grid of solutions
#' @param GRID The \code{GRID} object output from \code{grid_ITH_optim}.
#'
#' @export
vis_GRID = function(GRID){
	
	GRID$solu = factor(seq(nrow(GRID)))
	GRID$pBIC = 0.5 * GRID$BIC
	GRID$pBIC = exp(GRID$pBIC - logSumExp(GRID$pBIC))
	GRID$pAIC = 0.5 * GRID$AIC
	GRID$pAIC = exp(GRID$pAIC - logSumExp(GRID$pAIC))

	tmp_df = melt(GRID,id.vars = "solu",c("pBIC","pAIC"))
	# tmp_df
	
	solu = value = variable = NULL
	gg = ggplot(data = tmp_df,mapping = aes(x = solu,y = value,
		group = variable,color = variable)) +
		geom_line(size = 1.3) + geom_point(size = 3) +
		scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
		xlab("Solution") + ylab("Posterior Probability") +
		labs(color = "Information Criterion") +
		theme(legend.position = "bottom",
			text = element_text(size = 20))
	
	return(gg)
}

#' @importFrom stats dgeom rbinom rnbinom runif
#' @importFrom smarter smart_df logSumExp
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line geom_point
#'	scale_x_discrete guide_axis xlab ylab labs theme
#'	element_text
#' @importFrom Rcpp sourceCpp
#' @useDynLib SMASH
NULL

# Create package
# rm(list=ls()); library(smarter)
# smart_prepPack(pack_dir = "C:/Users/Admin/Desktop/github/SMASH",
#		pandoc = "C:/Program Files/RStudio/bin/pandoc",
#		make_vign = TRUE,cran = FALSE,build_dir = "C:/Users/Admin/Desktop")


###


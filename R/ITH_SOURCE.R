# Functions
smart_df = function(...){
	data.frame(...,stringsAsFactors = FALSE)
}
poss_mult = function(x,alloc){
	x2 = unique(c(1,x))
	output = x2[x2 > 0]
	if(alloc != 1){
		output = 1
	}
	output
}
calc_maf = function(purity,vec_q,mult,tCN,vec_alloc){
	purity * mult * sum(vec_q * vec_alloc) /
		(tCN * purity + 2*(1-purity))
}

#' @title gen_subj_truth
#' @description Simulates copy number states, multiplicities, allocations
#' @param mat_config A subclone configuration matrix pre-defined in R list \code{eS}
#' @param maxLOCI A positive integer number of simulated somatic variant calls
#' @export
#' @examples
#' data(eS)
#' gen_subj_truth(mat_config = eS[[2]][[1]],maxLOCI = 100)
gen_subj_truth = function(mat_config,maxLOCI,show_thres=FALSE){
	
	# Enumerate possible CN states
	all_CN_states = c()
	for(tCN in seq(5)){
		CN_1 = seq(0,floor(tCN/2))
		CN_2 = tCN - CN_1
		all_CN_states = rbind(all_CN_states,smart_df(tCN=tCN,CN_1=CN_1,CN_2=CN_2))
		rm(CN_1,CN_2)
	}
	all_CN_states$prob = dgeom(abs(all_CN_states[,2]-1)+abs(all_CN_states[,3]-1),prob=0.45)
	all_CN_states$prob = all_CN_states$prob / sum(all_CN_states$prob)
	tmp_all_CN_states = all_CN_states[which(all_CN_states$CN_1 == all_CN_states$CN_2),]
	tmp_all_CN_states$prob = tmp_all_CN_states$prob / sum(tmp_all_CN_states$prob)

	vartheta_to_eta = function(vv){
		xx = exp(vv)
		xx / (1 + sum(xx))
	}
	min_diff = function(xx){
		min(diff(sort(xx)))
	}
	
	while(TRUE){
		# Generate true_vartheta
		true_vartheta = runif(ncol(mat_config),-3,1)
		
		# Calculate eta and purity and q
		true_eta = vartheta_to_eta(true_vartheta)
		true_purity = sum(true_eta)
		true_q = true_eta / true_purity
		
		# Generate CN states
		num_CN_states = sample(seq(1,5),1)
		if( num_CN_states == 1 ){
			# sample either tCN = 2 or 4 with allelic balance
			poss_CN = tmp_all_CN_states[sample(seq(nrow(tmp_all_CN_states)),1,prob=tmp_all_CN_states$prob),
				c("CN_1","CN_2","prob")]
		} else if( num_CN_states > 1 ){
			poss_CN = all_CN_states[sample(seq(nrow(all_CN_states)),num_CN_states,
				replace=FALSE,prob=all_CN_states$prob),c("CN_1","CN_2","prob")]
		}
		poss_CN$prob = poss_CN$prob / sum(poss_CN$prob)

		# Generate multiplicity, allocation
		subj_truth = poss_CN[sample(seq(nrow(poss_CN)),maxLOCI,replace=TRUE,prob=poss_CN$prob),c("CN_1","CN_2")]
			rownames(subj_truth) = NULL
		subj_truth$tCN = rowSums(subj_truth)
		subj_truth$true_A = sample(seq(nrow(mat_config)),maxLOCI,replace=TRUE)
		subj_truth$true_M = apply(subj_truth[,c("CN_1","CN_2","true_A")],1,
			function(x) sample(poss_mult(x = x[1:2],alloc = x[3]),1))
		subj_truth$true_MAF = apply(subj_truth[,c("true_M","tCN","true_A")],1,
			function(x) calc_maf(purity = true_purity,vec_q = true_q,
				mult = x[1],tCN = x[2],vec_alloc = mat_config[x[3],]))
		subj_truth$true_CP = sapply(subj_truth$true_A,function(x) sum(mat_config[x,] * true_q))
		
		# Check at least 2 mutations per subclone
		tmp_tab = table(subj_truth$true_A)
		min_alloc = min(tmp_tab)
		
		# Check min_diff between CP per CN state
		min_CP_diff = 100
		uniq_data = unique(subj_truth[,c("CN_1","CN_2","true_CP")])
		uniq_data2 = unique(subj_truth[,c("CN_1","CN_2")])
		for(ii in seq(nrow(uniq_data2))){
			# ii = 2
			vec_CP = uniq_data$true_CP[which(uniq_data$CN_1 == uniq_data2$CN_1[ii] 
				& uniq_data$CN_2 == uniq_data2$CN_2[ii])]
			if( length(vec_CP) > 1 ){
				if( min_diff(vec_CP) < min_CP_diff ) min_CP_diff = min_diff(vec_CP)
			}
		}
		
		# Check min_diff(MAF) within CN state
		min_MAF_diff = 100
		# subj_truth = blah$subj_truth
		uniq_data = unique(subj_truth[,c("CN_1","CN_2","true_MAF")])
		uniq_data2 = unique(subj_truth[,c("CN_1","CN_2")])
		for(ii in seq(nrow(uniq_data2))){
			# ii = 1
			vec_MAF = uniq_data$true_MAF[which(uniq_data$CN_1 == uniq_data2$CN_1[ii] 
				& uniq_data$CN_2 == uniq_data2$CN_2[ii])]
			if( length(vec_MAF) > 1 ){
				min_MAF_diff = c(min_MAF_diff,min_diff(vec_MAF))
			}
		}
	
		if(show_thres){
			cat("min_alloc:\n\t"); cat(tmp_tab,"\n")
			cat("min_CP_diff:\n\t"); cat(min_CP_diff,"\n")
			cat("min_MAF_diff:\n\t"); cat(min_MAF_diff,"\n")
		}
		
		balance = length(which(uniq_data2$CN_1 == uniq_data2$CN_2 & uniq_data2$CN_2 %in% c(1,2))) > 0
		if( min(subj_truth$true_MAF) >= 0.025 # set lower bound on mean VAFs
			&& min(min_MAF_diff) > 0.01 # keep MAFs well separated
			&& balance # ensure some VCs exist in regions of allelic balance
			&& min(true_q) >= 0.05 # ensure the smallest subclone proportion among cancer cells to be >= 5%
			&& min_alloc >= 2 # ensure each subclone has at least 2 VCs classified to it
			&& length(tmp_tab) == ncol(mat_config) # ensure each subclone has VCs classified to them
			&& min_CP_diff > 0.05 # ensure mean CP of subclones well separated
			&& true_purity < 1 && true_purity > 0 # assuming bulk tissue, purity should be > 0 and < 1
			) break
	}
	
	# Output
	list(subj_truth = subj_truth,
		vartheta = true_vartheta,
		purity = true_purity,
		eta = true_eta,
		q = true_q)
}

#' @title gen_ITH_RD
#' @description Simulates observed alternate and reference read counts
#' @param DATA The output data.frame from \code{gen_subj_truth} 
#' @param RD A positive integer for the mean read depth generated from the negative binomial distribution
#' @export
#' @examples
#' data(eS)
#' truth = gen_subj_truth(mat_config = eS[[2]][[1]],maxLOCI = 100)
#' gen_ITH_RD(DATA = truth$subj_truth,RD = 100)
gen_ITH_RD = function(DATA,RD){
	
	B = nrow(DATA)
	while(TRUE){
		DATA$tDP = rnbinom(B,mu=RD,size=2) + 30
		DATA$tAD = rbinom(B,DATA$tDP,DATA$true_MAF)
		if(mean(DATA$tAD == 0) < 0.025) break
	}
	DATA$tRD = DATA$tDP - DATA$tAD
	as.matrix(DATA[,c("tAD","tRD")])
}


# Clustering read counts
#' @title ITH_optim
#' @description Performs EM algorithm for a given configuration matrix
#' @param my_data An R data frame containing columns
#' \enumerate{
#'		\item \code{tAD}, tumor alternate read counts
#'		\item \code{tRD}, tumor reference read counts
#'		\item \code{CN_1}, minor allele count
#'		\item \code{CN_2}, major allele count, where \code{CN_1} <= \code{CN_2}
#'		\item \code{tCN} = \code{CN_1 + CN_2}
#' }
#' @param my_purity A single numeric value of known/estimated purity
#' @param init_eS A subclone configuration matrix pre-defined in R list \code{eS}
#' @param my_unc_q An optimal initial vector for the unconstrained \code{q} vector, 
#'		useful after running \code{grid_ITH_optim}. If this variable is \code{NULL},
#'		then the subclone proportions, \code{q}, are randomly generated.
#' @param max_iter Positive integer, preferably 1000 or more, setting the maximum number of iterations
#' @param my_epsilon Convergence criterion threshold for changes in the log likelihood, preferably 1e-6 or smaller
#' @return If the EM algorithm converges, the output will be a list containing
#' \itemize{
#'		\item \code{iter}, number of iterations
#'		\item \code{converge}, convergence status
#'		\item \code{unc_q0}, initial unconstrained subclone proportions parameter
#'		\item \code{unc_q}, unconstrained estimate of \code{q}
#'		\item \code{q}, estimated subclone proportions among cancer cells
#'		\item \code{CN_MA_pi}, estimated mixture probabilities of multiplicities and allocations given copy number states
#'		\item \code{eta}, estimated subclone proportion among tumor cells
#'		\item \code{purity}, user-inputted tumor purity
#'		\item \code{entropy}, estimated entropy
#'		\item \code{infer}, R data frame of inferred variant allocations 
#'			\code{infer_A}, multiplicities \code{infer_M}, cellular prevalences \code{infer_CP}
#'		\item \code{ms}, model size, number of parameters within parameter space
#'		\item \code{LL}, observed log likelihood evaluated at maximal estimates
#'		\item \code{AIC = 2 * LL - 2 * ms}, Negative AIC, used for model selection
#'		\item \code{BIC = 2 * LL - ms * log(LOCI)}, Negative BIC, used for model selection
#'		\item \code{LOCI}, number of inputted somatic variants
#' }
#' @export
ITH_optim = function(my_data,my_purity,init_eS,my_unc_q=NULL,max_iter=4e3,my_epsilon=1e-6){

	if( !is.data.frame(my_data) ){
		stop("my_data is not a data.frame")
	}
	req_vars = c("tCN","CN_1","CN_2","tAD","tRD")
	miss_vars = req_vars[!req_vars %in% names(my_data)]
	if( length(miss_vars) > 0 ){
		stop(paste0(paste(miss_vars,collapse=", "),"are missing from my_data"))
	}
	
	# Initialize Params
	rownames(my_data) = NULL
	uniq_CN = unique(my_data[,c("tCN","CN_2","CN_1")])
	uniq_CN = uniq_CN[order(uniq_CN$tCN,uniq_CN$CN_2),]
	uniq_CN_MA = c()
	for(ii in seq(nrow(uniq_CN))){
		one_row = uniq_CN[ii,c("tCN","CN_2","CN_1")]
		for(kk in seq(nrow(init_eS))){
			tmp_mult = poss_mult(as.numeric(one_row[,c("CN_1","CN_2")]),alloc=kk)
		for(jj in seq(length(tmp_mult))){
			uniq_CN_MA = rbind(uniq_CN_MA,smart_df(one_row,M=tmp_mult[jj],A=kk))
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
	
	if( is.null(my_unc_q) ){
		init_unc_q = runif(ncol(init_eS)-1,-3,1)
	} else {
		init_unc_q = my_unc_q
	}
	unc_q0 = init_unc_q
	
	# Create B matrix: Stores 0s and 1s for what's feasible given copy number state
	num_LOCI = nrow(my_data)
	my_mat_B = matrix(0,num_LOCI,nrow(uniq_CN_MA))
	for(bb in seq(num_LOCI)){
		one_row = my_data[bb,]
		one_index = which(uniq_CN_MA$tCN == one_row$tCN & uniq_CN_MA$CN_2 == one_row$CN_2)
		my_mat_B[bb,one_index] = 1
	}

	# Optimize
	if( ncol(init_eS) == 1 ){
		ith_out = Rcpp_ITH_opt(RD = as.matrix(my_data[,c("tAD","tRD")]),
			BB = my_mat_B,uniq_BB = unique(my_mat_B),
			uniq_CN_MA = as.matrix(uniq_CN_MA),eS = init_eS,
			purity = my_purity,tCN = my_data$tCN,
			pi0 = init_pi_CN_MA,unc_qq0 = unc_q0,
			mstep = FALSE,max_iter = max_iter,
			eps = my_epsilon,show = FALSE)
	} else {
		ith_out = Rcpp_ITH_opt(RD = as.matrix(my_data[,c("tAD","tRD")]),
			BB = my_mat_B,uniq_BB = unique(my_mat_B),
			uniq_CN_MA = as.matrix(uniq_CN_MA),eS = init_eS,
			purity = my_purity,tCN = my_data$tCN,
			pi0 = init_pi_CN_MA,unc_qq0 = unc_q0,
			mstep = TRUE,max_iter = max_iter,
			eps = my_epsilon,show = FALSE)
	}
	
	# Output
	if( ith_out$converge == 1 ){
		uniq_CN_MA$infer_MAF = apply(uniq_CN_MA,1,function(x) 
			calc_maf(my_purity,ith_out$qq,x[4],x[1],init_eS[x[5],]))
		uniq_CN_MA$infer_CP = sapply(uniq_CN_MA$A,function(x) 
			sum(init_eS[x,] * ith_out$qq))
		infer_out = uniq_CN_MA[apply(ith_out$Z,1,which.max),c("M","A","infer_MAF","infer_CP")]
			rownames(infer_out) = NULL
			names(infer_out)[1:2] = paste0("infer_",names(infer_out)[1:2])
		list(converge=TRUE,unc_q0=unc_q0,iter=ith_out$iter,q=ith_out$qq,
			CN_MA_pi=smart_df(uniq_CN_MA,pi=ith_out$CN_MA_pi),
			eta=ith_out$qq*my_purity,purity=my_purity,
			entropy=round(ith_out$entropy,2),infer=infer_out,
			ms=ith_out$ms,LL=ith_out$LL,AIC=ith_out$AIC,
			BIC=ith_out$BIC,LOCI=num_LOCI)
	} else {
		list(converge=FALSE)
	}
	
}

#' @title grid_ITH_optim
#' @description This function performs a grid search over enumerated 
#' 	configurations within the pre-defined list \code{eS}
#' @inheritDotParams ITH_optim my_data my_purity max_iter my_epsilon
#' @param list_eS A nested list of subclone configuration matrices
#' @param trials Positive integer, number of random initializations of subclone proportions
#' @return A R list containing two objects. \code{GRID} is a data frame 
#'		where each row denotes a feasible subclone configuration with 
#'		corresponding subclone proportion estimates \code{q}
#'		and somatic variant allocations \code{alloc}. \code{INFER} 
#'		is a list where \code{INFER[[i]]} corresponds to the i-th 
#'		row or model of \code{GRID}.
#' @export
#' @examples
#' data(eS)
#' set.seed(1); truth = gen_subj_truth(mat_config = eS[[2]][[1]],maxLOCI = 100)
#' truth[c("purity","q","eta")] # the underlying parameters
#' obs_data = gen_ITH_RD(DATA = truth$subj_truth,RD = 200)
#' obs_data = data.frame(obs_data,truth$subj_truth,stringsAsFactors=FALSE)
#' opt_out = grid_ITH_optim(my_data = obs_data[,c("tAD","tRD","CN_1","CN_2","tCN")],
#'		my_purity = truth$purity,list_eS = eS,trials = 50)
#' opt_out$GRID
#' opt_out$GRID[order(-opt_out$GRID$BIC),] # sort from best to worst model
#' opt_out$INFER[[2]][1:5,] # inferred results for the first five somatic variants and second model presented in opt_out$GRID
grid_ITH_optim = function(my_data,my_purity,list_eS,trials=20,max_iter=4e3,my_epsilon=1e-6){
	
	# Run EM on all provided subclone configurations
	all_ck = c()
	for(cc in seq(length(list_eS))){
		cat(paste0("cc = ",cc,"; "))
		tmp_ck0 = c()
		for(kk in seq(length(list_eS[[cc]]))){
			cat(paste0("kk = ",kk," "))
			tmp_ck = c()
			for(count in seq(trials)){
				if(count %% 5 == 0) cat("*")
				ith_out = ITH_optim(my_data = my_data,my_purity = my_purity,
					init_eS = list_eS[[cc]][[kk]], my_unc_q = NULL, 
					max_iter = max_iter,my_epsilon = my_epsilon)
				
				if( ith_out$converge ){
					# Even though we have convergence, check all subclones have an allocation
					tab_A = table(ith_out$infer$infer_A)
					len_A = length(tab_A)
					if( len_A == cc && min(round(ith_out$q,2)) > 0 ){
						tmp_ck = rbind(tmp_ck,smart_df(cc=cc,kk=kk,ms=ith_out$ms,
							entropy=ith_out$entropy,
							LL=ith_out$LL,AIC=ith_out$AIC,BIC=ith_out$BIC,
							q=paste(round(ith_out$q,2),collapse=","),
							unc_q0=paste(round(ith_out$unc_q0,6),collapse=","),
							alloc=paste(paste0(names(tab_A),"|",tab_A),collapse=";")))
					}
				}
			}
			cat(" ")
			if( !is.null(tmp_ck) && nrow(tmp_ck) > 0 ){
				tmp_ck = tmp_ck[!duplicated(tmp_ck[,c("cc","kk","ms","entropy","LL")]),]
				tmp_ck0 = rbind(tmp_ck0,tmp_ck)
			}
		
		}
		
		if( !is.null(tmp_ck0) && nrow(tmp_ck0) > 0 ){
			all_ck = rbind(all_ck,tmp_ck0)
		}
		cat("\n")
	}
	all_ck = all_ck[order(all_ck$cc,all_ck$kk,all_ck$ms),]
	rownames(all_ck) = NULL
	# all_ck
	
	# Obtain details of inferred results
	infer_list = list()
	for(ii in seq(nrow(all_ck))){
		if( all_ck$cc[ii] == 1 ){
			ith_out = ITH_optim(my_data = my_data,
				my_purity = my_purity,
				init_eS = list_eS[[all_ck$cc[ii]]][[all_ck$kk[ii]]],
				my_unc_q = NULL,
				max_iter = max_iter,
				my_epsilon = my_epsilon)
		} else if( all_ck$cc[ii] > 1 ){
			ith_out = ITH_optim(my_data = my_data,
				my_purity = my_purity,
				init_eS = list_eS[[all_ck$cc[ii]]][[all_ck$kk[ii]]],
				my_unc_q = as.numeric(strsplit(all_ck$unc_q0[ii],",")[[1]]),
				max_iter = max_iter,
				my_epsilon = my_epsilon)
			if( !ith_out$converge ){
				while(TRUE){
					ith_out = ITH_optim(my_data = my_data,
						my_purity = my_purity,
						init_eS = list_eS[[all_ck$cc[ii]]][[all_ck$kk[ii]]],
						my_unc_q = NULL,
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
	
	list(GRID = all_ck,INFER = infer_list)
}



###


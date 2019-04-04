#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends("RcppArmadillo")]]

template<typename T>
void printR_obj(const T& obj){
	Rcpp::Rcout << obj << std::endl;	
}

// --------------------
// Intermediate Functions
// --------------------

// [[Rcpp::export]]
double Rcpp_norm(const arma::vec& a){
	return arma::norm(a);
}

// [[Rcpp::export]]
double Rcpp_logSumExp(const arma::vec& log_x){
	if( log_x.n_elem == 1 ){
		return log_x.at(0);
	} else {
		double max_val = arma::max(log_x);
		arma::vec log_x_2 = log_x - max_val;
		return std::log(arma::sum(arma::exp(log_x_2))) + max_val;
	}
}


// --------------------
// Minor ITH functions
// --------------------

// [[Rcpp::export]]
arma::vec Rcpp_vartheta_to_eta(const arma::vec& x){
	if(x.n_elem == 1){
		return 1.0 / (1.0 + arma::exp(-1.0 * x));
	} else {
		arma::vec exp_x = arma::exp(x);
		return exp_x / (1.0+arma::sum(exp_x));
	}
}

// [[Rcpp::export]]
arma::vec Rcpp_vartheta_to_qq(const arma::vec& xx){
	arma::uword num_q = xx.n_elem;
	arma::vec qq = arma::zeros<arma::vec>(num_q + 1);
	qq.subvec(0,num_q-1) = Rcpp_vartheta_to_eta(xx);
	qq.at(num_q) = 1.0 - arma::sum(qq.subvec(0,num_q-1));
	return qq;
}

double Rcpp_calc_entropy2(const arma::vec& qq){
	// return -arma::sum(qq % arma::log(qq));
	return -1.0 * arma::dot(qq, arma::log(qq));
}

arma::mat Rcpp_calc_jacobian2(const arma::vec& eta){
	arma::uword ii, jj, D = eta.n_elem;
	arma::mat out_jacob = arma::zeros<arma::mat>(D,D);
	
	for(ii = 0; ii < D; ii++){
	for(jj = 0; jj < D; jj++){
		out_jacob.at(ii,jj) = eta.at(ii)*((ii==jj) - eta.at(jj));
	}
	}
	return out_jacob;
}

double Rcpp_calc_maf(const double& purity,
  const arma::vec& q,const double& mult,
	const double& tCN,const arma::vec& u){
  
	return purity * mult * arma::dot(q,u) / 
		(tCN * purity + 2.0 * (1.0 - purity));
}

double Rcpp_round(const double& aa,
  const arma::uword& digits){
  
	double factor = std::pow(10,digits);
	double out = std::round(aa * factor) / factor;
	return out;
}


// --------------------
// Optimization Functions
// --------------------

double Rcpp_ITH_obsLL_Z(const arma::mat& RD,
  const arma::vec& log_DP,const arma::vec& LBC,arma::mat& ZZ,
  const arma::mat& uniq_CN_MA,const arma::mat& log_BB,
  const double& pi_eps,const arma::vec& pi,const arma::mat& eS,
  const double& purity,const arma::vec& qq){
  
  arma::uword vv, bb;
  double maf,tmp_num;
  
  for(vv = 0; vv < log_BB.n_cols + 1; vv++){
    if( vv == 0 ){
      ZZ.col(vv) = std::log(pi_eps) - log_DP;
    } else {
      maf = Rcpp_calc_maf(purity,
        qq,
        uniq_CN_MA.at(vv-1,3),
        uniq_CN_MA.at(vv-1,0),
        eS.row( uniq_CN_MA.at(vv-1,4)-1 ).t());
      ZZ.col(vv) = LBC + log_BB.col(vv-1) + std::log(1.0 - pi_eps) + std::log(pi.at(vv-1)) +
        RD.col(0) * std::log(maf) + RD.col(1) * std::log(1.0 - maf);
    }
  }
  
  // Calc obs_LL and E-step
  double out_LL = 0.0;
  for(bb = 0; bb < RD.n_rows; bb++){
    tmp_num = Rcpp_logSumExp(ZZ.row(bb).t());
    out_LL += tmp_num;
    ZZ.row(bb) = arma::exp(ZZ.row(bb).t() - tmp_num).t();
  }
  
  return out_LL;
}

double Rcpp_ITH_compLL(const arma::mat& RD,
  const arma::mat& ZZ,const arma::mat& uniq_CN_MA,
  const arma::vec& pi,const arma::mat& eS,
	const double& purity,const arma::vec& qq){
	
	arma::uword dd;
	double maf, out_LL = 0.0;
	
	for(dd = 0; dd < uniq_CN_MA.n_rows; dd++){
		maf = Rcpp_calc_maf(purity,
			qq,
			uniq_CN_MA.at(dd,3),
			uniq_CN_MA.at(dd,0),
			eS.row( uniq_CN_MA.at(dd,4)-1 ).t());
		// small change here for if pi.at(dd) == 0, don't sum
		if( pi.at(dd) > 0 ){
			out_LL += arma::sum(
				ZZ.col(dd + 1) % ( std::log(pi.at(dd)) + 
					RD.col(0) * std::log(maf) + RD.col(1) * std::log(1.0 - maf) )
				);
		}
	}
	
	return out_LL;
}

arma::vec Rcpp_ITH_compGRAD(const arma::mat& RD,
  const arma::mat& ZZ,const arma::mat& uniq_CN_MA,
  const arma::mat& eS,const double& purity,
	const arma::vec& qq,const arma::vec& tCN){
	// partial Q / partial unc_q = partial q* / partial unc_q * partial q / partial q* * partial Q / partial q
	
	// Note: This function will only apply to length(qq) > 1 aka 2 or more subclones
	arma::uword num_qq = qq.n_elem - 1;
	arma::mat jacob = Rcpp_calc_jacobian2(qq.subvec(0,num_qq - 1));
	arma::mat part_q_qstar = arma::zeros<arma::mat>(num_qq,num_qq+1);
		part_q_qstar.fill(-1.0);
		part_q_qstar.cols(0,num_qq-1) = arma::eye<arma::mat>(num_qq,num_qq);

	arma::uword bb,dd;
	double maf, part_maf;
	arma::vec out_GRAD = arma::zeros<arma::vec>(num_qq + 1);
	
	for(bb = 0; bb < RD.n_rows; bb++){
	for(dd = 0; dd < uniq_CN_MA.n_rows; dd++){
		maf = Rcpp_calc_maf(purity,
			qq,
			uniq_CN_MA.at(dd,3),
			uniq_CN_MA.at(dd,0),
			eS.row( uniq_CN_MA.at(dd,4)-1 ).t());
		part_maf = purity * uniq_CN_MA.at(dd,3) / (purity * tCN.at(bb) + 2.0 * (1.0 - purity));
		out_GRAD += ZZ.at(bb,dd + 1) * 
			( RD.at(bb,0) / maf - RD.at(bb,1) / (1.0 - maf) ) * 
			part_maf * eS.row( uniq_CN_MA.at(dd,4)-1 ).t();
	}
	}

	return jacob * part_q_qstar * out_GRAD;
}

void Rcpp_ITH_Mstep_BFGS(const arma::mat& RD,
  const arma::mat& ZZ,const arma::mat& uniq_CN_MA,
  const arma::vec& pi,const arma::mat& eS,
  const double& purity,const arma::vec& tCN,
  const arma::vec& params0,arma::vec& qq,
  const arma::uword& max_iter = 4e3,
  const double& eps = 1e-7,const bool& show = true){
  
  arma::uword np = params0.n_elem,
    iter = 0, jj, uu;
  
  arma::vec xk = params0;
  arma::mat I_np = arma::eye<arma::mat>(np,np);
  arma::mat inv_Bk = I_np;
  arma::vec curr_xk = arma::zeros<arma::vec>(np);
  arma::vec new_xk = curr_xk, gr_k = curr_xk,
    p_k = curr_xk, s_k = curr_xk, y_k = curr_xk;
  arma::mat ISYT = arma::zeros<arma::mat>(np,np);
  
  double old_LL,new_LL,inv_norm_p_k,tmp_alpha,ys;
  double fnscale = -1.0; // For maximization
  double curr_LL = 0.0;
  arma::vec tmp_qq = arma::zeros<arma::vec>(np+1);
  
  while(iter < max_iter){
    // Check any qq elements too small
    tmp_qq = Rcpp_vartheta_to_qq(xk);
    if( arma::any(tmp_qq < 1e-3) ){
      break;
    }		
    
    // Calculate Direction p_k
    gr_k = fnscale * Rcpp_ITH_compGRAD(RD,ZZ,uniq_CN_MA,
      eS,purity,tmp_qq,tCN);
    p_k = -1.0 * inv_Bk * gr_k;
    inv_norm_p_k = 1.0 / std::max(1.0,Rcpp_norm(p_k));
    
    // Line search for new xk
    uu = 0;
    old_LL = fnscale * Rcpp_ITH_compLL(RD,ZZ,uniq_CN_MA,
      pi,eS,purity,tmp_qq);
    for(jj = 0; jj <= 30; jj++){
      tmp_alpha = inv_norm_p_k / std::pow(4,jj);
      new_xk = xk + tmp_alpha * p_k;
      tmp_qq = Rcpp_vartheta_to_qq(new_xk);
      new_LL = fnscale * Rcpp_ITH_compLL(RD,ZZ,uniq_CN_MA,
        pi,eS,purity,tmp_qq);
      if(new_LL < old_LL){ // minimizing
        s_k = tmp_alpha * p_k;
        y_k = fnscale * Rcpp_ITH_compGRAD(RD,ZZ,uniq_CN_MA,
          eS,purity,tmp_qq,tCN) - gr_k;
        ys = arma::dot(y_k,s_k);
        if( ys > 0.0 ){
          if(show) Rcpp::Rcout << "Update xk and inv_Bk, ";
          ISYT = I_np - (s_k * y_k.t()) / ys;
          inv_Bk = ISYT * inv_Bk * ISYT.t() + s_k * s_k.t() / ys;
        } else {
          if(show) Rcpp::Rcout << "Update xk only, ";
        }
        xk = new_xk;
        old_LL = new_LL;
        uu = 1;
        break;
      }
    }
    
    if( uu == 0 ){ // aka no update
      if( Rcpp_norm(gr_k) > 1.0 ){
        if(show) Rcpp::Rcout << "Reset inv_Bk, ";
        inv_Bk = I_np;
      } else {
        if(show) Rcpp::Rcout << "Failed line search, ";
        break;
      }
    }
    
    // Check Convergence
    if( iter > 0 ){
      if( std::abs(curr_LL - old_LL) < eps &&
          Rcpp_norm(curr_xk - xk) < eps ){
        gr_k = Rcpp_ITH_compGRAD(RD,ZZ,uniq_CN_MA,
          eS,purity,Rcpp_vartheta_to_qq(xk),tCN);
        if( Rcpp_norm(gr_k) < eps ){
          break;
        }
      }
    }
    
    curr_xk = xk;
    curr_LL = old_LL;
    iter++;
  }
  
  if(show) Rcpp::Rcout << std::endl;
  qq = Rcpp_vartheta_to_qq(xk);
}

// [[Rcpp::export]]
Rcpp::List Rcpp_ITH_opt(const arma::mat& RD,
  const arma::vec& log_DP,const arma::vec& LBC,const arma::mat& BB,
  const arma::umat& uniq_BB,const arma::mat& uniq_CN_MA,
  const arma::mat& eS,const double& purity,
  const arma::vec& tCN,const double& pi_eps0,const arma::vec& pi0,
  const arma::vec& unc_qq0,const bool& mstep = true,
  const arma::uword& max_iter = 4e3,
  const double& eps = 1e-6,const bool& show = true){
  
  arma::uword LOCI = RD.n_rows, DD = uniq_CN_MA.n_rows + 1,
    num_uniq_BB = uniq_BB.n_rows, num_ss = eS.n_cols,
    converge = 0, bb;
  double old_LL, curr_LL;
  arma::mat ZZ = arma::zeros<arma::mat>(LOCI,DD);
  arma::mat log_BB = arma::log(BB);
  arma::vec colsum_ZZ = arma::zeros<arma::vec>(DD);
  arma::mat CLASS = arma::zeros<arma::mat>(LOCI,DD);
  arma::vec freq_class = arma::zeros<arma::vec>(DD);
  
  arma::vec unc_qq = unc_qq0;
  arma::vec tmp_qq = arma::ones<arma::vec>(num_ss);
  double pi_eps = pi_eps0;
  arma::uword reset_opt = 0;
  arma::vec pi = pi0;
  arma::vec curr_unc_qq = unc_qq;
  double curr_pi_eps = pi_eps;
  arma::vec curr_pi = pi;
  curr_LL = 0.0;
  
  arma::uword iter = 0;
  while( iter < max_iter ){
    // E-step
    if(num_ss == 1){
      old_LL = Rcpp_ITH_obsLL_Z(RD,log_DP,LBC,
        ZZ,uniq_CN_MA,log_BB,pi_eps,pi,eS,purity,tmp_qq);
    } else {
      old_LL = Rcpp_ITH_obsLL_Z(RD,log_DP,LBC,
        ZZ,uniq_CN_MA,log_BB,pi_eps,pi,eS,purity,Rcpp_vartheta_to_qq(unc_qq));
    }
    
    // Check convergence
    if( iter > 0 ){
      if( std::abs(old_LL - curr_LL) < eps ){
        if( std::abs(pi_eps - curr_pi_eps) < eps 
              && Rcpp_norm(pi - curr_pi) < eps 
              && Rcpp_norm(unc_qq - curr_unc_qq) < eps ){
              converge = 1;
          break;
        }
      }
    }
    
    // M-steps
    // Update pi
    colsum_ZZ = arma::sum(ZZ,0).t();
    pi_eps = colsum_ZZ.at(0) / LOCI;
    for(bb = 0; bb < num_uniq_BB; bb++){
      arma::uvec tmp_index = arma::find(uniq_BB.row(bb).t() == 1);
      pi(tmp_index) = colsum_ZZ(tmp_index + 1) / arma::sum(colsum_ZZ(tmp_index + 1));
    }
    
    // Update unc_qq: If num_ss > 1
    if( mstep && num_ss > 1 ){
      Rcpp_ITH_Mstep_BFGS(RD,ZZ,uniq_CN_MA,
        pi,eS,purity,tCN,unc_qq,tmp_qq,4e3,1e-7,show);
      
      if( arma::any(tmp_qq < 1e-3) ){
        break;
      }
      unc_qq = arma::log(tmp_qq.subvec(0,num_ss - 2) / tmp_qq.at(num_ss - 1));
    }
    
    // Code to set pi_eps and some pi's to zero if nothing classified
    CLASS.fill(0.0);
    for(bb = 0; bb < LOCI; bb++){
      CLASS.at(bb,arma::index_max(ZZ.row(bb).t())) = 1;
    }
    freq_class = arma::sum(CLASS,0).t();
    if( freq_class.at(0) == 0 ) pi_eps = 0.0;
    for(bb = 0; bb < num_uniq_BB; bb++){
      if( arma::any(freq_class(arma::find(uniq_BB.row(bb).t() == 1) + 1) == 0.0) ){
        arma::uvec one_index = arma::find(uniq_BB.row(bb).t() == 1);
        arma::vec tmp_prop = pi(one_index);
        tmp_prop(arma::find(freq_class(one_index + 1) == 0)).fill(0.0);
        if( arma::sum(tmp_prop) > 0.0 ){
          tmp_prop /= arma::sum(tmp_prop);
        } else {
          if(show) Rcpp::Rcout << "Reset optimization but set pi_eps = 0\n";
          iter = 0;
          pi_eps = 0.0;
          pi = pi0;
          unc_qq = unc_qq0;
          reset_opt = 1;
          break;
        }
        pi(one_index) = tmp_prop;
      }
    }
    
    if( reset_opt == 1 ){
      reset_opt = 0;
      continue;
    }
    
    if(show){
      Rcpp::Rcout << "old_LL = " << old_LL << "\n";
      Rcpp::Rcout << "pi_eps = " << pi_eps << "\n";
      Rcpp::Rcout << "pi = \n\t" << pi.t() << "\n";
    }
    
    curr_LL = old_LL;
    curr_unc_qq = unc_qq;
    curr_pi_eps = pi_eps;
    curr_pi = pi;
    iter++;
  }
  
  arma::vec qq = arma::ones<arma::vec>(num_ss);
  if( num_ss == 1 ){
    old_LL = Rcpp_ITH_obsLL_Z(RD,log_DP,LBC,
      ZZ,uniq_CN_MA,BB,pi_eps,pi,eS,purity,qq);
  } else {
    qq = Rcpp_vartheta_to_qq(unc_qq);
    old_LL = Rcpp_ITH_obsLL_Z(RD,log_DP,LBC,
      ZZ,uniq_CN_MA,BB,pi_eps,pi,eS,purity,qq);
  }
  
  // Calculate num_params
  arma::uword num_params = 1*(pi_eps > 0.0) + num_ss - 1;
  for(bb = 0; bb < num_uniq_BB; bb++){
    num_params += arma::sum(1*(uniq_BB.row(bb).t() % pi > 0.0)) - 1;
  }
  
  old_LL = Rcpp_round(old_LL,2);
  double BIC = 2.0 * old_LL - num_params * std::log(LOCI);
  double AIC = 2.0 * old_LL - 2.0 * num_params;
  double entropy = Rcpp_calc_entropy2(qq);
  
  return Rcpp::List::create(
    Rcpp::Named("converge",converge),
    Rcpp::Named("LL",old_LL),
    Rcpp::Named("AIC",AIC),
    Rcpp::Named("BIC",BIC),
    Rcpp::Named("iter",iter),
    Rcpp::Named("ms",num_params),
    Rcpp::Named("Z",ZZ),
    Rcpp::Named("entropy",entropy),
    Rcpp::Named("pi_eps",pi_eps),
    Rcpp::Named("CN_MA_pi",Rcpp::NumericVector(pi.begin(),pi.end())),
    Rcpp::Named("qq",Rcpp::NumericVector(qq.begin(),qq.end())));
}




// --------------------



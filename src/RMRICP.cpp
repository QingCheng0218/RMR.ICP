#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>

#include "GibbsAlpGamEtaW_ptr.hpp"

#include <ctime>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
umat blockfun(ivec NB){
  uword nblocks = NB.size();
  umat block_inf = ones<umat>(nblocks, 2);
  ivec NBinf = cumsum(NB);
  for(int l = 0; l < nblocks; l++){
    if(l==0){
      block_inf(l, 0) = 0;
    }else{
      block_inf(l, 0) = NBinf[l - 1];
    }
    block_inf(l, 1) = NBinf[l] - 1;
  }
  return block_inf;
}

//[[Rcpp::export]]
double normal_pdf(double x, double m, double s)
{
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

ObjRMRindep RMRindepObj(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2,
                        Options_RMRindep* opts){
  // ----------------------------------------------------------------------
  // check number of input arguments
  double agm = opts -> agm;
  double bgm = opts -> bgm;
  double aal = opts -> aal;
  double bal = opts -> bal;
  double a = opts -> a;
  double b = opts -> b;
  double v = opts -> v;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  // ----------------------------------------------------------------------
  // initial values
  int p = Gammah.n_elem;
  double sgga2 = 0.001; double sgal2 = 0.001; double beta1 = 0.01; double ome = 0.1;
  double beta2 = 0.01; double logw = log(ome  / (1 - ome));
  int numsave = maxIter / thin;
  vec Beta1res = ones(numsave, 1);
  vec Beta2res = ones(numsave, 1);
  vec Sgal2Res = ones(numsave, 1);
  vec Sgga2Res = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  
  vec mu = 0.01*ones(p, 1);
  vec muA = 0.01*ones(p, 1);
  vec W = ones(p, 1);
  ivec Eta = zeros<ivec>(p, 1);
  // ----------------------------------------------------------------------
  vec sG2 = se2%se2;
  vec sg2 = se1%se1;
  vec invsG2 = 1. / sG2;
  vec invsg2 = 1. / sg2;
  vec GinvsG2 = Gammah / sG2;
  vec ginvsg2 = gammah / sg2;
  
  double b12 = beta1*beta1;
  double b22 = beta2*beta2;
  int l = 0;
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    double invsgga2 = 1. / sgga2;
    double invsgal2 = 1. / sgal2;
    vec Winsgal2 = W * invsgal2;
    
    vec v21, v20, v2A = zeros(p, 1);
    v20 = 1. / (invsg2 + b12*invsG2 + invsgga2);
    v21 = 1. / (invsg2 + b22*invsG2 + invsgga2);
    v2A = 1. / (invsG2 + Winsgal2);
    
    // cout<<"v20:"<<sum(v20)<<"  v21:" << sum(v21) << "  v2A:"<< sum(v2A)<<endl;
    //  ------------------------------------------------------- 
    //  check variance of Gamma to decide which one is beta1.
    //  var2 should be larger than var1.
    // double var1, var2, betatmp;
    
    // if (sum(Eta) != 0 && sum(Eta) != p) {
    //   var1 = var(beta1*mu.elem(Eta == 0));
    //   var2 = var(beta2*mu.elem(Eta == 1) + muA.elem(Eta == 1));
    //   
    //   if (var1 > var2) {
    //     Rcpp::Rcout << "Reorder beta1 and beta2" << std::endl;
    //     betatmp = beta1;
    //     beta1 = beta2;
    //     beta2 = betatmp;
    //   }
    // }
    // ------------------------------------------------- #
    // for b0k and D0k
    
    vec dm0, invdm0, bm0, bm02, dm0Winvsgal2;
    vec Dm111, Dm122, Dm112;
    mat bm1;
    
    dm0 = invsgga2 + invsg2 + b12*invsG2;
    invdm0 = 1. / dm0;
    bm0 = ginvsg2 + beta1*GinvsG2;
    bm02 = bm0 % bm0;
    dm0Winvsgal2 = dm0 % Winsgal2;
    
    // for b1k and D1k
    bm1 = join_rows(ginvsg2 + beta2*GinvsG2, GinvsG2);
    Dm111 = invsgga2 + invsg2 + b22*invsG2;
    Dm122 = Winsgal2 + invsG2;
    Dm112 = beta2 * invsG2;
    mat Dm1j = zeros(2, 2);
    
    // ------------------------------------------------- #
    double mu1, mu0, muA0, wa, wb, lik0, lik1, prob0, prob;
    for(int j = 0; j < p; j++){
      // ------------------ //
      // update gamma
      // ------------------ //
      if (Eta[j] == 1) {
        mu1 = (ginvsg2[j] + beta2 * GinvsG2[j] - beta2 * muA[j] * invsG2[j]) * v21[j];
        mu[j] = mu1 + randn()*sqrt(v21[j]);
        // mu[j] = mu1;
      } else {
        mu0 = (ginvsg2[j] + beta1 * GinvsG2[j]) * v20[j];
        mu[j] = mu0 + randn()*sqrt(v20[j]);
        // mu[j] = mu0;
      }
      // ------------------ //
      // update alpha
      // ------------------ //
      
      if (Eta[j] == 1) {
        muA0 = (GinvsG2[j] - beta2 * invsG2[j] * mu[j]) * v2A[j];
        muA[j] = muA0 + randn()*sqrt(v2A[j]);
        // muA[j] = muA0;
      } else {
        muA[j] = randn()*sqrt(sgal2 / W[j]);
        // muA[j] = 0;
      }
      
      // ------------------ //
      // update W
      // ------------------ //
      wa = (v + 1) / 2;
      wb = (v + muA[j] * muA[j] * invsgal2) / 2;
      W[j] =  randg<double>(distr_param(wa, 1./wb));
      // W[j] = 1;
      // ------------------ //
      // update eta
      // ------------------ //
      rowvec bm1j = bm1.row(j);
      Dm1j(0, 0) = Dm111[j];
      Dm1j(0, 1) = Dm112[j];
      Dm1j(1, 0) = Dm112[j];
      Dm1j(1, 1) = Dm122[j];
      mat invDm1j = inv(Dm1j);
      
      lik1 = 0.5*as_scalar(bm1j*invDm1j*bm1j.t()) - 0.5*log(det(Dm1j));
      lik0 = 0.5*invdm0[j]*bm02[j] - 0.5*log(dm0Winvsgal2[j]);
      prob0 = logw  + lik1 - lik0;
      prob = 1. / ( 1 + exp(-prob0));
      Eta[j] = R::rbinom(1, prob);
      
    }
    
    // cout << "mu::" << mu.subvec(0, 4).t() << " sum::" << sum(mu) <<   endl;
    // cout <<  "muA::" << muA.subvec(0, 4).t() << " sum::"  << sum(muA) << endl;
    // ------------------ //
    // update beta1
    // ------------------ //
    vec muEta, muEta1, muAEta;
    muEta = mu%Eta;
    muEta1 = mu%(1 - Eta);
    muAEta = muA%Eta;
    
    // cout <<"muEta::" << sum(muEta) << "muEta1" << sum(muEta1) << "muAEta:" << sum(muAEta) << endl;
    
    double sig2b1, mub1;
    if(sum(Eta)==p){
      beta1 = 0;
    }else{
      sig2b1 = 1. / sum(muEta1%invsG2%muEta1);
      mub1 = sum(GinvsG2%muEta1)*sig2b1;
      beta1 = mub1 + randn()*sqrt(sig2b1);
      // beta1 = mub1;
      // cout <<"sib2b1: "<< sig2b1 << " sum::" << sum(GinvsG2%muEta1) << endl;
    }
    
    // ------------------ //
    // update beta2
    // ------------------ //
    double sig2b2, mub2;
    if(sum(Eta)==0){
      beta2 = 0;
    }else{
      sig2b2 = 1. / sum(muEta%invsG2%muEta);
      mub2 = (sum(GinvsG2%muEta) - sum(muEta%invsG2%muAEta))*sig2b2;
      beta2 = mub2 + randn()*sqrt(sig2b2);
      // beta2 = mub2;
      // cout <<" sig2b2: "<<  sig2b2 << " sum::" << (sum(GinvsG2%muEta) - sum(muEta%invsG2%muAEta)) << endl;
    }
    b12 = beta1*beta1;
    b22 = beta2*beta2;
    // cout << "beta1: " << beta1 <<  " beta2::" << beta2 << endl;
    // ------------------ //
    // update sgga2
    // ------------------ //
    double tagm, tbgm, taal, tbal;
    tagm = agm + p / 2;
    tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
    // sgga2 = tbgm;
    // ------------------ //
    // update sgal2
    // ------------------ //
    taal = aal + p / 2;
    tbal = as_scalar(muA.t()*muA) / 2 + bal;
    sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
    // sgal2 = tbal;
    // cout << " sgga2 : " << sgga2 << " sgal2::" << sgal2  << endl;
    // ------------------ //
    // update omega
    // ------------------ //
    // double wpa1, wpa2, w, logw;
    double wpa1, wpa2;
    wpa1 = a + sum(Eta);
    wpa2 = b + Eta.n_elem - sum(Eta);
    ome = R::rbeta(wpa1, wpa2);
    logw = log(ome  / (1 - ome));
    // cout << "logw::" << logw <<  endl;
    // 
    // cout << "iter:: " << iter << " sum::" << sum(Eta) << endl;
    // 
    // cout << "-----------------------------------" << endl;
    
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        Sgga2Res[l] = sgga2;
        Sgal2Res[l] = sgal2;
        Beta1res[l] = beta1;
        Beta2res[l] = beta2;
        EtaAll.col(l) = Eta;
        l += 1;
      }
    }
  }
  
  arma::colvec Etares = mean(arma::conv_to<arma::mat>::from(EtaAll), 1);
  
  
  ObjRMRindep obj;
  
  obj.Etares = Etares;
  obj.Beta1res = Beta1res;
  obj.Beta2res = Beta2res;
  obj.Sgga2Res = Sgga2Res;
  obj.Sgal2Res = Sgal2Res;
  
  return obj;
  
}

//[[Rcpp::export]]
List RMRICPindep(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1,
                  arma::vec &se2, SEXP opts = R_NilValue)
{
  
  Options_RMRindep* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_RMRindep(opt["agm"], opt["bgm"], opt["aal"], opt["bal"], opt["a"], opt["b"], opt["v"], 
                                  opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_RMRindep();
  }
  
  ObjRMRindep obj = RMRindepObj(gammah, Gammah, se1, se2, lp_opt);
  
  vec Beta1res = obj.Beta1res;
  double bhat = mean(Beta1res);
  double se = stddev(Beta1res);
  double pvalue = 2*(R::pnorm(abs(bhat / se), 0, 1, 0, 0));
  
  
  List output = List::create(
    Rcpp::Named("beta.hat") = bhat,
    Rcpp::Named("beta.se") = se,
    Rcpp::Named("beta.p.value") = pvalue,
    Rcpp::Named("Etares") = Rcpp::wrap(obj.Etares),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta2res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Sgal2Res") = Rcpp::wrap(obj.Sgal2Res)
    
  );
  
  return output;
}



ObjRMR RMRLDobj(arma::field<vec> F4gammah, arma::field<vec> F4Gammah,
                arma::field<vec> F4se1, arma::field<vec> F4se2,
                arma::field<mat> F4Rblock,  Options_RMR* opts)
{
  
  // ----------------------------------------------------------------------
  // check number of input arguments
  double agm = opts -> agm;
  double bgm = opts -> bgm;
  double aal = opts -> aal;
  double bal = opts -> bal;
  double a = opts -> a;
  double b = opts -> b;
  double v = opts -> v;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  uword coreNum = opts -> coreNum;
  
  // ----------------------------------------------------------------------
  // initial values
  double sgga2 = 1; double sgal2 = 1; 
  double beta1 = 0.01; double beta2 = 0.01; 
  double ome = 0.1; 
  
  uword nblocks = F4Rblock.size();
  ivec Eta = zeros<ivec>(nblocks, 1);
  // ivec Eta = ones<ivec>(nblocks, 1);
  
  ivec NB = zeros<ivec>(nblocks, 1);
  
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    NB[nn] = F4gammah(nn, 0).size();
  }
  int p = sum(NB);
  umat block_inf = blockfun(NB);
  // ------------------------------------------ //
  // save the result.
  int numsave = maxIter / thin;
  // imat EtaAll = ones<imat>(nblocks, numsave);
  imat EtaAll = zeros<imat>(nblocks, numsave);
  vec Beta1res = ones(numsave, 1);
  vec Beta2res = ones(numsave, 1);
  vec Sgal2Res = ones(numsave, 1);
  vec Sgga2Res = ones(numsave, 1);
  // ------------------------------------------ //
  // cout << "check error 0" << endl;
  field<vec> F4sg2(nblocks, 1), F4sG2(nblocks, 1), F4invsg2(nblocks, 1), F4invsG2(nblocks, 1), F4GinvsG2(nblocks, 1), F4ginvsg2(nblocks, 1);
  field<vec> F4diaginsGRinsG(nblocks, 1), F4diaginsgRinsg(nblocks, 1), F4RinsGmu(nblocks, 1), F4Rinsgmu(nblocks, 1), F4RinsGmuA(nblocks, 1);
  field<mat> F4insGRinsG(nblocks, 1), F4insgRinsg(nblocks, 1), F4insgRinsG(nblocks, 1), F4Rinsg(nblocks, 1), F4RinsG(nblocks, 1);
  field<vec> F4W(nblocks, 1), F4DinsGRinsG(nblocks, 1), F4mu(nblocks, 1), F4muA(nblocks, 1);
  
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    
    
    NB[nn] = block_inf(nn, 1) - block_inf(nn, 0) + 1;
    // cout << "NB[nn] ::" << NB[nn] << endl;
    vec se1_block = F4se1(nn, 0);
    vec se2_block = F4se2(nn, 0);
    vec sg2_block = pow(se1_block, 2);
    vec sG2_block = pow(se2_block, 2);
    vec bh1_block = F4gammah(nn, 0);
    vec bh2_block = F4Gammah(nn, 0);
    vec mu_block = 0.01*ones(NB[nn], 1);
    vec muA_block = 0.01*ones(NB[nn], 1);
    
    
    mat R_block =  symmatu(F4Rblock(nn, 0));
    
    F4mu(nn, 0) = mu_block;
    F4muA(nn, 0) = muA_block;
    F4DinsGRinsG(nn, 0) = diagvec(R_block) / sG2_block;
    F4GinvsG2(nn, 0) = bh2_block / sG2_block;
    F4ginvsg2(nn, 0) = bh1_block / sg2_block;
    F4invsg2(nn, 0) = 1. / sg2_block;
    F4invsG2(nn, 0) = 1. / sG2_block;
    
    F4insGRinsG(nn, 0) = diagmat(1. / se2_block)*R_block*diagmat(1. / se2_block);
    F4insgRinsg(nn, 0) = diagmat(1. / se1_block)*R_block*diagmat(1. / se1_block);
    
    F4Rinsg(nn, 0) = R_block*diagmat(1 / se1_block);
    F4RinsG(nn, 0) = R_block*diagmat(1 / se2_block);
    
    F4RinsGmu(nn, 0) = R_block*diagmat(1 / se2_block)*mu_block;
    F4Rinsgmu(nn, 0) = R_block*diagmat(1 / se1_block)*mu_block;
    F4RinsGmuA(nn, 0) = R_block*diagmat(1. / se2_block)*muA_block;
    
    F4W(nn, 0) = ones(NB[nn], 1);
  }
  
  
  int lt = 0;
  
  vec beta1mean = zeros(nblocks, 1); // mean for each block.
  vec beta1sig = zeros(nblocks, 1); // sigma for each block.
  vec beta2mean = zeros(nblocks, 1); // mean for each block.
  vec beta2sig = zeros(nblocks, 1); // sigma for each block.
  
  double invsgga2, invsgal2, logw;
  vec Mu2 = zeros(nblocks, 1);
  vec MuA2 = zeros(nblocks, 1);
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    logw = log(ome /( 1 - ome));
    invsgga2 = 1. / sgga2;
    invsgal2 = 1. / sgal2;
    
    // ------------------------------------------------------------------
    // set parallel computation for gamma, alpha, eta, W
    paraBlock_GamAlpEtaW parobj_GamAlpEtaW(nblocks, F4se1, F4se2, F4ginvsg2, F4GinvsG2, F4mu, F4muA,F4Rblock,
                                           F4DinsGRinsG, F4insgRinsg, F4insGRinsG,F4Rinsg, F4RinsG,
                                           F4Rinsgmu, F4RinsGmu, F4RinsGmuA, F4W, F4invsg2, F4invsG2,
                                           beta1sig, beta2sig, beta1mean, beta2mean, Eta, Mu2, MuA2,
                                           beta1, beta2, logw, invsgga2, invsgal2, v);
    
    
    const int n_thread = coreNum;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_GamAlpEtaW::update_by_thread_GamAlpEtaW, &parobj_GamAlpEtaW, i_thread);
    }
    
    // cout << "check error 1::" << endl;
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    // cout << "check error 0::" << endl;
    // save the parallel result
    beta1sig =  parobj_GamAlpEtaW.beta1sig;
    beta2sig =  parobj_GamAlpEtaW.beta2sig;
    
    beta1mean =  parobj_GamAlpEtaW.beta1mean;
    beta2mean =  parobj_GamAlpEtaW.beta2mean;
    
    Mu2 = parobj_GamAlpEtaW.Mu2;
    MuA2 = parobj_GamAlpEtaW.MuA2;
    
    F4mu = parobj_GamAlpEtaW.F4mu;
    F4muA = parobj_GamAlpEtaW.F4muA;
    F4Rinsgmu = parobj_GamAlpEtaW.F4Rinsgmu;
    F4RinsGmu = parobj_GamAlpEtaW.F4RinsGmu;
    F4RinsGmuA = parobj_GamAlpEtaW.F4RinsGmuA;
    F4W = parobj_GamAlpEtaW.F4W;
    Eta = parobj_GamAlpEtaW.Eta;
    
    // cout << "beta1sig :" << beta1sig .t() << endl;
    // cout << "MuA2:" << MuA2.t() << endl;
    // --------------------------------------------- #
    // update beta1
    // --------------------------------------------- #
    double sig2b1, mub1;
    if(sum(Eta)==(int)nblocks){
      beta1 = 0;
    }else{
      sig2b1 = 1. / sum(beta1sig);
      mub1 = sum(beta1mean) * sig2b1;
      beta1 = mub1 + randn()*sqrt(sig2b1);
      // beta1 = mub1;
    }
    // cout<<"Mb0:"<< beta0 << "--sigb0:" << sqrt(sig2b0) << endl;
    // --------------------------------------------- #
    // update beta2
    // --------------------------------------------- #
    double sig2b2, mub2;
    if(sum(Eta)==0){
      beta2 = 0;
    }else{
      sig2b2 = 1. / sum(beta2sig);
      mub2 = sum(beta2mean) * sig2b2;
      // cout << "sum(beta2mean): "<< sum(beta2mean)<<endl;
      // cout << "sig2b2: "<< sig2b2 << endl;
      beta2 = mub2 + randn()*sqrt(sig2b2);
      // beta2 = mub2;
    }
    
    // cout <<"sum(beta1mean): " << sum(beta1mean) << endl;
    // cout << "sum(beta1sig): " << sum(beta1sig) << endl;
    // cout <<"sum(beta2mean): " << sum(beta2mean) << endl;
    // cout << "sum(beta2sig): " << sum(beta2sig) << endl;

    // ------------------ //
    // update sgga2
    // ------------------ //
    double tagm, tbgm, taal, tbal;
    tagm = agm + p / 2;
    tbgm = sum(Mu2) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm));
    // sgga2 = tbgm;
    // ------------------ //
    // update sgal2
    // ------------------ //
    taal = aal + p / 2;
    tbal = sum(MuA2) / 2 + bal;
    sgal2 =  1 / randg<double>(distr_param(taal, 1/tbal));
    // sgal2 = tbal;
    
    // ------------------ //
    // update omega
    // ------------------ //
    double wpa1, wpa2;
    wpa1 = a + sum(Eta);
    wpa2 = b + Eta.n_elem - sum(Eta);
    ome = R::rbeta(wpa1, wpa2);
    
    // cout << "iter:" << iter << endl;
    // cout <<"beta1:" << beta1<< endl;
    // cout << "beta2:" << beta2 << endl;
    // cout << "sgga2:" <<sgga2 << endl;
    // cout << "sgal2:" << sgal2<< endl;
    // cout <<"wpa2:" << wpa2 << endl;
    // 
    // cout << "--------------------------------------" << endl;
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        Sgga2Res[lt] = sgga2;
        Sgal2Res[lt] = sgal2;
        Beta1res[lt] = beta1;
        Beta2res[lt] = beta2;
        EtaAll.col(lt) = Eta;
        lt += 1;
      }
    }
    
  }
  
  arma::colvec Etares = mean(arma::conv_to<arma::mat>::from(EtaAll), 1);
  
  
  ObjRMR obj;
  obj.Etares = Etares;
  obj.Beta1res = Beta1res;
  obj.Beta2res = Beta2res;
  obj.Sgga2Res = Sgga2Res;
  obj.Sgal2Res = Sgal2Res;
  
  
  return obj;
}


// [[Rcpp::export]]
Rcpp::List RMRICPSim(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2,
                  arma::mat &R, arma::umat block_inf, SEXP opts = R_NilValue){
  
  uword nblocks = block_inf.n_rows;
  field<mat> F4Rblock(nblocks, 1);
  field<vec> F4gammah(nblocks, 1);
  field<vec> F4Gammah(nblocks, 1);
  field<vec> F4se1(nblocks, 1);
  field<vec> F4se2(nblocks, 1);
  
  for(int i=0; i< (int)(nblocks); i++){
    F4Rblock(i, 0) = R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1));
    F4gammah(i, 0) = gammah.subvec(block_inf(i, 0), block_inf(i, 1));
    F4Gammah(i, 0) = Gammah.subvec(block_inf(i, 0), block_inf(i, 1));
    F4se1(i, 0) = se1.subvec(block_inf(i, 0), block_inf(i, 1));
    F4se2(i, 0) = se2.subvec(block_inf(i, 0), block_inf(i, 1));
  }
  
  
  Options_RMR* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_RMR(opt["agm"], opt["bgm"], opt["aal"], opt["bal"],
                             opt["a"], opt["b"], opt["v"], opt["maxIter"], opt["thin"], opt["burnin"], opt["coreNum"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_RMR();
  }
  
  
  
  ObjRMR obj = RMRLDobj(F4gammah, F4Gammah, F4se1, F4se2, F4Rblock,
                        lp_opt);
  
  vec Beta1res = obj.Beta1res;
  double bhat = mean(Beta1res);
  double se = stddev(Beta1res);
  double pvalue = 2*(R::pnorm(abs(bhat / se), 0, 1, 0, 0));
  
  
  List output = List::create(
    Rcpp::Named("beta.hat") = bhat,
    Rcpp::Named("beta.se") = se,
    Rcpp::Named("beta.p.value") = pvalue,
    Rcpp::Named("Etares") = Rcpp::wrap(obj.Etares),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta2res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Sgal2Res") = Rcpp::wrap(obj.Sgal2Res)
  );
  
  return output;
}


// [[Rcpp::export]]
Rcpp::List RMRICP(arma::field<vec> F4gammah, arma::field<vec> F4Gammah,
                   arma::field<vec> F4se1, arma::field<vec> F4se2,
                   arma::field<mat> F4Rblock,  SEXP opts = R_NilValue){
  
  
  Options_RMR* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_RMR(opt["agm"], opt["bgm"], opt["aal"], opt["bal"],
                             opt["a"], opt["b"], opt["v"], opt["maxIter"], opt["thin"], opt["burnin"], opt["coreNum"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_RMR();
  }
  
  
  
  ObjRMR obj = RMRLDobj(F4gammah, F4Gammah, F4se1, F4se2, F4Rblock,
                        lp_opt);
  
  vec Beta1res = obj.Beta1res;
  double bhat = mean(Beta1res);
  double se = stddev(Beta1res);
  double pvalue = 2*(R::pnorm(abs(bhat / se), 0, 1, 0, 0));
  
  
  List output = List::create(
    Rcpp::Named("beta.hat") = bhat,
    Rcpp::Named("beta.se") = se,
    Rcpp::Named("beta.p.value") = pvalue,
    Rcpp::Named("Etares") = Rcpp::wrap(obj.Etares),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta2res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Sgal2Res") = Rcpp::wrap(obj.Sgal2Res)
  );
  
  return output;
}
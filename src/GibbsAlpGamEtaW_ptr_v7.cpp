#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "GibbsAlpGamEtaW_ptr_v7.hpp"
#include <random>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


void paraBlock_GamAlpEtaW::loop_by_block_gibbs_GamAlpEtaW(int l){
  
  double vl = v[l];
  vec se1l = F4se1(l, 0);
  vec se2l = F4se2(l, 0);
  vec ginvsg2l = F4ginvsg2(l, 0);
  vec GinvsG2l = F4GinvsG2(l, 0);
  vec mul = F4mu(l, 0);
  vec muAl = F4muA(l, 0);
  mat Rl = F4Rblock(l, 0);
  vec Rdiagl = diagvec(Rl);
  vec DinsGRinsGl = F4DinsGRinsG(l, 0);
  
  mat insgRinsgl = F4insgRinsg(l, 0);
  mat insGRinsGl = F4insGRinsG(l, 0);
  mat Rinsgl = F4Rinsg(l, 0);
  mat RinsGl = F4RinsG(l, 0);
  
  vec Rinsgmul = F4Rinsgmu(l, 0);
  vec RinsGmul = F4RinsGmu(l, 0);
  vec RinsGmuAl = F4RinsGmuA(l, 0);
  vec Wl = F4W(l, 0);
  vec invsg2l = F4invsg2(l, 0);
  vec invsG2l = F4invsG2(l, 0);
  
  vec Winvsgal2l = Wl * invsgal2;
  int pl = Wl.n_elem;
  double b12 = beta1 * beta1;
  double b22 = beta2 * beta2;
  double wa, wb;
  // --------------------------------------------- #
  // update gamma
  // --------------------------------------------- #
  vec v20l= 1. / (invsg2l + b12*DinsGRinsGl + invsgga2);
  vec v21l = 1. / (invsg2l + b22*DinsGRinsGl+ invsgga2);
  
  if(Eta[l]==0){
    for(int j = 0; j < pl; j++){
      vec tmp1, tmp2;
      double RinSmujj1, RinSmujj2, mu0;
      tmp1 = Rinsgmul - Rinsgl.col(j)*mul[j];
      tmp2 = RinsGmul - RinsGl.col(j)*mul[j];
      RinSmujj1 = Rinsgmul[j] - Rdiagl[j]*mul[j] / se1l[j];
      RinSmujj2 = RinsGmul[j] - Rdiagl[j]*mul[j] / se2l[j];
      
      mu0 = (ginvsg2l[j] + beta1*GinvsG2l[j] - RinSmujj1 / se1l[j] - b12/se2l[j] * RinSmujj2)*v20l[j];
      mul[j] = mu0 + randn()*sqrt(v20l[j]);
      // mul[j] = mu0;
      Rinsgmul = tmp1 + Rinsgl.col(j)*mul[j];
      RinsGmul = tmp2 + RinsGl.col(j)*mul[j];
    }
  }else{
    for(int j = 0; j < pl; j++){
      vec tmp1, tmp2;
      double RinSmujj1, RinSmujj2, mu1;
      tmp1 = Rinsgmul - Rinsgl.col(j)*mul[j];
      tmp2 = RinsGmul - RinsGl.col(j)*mul[j];
      RinSmujj1 = Rinsgmul[j] - Rdiagl[j]*mul[j] / se1l[j];
      RinSmujj2 = RinsGmul[j] - Rdiagl[j]*mul[j] / se2l[j];
      
      mu1 = (ginvsg2l[j] + beta2*GinvsG2l[j] - RinSmujj1 / se1l[j] - b22/se2l[j] * RinSmujj2 - beta2/se2l[j]*RinsGmuAl[j])*v21l[j];
      
      mul[j] = mu1 + randn()*sqrt(v20l[j]);
      // mul[j] = mu1;
      Rinsgmul = tmp1 + Rinsgl.col(j)*mul[j];
      RinsGmul = tmp2 + RinsGl.col(j)*mul[j];
      
    }
  }
  
  // --------------------------------------------- #
  // update alpha
  // --------------------------------------------- #
  
  vec v2Al = 1. / (invsG2l + Winvsgal2l);
  
  
  if(Eta[l]==0){
    for(int k = 0; k < pl; k++){
      muAl[k] = randn()*sqrt(1. / Winvsgal2l[k]);
      // muAl[k] = 0.01;
    }
    RinsGmuAl = Rl * diagmat(1. /se2l) * muAl; // remember this update !!
    
  }else{
    for(int k = 0; k < pl; k++){
      vec tmp3;
      double RinSmuAkk, muA1;
      tmp3 = RinsGmuAl - RinsGl.col(k)*muAl[k];
      RinSmuAkk = RinsGmuAl[k] - Rdiagl[k]*muAl[k]/se2l[k];
      muA1 = (GinvsG2l[k] - beta2/se2l[k]*RinsGmul[k] - 1/se2l[k]*RinSmuAkk)*v2Al[k];
      muAl[k] = muA1 + randn()*sqrt(v2Al[k]);
      // muAl[k] = muA1;
      RinsGmuAl = tmp3 + RinsGl.col(k)*muAl[k];
    }
    
  }
  
  // --------------------------------------------- #
  // update W set Wlk = 1(normal distribution), if pl < 10
  // --------------------------------------------- #
  if(pl<10){
    Wl = ones(pl, 1);
  }else{
    for(int k = 0; k < pl; k++){
      wa = (vl + 1) / 2;
      wb = (vl + muAl[k] * muAl[k] * invsgal2) / 2;
      Wl[k] = randg<double>(distr_param(wa, 1. / wb)); // Wl[k] = 1;
    }
  }

  // --------------------------------------------- #
  // update Eta
  // --------------------------------------------- #
  int eta;
  mat D0l, D1l, b1l, invD1l, invD0l;
  vec b0l;
  double lik1, lik0, prob0, prob;
  vec muEta1, muEta0, muAEta1;
  
  D0l = invsgga2*diagmat(ones(pl, 1)) + insgRinsgl + b12*insGRinsGl;
  b0l = ginvsg2l + beta1 * GinvsG2l;
  
  D1l = join_cols(join_rows(invsgga2*diagmat(ones(pl, 1)) + insgRinsgl + b22*insGRinsGl, beta2*insGRinsGl),
                  join_rows(beta2*insGRinsGl, diagmat(Winvsgal2l) + insGRinsGl));
  b1l = join_cols(ginvsg2l + beta2*GinvsG2l, GinvsG2l);
  
  invD1l = inv(D1l);
  invD0l = inv(D0l);
  // lik1 = 0.5*as_scalar(b1l.t()*invD1l*b1l) - 0.5*log(det(D1l));
  lik1 = 0.5*as_scalar(b1l.t()*invD1l*b1l) - sum(log(diagvec(chol(D1l))));
  // lik0 = 0.5*as_scalar(b0l.t()*invD0l*b0l) - 0.5*log(det(D0l)) - 0.5*log(det(diagmat(Winvsgal2l)));
  // lik0 = 0.5*as_scalar(b0l.t()*invD0l*b0l) - 0.5*log(det(D0l)) - 0.5*log(prod(Winvsgal2l));
  lik0 = 0.5*as_scalar(b0l.t()*invD0l*b0l) - sum(log(diagvec(chol(D0l)))) - 0.5*log(prod(Winvsgal2l));
  
  
  // if(l==2){
  //   // cout<< D1l <<endl;
  //   cout <<"det(D1l): "<< 0.5*log(det(D1l))<< " sum(log(diagvec(chol(D1l)))):" << sum(log(diagvec(chol(D1l)))) <<"\n" << endl;
  //   cout << "diff1" << 0.5*log(det(D1l)) - sum(log(diagvec(chol(D1l)))) << endl;
  //   cout << "diff0" << 0.5*log(det(D0l)) - sum(log(diagvec(chol(D0l)))) << endl;
  //   cout << "---------------------"<<"\n"<<endl;
  // }
  prob0 = lik1 - lik0 + logw;
  prob = 1. / (1 + exp(-prob0));
  eta = R::rbinom(1, prob);
  Eta[l] = eta;
  // cout << "l:" << l<< log(prod(Winvsgal2l))  << "::" << log(det(diagmat(Winvsgal2l))) << "diff:" << log(prod(Winvsgal2l))- log(det(diagmat(Winvsgal2l))) <<endl;
  // --------------------------------------------- #
  // update terms in each blocks
  // --------------------------------------------- #
  F4mu(l, 0) = mul;
  F4muA(l, 0) = muAl;
  F4Rinsgmu(l, 0) = Rinsgmul;
  F4RinsGmu(l, 0) = RinsGmul;
  F4RinsGmuA(l, 0) = RinsGmuAl;
  F4W(l, 0) = Wl;
  // F4Winvsgal2(l, 0) =   Wl / sgal2;
  
  // --------------------------------------------- #
  // update terms for updating Beta
  // --------------------------------------------- #
  muEta1 = mul * eta;
  muEta0 = mul * (1 - eta);
  muAEta1 = muAl * eta;
  
  beta1sig[l] = as_scalar(muEta0.t()*insGRinsGl*muEta0);
  beta1mean[l] = as_scalar(GinvsG2l.t()*muEta0);
  
  beta2sig[l] = as_scalar(muEta1.t()*insGRinsGl*muEta1);
  beta2mean[l] = as_scalar(GinvsG2l.t()*muEta1 - muEta1.t()*insGRinsGl*muAEta1);

  // -----------------------------------------------------------------------
  // for sgga2 and sgal2 iteration
  Mu2[l] = sum(mul % mul);
  MuA2[l] = sum(muAl % muAl);
  // // -----------------------------------------------------
  se1l.reset();
  se2l.reset();
  Rinsgl.reset();
  RinsGl.reset();
  insGRinsGl.reset();
  insgRinsgl.reset();
  Rl.reset();
  ginvsg2l.reset();
  GinvsG2l.reset();
  
}

std::mutex _mtx0;
int paraBlock_GamAlpEtaW::next_GamAlpEtaW(){
  std::lock_guard<std::mutex> lockGuard(_mtx0);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_GamAlpEtaW::update_by_thread_GamAlpEtaW(int thread_id){
  while(true){
    int idx = next_GamAlpEtaW();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_GamAlpEtaW(idx);
  }
}
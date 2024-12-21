#ifndef GibbsAlpGamEtaW_ptr_v7_hpp
#define GibbsAlpGamEtaW_ptr_v7_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

class Options_RMRindep{
public:
  // Constructor definition
  // The complier deciedes which constructor to be called depending on 
  // the number of argument present with the object
  Options_RMRindep(){
    this -> agm = 0;
    this -> bgm = 0;
    this -> aal = 0;
    this -> bal = 0;
    this -> a = 1;
    this -> b = 1;
    this -> v = 100;
    this -> maxIter = 3000;
    this -> thin = 10;
    this -> burnin = 2000;
  }
  
  Options_RMRindep(double agm, double bgm, double aal, double bal, 
              double a, double b, double v, uword maxIter, uword thin, uword burnin){
    
    this -> agm = agm;
    this -> bgm = bgm;
    this -> aal = aal;
    this -> bal = bal;
    this -> a = a;
    this -> b = b;
    this -> v = v;
    this -> maxIter = maxIter;
    this -> thin = thin;
    this -> burnin = burnin;
  }
  
  double agm;
  double bgm;
  double aal;
  double bal;
  double a;
  double b;
  double v;
  uword maxIter;
  uword thin;
  uword burnin;
  
};
struct ObjRMRindep{
  vec Beta1res;
  vec Beta2res;
  vec Sgga2Res;
  vec Sgal2Res;
  vec Etares;
};


class Options_RMR{
public:
  // Constructor definition
  // The complier deciedes which constructor to be called depending on 
  // the number of argument present with the object
  Options_RMR(uword nblocks){
    this -> agm = 0;
    this -> bgm = 0;
    this -> aal = 0;
    this -> bal = 0;
    this -> a = 1;
    this -> b = 1;
    this -> v = 100*ones(nblocks, 1);
    this -> maxIter = 3000;
    this -> thin = 10;
    this -> burnin = 2000;
    this -> coreNum = 20;
  }
  
  Options_RMR(double agm, double bgm, double aal, double bal, 
                 double a, double b, vec v, uword maxIter, uword thin, uword burnin, uword coreNum){
    
    this -> agm = agm;
    this -> bgm = bgm;
    this -> aal = aal;
    this -> bal = bal;
    this -> a = a;
    this -> b = b;
    this -> v = v;
    this -> maxIter = maxIter;
    this -> thin = thin;
    this -> burnin = burnin;
    this -> coreNum = coreNum;
    
  }
  
  double agm;
  double bgm;
  double aal;
  double bal;
  double a;
  double b;
  vec v;
  uword maxIter;
  uword thin;
  uword burnin;
  uword coreNum;
  
};
struct ObjRMR{
  vec Beta1res;
  vec Beta2res;
  vec Sgga2Res;
  vec Sgal2Res;
  vec Etares;
};



class paraBlock_GamAlpEtaW{
public:
  int current_idx=0;
  // uword Ngene_active;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  
  int conspar;
  ivec Eta;
  vec beta1sig, beta2sig, beta1mean, beta2mean, Mu2, MuA2, v;
  double beta1, beta2, logw, invsgga2, invsgal2;
  
  field<vec> F4se1, F4se2, F4ginvsg2, F4GinvsG2, F4DinsGRinsG;
  field<vec> F4mu, F4muA, F4Rinsgmu,  F4RinsGmu, F4RinsGmuA;
  field<mat> F4Rins, F4Rins2, F4Rblock, F4insGRinsG, F4insgRinsg, F4Rinsg, F4RinsG;
  field<vec> F4W, F4invsg2, F4invsG2;
  
  // field<mat> F4sgRinvsg, F4sGRinvsG, F4sgRsg, F4sGRsG;
  
  
  paraBlock_GamAlpEtaW(uword &nblocks, field<vec> &F4se1, field<vec> &F4se2, field<vec> &F4ginvsg2,
                       field<vec> &F4GinvsG2, field<vec> &F4mu, field<vec> &F4muA,field<mat> &F4Rblock,
                       field<vec> &F4DinsGRinsG,field<mat> &F4insgRinsg, field<mat> &F4insGRinsG,
                       field<mat> &F4Rinsg, field<mat> &F4RinsG,
                       field<vec> &F4Rinsgmu, field<vec> &F4RinsGmu, field<vec> &F4RinsGmuA,
                       field<vec> &F4W, field<vec> &F4invsg2, field<vec> &F4invsG2,
                       arma::vec &beta1sig, arma::vec &beta2sig, arma::vec &beta1mean, arma::vec &beta2mean,
                       arma::ivec Eta, arma::vec &Mu2, arma::vec &MuA2,
                       double &beta1, double &beta2, double &logw, double &invsgga2, double &invsgal2, vec &v){ 
    this -> nblocks = nblocks;
    this -> F4se1 = F4se1;
    this -> F4se2 = F4se2;
    this -> F4ginvsg2 = F4ginvsg2;
    this -> F4GinvsG2 = F4GinvsG2;
    this -> F4mu = F4mu;
    this -> F4muA = F4muA;
    this -> F4Rblock = F4Rblock;
    this -> F4DinsGRinsG = F4DinsGRinsG;
    this -> F4insgRinsg = F4insgRinsg;
    this -> F4insGRinsG = F4insGRinsG;
    this -> F4RinsGmu = F4RinsGmu;
    this -> F4RinsGmuA = F4RinsGmuA;
    this -> F4Rinsg = F4Rinsg;
    this -> F4RinsG = F4RinsG;
    this -> F4Rinsgmu = F4Rinsgmu;
    this -> F4RinsGmu = F4RinsGmu;
    this -> F4RinsGmuA = F4RinsGmuA;
    this -> F4W = F4W;
    this -> F4invsg2 = F4invsg2;
    this -> F4invsG2 = F4invsG2;
    this -> beta1sig = beta1sig;
    this -> beta1mean = beta1mean;
    this -> beta2sig = beta2sig;
    this -> beta2mean = beta2mean;
    this -> Mu2 = Mu2;
    this -> MuA2 = MuA2;
    this -> Eta = Eta;
    this -> beta1 = beta1;
    this -> beta2 = beta2;
    this -> logw = logw;
    this -> invsgga2 = invsgga2;
    this -> invsgal2 = invsgal2;
    this -> v = v;

  }
  
  
  int  next_GamAlpEtaW();
  void loop_by_block_gibbs_GamAlpEtaW(int l);
  void update_by_thread_GamAlpEtaW(int thread_id);
  
};


#endif /* GibbsAlpGamEtaW_ptr_v7_hpp */

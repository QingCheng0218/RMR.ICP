#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::ivec std_setdiff(arma::ivec& x, arma::ivec& y) {

  // std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  // std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;

  std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
                      std::inserter(out, out.end()));

  // return arma::conv_to<arma::ivec>::from(out);
  return out;
}

// [[Rcpp::export]]
ivec LDclump(arma::mat &R, double ld_r2_thresh){
  int nsnp = R.n_rows;

  R.diag() = zeros(nsnp);
  R = trimatl(R);

  ivec a = regspace<ivec>(0, 1, nsnp - 1);
  ivec idx1 = zeros<ivec>(nsnp*nsnp, 1);
  ivec idx2 = zeros<ivec>(nsnp*nsnp, 1);
  for(int i=0; i< nsnp; i++){
    uvec ii = regspace<uvec>(i*nsnp, 1, ((i+1)*nsnp - 1));
    idx1.elem(ii) = a;
    idx2.elem(ii) = i*ones<ivec>(nsnp, 1);
  }

  uvec q1 = find(R%R > ld_r2_thresh);
  ivec id1, id2;

  ivec IDsave;

  if(q1.n_elem){
    id1 = idx2.elem(q1);
    id2 = idx1.elem(q1);

    idx1.reset();
    idx2.reset();

    ivec slct_indx = join_cols(id1, id2);
    ivec indxbuf = arma::unique(slct_indx);

    // count how many SNPs in high LD
    int nproc = indxbuf.n_elem;
    ivec n_slct_snp = zeros<ivec>(nproc, 1);

    for( int i = 0; i < nproc; i = i + 1 ){
      n_slct_snp[i] = sum(slct_indx==indxbuf[i]);
    }

    // decide the index to remove
    nproc = id1.n_elem;
    int n1, n2;
    for(int i = 0; i < nproc; i = i + 1){
      n1 = n_slct_snp[as_scalar(find(indxbuf==id1[i]))];
      n2 = n_slct_snp[as_scalar(find(indxbuf==id2[i]))];
      if(n1 < n2){
        int t;
        t = id1[i];
        id1[i] = id2[i];
        id2[i] = t;
      }

    }

    id1 = arma::unique(id1);
    // return the index to save(C indicator start from 0)
    IDsave = std_setdiff(a, id1);
  }else{
    IDsave = a;

  }

  return IDsave;

}

// [[Rcpp::export]]
List blockinffun(arma::field<vec> F4gammah){
  int L = F4gammah.size();
  ivec blinf = zeros<ivec>(L, 1);
  imat block_inf = zeros<imat>(L, 2);
  for(int j = 0; j < L ; j++){
    blinf[j] = F4gammah[j].n_elem;
  }
  ivec blsum = cumsum(blinf);
  int totalp = sum(blinf);
  ivec col0 = zeros<ivec>(L, 1);
  col0[0] = 0;
  col0.subvec(1, L - 1) = blsum.subvec(0, L - 2);
  block_inf.col(0) = col0;
  block_inf.col(1) = blsum - 1;

  List result = List::create(Named("blinf") = blinf,
                             Named("block_inf") = block_inf,
                             Named("totalp") = totalp);

  return result;
}


// [[Rcpp::export]]
List LDclumpfun(arma::field<vec> F4gammah, arma::field<vec> F4Gammah,
                 arma::field<vec> F4se1, arma::field<vec> F4se2,
                 arma::field<mat> F4Rblock, double ld_r2_thresh) {

  vec gammahind, Gammahind, se1ind, se2ind;
  List resultFromblockfun = blockinffun(F4gammah);
  ivec blinf = resultFromblockfun["blinf"];
  imat block_inf = resultFromblockfun["block_inf"];
  int totalp = resultFromblockfun["totalp"];
  mat R = zeros(totalp, totalp);
  int L = F4gammah.size();
  ivec id4ld;
  for(int l = 0; l < L; l++){
      vec F4gammah_i = F4gammah[l];
      vec F4Gammah_i = F4Gammah[l];
      vec F4se1_i = F4se1[l];
      vec F4se2_i = F4se2[l];
      mat ldrho = F4Rblock[l];
      R.submat(block_inf(l, 0), block_inf(l, 0), block_inf(l, 1), block_inf(l, 1)) = ldrho;

      ivec index = LDclump(ldrho, ld_r2_thresh);
      ivec blockid = regspace<ivec>(block_inf(l, 0), block_inf(l, 1));
      ivec id4ldtmp = zeros<ivec>(index.size(), 1);
      vec gammahtmp = zeros(index.size(), 1);
      vec Gammahtmp = zeros(index.size(), 1);
      vec se1tmp = zeros(index.size(), 1);
      vec se2tmp = zeros(index.size(), 1);
      for(int i = 0; i < index.size(); ++i) {
        id4ldtmp[i] = blockid[index[i]];
        gammahtmp[i] = F4gammah_i[index[i]];
        Gammahtmp[i] = F4Gammah_i[index[i]];
        se1tmp[i] = F4se1_i[index[i]];
        se2tmp[i] = F4se2_i[index[i]];
      }
      id4ld = join_cols(id4ld, id4ldtmp);
      gammahind = join_cols(gammahind, gammahtmp);
      Gammahind = join_cols(Gammahind, Gammahtmp);
      se1ind = join_cols(se1ind, se1tmp);
      se2ind = join_cols(se2ind, se2tmp);
  }
  
  List result = List::create(
    Named("block_inf") = block_inf,
    Named("id4ld") = id4ld + 1, // the index start from 1 (R environment).
    Named("R") = R,
    Named("gammahind") = gammahind,
    Named("Gammahind") = Gammahind,
    Named("se1ind") = se1ind,
    Named("se2ind") = se2ind);

  return result;
}





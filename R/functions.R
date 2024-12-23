genRawGeno <- function(maf, L, M, rho, n){
  SIGMA = matrix(nrow=M,ncol=M)
  for (i in 1:M){
    for (j in 1:M){
      SIGMA[i,j] = rho^(abs(i-j));
    }
  }
  
  nsnp = L*M;
  X = NULL;
  for ( l in 1:L ){
    
    index = (M*(l-1)+1): (M*l);
    AAprob = maf[index]^2.;
    Aaprob = 2*maf[index]*(1-maf[index]);
    quanti = matrix(c(1-Aaprob-AAprob, 1- AAprob),M,2);
    Xt = rmvnorm(n, mean=rep(0,M), sigma=SIGMA, method="chol")
    Xt2 = matrix(0,n,M);
    for (j in 1:M){
      cutoff = qnorm(quanti[j,]);
      Xt2[Xt[,j] < cutoff[1],j] = 0;
      Xt2[Xt[,j] >= cutoff[1] & Xt[,j] < cutoff[2],j] = 1;  ## attention
      Xt2[Xt[,j] >= cutoff[2],j] = 2;
    }
    X <- cbind(X,Xt2);
  }
  return(X)
}

genSumStat <- function(x12, n1, n2, M, L, b1, Alrate, sigma2g, dfA, delta, h2a, h2g){
  
  p <- M*L;
  block_inf <- cbind(seq(1, p, M), seq(M, p, M));
  block_inf1 <- block_inf - 1;
  # ------------------------------------------------------------------------
  # The correlated horizontal pleiotropy (alpha).
  # First:generate alpha titlde (came from t-dsitribution) and delta to obtain alpha.
  gamma = rnorm(p)*sqrt(sigma2g);
  alphaT = rt(p, df = dfA);
  alpha = delta*gamma + alphaT;
  # ------------------------------------------------------------------------
  q = 50
  u = matrix(rnorm( (n1+n2) * q),ncol=q);
  Su = matrix(c(1,0.8,0.8,1),nrow=2)
  bu = rmvnorm(q,mean=rep(0,2), sigma = Su,method="chol")
  by = bu[,1]; bz = bu[,2];
  uby = u%*%by; ubz = u%*%bz;
  uby = uby/sqrt(as.numeric(var(uby)/0.6));
  ubz = ubz/sqrt(as.numeric(var(ubz)/0.2));
  
  x12g = x12%*%gamma;
  
  if(b1!=0){
    h2ga = (h2g *( 1 + b1^2))/(b1^2 * (1 - h2g));
    gamma0 = gamma/sqrt(as.numeric(var(x12g)/h2ga));
    x12g = x12%*%gamma0;
  }
  
  yall = x12g + uby + rnorm(n1+n2)*as.numeric(sqrt(1-var(uby)));
  # ------------------------------------------------------------------------
  # The direct effects on Z
  h2yb = var(b1*yall);
  h2al = (h2a + h2a*h2yb)/(1 - h2a);
  
  if(h2a==0){
    alpha = rep(0, p);
    x12a = x12%*%alpha;
  }else{
    if(Alrate!=1){
      alno = floor(L*Alrate);
      ind = sample(1:L,alno);
      if(length(ind)==1){
        indxAL = block_inf[ind, 1]:block_inf[ind, 2]
        alpha[-indxAL] = 0;
      }else{
        indxAL = NULL;
        for(i in 1:length(ind)){
          tmp = block_inf[ind[i], 1]:block_inf[ind[i], 2]
          indxAL = append(indxAL, tmp);
        }
        alpha[-indxAL] = 0;
        x12a = x12%*%alpha;
        alpha0 = alpha/sqrt(as.numeric(var(x12a)/(h2al)));
        x12a = x12%*%alpha0;
      }
    }
  }
  
  # ------------------------------------------------------------------------
  resz = ubz + rnorm(n1 + n2)*as.numeric(sqrt(1-var(ubz)));
  zall = b1*yall  + x12a +  resz;
  var(x12a)/var(zall);
  var(b1*x12g)/var(zall);
  
  
  y = yall[1:n1];
  z = zall[(n1+1):(n1+n2)];
  
  x1 = x12[1:n1, ];
  x2 = x12[(n1+1):(n1+n2), ]
  
  
  # create summary statistics
  gammaSS = fastSigLm(y, x1);
  GammaSS = fastSigLm(z, x2);
  gammah = gammaSS$coef;
  se1 = gammaSS$std;
  Gammah = GammaSS$coef;
  se2 = GammaSS$std;
  
  return(list(gammah = gammah, se1 = se1, se2 = se2, 
              Gammah = Gammah, CHPindex = sort(ind)))
  
}


traceplot <- function(bhatpoint){
  y <- bhatpoint
  x <- 1:length(bhatpoint);
  da <- cbind(x, y);
  dat <- data.frame(da);
  p1 <- ggplot(data = dat, aes(x= x, y = y))+  geom_line()  +
    labs( x = paste0("GibbsSampleIndex"), y =  expression(hat(beta[1])));
  p1 = p1 + theme(axis.title.x = element_text(size=10,face = "bold"),
                  axis.text.x = element_text(size=12,face = "bold"),
                  axis.title.y = element_text(size=10,face = "bold"),
                  axis.text.y = element_text(size=12,face = "bold"));
  return(p1);
}

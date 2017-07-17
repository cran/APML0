

###############################
#####  Linear Regression  #####
###############################


######################
###  Enet (L1+L2)  ###
######################

EnetLm=function(x, y, alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, wbeta=rep(1,ncol(x)), isd=FALSE, keep.beta=FALSE, thresh=1e-6, maxit=1e+5) {
  # alpha=alphai; rlambda=NULL; inzero=TRUE; isd=FALSE; keep.beta=FALSE; thresh=1e-5; maxit=1e+5; lambda=lambdai

  penalty=ifelse(alpha==1, "Lasso", "Enet")

  N0=nrow(x); p=ncol(x)

  ### Adaptive
  adaptive=ifelse(any(wbeta!=1), TRUE, FALSE)

  aPen=ifelse(all(wbeta>0), TRUE, FALSE)

  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }

  #####  Run  #####
  out=EnetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, as.integer(isd), p, N0, thresh, maxit, 1e-5)
  nlambdai=out$nlambda ## number of lambdas
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]

  out$Beta=Matrix(out$Beta[, 1:nlambdai], sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]

  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
  } else {

    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))


    #####  Cross-validation estimates  #####
    # ypred=matrix(0,nrow=N0,ncol=nlambdai)

    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      outi[[i]]=cvEnetLmC(x[!temid, ], y[!temid], alpha, lambdai, nlambdai, wbeta, N0i[i],p, thresh, maxit, x[temid, ], y[temid], Nf[i], 1e-5)

      # ypred[temid,]=outi[[i]]$pred
      cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda] ## for ith fold
    }


    cvRSS=matrix(cvRSS[, 1:nlambdai], ncol=nlambdai)
    cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))


    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)

    if (!inzero) {
      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }

    #####  Cross-validation hard threshold  #####
    il0=indexi; cvm=list(); cv.min=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=sapply(outi, function(x){x$Beta[, il0]})

      Betao=apply(Betai!=0, 2, sum)
      numi2=min(max(Betao), numi)

      a0i=sapply(outi, function(x){x$a0S[il0]})

      if (numi2>0) {
        cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; temid=foldid==i
          numj=min(Betao[i], numi)
          if (numj==0) {
            cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
          } else {
            Betajj=Betaj
            Betajj[wbeta==0]=max(abs(Betaj))+1
            temo=rank(-abs(Betajj), ties.method="min")

            temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
            temo=temo[order(temo[, 1]), ]
            cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i], a0i[i])

          }
        }
      } else {
        cvRSS=matrix(NA, nrow=nfolds, ncol=1)
        for (i in 1:nfolds) {
          temid=foldid==i
          cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
        }
      }

      cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      temi=cvm[[il0]]

      if (aPen) {
        cv.min[il0]=min(temi)
      } else {
        cv.min[il0]=ifelse(length(temi)>sum(wbeta==0),min(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
      }
      # cv.min[il0]=ifelse(length(temi)>sum(wbeta==0),min(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
      # cv.min[il0]=min(cvm[[il0]])

      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.min[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=sapply(outi, function(x){x$Beta[, il1[j]]})
            Betao=apply(Betai!=0, 2, sum)
            numi2=min(max(Betao), numi)

            if (numi2>0) {
              cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
              for (i in 1:nfolds ){
                Betaj=Betai[, i]; temid=foldid==i
                numj=min(Betao[i], numi)
                if (numj==0) {
                  cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
                } else {
                  Betajj=Betaj
                  Betajj[wbeta==0]=max(abs(Betaj))+1
                  temo=rank(-abs(Betajj), ties.method="min")

                  temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
                  temo=temo[order(temo[, 1]), ]
                  cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i], a0i[i])
                }
              }
            } else {
              cvRSS=matrix(NA, nrow=nfolds, ncol=1)
              for(i in 1:nfolds) {
                temid=foldid==i
                cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
              }
            }
            cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum)
            rm(cvRSS)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)

            temi=cvm[[il1[j]]]
            if (aPen) {
              cv.min[il1[j]]=min(temi)
            } else {
              cv.min[il1[j]]=ifelse(length(temi)>sum(wbeta==0),min(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
            }
            # cv.min[il1[j]]=ifelse(length(temi)>sum(wbeta==0),min(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
            # cv.min[il1[j]]=min(cvm[[il1[j]]])

          }
        } else {
          break
        }
      }
      if (il1[j]==1 | il1[j]==nlambdai)
        break
      if (il0==which.min(cv.min)) {
        break
      } else {
        il0=which.min(cv.min)
      }
    }
    index0=which.min(cv.min)

    Beta0=out$Beta[,index0]

    temi=cvm[[index0]]
    if (aPen) {
      cuti=which.min(temi)
    } else {
      cuti=ifelse(length(temi)>sum(wbeta==0),which.min(temi[-c(1:sum(wbeta==0))])+sum(wbeta==0),sum(wbeta==0))
    }
    # cuti=which.min(cvm[[index0]])

    Beta0j=Beta0
    Beta0j[which(wbeta==0)]=max(abs(Beta0j))+1

    Beta0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0

    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.min[index0],nzero=cuti)

    if (!keep.beta) {
      # cv.nzero=cvm[[index0]]
      return(list(Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
  }
}



############################
###  Net (L1+Laplacian)  ###
############################

NetLm=function(x, y, Omega=NULL, alpha=1.0, lambda=NULL, nlambda=50, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), isd=FALSE, keep.beta=FALSE, thresh=1e-6, maxit=1e+5){
  # Omega=Omegai;alpha=alphai; rlambda=NULL; inzero=TRUE; isd=FALSE; keep.beta=FALSE; thresh=1e-5; maxit=1e+5; sgn=rep(1,ncol(x)); nlambda=20; lambda=lambdai;

  penalty=ifelse(alpha==1, "Lasso", "Net")

  N0=nrow(x);p=ncol(x)

  ### Adaptive
  adaptive=c(ifelse(any(wbeta!=1), TRUE, FALSE),ifelse(any(sgn!=1), TRUE, FALSE))

  aPen=ifelse(all(wbeta>0), TRUE, FALSE)

  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }

  ## Check Omega: positive off-diagonal, zero diagonal
  if (any(diag(Omega)!=0)) {
    diag(Omega)=0
    # cat("Diagonal of Omega was set to all zeros\n")
  }
  if (any(Omega<0)) {
    Omega=abs(Omega)
    # cat("Off-diagonal of Omega was foced to non-negative values\n")
  }

  ### Correlation/Adjacency matrix
  if (inherits(Omega, "dgCMatrix")) {
    W=OmegaSC(Omega, sgn); W$loc=W$loc+1
  } else {
    W=OmegaC(Omega, sgn); W$loc=W$loc+1
  }
  rm(Omega)

  #####  Run  #####
  out=NetLmC(x, y, alpha, lambda, nlambda, ilambda, wbeta, W$Omega, W$loc, W$nadj, as.integer(isd), p, N0, thresh, maxit, 1e-5)
  nlambdai=out$nlambda
  if (nlambdai==0)
    return(NULL)
  lambdai=out$lambda[1:nlambdai]

  out$Beta=Matrix(out$Beta[, 1:nlambdai], sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]
  out$rsq=out$rsq[1:nlambdai]

  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, rsq=out$rsq, nzero=out$nzero)
    return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
  } else {

    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid)
    N0i=numeric(nfolds); Nf=numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
    }
    weighti=as.vector(tapply(rep(1,N0), foldid, sum))

    ###  Cross-validation estimates  ###
    outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      temid=(foldid==i)
      outi[[i]]=cvNetLmC(x[!temid, ], y[!temid], alpha, lambdai, nlambdai, wbeta, W$Omega, W$loc, W$nadj, N0i[i],  p,  thresh, maxit, x[temid, ], y[temid], Nf[i], 1e-5)
      cvRSS[i, 1:outi[[i]]$nlambda]=outi[[i]]$RSSp[1:outi[[i]]$nlambda]
    }

    cvRSS=matrix(cvRSS[, 1:nlambdai], ncol=nlambdai)
    cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvRSS) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))

    indexi=which.min(cvm)
    indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai, rsq=out$rsq, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi, stringsAsFactors=FALSE)

    if (!inzero) {
      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag))
      }
    }

    #####  Cross-validation hard threshold  #####
    il0=indexi; cvm=list(); cv.min=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=sapply(outi, function(x){x$Beta[, il0]})
      Betao=apply(Betai!=0, 2, sum)
      numi2=min(max(Betao), numi)

      a0i=sapply(outi, function(x){x$a0S[il0]})

      if (numi2>0) {
        cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; temid=foldid==i
          numj=min(Betao[i], numi)
          if (numj==0) {
            cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
          } else {

            Betajj=Betaj
            Betajj[wbeta==0]=max(abs(Betaj))+1
            temo=rank(-abs(Betajj), ties.method="min")

            temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
            temo=temo[order(temo[, 1]), ]
            cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i], a0i[i])
          }
        }
      } else {
        cvRSS=matrix(NA, nrow=nfolds, ncol=1)
        for (i in 1:nfolds) {
          temid=foldid==i
          cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
        }
      }

      cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvRSS) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      temi=cvm[[il0]]
      if (aPen) {
        cv.min[il0]=min(temi)
      } else {
        cv.min[il0]=ifelse(length(temi)>sum(wbeta==0),min(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
      }
      # cv.min[il0]=min(cvm[[il0]])

      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.min[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=sapply(outi, function(x){x$Beta[, il1[j]]})
            Betao=apply(Betai!=0, 2, sum)
            numi2=min(max(Betao), numi)

            if (numi2>0) {
              cvRSS=matrix(NA, nrow=nfolds, ncol=numi2)
              for (i in 1:nfolds) {
                Betaj=Betai[, i]; temid=foldid==i
                numj=min(Betao[i], numi)
                if (numj==0) {
                  cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), numj, numi2, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
                } else {

                  Betajj=Betaj
                  Betajj[wbeta==0]=max(abs(Betaj))+1
                  temo=rank(-abs(Betajj), ties.method="min")

                  temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
                  temo=temo[order(temo[, 1]), ]
                  cvRSS[i, ]=cvTrimLmC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, x[temid,], y[temid], Nf[i], a0i[i])
                }
              }
            } else {
              cvRSS=matrix(NA, nrow=nfolds, ncol=1)
              for(i in 1:nfolds) {
                temid=foldid==i
                cvRSS[i, ]=cvTrimLmC(c(0.0, 0.0), 0, 0, c(0, 0), x[temid,], y[temid], Nf[i], a0i[i])
              }
            }

            cvraw=cvRSS/weighti;nfoldi=apply(!is.na(cvraw), 2, sum)
            rm(cvRSS)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)

            temi=cvm[[il1[j]]]
            if (aPen) {
              cv.min[il1[j]]=min(temi)
            } else {
              cv.min[il1[j]]=ifelse(length(temi)>sum(wbeta==0),min(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
            }
            # cv.min[il1[j]]=min(cvm[[il1[j]]])
          }
        } else {
          break
        }
      }
      if(il1[j]==1 | il1[j]==nlambdai)
        break
      if (il0==which.min(cv.min)) {
        break
      } else {
        il0=which.min(cv.min)
      }
    }
    index0=which.min(cv.min)

    Beta0=out$Beta[,index0]

    temi=cvm[[index0]]
    if (aPen) {
      cuti=which.min(temi)
    } else {
      cuti=ifelse(length(temi)>sum(wbeta==0),which.min(temi[-c(1:sum(wbeta==0))])+sum(wbeta==0),sum(wbeta==0))
    }
    # cuti=which.min(cvm[[index0]])

    Beta0j=Beta0
    Beta0j[which(wbeta==0)]=max(abs(Beta0j))+1

    Beta0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0

    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.min[index0],nzero=cuti)

    if (!keep.beta) {
      # lambda.1se=lambdai[indexij], cv.nzero=cvm[[index0]]
      return(list(Beta=out$Beta[, c(indexij, indexi)], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$Beta, Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
  }
}





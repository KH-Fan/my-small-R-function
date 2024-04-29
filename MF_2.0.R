
rm(list=ls())

"%nin%" <- function(x,y){!x%in%y}

"w/o" <- function(x,y){x[!x%in%y]}

comb_dist <- function(n){  
  if(n>9){
    x <- comb_list(n)
    z <- NULL
    for(i in 1:length(x)){
      y <- letters[1:n]
      y[x[[i]]] <- ""
      z <- c(z,paste(y,collapse=""))
    }
  }else{
    x <- rev(comb_list(n))
    z <- NULL
    for(i in 1:length(x)){
      y <- 1:n
      y[x[[i]]] <- 0
      z <- c(z,as.integer(paste(y,collapse="")))
    }
  }
  return(z)
}

comb_list <- function(n, i) {
  if(missing(i)) i <- n
  if(i==0) combn(n,0,simplify=FALSE)
  else append(combn(n,i,simplify=FALSE), comb_list(n, i-1))
}


cons <- function(fits,convt,digits=1){
  
  est <- coefficients(fits)
  varcov <- vcov(fits)
  
  est.sum <- t(convt) %*% est
  se.sum <- sqrt(t(convt) %*% varcov %*% convt)
  
  est <- sapply(exp(est.sum), formatC, format='f', digit=1)
  ucl <- sapply(exp(est.sum + qnorm(.975) * se.sum), formatC, format='f', digit=digits)
  lcl <- sapply(exp(est.sum - qnorm(.975) * se.sum), formatC, format='f', digit=digits)
  haz <- paste(est," ", "(",lcl,", ", ucl, ")", sep="")
  
  return(haz)
}

stderr <-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

cpct <- function(data.matrix,digits=0,latex=FALSE){
  ds <- ifelse(latex,"//","")
  col.perc <- matrix(formatC(t(data.matrix)/apply(data.matrix,2,sum)*100, format="f", digit=digits),
                     nrow=nrow(data.matrix),ncol=ncol(data.matrix),byrow=TRUE)  
  cpct <- matrix(paste(data.matrix," (",col.perc,ds,"%)",sep=""),
                 nrow=nrow(data.matrix),ncol=ncol(data.matrix),byrow=FALSE)
  rownames(cpct) <- rownames(data.matrix)
  colnames(cpct) <- colnames(data.matrix)
  return(cpct)
}

rpct <- function(data.matrix,digits=0,latex=FALSE){
  ds <- ifelse(latex,"//","")
  row.perc <- matrix(formatC(data.matrix/apply(data.matrix,1,sum)*100, format="f", digit=digits),
                     nrow=nrow(data.matrix),ncol=ncol(data.matrix),byrow=FALSE)  
  rpct <- matrix(paste(data.matrix," (",row.perc,ds,"%)",sep=""),
                 nrow=nrow(data.matrix),ncol=ncol(data.matrix),byrow=FALSE)
  rownames(rpct) <- rownames(data.matrix)
  colnames(rpct) <- colnames(data.matrix)
  return(rpct)
}

adj.p <- function(x,r=3,latex=FALSE){
  ds <- ifelse(latex,"$","")
  um <- ifelse(latex,"^","")
  sapply(x,function(x){
  y <- round(x,r)
  z <- ifelse(y<1/10^r,paste(ds,um,"*<",1/10^r,ds,sep=""),
              ifelse(y<0.05,paste(ds,um,"*",add0(x,r),ds,sep=""),add0(x,r)))
  return(z)
  })
}


add0 <- function(x,r=3){
  sapply(x,function(x){
    if(abs(x)<10^(-r-1)){
      return(paste("0",paste(rep("0",r),collapse=""),sep="."))
    }else if(!is.na(x)){
      y <- as.character(round(x,r))
      n <- nchar(y)
      if(length(grep("\\.",y))==0){
        if(r==0){
          return(y)
        }else{
          return(paste(y,substr(10^(r),2,(r+1)),sep="."))
        }
      }else{
        add.n <- n-which(unlist(strsplit(y,split=""))==".")
        if(add.n<r){
          return(paste(y,substr(10^(r-add.n),2,(r-add.n+1)),sep=""))
        }else{
          return(y)
        }
      }
    }else{
      return(NA)
    }})
}

hocoef <- function(fits,digits=2,pv.digits=3,def.ho="Odds Ratio",latex=FALSE){
  temp <- summary(fits)$coeff
  hor <- add0(exp(temp[-1,1]),r=digits)
  lcl <- add0(exp(temp[-1,1] - qnorm(.975) * temp[-1,2]),r=digits)
  ucl <- add0(exp(temp[-1,1] + qnorm(.975) * temp[-1,2]),r=digits)
  cis <- paste(lcl,ucl,sep=" to ")
  pvs <- adj.p(temp[-1,4],r=pv.digits,latex)
  hocoef <- cbind(hor,cis,pvs)
  ds <- ifelse(latex,"\\","")
  cin <- paste("95",ds,"% CI",sep="")
  colnames(hocoef) <- c(def.ho,cin,"P-value")
  return(hocoef)
}

setwd("~/Desktop/")

save.image("MF_2.0.RData")


####
# These functions are from SPP package https://github.com/hms-dbmi/spp commit 847275a under R/zroutines.R file
# They had been removed from the SPP namespace and we explicitly put them here
####
wtd <- function(x,y,s,e,whs=200,return.peaks=T,min.thr=5,min.dist=200,step=1,direct.count=F,tag.weight=1,bg.x=NULL,bg.y=NULL,bg.weight=1,mask.x=NULL,mask.y=NULL,ignore.masking=F, bg.whs=whs, round.up=F, ...) {
  ignore.masking <- ignore.masking | (is.null(mask.x) & is.null(mask.y));
  if(step>1) {
    x <- floor(x/step+0.5); y <- floor(y/step+0.5)

    if(!is.null(bg.x)) {
      bg.x <- floor(bg.x/step+0.5); bg.y <- floor(bg.y/step+0.5)
    }

    if(!is.null(mask.x)) {
      mask.x <- floor(mask.x/step+0.5); mask.y <- floor(mask.y/step+0.5)
    }


    whs <- floor(whs/step+0.5);
    bg.whs <- floor(bg.whs/step+0.5);
    min.dist <- floor(min.dist/step +0.5);
    s <- floor(s/step+0.5)
    e <- floor(e/step+0.5)
  }

  # scale bg.weight, since within calculation they are considered independent
  bg.weight <- bg.weight*tag.weight;

  rx <- c(s-whs,e+whs);

  # compile tag vectors
  xt <- table(x);
  xh <- integer(diff(rx)+1);
  xh[as.integer(names(xt))-rx[1]+1] <- as.integer(xt);

  yt <- table(y);
  yh <- integer(diff(rx)+1);
  yh[as.integer(names(yt))-rx[1]+1] <- as.integer(yt);

  # compile background vectors
  if(!is.null(bg.x) && length(bg.x)>0) {
    bg.subtract <- 1;

    bg.xt <- table(bg.x);
    bg.xh <- integer(diff(rx)+1);
    bg.xh[as.integer(names(bg.xt))-rx[1]+1] <- as.integer(bg.xt);
    rm(bg.xt);

    bg.yt <- table(bg.y);
    bg.yh <- integer(diff(rx)+1);
    bg.yh[as.integer(names(bg.yt))-rx[1]+1] <- as.integer(bg.yt);
    rm(bg.yt);

    # adjust bg.weight according to bg.whs
    if(bg.whs!=whs) {
      bg.weight <- bg.weight*whs/bg.whs;
    }
  } else {
    bg.subtract <- 0;
    bg.xh <- bg.yh <- c();
  }

  # record masked positions
  if(!ignore.masking) {
    if(!is.null(mask.x) && length(mask.x)>0) {
      mvx <- unique(mask.x); mvx <- setdiff(mvx,as.numeric(names(xt)));
      mvx <- mvx[mvx>=rx[1] & mvx<=rx[2]];
      xh[mvx-rx[1]+1] <- -1;
    }

    if(!is.null(mask.y) && length(mask.y)>0) {
      mvy <- unique(mask.y); mvy <- setdiff(mvy,as.numeric(names(yt)));
      mvy <- mvy[mvy>=rx[1] & mvy<=rx[2]];
      yh[mvy-rx[1]+1] <- -1;
    }
  }

  rm(xt,yt);

  if(round.up) { round.up <- 1; } else { round.up <- 0; }

  storage.mode(xh) <- storage.mode(yh) <- "integer";
  storage.mode(bg.xh) <- storage.mode(bg.yh) <- "integer";
  nx <- length(xh);   storage.mode(nx) <- storage.mode(whs) <- storage.mode(bg.whs) <- "integer";
  rp <- as.integer(return.peaks);
  dcon <- as.integer(direct.count);
  storage.mode(rp) <- storage.mode(min.dist) <- "integer";
  storage.mode(min.thr) <- "double";
  storage.mode(dcon) <- "integer";
  storage.mode(tag.weight) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(bg.subtract) <- "integer";
  storage.mode(round.up) <- "integer";
  im <- as.integer(ignore.masking);
  storage.mode(im) <- "integer";
  z <- .Call("spp_wtd",xh,yh,whs,rp,min.dist,min.thr,dcon,tag.weight,im,bg.subtract,bg.xh,bg.yh,bg.whs,bg.weight,round.up,PACKAGE="spp");
  if(return.peaks) {
    return(data.frame(x=(z$x+rx[1])*step,y=z$v));
  } else {
    return(list(x=rx*step,y=z));
  }
}

#' @export
tag.wtd <- function(ctv,s,e,return.peaks=T, bg.ctv=NULL,  mask.ctv=NULL, ...) {
  x <- ctv[ctv>=s & ctv<=e];
  y <- (-1)*ctv[ctv<=-s & ctv>=-e];

  if(!is.null(bg.ctv)) {
    bg.x <- bg.ctv[bg.ctv>=s & bg.ctv<=e];
    bg.y <- (-1)*bg.ctv[bg.ctv<=-s & bg.ctv>=-e];
  } else {
    bg.x <- bg.y <- NULL;
  }

  if(!is.null(mask.ctv)) {
    mask.x <- mask.ctv[mask.ctv>=s & mask.ctv<=e];
    mask.y <- (-1)*mask.ctv[mask.ctv<=-s & mask.ctv>=-e];
  } else {
    mask.x <- mask.y <- NULL;
  }

  if(length(x)==0 | length(y) ==0) {
    if(return.peaks) {
      return(data.frame(x=c(),y=c()));
    } else {
      rx <- range(c(x,y));
      return(list(x=rx,y=numeric(diff(rx)+1)));
    }
  } else {
    return(wtd(x,y,s,e,return.peaks=return.peaks,  bg.x=bg.x,bg.y=bg.y, mask.x=mask.x,mask.y=mask.y, ...))
  }
}

# calculate window cross-correlation
lwcc <- function(x,y,s,e,whs=100,isize=20,return.peaks=T,min.thr=1,min.dist=100,step=1,tag.weight=1,bg.x=NULL,bg.y=NULL,bg.weight=NULL,mask.x=NULL,mask.y=NULL,bg.whs=whs,round.up=F) {
  if(step>1) {
    x <- floor(x/step+0.5); y <- floor(y/step+0.5)

    if(!is.null(bg.x)) {
      bg.x <- floor(bg.x/step+0.5); bg.y <- floor(bg.y/step+0.5)
    }

    if(!is.null(mask.x)) {
      mask.x <- floor(mask.x/step+0.5); mask.y <- floor(mask.y/step+0.5)
    }

    whs <- floor(whs/step+0.5);
    bg.whs <- floor(bg.whs/step+0.5);
    isize <- floor(isize/step+0.5);
    min.dist <- floor(min.dist/step +0.5);
    s <- floor(s/step+0.5)
    e <- floor(e/step+0.5)
  }

  # scale bg.weight, since within calculation they are considered independent
  bg.weight <- bg.weight*tag.weight;


  rx <- c(s-whs,e+whs);
  xt <- table(x);
  xh <- integer(diff(rx)+1);
  xh[as.integer(names(xt))-rx[1]+1] <- as.integer(xt);

  yt <- table(y);

  yh <- integer(diff(rx)+1);
  yh[as.integer(names(yt))-rx[1]+1] <- as.integer(yt);

  # compile background vectors
  if(!is.null(bg.x) && length(bg.x)>0) {
    bg.subtract <- 1;

    bg.xt <- table(bg.x);
    bg.xh <- integer(diff(rx)+1);
    bg.xh[as.integer(names(bg.xt))-rx[1]+1] <- as.integer(bg.xt);
    rm(bg.xt);

    bg.yt <- table(bg.y);
    bg.yh <- integer(diff(rx)+1);
    bg.yh[as.integer(names(bg.yt))-rx[1]+1] <- as.integer(bg.yt);
    rm(bg.yt);

    # adjust bg.weight according to bg.whs
    bg.weight <- bg.weight*(whs-isize)/bg.whs;
  } else {
    bg.subtract <- 0;
    bg.xh <- bg.yh <- c();
  }

  # record masked positions
  if(!is.null(mask.x) && length(mask.x)>0) {
    mvx <- unique(mask.x); mvx <- setdiff(mvx,as.numeric(names(xt)));
    mvx <- mvx[mvx>=rx[1] & mvx<=rx[2]];

    xh[mvx-rx[1]+1] <- -1;
  }

  if(!is.null(mask.y) && length(mask.y)>0) {
    mvy <- unique(mask.y); mvy <- setdiff(mvy,as.numeric(names(yt)));
    mvy <- mvy[mvy>=rx[1] & mvy<=rx[2]];
    yh[mvy-rx[1]+1] <- -1;
  }

  rm(xt,yt);
  if(round.up) { round.up <- 1; } else { round.up <- 0; }

  storage.mode(xh) <- storage.mode(yh) <- "integer";
  storage.mode(bg.xh) <- storage.mode(bg.yh) <- "integer";
  nx <- length(xh);   storage.mode(nx) <- storage.mode(whs) <- storage.mode(isize) <- storage.mode(bg.whs) <- "integer";
  rp <- as.integer(return.peaks);
  storage.mode(rp) <- storage.mode(min.dist) <- "integer";
  storage.mode(min.thr) <- "double";
  storage.mode(tag.weight) <- "double";
  storage.mode(bg.weight) <- "double";
  storage.mode(bg.subtract) <- "integer";
  storage.mode(round.up) <- "integer";

  # allocate return arrays
  #cc <- numeric(nx); storage.mode(cc) <- "double";
  z <- .Call("spp_lwcc",xh,yh,whs,isize,rp,min.dist,min.thr,tag.weight,bg.subtract,bg.xh,bg.yh,bg.whs,bg.weight,round.up,PACKAGE="spp");
  if(return.peaks) {
    return(data.frame(x=(z$x+rx[1])*step,y=z$v));
  } else {
    return(list(x=rx*step,y=z));
  }
}

#' @export
tag.lwcc <- function(ctv,s,e,return.peaks=T, bg.ctv=NULL, mask.ctv=NULL, ...) {
  x <- ctv[ctv>=s & ctv<=e];
  y <- (-1)*ctv[ctv<=-s & ctv>=-e];

  if(!is.null(bg.ctv)) {
    bg.x <- bg.ctv[bg.ctv>=s & bg.ctv<=e];
    bg.y <- (-1)*bg.ctv[bg.ctv<=-s & bg.ctv>=-e];
  } else {
    bg.x <- bg.y <- NULL;
  }

  if(!is.null(mask.ctv)) {
    mask.x <- mask.ctv[mask.ctv>=s & mask.ctv<=e];
    mask.y <- (-1)*mask.ctv[mask.ctv<=-s & mask.ctv>=-e];
  } else {
    mask.x <- mask.y <- NULL;
  }

  if(length(x)==0 | length(y) ==0) {
    if(return.peaks) {
      return(data.frame(x=c(),y=c()));
    } else {
      rx <- range(c(x,y));
      return(list(x=rx,y=numeric(diff(rx)+1)));
    }
  } else {
    return(lwcc(x,y, s,e,return.peaks=return.peaks, bg.x=bg.x,bg.y=bg.y,  mask.x=mask.x,mask.y=mask.y, ...))
  }
}
####
# End of SPP section
####

####
# Bay2Ctrls functions
####
#' @export
shiftTag <- function (npl, tagShift = 0, tagLen = 50, outfile = 'test') {
    chrl <- names(npl)
    names(chrl) <- chrl
    cat("" , file = outfile, sep = "\t", fill = FALSE, labels = NULL,
         append = FALSE)
    npl <- lapply(chrl, function(chr) {
        rp.chr <- npl[[chr]]
    rchr <- rep(chr,length(rp.chr))
    rstrand <- rep(1,length(rp.chr))
    rstrand[rp.chr<0] <- 0
    rp <- abs(rp.chr + tagShift)
    rs <- rp
    re <- rp
    rs[rstrand == 0] <- rp[rstrand == 0] - tagLen
    re[rstrand == 1] <- rp[rstrand == 1] + tagLen
    rs[rs<0] <- 1
      utils::write.table(cbind(chr,rs,re), file = outfile, append = TRUE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = FALSE, qmethod = c("escape", "double"),
                 fileEncoding = "")
        return(rp.chr)
    })
    return(npl)
}

#' @export
mcbin <- function (n.a, n.ap, n.b, n.bp, I = 1e+06) {
    theta1 <- stats::rbeta(I, n.a + 1, n.ap + 1)
    theta2 <- stats::rbeta(I, n.b + 1, n.bp + 1)
    diff <- theta1 - theta2
    diff.hit <- diff[diff > 0 & theta1 > 0.5]
    p <- log(I - length(diff.hit) + 1) - log(I)
    return(p)
}

#' @export
mcbin.a2ap <- function (n.a, n.ap, n.b, n.bp, I = 1e+06) {
    theta1 <- stats::rbeta(I, n.a + 1, n.ap + 1)
#    theta2 <- stats::rbeta(I, n.b + 1, n.bp + 1)
#    diff <- theta1 - theta2
#    diff.hit <- diff[diff > 0 & theta1 > 0.5]
    diff.hit <- theta1[theta1 > 0.5]
    p <- log(I - length(diff.hit) + 1) - log(I)
    return(p)
}

#' @export
mcbin.a2b <- function (n.a, n.ap, n.b, n.bp, I = 1e+06) {
#    theta1 <- stats::rbeta(I, n.a + 1, n.ap + 1)
    theta1 <- stats::rbeta(I, n.a + 1, n.b + 1)
#    theta2 <- stats::rbeta(I, n.b + 1, n.bp + 1)
#    diff <- theta1 - theta2
#    diff.hit <- diff[diff > 0 & theta1 > 0.5]
    diff.hit <- theta1[theta1 > 0.5]
    p <- log(I - length(diff.hit) + 1) - log(I)
    return(p)
}

#' @export
mcbin.test <- function (mx, I = 1e+06) {
    apply(mx, 1, function(mmx) {
        return(mcbin(mmx[1], mmx[2], mmx[3], mmx[4], I))
    })
}

#' @export
mcbin.test.a2ap <- function (mx, I = 1e+06) {
    apply(mx, 1, function(mmx) {
        return(mcbin.a2ap(mmx[1], mmx[2], mmx[3], mmx[4], I))
    })
}

#' @export
mcbin.test.a2b <- function (mx, I = 1e+06) {
    apply(mx, 1, function(mmx) {
        return(mcbin.a2b(mmx[1], mmx[2], mmx[3], mmx[4], I))
    })
}

#' @export
add.mcbin2matrix <- function (nta, ntb, ntap, ntbp, I = 1e+06) {
    nta <- as.numeric(nta)
    ntb <- as.numeric(ntb)
    ntap <- as.numeric(ntap)
    ntbp <- as.numeric(ntbp)
    mcbin.pv <- 0 - mcbin.test(cbind(nta, ntap, ntb, ntbp),I)
    return(as.numeric(mcbin.pv))
}

#' @export
add.mcbin2matrix.a2ap <- function (nta, ntb, ntap, ntbp, I = 1e+06) {
    nta <- as.numeric(nta)
    ntb <- as.numeric(ntb)
    ntap <- as.numeric(ntap)
    ntbp <- as.numeric(ntbp)
    mcbin.pv <- 0 - mcbin.test.a2ap(cbind(nta, ntap, ntb, ntbp),I)
    return(as.numeric(mcbin.pv))
}

#' @export
add.mcbin2matrix.a2b <- function (nta, ntb, ntap, ntbp, I = 1e+06) {
    nta <- as.numeric(nta)
    ntb <- as.numeric(ntb)
    ntap <- as.numeric(ntap)
    ntbp <- as.numeric(ntbp)
    mcbin.pv <- 0 - mcbin.test.a2b(cbind(nta, ntap, ntb, ntbp),I)
    return(as.numeric(mcbin.pv))
}

#' @export
add.count2 <- function (rp, r.a, r.b, readSize=10000000, tagShift = 0, wdfold = 0) {
    chrl <- names(rp$npl)
    names(chrl) <- chrl
    a.sum <- sum(unlist(lapply(r.a, length)))
    b.sum <- sum(unlist(lapply(r.b, length)))
    max.sum <- readSize
    rp$npl <- lapply(chrl, function(chr) {
        rp.chr <- rp$npl[[chr]]
        if(wdfold != 0 & tagShift !=0){
            xs <- rp.chr$x - tagShift * wdfold
            xe <- rp.chr$x + tagShift * wdfold
        }
        if(wdfold == 0){
            xs <- rp.chr$rs
            xe <- rp.chr$re
        }
        xs[xs <= 0] <- 1
        rp.chr$nt.a <- spp::points_within(abs(r.a[[chr]] + tagShift), xs, xe, return.point.counts = T)
        rp.chr$nt.a <- round(rp.chr$nt.a * (max.sum/a.sum), digits = 0) + 1
        rp.chr$nt.a[is.na(rp.chr$nt.a)] <- 1

        rp.chr$nt.b <- spp::points_within(abs(r.b[[chr]] + tagShift), xs, xe, return.point.counts = T)
        rp.chr$nt.b <- round(rp.chr$nt.b * (max.sum/b.sum), digits = 0) + 1
        rp.chr$nt.b[is.na(rp.chr$nt.b)] <- 1
    rp.chr$nt.ap <- rp.chr$nt.b
    rp.chr$nt.bp <- rp.chr$nt.b

        rp.chr$nt.diff <- rp.chr$nt.a - rp.chr$nt.b
        rp.chr$nt.diffaap <- rp.chr$nt.a - rp.chr$nt.b
        rp.chr$nt.diffab <- rp.chr$nt.a - rp.chr$nt.b
        rp.chr$nt.ratioaap <- log(rp.chr$nt.a) - log(rp.chr$nt.b)
        rp.chr$nt.ratioab <- log(rp.chr$nt.a) - log(rp.chr$nt.b)
        rp.chr$nt.ratio <- log(rp.chr$nt.a) - log(rp.chr$nt.b)
        rp.chr$nt.diffRatio <- (rp.chr$nt.a - rp.chr$nt.b)/rp.chr$nt.b
        return(rp.chr)
    })
    return(rp$npl)
}

#' @export
add.count4 <- function (rp, r.a, r.b, r.ap, r.bp, readSize=10000000, tagShift = 0, wdfold = 0) {
    chrl <- names(rp$npl)
    names(chrl) <- chrl
    a.sum <- sum(unlist(lapply(r.a, length)))
    ap.sum <- sum(unlist(lapply(r.ap, length)))
    b.sum <- sum(unlist(lapply(r.b, length)))
    bp.sum <- sum(unlist(lapply(r.bp, length)))
    max.sum <- readSize
    rp$npl <- lapply(chrl, function(chr) {
        rp.chr <- rp$npl[[chr]]
    if(wdfold != 0 & tagShift !=0){
      xs <- rp.chr$x - tagShift * wdfold
      xe <- rp.chr$x + tagShift * wdfold
    }
    if(wdfold == 0){
      xs <- rp.chr$rs
      xe <- rp.chr$re
    }
        xs[xs <= 0] <- 1
        rp.chr$nt.a <- spp::points_within(abs(r.a[[chr]] + tagShift), xs, xe, return.point.counts = T)
        rp.chr$nt.ap <- spp::points_within(abs(r.ap[[chr]] + tagShift), xs, xe, return.point.counts = T)
        rp.chr$nt.b <- spp::points_within(abs(r.b[[chr]] + tagShift), xs, xe, return.point.counts = T)
        rp.chr$nt.bp <- spp::points_within(abs(r.bp[[chr]] + tagShift), xs, xe, return.point.counts = T)
        rp.chr$nt.a <- round(rp.chr$nt.a * (max.sum/a.sum), digits = 0) + 1
        rp.chr$nt.a[is.na(rp.chr$nt.a)] <- 1
        rp.chr$nt.ap <- round(rp.chr$nt.ap * (max.sum/ap.sum), digits = 0) + 1
        rp.chr$nt.ap[is.na(rp.chr$nt.ap)] <- 1
        rp.chr$nt.b <- round(rp.chr$nt.b * (max.sum/b.sum), digits = 0) + 1
        rp.chr$nt.b[is.na(rp.chr$nt.b)] <- 1
        rp.chr$nt.bp <- round(rp.chr$nt.bp * (max.sum/bp.sum), digits = 0) + 1
        rp.chr$nt.bp[is.na(rp.chr$nt.bp)] <- 1
        rp.chr$nt.diff <- (rp.chr$nt.a - rp.chr$nt.ap) - (rp.chr$nt.b - rp.chr$nt.bp)
        rp.chr$nt.diffaap <- (rp.chr$nt.a - rp.chr$nt.ap)
        rp.chr$nt.diffab <- (rp.chr$nt.a - rp.chr$nt.b)
        rp.chr$nt.ratioaap <- log(rp.chr$nt.a) - log(rp.chr$nt.ap)
        rp.chr$nt.ratioab <- log(rp.chr$nt.a) - log(rp.chr$nt.b)
        rp.chr$nt.ratio <- (log(rp.chr$nt.a) - log(rp.chr$nt.ap)) - (log(rp.chr$nt.b) - log(rp.chr$nt.bp))
        rp.chr$nt.diffRatio <- (rp.chr$nt.a - rp.chr$nt.b)/rp.chr$nt.ap
        return(rp.chr)
    })
    return(rp$npl)
}

#' @export
add.rsre <- function (bd, fname, margin = bd$whs) {
    if (is.null(margin)) {
        margin <- 50
    }
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    npl <- lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        x <- df$x
        rs <- df$rs
        if (is.null(rs)) {
            rs <- rep(NA, length(x))
        }
        re <- df$re
        if (is.null(re)) {
            re <- rep(NA, length(x))
        }
        ivi <- which(is.na(rs))
        if (any(ivi)) {
            rs[ivi] <- x[ivi] - margin
        }
        ivi <- which(is.na(re))
        if (any(ivi)) {
            re[ivi] <- x[ivi] + margin
        }
        df$rs <- rs
        df$re <- re
        return(df)
    })
    return(npl)
}

#' @export
pois.test <- function (mx) {
    apply(mx, 1, function(mmx) {
        return(stats::poisson.test(mmx, alternative = "g")$p.value)
    })
}

#' @export
add.pois <- function (rp, flg = "aap"){
    chrl <- names(rp$npl)
    names(chrl) <- chrl
    rp$npl <- lapply(chrl, function(chr) {
        rp.chr <- rp$npl[[chr]]
        if (!is.null(rp.chr)) {
            if (dim(rp.chr)[1] > 0) {
                if(flg == "aap"){
                    rp.chr$pvaap <- pois.test(cbind(rp.chr$nt.a, rp.chr$nt.ap))
                    rp.chr$pvaap <- rp.chr$pvaap + 2.225074e-307
                    rp.chr$pvaap <- 0 - log(rp.chr$pvaap)
                    rp.chr$pvab <- rp.chr$pvaap

                }else if(flg == "ab"){
                    rp.chr$pvab <- pois.test(cbind(rp.chr$nt.a, rp.chr$nt.b))
                    rp.chr$pvab <- rp.chr$pvab + 2.225074e-307
                    rp.chr$pvab <- 0 - log(rp.chr$pvab)
                    rp.chr$pvaap <- rp.chr$pvab

                }else if(flg == "both"){
                    rp.chr$pvab <- pois.test(cbind(rp.chr$nt.a, rp.chr$nt.b))
                    rp.chr$pvab <- rp.chr$pvab + 2.225074e-307
                    rp.chr$pvab <- 0 - log(rp.chr$pvab)
                    rp.chr$pvaap <- pois.test(cbind(rp.chr$nt.a, rp.chr$nt.ap))
                    rp.chr$pvaap <- rp.chr$pvaap + 2.225074e-307
                    rp.chr$pvaap <- 0 - log(rp.chr$pvaap)
                }
                return(rp.chr)
            }
        }
    })
    return(rp$npl)
}

#' @export
write.narrowpeak.mcbin <- function (bd, fname, rby = "spp", mcstep = 1e+06, coe.mc = 1, coe.spp = 1, npeaks = NA) {
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    mcstep <- as.numeric(mcstep)
    npeaks <- as.numeric(npeaks)
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        if (!is.null(df)) {
            if (dim(df)[1] > 0) {
            x <- df$x; rs <- df$rs; re <- df$re; pv <- df$pv
                cbind(chr, rs, re, ".", "0", ".", as.numeric(df$y), 0, df$nt.ratio, x - rs, df$nt.diff, df$pvab, df$pvaap, df$nt.a, df$nt.b, df$nt.ap, df$nt.bp, "0")
                #      1   2   3    4    5    6     7               8   9           10       11          12       13        14       15        16       17        18
            }
        }
    }))

    if (!is.na(npeaks)) {
        npeaks <- min(nrow(md), npeaks)
    }

    if(rby == "pvmc"){
        mcbin.pv <- add.mcbin2matrix(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
        md <- md[order(as.numeric(md[, 8]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
    }else if(rby == "spp"){
        r1.col <- 7
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 12]), as.numeric(md[, 13]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }else{
        if(rby == "pvab"){
            r1.col <- 12
        }else if(rby == "pvaap"){
            r1.col <- 13
        }else{
            stop("-rankby is missing")
        }
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }

    val.mcbinPV <- as.numeric(md[, 8])
    val.spp <- as.numeric(md[, 7])
    val.pois4ab.pv <- as.numeric(md[, 12])
    val.pois4aap.pv <- as.numeric(md[, 13])
#    rank.mcbinPV <- rank(val.mcbinPV, ties.method = c("max"))
#    rank.spp <- rank(val.spp, ties.method = c("max"))
#    rank.pois4ab.pv <- rank(val.pois4ab.pv, ties.method = c("max"))
#    rank.pois4aap.pv <- rank(val.pois4aap.pv, ties.method = c("max"))
#    rank.sum <- rank.mcbinPV*coe.mc + rank.spp*coe.spp
#    rank.mc <- order(as.numeric(val.mcbinPV), as.numeric(val.spp), decreasing = T)

#    rank.mcbinPV <- order(val.mcbinPV,val.spp,val.pois4ab.pv,val.pois4aap.pv,decreasing = F)
#    rank.spp <- order(val.spp,val.mcbinPV,val.pois4ab.pv,val.pois4aap.pv,decreasing = F)
#    rank.pois4ab.pv <- order(val.pois4ab.pv,val.spp,decreasing = F)
#    rank.pois4aap.pv <- order(val.pois4aap.pv,val.spp,decreasing = F)
#    rank.sum <- rank.mcbinPV*coe.mc + rank.spp*coe.spp

    rank.mcbinPV <- rank(as.numeric(val.mcbinPV))
    rank.spp <- rank(as.numeric(val.spp))
    rank.pois4ab.pv <- rank(as.numeric(val.pois4ab.pv))
    rank.pois4aap.pv <- rank(as.numeric(val.pois4aap.pv))
    rank.sum <- rank.mcbinPV*coe.mc + rank.spp*coe.spp

    md[,9] <- md[,8]
    md[,8] <- rank.sum

    md <- md[order(as.numeric(md[, 8]), decreasing = T), ]
    md <- md[, 1:10]
    utils::write.table(md, file = fname, col.names = F, row.names = F, quote = F, sep = "\t", append = F)
}


#' @export
write.narrowpeak.mcbin.alone <- function (bd, fname, rby = "spp", mcstep = 1e+06, coe.mc = 1, coe.spp = 1, npeaks = NA) {
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    mcstep <- as.numeric(mcstep)
    npeaks <- as.numeric(npeaks)
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        if (!is.null(df)) {
            if (dim(df)[1] > 0) {
            x <- df$x; rs <- df$rs; re <- df$re; pv <- df$pv
                cbind(chr, rs, re, ".", "0", ".", as.numeric(df$y), 0, df$nt.ratio, x - rs, df$nt.diff, df$pvab, df$pvaap, df$nt.a, df$nt.b, df$nt.ap, df$nt.bp, "0")
                #      1   2   3    4    5    6     7               8   9           10       11          12       13        14       15        16       17        18
            }
        }
    }))

    if (!is.na(npeaks)) {
        npeaks <- min(nrow(md), npeaks)
    }

    if(rby == "pvmc"){
        mcbin.pv <- add.mcbin2matrix(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
        md <- md[order(as.numeric(md[, 8]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
    }else if(rby == "spp"){
        r1.col <- 7
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 12]), as.numeric(md[, 13]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }else{
        if(rby == "pvab"){
            r1.col <- 12
        }else if(rby == "pvaap"){
            r1.col <- 13
        }else{
            stop("-rankby is missing")
        }
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }

    val.mcbinPV <- as.numeric(md[, 8])
    val.spp <- as.numeric(md[, 7])
    val.pois4ab.pv <- as.numeric(md[, 12])
    val.pois4aap.pv <- as.numeric(md[, 13])
    rank.mcbinPV <- rank(as.numeric(val.mcbinPV))
    rank.spp <- rank(as.numeric(val.spp))
    rank.pois4ab.pv <- rank(as.numeric(val.pois4ab.pv))
    rank.pois4aap.pv <- rank(as.numeric(val.pois4aap.pv))
    rank.sum <- rank.mcbinPV*coe.mc + rank.spp*coe.spp

    md[,9] <- md[,8]
    md[,8] <- rank.mcbinPV

    md <- md[order(as.numeric(md[, 8]), decreasing = T), ]
    md <- md[, 1:10]
    utils::write.table(md, file = fname, col.names = F, row.names = F, quote = F, sep = "\t", append = F)
}

#' @export
write.narrowpeak.mcbin.alone.a2ap <- function (bd, fname, rby = "spp", mcstep = 1e+06, coe.mc = 1, coe.spp = 1, npeaks = NA) {
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    mcstep <- as.numeric(mcstep)
    npeaks <- as.numeric(npeaks)
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        if (!is.null(df)) {
            if (dim(df)[1] > 0) {
            x <- df$x; rs <- df$rs; re <- df$re; pv <- df$pv
                cbind(chr, rs, re, ".", "0", ".", as.numeric(df$y), 0, df$nt.ratio, x - rs, df$nt.diff, df$pvab, df$pvaap, df$nt.a, df$nt.b, df$nt.ap, df$nt.bp, "0")
                #      1   2   3    4    5    6     7               8   9           10       11          12       13        14       15        16       17        18
            }
        }
    }))

    if (!is.na(npeaks)) {
        npeaks <- min(nrow(md), npeaks)
    }

    if(rby == "pvmc"){
        mcbin.pv <- add.mcbin2matrix.a2ap(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
        md <- md[order(as.numeric(md[, 8]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
    }else if(rby == "spp"){
        r1.col <- 7
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 12]), as.numeric(md[, 13]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix.a2ap(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }else{
        if(rby == "pvab"){
            r1.col <- 12
        }else if(rby == "pvaap"){
            r1.col <- 13
        }else{
            stop("-rankby is missing")
        }
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix.a2ap(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }

    val.mcbinPV <- as.numeric(md[, 8])
    val.spp <- as.numeric(md[, 7])
    val.pois4ab.pv <- as.numeric(md[, 12])
    val.pois4aap.pv <- as.numeric(md[, 13])
    rank.mcbinPV <- rank(as.numeric(val.mcbinPV))
    rank.spp <- rank(as.numeric(val.spp))
    rank.pois4ab.pv <- rank(as.numeric(val.pois4ab.pv))
    rank.pois4aap.pv <- rank(as.numeric(val.pois4aap.pv))
    rank.sum <- rank.mcbinPV*coe.mc + rank.spp*coe.spp

    md[,9] <- md[,8]
    md[,8] <- rank.mcbinPV

    md <- md[order(as.numeric(md[, 8]), decreasing = T), ]
    md <- md[, 1:10]
    utils::write.table(md, file = fname, col.names = F, row.names = F, quote = F, sep = "\t", append = F)
}

#' @export
write.narrowpeak.mcbin.alone.a2b <- function (bd, fname, rby = "spp", mcstep = 1e+06, coe.mc = 1, coe.spp = 1, npeaks = NA) {
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    mcstep <- as.numeric(mcstep)
    npeaks <- as.numeric(npeaks)
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        if (!is.null(df)) {
            if (dim(df)[1] > 0) {
            x <- df$x; rs <- df$rs; re <- df$re; pv <- df$pv
                cbind(chr, rs, re, ".", "0", ".", as.numeric(df$y), 0, df$nt.ratio, x - rs, df$nt.diff, df$pvab, df$pvaap, df$nt.a, df$nt.b, df$nt.ap, df$nt.bp, "0")
                #      1   2   3    4    5    6     7               8   9           10       11          12       13        14       15        16       17        18
            }
        }
    }))

    if (!is.na(npeaks)) {
        npeaks <- min(nrow(md), npeaks)
    }

    if(rby == "pvmc"){
        mcbin.pv <- add.mcbin2matrix.a2b(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
        md <- md[order(as.numeric(md[, 8]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
    }else if(rby == "spp"){
        r1.col <- 7
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 12]), as.numeric(md[, 13]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix.a2b(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }else{
        if(rby == "pvab"){
            r1.col <- 12
        }else if(rby == "pvaap"){
            r1.col <- 13
        }else{
            stop("-rankby is missing")
        }
        md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 7]), decreasing = T), ]
        md <- md[1:npeaks, ]
        mcbin.pv <- add.mcbin2matrix.a2b(md[, 14], md[, 15], md[, 16], md[, 17], I = mcstep)
        md[,8] <- mcbin.pv
    }

    val.mcbinPV <- as.numeric(md[, 8])
    val.spp <- as.numeric(md[, 7])
    val.pois4ab.pv <- as.numeric(md[, 12])
    val.pois4aap.pv <- as.numeric(md[, 13])
    rank.mcbinPV <- rank(as.numeric(val.mcbinPV))
    rank.spp <- rank(as.numeric(val.spp))
    rank.pois4ab.pv <- rank(as.numeric(val.pois4ab.pv))
    rank.pois4aap.pv <- rank(as.numeric(val.pois4aap.pv))
    rank.sum <- rank.mcbinPV*coe.mc + rank.spp*coe.spp

    md[,9] <- md[,8]
    md[,8] <- rank.mcbinPV

    md <- md[order(as.numeric(md[, 8]), decreasing = T), ]
    md <- md[, 1:10]
    utils::write.table(md, file = fname, col.names = F, row.names = F, quote = F, sep = "\t", append = F)
}

#' @export
write.narrowpeak.pois <- function (bd, fname, rby = "pvalue", npeaks = NA) {
    chrl <- names(bd$npl)
    names(chrl) <- chrl
    npeaks <- as.numeric(npeaks)
    md <- do.call(rbind, lapply(chrl, function(chr) {
        df <- bd$npl[[chr]]
        if (!is.null(df)) {
            if (dim(df)[1] > 0) {
                x <- df$x; rs <- df$rs; re <- df$re; pv <- df$pv
                cbind(chr, rs, re, ".", "0", ".", as.numeric(df$y), 0, df$nt.ratio, x - rs, df$nt.diff, df$pvab, df$pvaap, df$nt.a, df$nt.b, df$nt.ap, df$nt.bp, "0")
                #      1   2   3    4    5    6     7               8   9           10       11          12       13        14       15        16       17        18
            }
        }
    }))

    if (!is.na(npeaks)) {
        npeaks <- min(nrow(md), npeaks)
    }

    r2.col <- 8
    if(rby == "pvalue"){
        r1.col <- 13
    }else if(rby == "ratio"){
        r1.col <- 9
        r2.col <- 7
    }else{
        stop("-rankby is missing")
    }
    md <- md[order(as.numeric(md[, r1.col]), as.numeric(md[, 7]), decreasing = T), ]
    md <- md[1:npeaks, ]

    md.tmp <- md[,9]; md[,9] <- md[,7]; md[,7] <- md.tmp
    md[,8] <- md[,12]
    md <- md[, 1:10]
    md <- md[order(as.numeric(md[, r2.col]), decreasing = T), ]
    utils::write.table(md, file = fname, col.names = F, row.names = F, quote = F, sep = "\t", append = F)
}

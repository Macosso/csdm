# Run from repository root: Rscript --vanilla review/reproduce.R
# Sources the reviewed checkout, without installing or modifying the package.
for (f in list.files('R', pattern = '\\.R$', full.names = TRUE)) source(f)
check <- function(label, expr) {
  cat('\n--- ', label, ' ---\n', sep = '')
  tryCatch(print(suppressMessages(suppressWarnings(force(expr)))),
           error = function(e) cat('ERROR:', conditionMessage(e), '\n'))
}
set.seed(92026)
d <- expand.grid(time = 1:40, id = 1:8)
d$x <- rnorm(nrow(d)); d$z <- rnorm(nrow(d))
d$y <- 1 + .6*d$x + .2*d$z + rnorm(nrow(d))
fit <- function(dat=d, formula=y~x+z, model='mg', ...) {
  suppressMessages(suppressWarnings(csdm(formula, dat, 'id', 'time', model=model, ...)))
}
m <- fit()
check('Baseline MG vs independent lm means', {
  b <- sapply(split(d,d$id), function(s) coef(lm(y~x+z,s)))
  c(coef_max_error=max(abs(coef(m)-rowMeans(b))),
    vcov_max_error=max(abs(vcov(m)-cov(t(b))/ncol(b))))
})
check('HC0 covariance: nonnegative expected', {
  X <- matrix(1,4,1); u <- c(-1,-1,-1,-1)
  c(actual=sandwich_vcov(X,u)[1,1], expected=1/4)
})
check('MG utility variance scaling', {
  B <- matrix(1:4,ncol=1)
  c(actual=pooled_vcov(B)[1,1], expected=var(1:4)/4)
})
check('Named CSA lags', csdm_csa(lags=c(y=2,x=1))$lags)
check('Named CSA lags fitting', fit(model='dcce',csa=csdm_csa(vars=c('y','x'),lags=c(y=2,x=1))))
check('Ignored covariance choices', sapply(c('mg','np','nw','wpn','ols'),function(v) max(abs(vcov(fit(vcov=csdm_vcov(v)))-vcov(m)))))
check('Ignored subset, weights, and MG lag request', c(
  subset_identical=identical(coef(fit(subset=id<4)),coef(m)),
  weights_identical=identical(coef(fit(weights=seq_len(nrow(d)))),coef(m)),
  lr_identical=identical(coef(fit(lr=csdm_lr(type='ardl',ylags=1))),coef(m))))
check('Transformed response', coef(fit(formula=I(y^2)~x+z)))
check('Dot formula includes bookkeeping predictor', names(coef(fit(formula=y~.))))
check('Duplicate indexes silently accepted / observations overwritten', {
  dd <- rbind(d,d[1,]); md <- fit(dd)
  c(input_rows=nrow(dd),reported_nobs=md$stats$nobs)
})
check('Time gap bridges lag', {
  dd <- d[!(d$id==1 & d$time==20),]
  mm <- fit(dd,model='dcce',csa=csdm_csa('_none'),lr=csdm_lr(type='ardl',ylags=1))
  c(residual_at_21=mm$residuals_e['1','21'], should_be_missing=TRUE)
})
check('Partial identification means/covariance use different units', {
  dd <- d; dd$z[dd$id==1] <- 0
  mm <- fit(dd)
  c(unit1_x=mm$coef_i['1','x'],unit1_z=mm$coef_i['1','z'],
    reported_x=coef(mm)['x'],mean_complete_units=mean(mm$coef_i[-1,'x']),
    actual_x_var=vcov(mm)['x','x'],all_x_var=var(mm$coef_i[,'x'])/8,
    dropped_units=length(mm$meta$dropped_units))
})
check('All units unidentified still return fit', {
  dd <- d; dd$z <- 0
  mm <- fit(dd); list(coef=coef(mm),vcov=vcov(mm))
})
ma <- fit(model='cs_ardl',csa=csdm_csa('_none'),lr=csdm_lr(type='ardl',ylags=1,xdlags=1))
check('CS-ARDL coefficient/covariance mismatch', list(coef=names(coef(ma)),vcov=colnames(vcov(ma)),confint=tryCatch(confint(ma),error=conditionMessage)))
check('Base modelling methods', lapply(c('nobs','fitted','df.residual','formula','model.frame'),function(g) {
  ans <- tryCatch(do.call(g,list(m)),error=conditionMessage)
  list(method=g,class=class(ans),length=length(ans),value=if(length(ans)<5) ans else head(as.vector(ans)))
}))
check('Fitting consumes global RNG', {
  set.seed(99); before <- .Random.seed; invisible(fit()); !identical(before,.Random.seed)
})
check('Leave-one-out mean at missing own observation', {
  a <- data.frame(id=1:3,time=1,x=c(NA,2,4))
  cross_sectional_avg(a,'id','time','x',leave_out=TRUE)$csa_x
})
check('na.rm FALSE still omits missing observations', {
  a <- data.frame(id=1:3,time=1,x=c(NA,2,4))
  cross_sectional_avg(a,'id','time','x',na.rm=FALSE)
})
check('Classic CD drops all times for one all-missing unit', {
  E <- m$residuals_e; E[1,] <- NA
  cd_test(E,type='CD')
})
check('CDstar accepts impossible PCA count', cd_test(matrix(rnorm(12),4,3),type='CDstar'))
check('CDw+ empirical null rejection rate (200 IID panels, N=20,T=100)', {
  set.seed(27182)
  p <- replicate(200, {
    E <- matrix(rnorm(20*100),20)
    a <- cd_test(E,type='CDw+')
    c(CDw=a$tests$CDw$p.value,CDw_plus=a$tests$CDw_plus$p.value)
  })
  rowMeans(p<.05)
})
check('Static CCE vs independently augmented lm', {
  av <- aggregate(cbind(y,x,z)~time,d,mean)
  names(av)[-1] <- c('cy','cx','cz')
  dd <- merge(d,av,by='time')
  b <- sapply(split(dd,dd$id),function(s) coef(lm(y~x+z+cy+cx+cz,s))[c('(Intercept)','x','z')])
  mm <- fit(model='cce')
  c(coef_max_error=max(abs(coef(mm)-rowMeans(b))),vcov_max_error=max(abs(vcov(mm)-cov(t(b))/8)))
})
check('CCE retains structural regressor fully spanned by CSA', {
  dd <- d; dd$x <- sin(dd$time)
  mm <- fit(dd,model='cce')
  c(reported_x=coef(mm)['x'],se_x=sqrt(vcov(mm)['x','x']),
    x_minus_own_csa=max(abs(dd$x-ave(dd$x,dd$time))))
})
check('Dynamic and long-run coefficients vs independent ARDL regressions', {
  av <- aggregate(cbind(y,x,z)~time,d,mean)
  names(av)[-1] <- c('cy','cx','cz')
  av$lcy <- c(NA,head(av$cy,-1)); av$lcx <- c(NA,head(av$cx,-1)); av$lcz <- c(NA,head(av$cz,-1))
  dd <- merge(d,av,by='time')
  bb <- lapply(split(dd,dd$id),function(s) {
    s <- s[order(s$time),]
    s$ly <- c(NA,head(s$y,-1)); s$lx <- c(NA,head(s$x,-1)); s$lz <- c(NA,head(s$z,-1))
    coef(lm(y~x+z+ly+lx+lz+cy+cx+cz+lcy+lcx+lcz,s))[c('(Intercept)','x','z','ly','lx','lz')]
  })
  b <- do.call(cbind,bb)
  mm <- fit(model='cs_ardl',csa=csdm_csa(lags=1),lr=csdm_lr(type='ardl',ylags=1,xdlags=1))
  expected_lr <- rowMeans(rbind(lr_y=b['ly',]-1,lr_x=(b['x',]+b['lx',])/(1-b['ly',]),lr_z=(b['z',]+b['lz',])/(1-b['ly',])))
  c(level_coef_max_error=max(abs(mm$coef_mg-rowMeans(b))),
    long_run_max_error=max(abs(tail(coef(mm),3)-expected_lr)))
})
check('CDstar vs unit-specific residual scales (Mata mean(matrix) is columnwise)', {
  set.seed(482)
  E <- outer(seq(.1,3,length.out=20),rnorm(100)) + matrix(rnorm(2000),20)
  A <- scale(t(E)); N <- ncol(A); TT <- nrow(A)
  F <- cbind(1,eigen(tcrossprod(A),symmetric=TRUE)$vectors[,1,drop=FALSE])
  B <- solve(crossprod(F),crossprod(F,A)); U <- A-F%*%B
  G <- sweep(B[-1,,drop=FALSE],1,sqrt(diag(tcrossprod(B[-1,,drop=FALSE])/N)),'/')
  s <- sqrt(colMeans(U^2))
  phi <- rowMeans(sweep(G,2,s,'/'))
  a <- as.numeric(1-t(sweep(G,2,s,'*'))%*%phi)/sqrt(N)
  ewt <- sum(a^2)
  cd <- sqrt(2*TT/(N*(N-1)))*sum(cor(U)[upper.tri(cor(U))])
  c(actual=cd_test(E,type='CDstar',n_pc=1)$tests$CDstar$statistic,
    unit_scale_reference=(cd+sqrt(TT/2)*(1-ewt))/ewt)
})

dvt<-function(betai,betaBi=NULL,treat,boot=TRUE,rand=1000,prefixi)
{  
  if(is.null(betaBi)) boot=FALSE
  dv.lmm<-function(betai,betaBi=NULL,treat,boot=FALSE,rand=1000,prefixi,save.data=FALSE)
  {
    warm.lev=unique(treat$treatment)
    betai3=dist.3col(betai)
    if(!is.null(betaBi)) betai3.B=dist.3col(betaBi)
    id1=match(betai3[,1],rownames(treat))
    id2=match(betai3[,2],rownames(treat))
    w1=treat$treatment[id1]
    w2=treat$treatment[id2]
    y1=treat$year[id1]
    y2=treat$year[id2]
    p1=treat$plot[id1]
    p2=treat$plot[id2]
    # consider block
    bys=paste0(treat$block,".",treat$year)
    by1=bys[id1]
    by2=bys[id2]
    id.use=which((by1==by2)&(w1!=w2))
    d.use=betai3[id.use,3]
    if(!is.null(betaBi)) d.use.B=betai3.B[id.use,3]
    t.use=y1[id.use]
    pp.use=sapply(id.use,function(x){paste(sort(c(p1[x],p2[x])),collapse = "")})
    if(!is.null(betaBi))
    {
      data.out=data.frame(Prefix=prefixi,betai3[id.use,1:2,drop=FALSE],t=t.use,D.A=d.use,
                          D.B=d.use.B,rand.eff=pp.use,stringsAsFactors = FALSE)
    }else{
      data.out=data.frame(Prefix=prefixi,betai3[id.use,1:2,drop=FALSE],t=t.use,D.A=d.use,rand.eff=pp.use,stringsAsFactors = FALSE)
    }
    if(save.data)  save.file(data.out,prefix=prefixi,filename = "Blk.Data.forFig")
    
    lmi=list()
    lmi[[1]]=lmer(D.A~t+((1+t)|rand.eff),data = data.out)
    if(!is.null(betaBi)) lmi[[2]]=lmer(D.B~t+((1+t)|rand.eff),data = data.out)
    
    out1=sapply(1:length(lmi),
                function(k)
                {
                  lmism=summary(lmi[[k]])
                  AIC1=AIC(lmi[[k]])
                  r2i=rsquared.glmm(lmi[[k]])
                  lmiCS=Anova(lmi[[k]],type = "II")
                  c(slop.fix=lmism$coefficients[2,1],slop.sd=attr(lmism$varcor[[1]],"stddev")[[2]],
                    R2M=r2i$Marginal,R2C=r2i$Conditional,AIC1=AIC1,AIC2=r2i$AIC,
                    P.typeII=lmiCS[[3]],Chisq=lmiCS[[1]])
                })
    if(boot)
    {
      output=list()
      pp.lev=unique(pp.use)
      t.lev=unique(t.use)
      ppt=paste0(pp.use,".",t.use)
      bt.ids=lapply(1:rand,
                    function(k)
                    {
                      pp.r=sample(pp.lev,length(pp.lev),replace = TRUE)
                      t.r=sample(t.lev,length(t.lev),replace = TRUE)
                      xx=1
                      while(length(unique(t.r))<3 & xx<100)
                      {
                        t.r=sample(t.lev,length(t.lev),replace = TRUE)
                      }
                      idr=expand.grid(1:length(pp.r),1:length(t.r))
                      ppt.r=sapply(1:nrow(idr),function(x){paste0(pp.r[idr[x,1]],".",t.r[idr[x,2]])})
                      match(ppt.r,ppt)
                    })
      trace.seq=seq(from=1,to=rand,by = 100)
      bt.s=t(sapply(1:rand,
                    function(k)
                    {
                      if(k %in% trace.seq) message("Now bootstrap k=",k,". ",date())
                      idk=bt.ids[[k]]
                      datak=data.frame(d=d.use[idk],dB=d.use.B[idk],t=t.use[idk],pp=pp.use[idk])
                      if(length(unique(pp.use[idk]))==1)
                      {
                        lmi1=lm(d~t,data=datak)
                        lmi2=lm(d~t,data=datak)
                      }else{
                        lmi1=lmer(d~t+((1+t)|pp),data = datak)
                        lmi2=lmer(dB~t+((1+t)|pp),data = datak)
                      }
                      s1=summary(lmi1)$coefficients[2,1]
                      s2=summary(lmi2)$coefficients[2,1]
                      c(s1,s2,s1-s2)
                    }))
      bt.s=rbind(c(out1[1,],out1[1,1]-out1[1,2]),bt.s)
      colnames(bt.s)=c("WA.slope","WB.slope","Slope.dif")
      rownames(bt.s)=c("obs",paste0("bt",1:rand))
      output$blk.bt=bt.s
      if(out1[1,1]>=out1[1,2]){pds1=sum(bt.s[,3]<=0)/nrow(bt.s)}else{pds1=sum(bt.s[,3]>=0)/nrow(bt.s)}
      output$p.ds.bt.blk=pds1
    }
    # not consider block
    id.use=which((y1==y2)&(w1!=w2))
    d.use=betai3[id.use,3]
    if(!is.null(betaBi)) d.use.B=betai3.B[id.use,3]
    t.use=y1[id.use]
    pp.use=sapply(id.use,function(x){paste(sort(c(p1[x],p2[x])),collapse = "")})
    if(!is.null(betaBi))
    {
      data.out=data.frame(Prefix=prefixi,betai3[id.use,1:2,drop=FALSE],t=t.use,D.A=d.use,
                          D.B=d.use.B,rand.eff=pp.use,stringsAsFactors = FALSE)
    }else{
      data.out=data.frame(Prefix=prefixi,betai3[id.use,1:2,drop=FALSE],t=t.use,D.A=d.use,rand.eff=pp.use,stringsAsFactors = FALSE)
    }
    if(save.data)  save.file(data.out,prefix=prefixi,filename = "noBlk.Data.forFig")
    
    lmi=list()
    lmi[[1]]=lmer(D.A~t+((1+t)|rand.eff),data = data.out)
    if(!is.null(betaBi)) lmi[[2]]=lmer(D.B~t+((1+t)|rand.eff),data = data.out)
    out2=sapply(1:length(lmi),
                function(k)
                {
                  lmism=summary(lmi[[k]])
                  AIC1=AIC(lmi[[k]])
                  r2i=rsquared.glmm(lmi[[k]])
                  lmiCS=Anova(lmi[[k]],type = "II")
                  c(slop.fix=lmism$coefficients[2,1],slop.sd=attr(lmism$varcor[[1]],"stddev")[[2]],
                    R2M=r2i$Marginal,R2C=r2i$Conditional,AIC1=AIC1,AIC2=r2i$AIC,
                    P.typeII=lmiCS[[3]],Chisq=lmiCS[[1]])
                })
    
    if(boot)
    {
      pp.lev=unique(pp.use)
      t.lev=unique(t.use)
      ppt=paste0(pp.use,".",t.use)
      wp=unique(treat$plot[which(treat$treatment==warm.lev[1])])
      cp=unique(treat$plot[which(treat$treatment==warm.lev[2])])
      
      bt.ids=lapply(1:rand,
                    function(k)
                    {
                      wp.r=sample(wp,length(wp),replace = TRUE)
                      cp.r=sample(cp,length(cp),replace = TRUE)
                      pp.r=sapply(1:length(wp),function(x){paste(sort(c(wp.r[x],cp.r[x])),collapse = "")})
                      t.r=sample(t.lev,length(t.lev),replace = TRUE)
                      xx=1
                      while(length(unique(t.r))<3 & xx<100)
                      {
                        t.r=sample(t.lev,length(t.lev),replace = TRUE)
                      }
                      idr=expand.grid(1:length(pp.r),1:length(t.r))
                      ppt.r=sapply(1:nrow(idr),function(x){paste0(pp.r[idr[x,1]],".",t.r[idr[x,2]])})
                      match(ppt.r,ppt)
                    })
      trace.seq=seq(from=1,to=rand,by = 100)
      bt.s=t(sapply(1:rand,
                    function(k)
                    {
                      if(k %in% trace.seq) message("Now bootstrap k=",k,". ",date())
                      idk=bt.ids[[k]]
                      datak=data.frame(d=d.use[idk],dB=d.use.B[idk],t=t.use[idk],pp=pp.use[idk])
                      if(length(unique(pp.use[idk]))==1)
                      {
                        lmi1=lm(d~t,data=datak)
                        lmi2=lm(d~t,data=datak)
                      }else{
                        lmi1=lmer(d~t+((1+t)|pp),data = datak)
                        lmi2=lmer(dB~t+((1+t)|pp),data = datak)
                      }
                      s1=summary(lmi1)$coefficients[2,1]
                      s2=summary(lmi2)$coefficients[2,1]
                      c(s1,s2,s1-s2)
                    }))
      bt.s=rbind(c(out2[1,],out2[1,1]-out2[1,2]),bt.s)
      colnames(bt.s)=c("WA.slope","WB.slope","Slope.dif")
      rownames(bt.s)=c("obs",paste0("bt",1:rand))
      output$nblk.bt=bt.s
      if(out2[1,1]>=out2[1,2]){pds2=sum(bt.s[,3]<=0)/nrow(bt.s)}else{pds2=sum(bt.s[,3]>=0)/nrow(bt.s)}
      output$p.ds.bt.nblk=pds2
    }
    outb=data.frame(out1,out2,stringsAsFactors = FALSE)
    if(!is.null(betaBi)){colnames(outb)=c("WA.block","WB.block","WA.noblock","WB.noblock")}else{colnames(outb)=c("block","noblock")}
    if(boot){output$summary=outb}else{output=outb}
    output
  }
  
  dvio=dv.lmm(betai=betai,betaBi = betaBi,treat = treat,boot = boot,rand = rand, prefixi = prefixi,save.data = TRUE)
  if(boot){dvi=dvio$summary}else{dvi=dvio}
  r2.obs=as.vector(as.matrix(dvi[3:4,]))
  aic.obs=as.vector(as.matrix(dvi[5:6,]))
  if(!is.null(betaBi)){ds.obs=c(dvi[1,1]-dvi[1,2],dvi[1,3]-dvi[1,4])}
  
  year.lev=unique(treat$year)
  year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
  trace.seq=seq(from=1,to=rand,by = 100)
  t1=Sys.time()
  ind.rand=lapply(1:nrow(year.perm),
                  function(k)
                  {
                    if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                    out=list()
                    idi=year.perm[k,]
                    perm.treat=treat
                    perm.treat[,"year"]=year.lev[idi[match(treat$year,year.lev)]]
                    dvr=dv.lmm(betai = betai, betaBi=betaBi, treat = perm.treat, boot = FALSE, prefixi = prefixi)
                    out$r2=as.vector(as.matrix(dvr[3:4,]))
                    out$aic=as.vector(as.matrix(dvr[5:6,]))
                    if(!is.null(betaBi)){out$ds=c(dvr[1,1]-dvr[1,2],dvr[1,3]-dvr[1,4])}
                    out
                  })
  format(Sys.time()-t1)
  r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
  aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
  if(!is.null(betaBi)) ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
  EPS <- sqrt(.Machine$double.eps)
  p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
  p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
  if(!is.null(betaBi))
  {
    p.ds1=(rowSums(ds.ran>=(matrix(ds.obs,nr=nrow(ds.ran),nc=ncol(ds.ran))-EPS))+1)/(ncol(ds.ran)+1)
    p.ds2=(rowSums(ds.ran<=(matrix(ds.obs,nr=nrow(ds.ran),nc=ncol(ds.ran))+EPS))+1)/(ncol(ds.ran)+1)
    p.dsp=p.ds1*(ds.obs>=0)+p.ds2*(ds.obs<0)
  }
  if(boot)
  {
    out=list()
    outp=rbind(matrix(p.r2,2,4),matrix(p.aic,2,4),
               c(dvio$p.ds.bt.blk,NA,dvio$p.ds.bt.nblk,NA),
               c(p.dsp[1],NA,p.dsp[2],NA))
    rownames(outp)=c("p.r2m","p.r2c","p.AIC1","P.AIC2","P.ds.bt","P.ds.perm")
    out$summary=rbind(as.matrix(dvi),outp)
    out$boot.blk=dvio$blk.bt
    out$boot.nblk=dvio$nblk.bt
  }else{
    if(!is.null(betaBi))
    {
      outp=rbind(matrix(p.r2,2,4),matrix(p.aic,2,4),
                 c(p.dsp[1],NA,p.dsp[2],NA))
      rownames(outp)=c("p.r2m","p.r2c","p.AIC1","P.AIC2","P.ds.perm")
    }else{
      outp=rbind(matrix(p.r2,2,2),matrix(p.aic,2,2))
      rownames(outp)=c("p.r2m","p.r2c","p.AIC1","P.AIC2")
    }
    out=rbind(as.matrix(dvi),outp)
  }
  out
}
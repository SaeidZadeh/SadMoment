library("copula")
library("BB")
library("stats4")
library("VGAM")
library("sads")
library("bbmle")
library("poilog")
library("splines")
library("BBmisc")
###########################################################################
SAD_Moment<-function(xpos=NULL,ypos=NULL,splist=NULL,ext.rate=1,dplot=T)
{
	prm<-proc.time()
	Y<-data.frame(splist,xpos,ypos);
	spp<-sort(intersect(splist,splist));
	W<-fg(splist,spp);
	TW<-W*ext.rate
	W<-log2(W[][W[]!=0])
	TW<-log2(TW[][TW[]!=0])
	larg.i2 <- round(max(TW))
	x<-round(max(xpos))
	y<-round(max(ypos))
	Area<-x*y
	Areap<-x*y*ext.rate
	ratio<-x/y;
	Ratio<-ratio
	a<-1:(ext.rate*x*y*1e-4);
	a<-a*1e4
	a<-sqrt(a/ratio)
	aa<-a[1:round(length(a)/ext.rate)];
	momdist<-matrix(data=0,nrow=20,ncol=length(aa)) 
	# nrow=20 corresponds to the case where maximum number of individuals of all species do not reach 741455. 
	momraw<-matrix(data=0,nrow=20,ncol=length(aa))
	momlin<-matrix(data=0,nrow=3,ncol=length(aa))
	if(ext.rate>1)
	{
		for(i in 1:length(aa))
		{
			n<-ceiling(2*y/aa[i]);
			tmp<-momfun(Y,aa[i],spp,n,20,ratio);
			momraw[,i]<-tmp[1,]
			momdist[,i]<-tmp[2,]
			momlin[,i]<-tmp[3,1:3]
		}
		Mom<-momdist[,length(aa)]
		ba<-ratio*a^2
		bba<-ratio*aa^2;
		momr<-extrapolation.raw(momraw,ratio,ba,bba)
		momd<-extrapolation.dist(momraw,momdist,ext.rate,momr)
		Momp<-momd
		S<-length(spp)
		larg.i <- round(max(W))
		Nbins<-larg.i+1
		Nbinsp<-larg.i2+1
		br <- seq(-0.5, larg.i2+0.5, 1)
		no.bins<-larg.i2+1
		no.m<-no.bins/2+1
		mom<-momd[1:no.m]
		no.bin<-no.bins+2
		C<-Coef.Tcheb.pol.2nd(no.bin,no.m)
		TchCo<-C
		Tp4<-C%*%mom
		TchMo<-Tp4
		tp4<-Tcheb.pol.2nd(no.bin,no.m)
		TchPo<-tp4
		fx4 <- rep(0,no.bins)
		for(i in 1:no.bins)
			fx4[i] <- sum(Tp4*tp4[,i])
		fx4 <- fx4*S
		no.binsp<-larg.i2+1
		fx4<-fx4[1:no.binsp]
		fx4[which(fx4<0)]<-0
	}
	if(ext.rate==1)
	{
		tmp<-momfun(Y,max(aa),spp,1,20,ratio)
		momd<-tmp[2,]
		S<-length(spp)
		larg.i <- round(max(W))
		br <- seq(-0.5, larg.i+0.5, 1)
		Mom<-momd
		Momp<-momd
		no.bins<-larg.i+1
		Nbins<-larg.i+1
		Nbinsp<-larg.i2+1
		no.m<-no.bins
		mom<-momd[1:no.m]
		no.bin<-no.bins
		C<-Coef.Tcheb.pol.2nd(no.bin,no.m)
		TchCo<-C
		Tp4<-C%*%mom
		TchMo<-Tp4
		tp4<-Tcheb.pol.2nd(no.bin,no.m)
		TchPo<-tp4
		fx4 <- rep(0,no.bins)
		for(i in 1:no.bins)
			fx4[i] <- sum(Tp4*tp4[,i])
		fx4 <- fx4*S
		fx4<-fx4[1:no.bins]
	}
	res<-hist(W, breaks=br,plot = F)
	lp<-res$counts
	if(dplot==T){
	  hist(W, breaks=br,main = paste0("Extrapolation to ",ext.rate," times of area"),ylab="Number of Species",ylim=c(0,max(c(fx4,lp))),xlab=substitute(paste(log[2],"(Number of Individuals)")))
	  lines(0:(larg.i2),fx4,col=4,lty=1)
	}
	INIT<-list(Area,Mom,Ratio,Nbins)
	names(INIT)<-c('Area','Moments','XY.Ratio','Number.of.bins');
	PRED<-list(Areap,Momp,Nbinsp,TchCo,TchPo,TchMo,fx4)
	names(PRED)<-c('Area','Moments','Number.of.bins','Coefficients.of.Tchebychev.Polynomials','Tchebychev.Polynomials.Values','Tchebychev.Moments.Values','Extrapolate.SAD')
	pttime<-proc.time()-prm
	dat<-list(spp,xpos,ypos,splist,length(spp),lp)
	names(dat)<-c('Species.Name','X.Position','Y.position','Species.In.Position','Number.Species','SAD')
	resu<-list(dat,INIT,PRED,pttime)
	names(resu)<-c('Data','Initial','Extrapolate','Time')
	return(invisible(resu))
}
###########################################################################
Dist_Moment<-function(xpos,ypos,splist,DIST="LS")
{
	Y<-data.frame(splist,xpos,ypos);
	spp<-sort(intersect(splist,splist));
	W<-fg(splist,spp);
	W<-log2(W[][W[]!=0])
	x<-round(max(xpos))
	y<-round(max(ypos))
	ratio<-x/y;
	aa<-1:(x*y*1e-4);
	aa<-aa*1e4
	aa<-sqrt(aa/ratio)
	momdist<-matrix(data=0,nrow=3,ncol=length(aa))
	momraw<-matrix(data=0,nrow=3,ncol=length(aa))
	momlin<-matrix(data=0,nrow=3,ncol=length(aa))
	No.Id<-rep(0,length(aa))
	No.Sp<-No.Id
	for(i in 1:length(aa))
	{
		n<-ceiling(2*y/aa[i]);
		tmp<-momfun(Y,aa[i],spp,n,3,ratio);
		momraw[,i]<-tmp[1,]
		momdist[,i]<-tmp[2,]
		momlin[,i]<-tmp[3,]
		temp<-datfun(Y,aa[i],spp,n,ratio)
		No.Sp[i]<-as.numeric(temp[1])
		No.Id[i]<-as.numeric(temp[2])
	}
	bba<-ratio*aa^2;
	tmp<-ro.c.z(No.Sp,No.Id,bba)
	ro<-tmp[1]
	c<-tmp[2]
	z<-tmp[3]
	x1<-aa*0
  	s1<-c*bba^z;
	c1<-lm(log(momlin[3,])~log(bba))
  	p2<-as.double(c1$coefficients)
  	c1<-lm(log(momdist[2,])~log(bba))
  	pl1<-as.double(c1$coefficients)
  	c1<-lm(log(momdist[3,])~log(bba))
  	pl2<-as.double(c1$coefficients)
  	mom<-matrix(data=0,nrow=20,ncol=length(bba));
  	if(DIST=="LS")
  	{
  		alpha1<-aa*0;
  		t<-1:length(aa)
  		for(i in t)
    		alpha1[i]<-calcalpha(ro*bba[i],c*bba[i]^z);
  		x1<-1-exp(-s1/alpha1)
  		for(i in t)
    		for(j in 1:20)
      			mom[j,i]<-LSmom(x1[i],j);
      	return(mom)
  	}
  	if(DIST=="PG")
  	{
  		alf<-1:length(bba)*0;
  		bet<-1:length(bba)*0;
  		for(i in 1:length(bba))
  		{
    		alf[i]<-alffinder(p2[1],p2[2],ro,c,z,bba[i])
    		bet[i]<-betfinder(p2[1],p2[2],ro,c,z,bba[i])
 	 	}
  		for(i in 1:length(bba))
   			for(j in 1:20)
      			mom[j,i]<-PGmom(alf[i]+1,bet[i],j);
      	return(mom)
	}
    if(DIST=="LN")
  	{
  		sigLN<-sigfinderLN(exp(pl2[1]),pl2[2],ro,c,z,bba);
  		muLN<-mufinderLN(exp(pl2[1]),pl2[2],bba)
  		for(i in 1:length(bba))
    		for(j in 1:20)
      			mom[j,i]<-LNmom(muLN[i],sigLN[i],j);
      	return(mom)
    }
    if(DIST=="PLN")
    {
    	sigPLN<-sigfinderPLN(exp(pl1[1]),pl1[2],exp(pl2[1]),pl2[2],bba);
  		muPLN<-mufinderPLN(exp(pl1[1]),pl1[2],bba)
  		for(i in 1:length(bba))
    		for(j in 1:20)
     			mom[j,i]<-PLNmom(muPLN[i],sigPLN[i],j);
		return(mom)
    }
}
###########################################################################
ro.c.z<-function(No.Sp,No.Id,bba)
{
	p<-lm(No.Id~bba+0)
	ro<-as.numeric(p$coefficients)
	p<-lm(log(No.Sp)~log(bba))
	tmp<-as.numeric(p$coefficients)
	return(c(ro,exp(tmp[1]),tmp[2]))
}
###########################################################################
extrapolation.raw<-function(momraw,ratio,ba,bba)
{
	A<-max(ba);
	t<-1:length(bba)
	momr<-1:20*0
	l<-1:20
	momr[1]<-1
	for(k in 2:20)
	{
		TMP<-log(momraw[k,t])
		X<-matrix(data=0,nrow=length(t),ncol=3)
		j<-1
		X[,1]<-log(bba[t])
		kk<-(log(bba)-min(log(bba)))/(max(log(bba))-min(log(bba)))*pi
		X[,2]<-sin(kk)
		X[,3]<-cos(kk)
		Y<-data.frame(X[,1:(3)])
		reslm1<-lm(TMP[t]~.,data=Y)
		l[k]<-summary(reslm1)$r.squared
		momr[k]<-exp(relpwave(reslm1$coefficients,A,bba))
	}
	return(momr)
}
###########################################################################
extrapolation.dist<-function(momraw,momdist,ext.rate,momr)
{
	p<-1:length(momraw[1,])*0
	for(k in 1:length(momraw[1,]))
	{
		c<-lm(log(momdist[,k])~log(momraw[,k])+0)
		p[k]<-as.double(c$coefficients)
	}
	momd<-momr^extpfun(p,ext.rate*length(momraw[1,]))[ext.rate*length(momraw[1,])]
	return(momd)
}
###########################################################################
manage_data<-function(df,spp,n,x,ratio)
{
	xc<-runif(n,0,round(ratio*max(df[,3]))-ratio*x);
	yc<-runif(n,0,round(max(df[,3]))-x);
	W<-rep(0,length(spp))
	for(i in 1:n)
	{
		g<-f(df,xc[i],yc[i],x,ratio);
		W<-W+fg(g,spp);
	}
	W<-W/n;
	return(W);
}
###########################################################################
provide_info<-function(L)
{
	i<-which.max(L[,1]);
	Llist<-c("Poisson Gamma","Log-Normal","Log-Series","Poisson Log-Normal")
	x<-(paste0("For the sub-area ( ",L[5,4]," ) the best fitting distribution is : ", Llist[i]))
	if(i == 1)
		y<-(paste0("alpha= ",L[1,2]," and beta= ",L[1,3]))
	if(i == 2)
		y<-(paste0("mu= ",L[2,2]," and sigma= ",L[2,3]))
	if(i == 3)
		y<-(paste0("x= ",L[3,2]))
	if(i == 4)
		y<-(paste0("mu= ",L[4,2]," , sigma= ",L[4,3]," and p= ",L[4,4]))
	return(rbind(unlist(x),unlist(y)))
}
###########################################################################
LLfun<-function(Y,spp)
{
	sa<-data_read();
	n<-ceiling(max(Y[,2])*max(Y[,3])/sa*5)
	ratio<-max(Y[,2])/max(Y[,3])
	y<-round(manage_data(Y,spp,n,sqrt(sa/ratio),ratio))
	y<-y[log(y)>=0];
	yp<-y;
	y<-log(y,2);
	L<-matrix(data=0,nrow=5,ncol=4)
	options(warn = -1)
	g<-function(x) -length(y)*digamma(x)-length(y)*log((x+mean(yp))/x)+sum(digamma(yp+x))
	t1<-try(uniroot(g,c(0.0001,10),tol = .Machine$double.eps))
	if(!is.error(t1))
	{
		t1<-uniroot(g,c(0.0001,10),tol = .Machine$double.eps)
		bet<-t1$root
		al<-bet/mean(yp,2);
		pg<-function(al,bet) -sum(lgamma(yp+bet)+bet*log(al)-lgamma(bet)-(bet+yp)*log(1+al)-lfactorial(yp))
		L[1,1]<-2*pg(al,bet)+4;
		L[1,2]<-al;
		L[1,3]<-bet;
	}
	mu<-mean(log(yp,2));
	sig<-sqrt(sum((log(yp,2)-mean(log(yp,2)))^2)/length(log(yp,2)));
	AICln<-function(mu,sig)
	{
		-2*(-length(log(yp,2))/2*log(2*pi*sig^2)-sum(log(yp,2))-sum((log(yp,2)-mu)^2/(2*sig^2)))+4
	}
	L[2,1]<-AICln(mu,sig)
	L[2,2]<-mu;
	L[2,3]<-sig;
	pmfls<-function(x,y)
	{
		-x^y/(y*log(1-x))
	}
	ls<-function(x) -sum(-log(-log(1-x))+yp*log(x)-log(yp,2))
	fit0 <- try(mle2(ls, start = list(x=0.9),method="Nelder-Mead",skip.hessian = TRUE))
	if(!is.error(fit0))
	{
		fit0 <- mle2(ls, start = list(x=0.9),method="Nelder-Mead",skip.hessian = TRUE)
		L[3,1]<-AIC(fit0);
		L[3,2]<-as.numeric(coef(fit0))
	}
	U<-poilogMLE(round(yp));
	L[4,1]<--2*U$logLval+4;
	L[4,2:3]<-as.numeric(U$par)
	L[4,4]<-as.numeric(U$p)
	L[5,4]<-sa
	options(warn = 1)
	provide_info(L)
}
###########################################################################
data_read<-function()
{
	if(interactive()) sa<-readline("Enter sub-area size: ")
	return(as.numeric(unlist(strsplit(sa,","))))
}
###########################################################################
momfun<-function(df,x,spp,n,no.m,ratio)
{
	xc<-runif(n,0,round(ratio*max(df[,3]))-ratio*x);
	yc<-runif(n,0,round(max(df[,3]))-x);
	mom<-rep(0,no.m)
	mom1<-mom
	mom2<-mom
	for(i in 1:n)
	{
		g<-f(df,xc[i],yc[i],x,ratio);
		W<-fg(g,spp);
		mom<-mom+r.mom(log2(W[][W[]!=0]),no.m)
		mom1<-mom1+r.mom(round(log2(W[][W[]!=0])),no.m)
		mom2[1:3]<-mom2[1:3]+r.mom(round((W[][W[]!=0])),3)
	}
	mom<-mom/n;
	mom1<-mom1/n;
	mom2<-mom2/n
	return(rbind(mom,mom1,mom2));
}
###########################################################################
datfun<-function(df,x,spp,n,ratio)
{
	xc<-runif(3*n,0,round(ratio*max(df[,3]))-ratio*x);
	yc<-runif(3*n,0,round(max(df[,3]))-x);
	S<-0
	N<-0
	for(i in 1:(3*n))
	{
		g<-f(df,xc[i],yc[i],x,ratio);
		W<-fg(g,spp);
		S<-S+length(intersect(g,g))
		N<-N+sum(W)
	}
	S<-S/n/3
	N<-N/n/3
	return(c(S,N));
}
###########################################################################
r.mom<-function(dat,no.m)
{
	S <- length(dat)
	mom <- rep(0,no.m)
	mom[1] <- 1
	for(i in 1:(no.m-1))
		mom[i+1] <- sum(dat^i)/S
	return(mom)
}
###########################################################################
f<-function(df,xcoor,ycoor,x,ratio) # list of species in an specific area
{
	h<-df[,1][df[,2]>=xcoor & df[,2]<=xcoor+ratio*x & df[,3]>=ycoor & df[,3]<=ycoor+x];
	return(h);
}
###########################################################################
fg<-function(g,sppnames)
{
	noin<-1:length(sppnames)*0;
	for(i in 1:length(sppnames))
		noin[i]<-length(g[][g[]==sppnames[i]]);
	return(noin);
}
###########################################################################
relpwave<-function(l,t,ba)
{
	tmp<-(log(t)-min(log(ba)))/(max(log(ba))-min(log(ba)))*pi
	for(i in 1:length(l))
		if(is.na(l[i])) 
			l[i]<-0
	k<-0
	X<-1:length(l)
	j<-1
	X[2]<-log(t)
	X[3]<-sin(tmp)
	X[4]<-cos(tmp)
	for(i in 1:length(l))
		k<-k+X[i]*l[i];
	return(k)
}
###########################################################################
extpfun<-function(dat,m)
{
	if(length(dat)>=m)
		return(dat[1:m])
	bba<-sqrt(1:length(dat)*10000)
	tmp<-length(bba)
	c<-lm(dat~bba)
	Cp<-as.double(residuals(c))
	kp<-as.double(c$coefficients)
	k1<-kp[1]+kp[2]*(sqrt(1+length(dat)))
	deg<-atan(Cp)
	dat<-c(dat,k1+(deg[tmp]))
	extpfun(dat,m)
}
###########################################################################
Coef.Tcheb.pol.2nd<-function(N,no.m)
{
	s<-matrix(data=0,nrow=no.m,ncol=no.m)
	u1<-matrix(data=0,nrow=2,ncol=1)
	u1[1,1]<-1/sqrt(N);
	s[1,1]<-u1[1,1];
	if(no.m==1) 
		return(s)
	u2<-matrix(data=0,nrow=2,ncol=2)
	u2[2,2]<-1
	u2[1,1]<-(1-N)/sqrt(N*(N^2-1)/3);
	u2[1,2]<-2/sqrt(N*(N^2-1)/3);
	if(no.m==2)
	{
		s[2,]<-u2[1,]
		return(s)
	}
	u3<-expnd(u2,no.m-1)
	s[2,]<-u3[1,]
	for(i in 2:(no.m-1))
	{
		al1<-(2/i)*sqrt((4*i^2-1)/(N^2-i^2));
		al2<-((1-N)/i)*sqrt((4*i^2-1)/(N^2-i^2));
		al3<-((i-1)/i)*sqrt((2*i+1)/(2*i-3))*sqrt((N^2-(i-1)^2)/(N^2-i^2));
		u3<-summ(summ(xmul(cmul(u2,al1)),cmul(u2,al2)),cmul(u1,-1*al3));
		u1<-u2;
		u2<-u3
		u3<-expnd(u3,no.m-1)
		s[i+1,]<-u3[1,];
	}
	return(s);
}
###########################################################################
expnd<-function(u,n)
{
	if(n==0)
	{
		m<-max(u[2,]);
		t<-0:m
		s<-matrix(data=0,nrow=2,ncol=(m+1));
		s[2,]<-t
		s[1,(u[2,]+1)]<-u[1,]
	}
	if(n>0)
	{
		m<-n;
		t<-0:m
		s<-matrix(data=0,nrow=2,ncol=(m+1));
		s[2,]<-t
		s[1,(u[2,]+1)]<-u[1,]
	}
	return(s);
}
###########################################################################
cmul<-function(u,c)
{
	u[1,]<-u[1,]*c;
	return(u)
}
###########################################################################
xmul<-function(u)
{
	u[2,]<-u[2,]+1
	return(u)
}
###########################################################################
summ<-function(u1,u2)
{
	m<-max(max(u1[2,]),max(u2[2,]))
	up1<-expnd(u1,m);
	up2<-expnd(u2,m)
	up<-up1+up2;
	up[2,]<-up1[2,]
	up<-shrt(up);
	return(up)
}
###########################################################################
shrt<-function(u)
{
	a<-which(u[1,]!=0)
	s<-u[,a];
	return(s);
}
###########################################################################
Tcheb.pol.2nd<-function(N,no.m)
{
	t<-matrix(data=0,nrow=no.m,ncol=N);
	t[1,1]<-tn(0,N)
	for(i in 2:(no.m))
		t[i,1]<-tn(i-1,N)
	for(i in 1:no.m)
		t[i,2]<-(1+((i-1)*i)/(1-N))*t[i,1]
	for(x in 2:floor(N/2))
		for(n in 1:no.m)
		{
			g1<-(-1*(n-1)*n-(2*x-1)*(x-N-1)-x)/(x*(N-x))
			g2<-((x-1)*(x-N-1))/(x*(N-x))
			t[n,x+1]<-g1*t[n,x]+g2*t[n,x-1]
		}
	for(x in floor(N/2):(N-1))
		for(n in 1:no.m)
			t[n,x+1]<-(-1)^(n-1)*t[n,N-x]
	return(t)
}
###########################################################################
tn<-function(n,N)
{
	if(n==0)
		return(1/bet.f(n,N))
	return(-1*sqrt((N-n)/(N+n))*sqrt((2*n+1)/(2*n-1))*tn(n-1,N))
}
###########################################################################
bet.f<-function(n,N)
{
	tmp<-N
	if(n!=0)
		for(i in 1:n)
			tmp<-tmp*(N^2-n^2)
	tmp<-tmp/(2*n+1)
	return(sqrt(tmp))
}
###########################################################################
calcalpha <- function(n.orig, s.orig)
{
	bads <- (n.orig<=0 | s.orig<=0 | n.orig<=s.orig)
	n<-n.orig
	s<-s.orig
	n[bads]<-100
	s[bads]<-50
	a <- rep(20,length(n))
	poorest <- rep(T,length(n))
	a<-3;
	while(length(a[poorest])>0)
	{
		a <- a-fo(a,n,s)/fprime(a,n,s)
		poorest <- abs(fo(a,n,s))>1e-10
		a[a<=0|is.nan(a)] <- 1
	}
	a[bads]<-(-1)
	return(a)
}
###########################################################################
fo <- function(a, n, s)
{
	return(a*log(1+(n/a))-s)
}
###########################################################################
fprime <- function(a, n, s)
{
	return(log(1+(n/a))-(n/(a+n)))
}
###########################################################################
moment1<-function(x,n) # moments calculator for Xi
{
	li<-polylog(x,1-n);
	return(-li/log(1-x));
}
###########################################################################
ptransform<-function(x,k)
{
	return(polylog(x,1,"sum",n.sum =2^(k+1)-1)-polylog(x,1,"sum",n.sum = 2^(k)-1));
}
###########################################################################
numsum<-function(n,x,k)
{
	s<-k^n*ptransform(x,k);
	return(s);
}
###########################################################################
LSmom<-function(x,n) # moments calculator for log2(Xi)
{
	tmp<-0;
	y<-7;
	s<-numsum(n,x,y);
	while(-tmp/log(1-x)!=-s/log(1-x)&y<=23) # recursively caculate moments if two consecutive values are equal we will consider it as corresponding moment
	{
		tmp<-s;
		y<-y+1;
		s<-s+numsum(n,x,y);
	}
	return(-s/log(1-x));
}
###########################################################################
BellB <- function(n, k, x)
{
	if(n == 0) return(rep(1, ifelse(missing(k), 1, length(k))))
	if(missing(k))
	{ # complete
		k <- 1:n
		v <- sapply(k, BellB, n=n, x=x)
		sum(v)
	}
	else
	{ # partial
		#' recursion function
		rf <- function(m, l)
		{
			if(l < 1)
			{
				x[m]
			} 
			else
      		{
				i <- l:(m-1)
				sum( choose(m, i) * x[m-i] * sapply(i, rf, l-1))
			}
		}
		K <- k-1
		v <- rf(n, K)
	    v/factorial(k)
	}
}
###########################################################################
st2k<-function(n,k)
{
	s<-factorial(k);
	z<-0;
	for(j in 0:k)
		z<-z+(-1)^(k-j)*choose(k,j)*j^n;
	return(z/s);
}
###########################################################################
M1<-function(ro,c,z,A)
{
	return(ro*A/(c*A^z));
#  return(ro*A/(c+z*log(A)))
}
###########################################################################
mo2nd<-function(res,slope,A)
{
	return(exp(slope*log(A)+res));
}
###########################################################################
betfinder<-function(res,slope,ro,c,z,A)
{
	q<-(mo2nd(res,slope,A)/M1(ro,c,z,A))-1-M1(ro,c,z,A);
	return(1/q);
}
###########################################################################
alffinder<-function(res,slope,ro,c,z,A)
{
	return(M1(ro,c,z,A)*betfinder(res,slope,ro,c,z,A));
}
###########################################################################
mmoment<-function(al,bet,k)
{
	x<-digamma(al)-log(bet);
	for(i in 2:k)
		x<-c(x,psigamma(al,i-1));
	s<-0;
	for(i in 1:k)
	{
		s<-s+BellB(k,i,x[1:(k-i+1)]);
	}
	return(s);
}
###########################################################################
PGmom<-function(al,bet,k)
{
	s<-0;
	for(i in 1:k)
		s<-s+st2k(k,i)*mmoment(al,bet,i);
	s<-s*((log(2))^k);
	return(s);
}
###########################################################################
mufinderLN<-function(res,slope,A)
{
	return(log(2)*exp(slope*log(A)+res));
}
###########################################################################
sigfinderLN<-function(res,slope,ro,c,z,A)
{
	return(sqrt(abs(2*(log(M1(ro,c,z,A))-mufinderLN(res,slope,A)))))
}
###########################################################################
dbfact<-function(i)
{
	s<-1;
	if(i%%2 == 1)
		for(j in 1:((i+1)/2))
			s<-s*(2*j-1);
	if(i%%2==0)
		for(j in 1:((i)/2))
			s<-s*(2*j);
	return(s);
}
###########################################################################
smom<-function(mu,sig,k)
{
	s<-0;
	if(k>=2)
		for(i in seq(2,k,2))
			s<-s+choose(k,i)*(mu^(k-i))*(sig^i)*dbfact(i-1);
	return(s)
}
###########################################################################
cmoment<-function(mu,sig,k)
{
	return(mu^k+smom(mu,sig,k))
}
###########################################################################
LNmom<-function(mu,sig,k)
{
	return(cmoment(mu,sig,k)/(log(2)^k))
}
###########################################################################
mufinderPLN<-function(res,slope,A)
{
	return(exp(slope*log(A)+res));
}
###########################################################################
sigfinderPLN<-function(res1,slope1,res2,slope2,A)
{
	return(mo2nd(res2,slope2,A)-(mufinderPLN(res1,slope1,A))^2-mufinderPLN(res1,slope1,A))
}
###########################################################################
PLNmom<-function(mu,sig,k)
{
	s<-0;
	for(i in 1:k)
		s<-s+st2k(k,i)*cmoment(mu,sig,i);
	return(s);
}
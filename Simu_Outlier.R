#!/u/local/apps/R/3.0.1/gcc-4.4.7/bin/Rscript

set.seed(854920324);

###############
#Calculate Posterior Probability
###############
myfun=function(fold_val,inter_val,a1,b1,a2,b2){
	text='dbeta(x,aa1,bb1)*pbeta((x-inter)/fold,aa2,bb2)';
	fun=as.function(alist(x=,fold=fold_val,inter=inter_val,aa1=a1,bb1=b1,aa2=a2,bb2=b2,eval(parse("",NULL,text))));
	return(fun);
}

#x>a*y+b, a is the fold_val, b is the inter_val
#a1,b1 is the Beta(a1,b1)
myintegrate=function(fold_val,inter_val,a1,b1,a2,b2){
	fun=myfun(fold_val,inter_val,a1,b1,a2,b2);
	#print(a1);print(b1);print(a2);print(b2);
	val=try(integrate(fun,lower=inter_val,upper=1,subdivisions=100),silent=TRUE);
	if (attr(val,'class')=='integrate'){return(val);}
	for (i in seq(200,1000,100)){
		if (attr(val,'class')=='integrate'){return(val);}
		val=try(integrate(fun,inter_val,1,i),silent=TRUE);
	}
	return('error');
}

###############
#END Calculate Posterior Probability
###############

#reformat the input of integration function
myintegrate.reformat=function(data_row,diff_cutoff){
	p1=myintegrate(1,diff_cutoff,data_row[1]+1,data_row[2]+1,data_row[3]+1,data_row[4]+1);
	p2=myintegrate(1,diff_cutoff,data_row[3]+1,data_row[4]+1,data_row[1]+1,data_row[2]+1);
	p=max(0,1-(p1$value+p2$value));return(p);
}


#This function simulates the posterior distribution functions based on the prior parameters
#In this simulation, the prior function is uniform. The likelihood function is Binomial. The posterior funciton is Beta.
#input: n1=exon_inclusion_count1/2+exon_skipping_count1; n2=exon_inclusion_count2/2+exon_skipping_count2;
mysimu=function(n1,n2,diff_cutoff,N=1000){
	theta1=runif(N,min=0,max=1);theta2=runif(N,min=0,max=1);
	myrbinom=function(prob=prob,size=size){rbinom(n=1,size=size,prob=prob)}
	exon_inc1=sapply(theta1,myrbinom,size=n1);
	exon_inc2=sapply(theta2,myrbinom,size=n2);
	data=cbind(exon_inc1+1,n1-exon_inc1+1,exon_inc2+1,n2-exon_inc2+1);
	p=apply(data,1,myintegrate.reformat,diff_cutoff=diff_cutoff)
}

mysimu.data.unif=function(n1,n2){
	N=1;
	theta1=runif(N,min=0,max=1);theta2=runif(N,min=0,max=1);
	exon_inc1=rbinom(n=1,size=n1,prob=theta1);
	exon_inc2=rbinom(n=1,size=n2,prob=theta2);
	data=c(exon_inc1,n1-exon_inc1,exon_inc2,n2-exon_inc2)
}

mysimu.data.unifplusnormal=function(n1,n2,sigma=0,min=0,max=0.1,sigma2=0.1){
	N=1;K=5;
	theta1=runif(N,min=0.001,max=0.999);theta2=-1;
	while((theta2>=1)|(theta2<=0)){theta2=theta1+sample(c(-1,1),1)*runif(N,min,max);}
	theta=sample(c(theta1,theta2));theta1=theta[1];theta2=theta[2];
	#exon_inc1=rbinom(n=1,size=n1,prob=theta1);
	#exon_inc2=rbinom(n=1,size=n2,prob=theta2);
	theta1_list=rnorm(K,log(theta1/(1-theta1)),sigma/(1-min(0.999,theta1))/max(0.001,theta1));
	theta2_list=rnorm(K,log(theta2/(1-theta2)),sigma/(1-min(0.999,theta2))/max(0.001,theta2));
	if (sigma2 > 0){ 
		theta1_list=c(theta1_list[1:(K-1)],log(theta1/(1-theta1))+1.96*sigma2/(1-min(0.999,theta1))/max(0.001,theta1))
		theta2_list=c(theta2_list[1:(K-1)],log(theta2/(1-theta2))-1.96*sigma2/(1-min(0.999,theta2))/max(0.001,theta2))
		#theta2_list=c(theta2_list[1:(K-1)],rnorm(1,log(theta2/(1-theta2)),sigma2/(1-min(0.999,theta2))/max(0.001,theta2)))
	}
	exon_inc1=NULL;exon_inc2=NULL;
	exon_inc1=round(n1*(1/(1+exp(-1*theta1_list))));
	exon_inc2=round(n2*(1/(1+exp(-1*theta2_list))));
	data=c(exon_inc1,n1-exon_inc1,exon_inc2,n2-exon_inc2)
}

mysimu.data.multiunif=function(n1,n2,sigma){
	N=1;
	theta=rmvnorm(N,mean=c(0,0),sigma=matrix(c(1,sigma,sigma,1),byrow=T,nrow=2));
	exon_inc1=rbinom(n=1,size=n1,pnorm(theta[1]));
	exon_inc2=rbinom(n=1,size=n2,pnorm(theta[2]));
	data=c(exon_inc1,n1-exon_inc1,exon_inc2,n2-exon_inc2)
}

mysimu.data.multiunif.FakeNormalize=function(n1,n2,sigma,diff=0,N=1){
	theta=rmvnorm(N,mean=c(0,0),sigma=matrix(c(1,sigma,sigma,1),byrow=T,nrow=2));
	exon_inc1=rbinom(N,size=n1,pnorm(theta[,1]));
	exon_inc2=rbinom(N,size=n2,pnorm(theta[,2]));
	data=cbind(exon_inc1,n1-exon_inc1,exon_inc2,n2-exon_inc2)
	now_diff=abs(data[,1]/(data[,1]+data[,2])-data[,3]/(data[,3]+data[,4]));
	return(mean(now_diff>=diff));
}


#This funciton normalizes the posterior probability
#input: n1=exon_inclusion_count1/2+exon_skipping_count1; n2=exon_inclusion_count2/2+exon_skipping_count2;
#input: p is the observed posterior probability
myNormalize=function(p,n1,n2,diff_cutoff){
	p_prior=mysimu(n1,n2,diff_cutoff);
	return(mean(round(p,5)>round(p_prior,5)));#round(,4) because of precision problem
}

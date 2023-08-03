#!/usr/bin/env R
# author: Tim Barraclough (t.barraclough@imperial.ac.uk)
# version: 1.0 

#################################################################################################
##Four locus model of habitat-based selection and assortative mating
##after Diehl and Bush 1989. The role of habitat preference in adaptation and speciation. In
##D. Otte and J.A. Endler, editors. Speciation and its consequences, pp 345-365. Sinauer, Sunderland.
#################################################################################################

##Used to plot figure 1 in Does selection favour the maintenance of porous species boundaries? 

##Locus a,A and b,B determine fitness in a-niche, i.e habitat type
##and probability of entering mating group 1 and 2
##Locus c,C and d,D determine fitness w.r.t environment, which can be shared or differ between a-niches

output<-NULL

##set up basic parameters
##mutation rate
	mut.vec<-c(0.0001)
##periodicity of environmental change
	period.vec<-c(50)
##selection coefficient on the niche loci
	alpha.selec.vec<-0.1 #seq(0,0.5,0.1)
##selection coefficient on the environment loci
	beta.selec.vec<-0.2 #seq(0,0.5,0.05)
##set up parameter table
	param.table<-expand.grid(mut.vec, period.vec, alpha.selec.vec, beta.selec.vec)
	nreps<-nrow(param.table)
##empty list to set environmental switch points later
	switch.points<-list()
##even periodicity or random shifts with average frequency
##TRUE = Figure 1A in paper, FALSE = a random trial of figure 1B in paper 
	even.periodicity<-TRUE

for (reps in (1:nreps)) {

##parameters
	num.gen<-1000
	mut<-param.table[reps,1]
	ri<-1.0
	selec<-c(param.table[reps,3],param.table[reps,3],param.table[reps,4],param.table[reps,4])		selec.div<-0.01  ##division in assortative mating parameter for evaluating selection gradient
	mod.mut<-2		 ##mutation parameter to scale fitness gradient to modify ri
##does environment fluctuate?
	fluct<-TRUE
	period<-param.table[reps,2]	
##do the a-niches start with same environment?
	same<-FALSE
##Random settlement as in Felsenstein or within niche as in Diehl and Bush
	random.settlement<-FALSE
##environment changes at these time points
##a) even periodicity
if (even.periodicity) {
	#switch.points[[1]]<-seq(100,num.gen,period)
	#switch.points[[2]]<-seq(100,num.gen,period) } else {	
##b) random with average period
	switch.points[[1]]<-sort(sample(100:num.gen, round((num.gen-99)*1/period,0),replace=F))
	switch.points[[2]]<-sort(sample(100:num.gen, round((num.gen-99)*1/period,0),replace=F))
	}
	
##set up genotype identities
	geno.combos<-expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
	geno.ids<-apply(expand.grid(c("a","A"),c("b","B"),c("c","C"),c("d","D")),1,function(x) paste(x,collapse=""))

##set up mutation matrix
	mut.list<-list()
	for (i in (1:4)) {
		mut.list[[i]]<-outer(substring(geno.ids,i,i),substring(geno.ids,i,i),"==")
		mut.list[[i]][mut.list[[i]]]<-(1-mut)
		mut.list[[i]][!mut.list[[i]]]<-mut
		}
	mut.matrix<-mut.list[[1]]*mut.list[[2]]*mut.list[[3]]*mut.list[[4]]

##Calculate random mating table
	alleles<-list()
	for (i in (1:4)) {
	alleles[[i]]<-strsplit(outer(substring(geno.ids,i,i),substring(geno.ids,i,i),paste)," ")}

	count.genotypes<-function(x) {
		genos<-outer(alleles[[1]][[x]],alleles[[2]][[x]],paste)
			genos<-gsub(" ","",as.vector(outer(genos,alleles[[3]][[x]],paste)))
			genos<-gsub(" ","",as.vector(outer(genos,alleles[[4]][[x]],paste)))
		return(geno.ids%in%genos/sum(geno.ids%in%genos))}

	RM.table<-t(simplify2array(lapply(1:length(alleles[[1]]),count.genotypes)))
	colnames(RM.table)<-geno.ids
	
allele.array<-array(unlist(alleles),dim=c(2,length(alleles[[1]]),length(alleles)))
	rownames(RM.table)<-paste(apply(allele.array[1,,],1,function(x) paste(x,collapse="")),".",
							  apply(allele.array[2,,],1,function(x) paste(x,collapse="")),sep="")

##starting environments determined here
if (same) {
	opt.1<-c(0,0,0,0)
	opt.2<-c(1,1,0,0)} else {
	opt.1<-c(0,0,0,0)
	opt.2<-c(1,1,1,1)}

##population genotypes with 50:50 optimal genotypes to start
	geno.freqs<-array(0,length(geno.ids)) 
	geno.freqs[apply(geno.combos,1,function(x) identical(as.vector(x),opt.1))]<-0.5
	geno.freqs[apply(geno.combos,1,function(x) identical(as.vector(x),opt.2))]<-0.5
	names(geno.freqs)<-geno.ids

##Calculate fitness at start
	optim1<-sweep(geno.combos,2,opt.1,"==")
	fitness1<-matrix(1,nrow=nrow(optim1),ncol=ncol(optim1))
	fitness1[optim1]<-fitness1[optim1]+sweep(optim1,2,selec,"*")[optim1]
	fitness1<-apply(fitness1,1,prod)

	optim2<-sweep(geno.combos,2,opt.2,"==")
	fitness2<-matrix(1,nrow=nrow(optim2),ncol=ncol(optim2))
	fitness2[optim2]<-fitness2[optim2]+sweep(optim2,2,selec,"*")[optim2]
	fitness2<-apply(fitness2,1,prod)

##and genotype frequencies in the two niches
	geno.freqs1<-geno.freqs
	geno.freqs2<-geno.freqs
	
##store initial results
	results<-c(geno.freqs,recomb=ri,fitness=prod(1+selec))

##START GENERATIONS
for (i in (1:num.gen)) {
		
	if ((fluct)&(i%in%switch.points[[1]])) {
	opt.1[3:4]<-!opt.1[3:4]}
	if ((fluct)&(i%in%switch.points[[2]])) {
	opt.2[3:4]<-!opt.2[3:4]}
	
	##mutation
	geno.freqs<-geno.freqs%*%mut.matrix

	check.gradient<-list()
	
	for (j in (1:3)) {
	
	if (j==1) ri.tmp<-ri
	if (j==2) {ri.tmp<-ri-selec.div}
	if (j==3) {ri.tmp<-ri+selec.div}
	
	if (ri.tmp>1)	ri.tmp<-1
	if (ri.tmp<0) 	ri.tmp<-0

	##Calculate probabilities of joining mating pools under assortative mating
	#ri<-x
	group1<-which(rowSums(geno.combos[,1:2])==0)
	group2<-which(rowSums(geno.combos[,1:2])==2)
	p.group<-matrix(0.5,ncol=length(geno.ids),nrow=2)
	p.group[1,group1]<-(1+ri.tmp)/2
	p.group[1,group2]<-(1-ri.tmp)/2
	p.group[2,group1]<-(1-ri.tmp)/2
	p.group[2,group2]<-(1+ri.tmp)/2
	colnames(p.group)<-geno.ids
	rownames(p.group)<-c("mating1","mating2")

	##mating pool freqs
	mating.pools<-sweep(p.group,2,geno.freqs,"*")
	
	##random mating in group 1
	mating.freqs1<-(mating.pools[1,]%o%mating.pools[1,])
	mating.freqs1<-mating.freqs1/sum(mating.freqs1)
	geno.freqs1<-colSums(RM.table*as.vector(mating.freqs1))

	##random mating in group 2
	mating.freqs2<-(mating.pools[2,]%o%mating.pools[2,])
	mating.freqs2<-mating.freqs2/sum(mating.freqs2)
	geno.freqs2<-colSums(RM.table*as.vector(mating.freqs2))

	##Random settlement or stay in their own niche
	if (random.settlement) {
	geno.freqs.tmp<-(geno.freqs1+geno.freqs2)/2
	geno.freqs1<-geno.freqs.tmp
	geno.freqs2<-geno.freqs.tmp }

	##calculate fitness, selection operates separately in each group
		optim1<-sweep(geno.combos,2,opt.1,"==")
		fitness1<-matrix(1,nrow=nrow(optim1),ncol=ncol(optim1))
		fitness1[optim1]<-fitness1[optim1]+sweep(optim1,2,selec,"*")[optim1]
		fitness1<-apply(fitness1,1,prod)

		optim2<-sweep(geno.combos,2,opt.2,"==")
		fitness2<-matrix(1,nrow=nrow(optim2),ncol=ncol(optim2))
		fitness2[optim2]<-fitness2[optim2]+sweep(optim2,2,selec,"*")[optim2]
		fitness2<-apply(fitness2,1,prod)

	mean.fitness<-mean(c(weighted.mean(fitness1,geno.freqs1),weighted.mean(fitness2,geno.freqs2)))

	check.gradient[[j]]<-list(fitness1,fitness2,mean.fitness,ri.tmp)
	
	}
	
			fitness1<-check.gradient[[1]][[1]] 
			fitness2<-check.gradient[[1]][[2]] 
			mean.fitness<-check.gradient[[1]][[3]]  

	##store gen results
	results<-rbind(results,c(geno.freqs,recomb=ri,fitness=mean(c(weighted.mean(fitness1,geno.freqs1),weighted.mean(fitness2,geno.freqs2)))))

	##recombination param for next iteration
	best<-which.max(c(check.gradient[[1]][[3]],check.gradient[[2]][[3]],check.gradient[[3]][[3]]))
	fitness.grad<-check.gradient[[best]][[3]]-check.gradient[[1]][[3]]
	ri<-ri+sign(check.gradient[[best]][[4]]-ri)*fitness.grad*mod.mut
	if (ri>1) ri<-1
	if (ri<0) ri<-0

	##apply selection
	geno.freqs1<-geno.freqs1*fitness1
	geno.freqs1<-geno.freqs1/sum(geno.freqs1)

	geno.freqs2<-geno.freqs2*fitness2
	geno.freqs2<-geno.freqs2/sum(geno.freqs2)
	
	##average to get overall frequences
	geno.freqs<-(geno.freqs1+geno.freqs2)/2
		
}	##end of gen loop


##output results
	#names(selec)<-paste("selec",c("a","b","c","d"),sep=".")
gm_mean = function(a){prod(a)^(1/length(a))}
output<-rbind(output,t(as.matrix(c(reps=reps,num.gen=num.gen,mut=mut,selec.alpha=selec[1],selec.beta=selec[3],fluct=fluct,period=period,same=same,mean.ri=mean(results[501:1000,ncol(results)-1]),max.ri.change=max(diff(results[,ncol(results)-1])),mean.fit=mean(results[501:1000,ncol(results)-1]),geo.mean.fit=gm_mean(results[501:1000,ncol(results)]),colMeans(results[501:1000,1:16])))))

}  ##end of reps loop

#matplot(results,type="l",lty=1,xlab="Time (generations)",ylab="Frequency of genotype",ylim=c(0,1))


par(xpd = NA, mar = c(10, 5, 8, 2) + 0.1)

plot(results[,ncol(results)-1],type="l",lty=1,xlab="Time (generations)",ylab="P(Assortative Mating)",ylim=c(0.97,1),las=1)

lines(c(0,1000),rep(mean(results[,ncol(results)-1]),2),lty=2)

# Note that the rectangle we make here has corner coordinates outside of
# our plotting device
plot.rect.1<-unique(c(0,switch.points[[1]],1000))
plot.rect.2<-unique(c(0,switch.points[[2]],1000))

rect(plot.rect.1[-length(plot.rect.1)], 1.005, plot.rect.1[-1], 1.0075, col=c("white","grey"))
rect(plot.rect.2[-length(plot.rect.2)], 1.0075, plot.rect.2[-1], 1.01, col=c("grey","white"))

text(x=-70,y=1.006,"patch 1")
text(x=-70,y=1.0085,"patch 2")


#!/usr/bin/env R
# author: Tim Barraclough (tim.barraclough@biology.ox.ac.uk)
# version: 1.0 

#################################################################################################
##Four locus model of habitat-based selection and assortative mating
##after Diehl and Bush 1989. The role of habitat preference in adaptation and speciation. In
##D. Otte and J.A. Endler, editors. Speciation and its consequences, pp 345-365. Sinauer, Sunderland.
#################################################################################################

##Used to plot figure 1 and 2 in Does selection favour the maintenance of porous species boundaries? Barraclough T.G. (2024) J. Evol. Biol.

##To plot figure 1 versus figure 2, change line 38 even.periodicity<-TRUE/FALSE

##Locus a,A and b,B determine fitness in a-niche, i.e habitat type
##and probability of entering mating group 1 and 2
##Locus c,C and d,D determine fitness w.r.t environment, which can be shared or differ between a-niches
##Locus E and e determine the probability that ab and AB genotypes enter the early or late mating pools. E means ab enters early mating pool with probability=1, AB to late mating pool with p=1
##e means ab and AB enter the mating pools at random


##set up basic parameters
##mutation rate
	mut.vec<-c(0.0001)
##periodicity of environmental change
	period.vec<-c(100)
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
##TRUE = Figure 1 in paper, FALSE = a random trial of figure 2 in paper 
	even.periodicity<-FALSE

for (reps in (1:nreps)) {

##parameters
	num.gen<-1000
	mut<-param.table[reps,1]
	ri<-1.0
	selec<-c(param.table[reps,3],param.table[reps,3],param.table[reps,4],param.table[reps,4])		
	selec.div<-0.01  ##division in assortative mating parameter for evaluating selection gradient
	mod.mut<-2		 ##mutation parameter to scale fitness gradient to modify ri
##does environment fluctuate?
	fluct<-TRUE
	period<-param.table[reps,2]	
##do the a-niches start with same environment?
	same<-FALSE
##Random settlement as in Felsenstein (TRUE) or within niche as in Diehl and Bush (FALSE)
	random.settlement<-FALSE
##environment changes at these time points
##a) even periodicity
if (even.periodicity) {
	switch.points[[1]]<-seq(200,num.gen,period)
	switch.points[[2]]<-seq(200,num.gen,period) } else {	
##b) random with average period
	switch.points[[1]]<-sort(sample(100:num.gen, round((num.gen-99)*1/period,0),replace=F))
	switch.points[[2]]<-sort(sample(100:num.gen, round((num.gen-99)*1/period,0),replace=F))
	}
	
##set up genotype identities
	geno.combos<-expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
	geno.ids<-apply(expand.grid(c("a","A"),c("b","B"),c("c","C"),c("d","D"),c("e","E")),1,function(x) paste(x,collapse=""))

##set up mutation matrix
	mut.list<-list()
	for (i in (1:5)) {
		mut.list[[i]]<-outer(substring(geno.ids,i,i),substring(geno.ids,i,i),"==")
		mut.list[[i]][mut.list[[i]]]<-(1-mut)
		mut.list[[i]][!mut.list[[i]]]<-mut
		}
	mut.matrix<-mut.list[[1]]*mut.list[[2]]*mut.list[[3]]*mut.list[[4]]*mut.list[[5]]

##Calculate random mating table
	alleles<-list()
	for (i in (1:5)) {
	alleles[[i]]<-strsplit(outer(substring(geno.ids,i,i),substring(geno.ids,i,i),paste)," ")}

	count.genotypes<-function(x) {
		genos<-outer(alleles[[1]][[x]],alleles[[2]][[x]],paste)
			genos<-gsub(" ","",as.vector(outer(genos,alleles[[3]][[x]],paste)))
			genos<-gsub(" ","",as.vector(outer(genos,alleles[[4]][[x]],paste)))
			genos<-gsub(" ","",as.vector(outer(genos,alleles[[5]][[x]],paste)))
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

##population genotypes with 50:50 optimal genotypes to start, but 100% E
	geno.freqs<-array(0,length(geno.ids)) 
	geno.freqs[apply(geno.combos,1,function(x) identical(as.vector(x[1:4]),opt.1))]<-c(0,0.5)
	geno.freqs[apply(geno.combos,1,function(x) identical(as.vector(x[1:4]),opt.2))]<-c(0,0.5)
	names(geno.freqs)<-geno.ids

##Calculate fitness at start
	optim1<-sweep(geno.combos[,1:4],2,opt.1,"==")
	fitness1<-matrix(1,nrow=nrow(optim1),ncol=ncol(optim1))
	fitness1[optim1]<-fitness1[optim1]+sweep(optim1,2,selec,"*")[optim1]
	fitness1<-apply(fitness1,1,prod)

	optim2<-sweep(geno.combos[,1:4],2,opt.2,"==")
	fitness2<-matrix(1,nrow=nrow(optim2),ncol=ncol(optim2))
	fitness2[optim2]<-fitness2[optim2]+sweep(optim2,2,selec,"*")[optim2]
	fitness2<-apply(fitness2,1,prod)

##and genotype frequencies in the two niches
	geno.freqs1<-geno.freqs
	geno.freqs2<-geno.freqs
	
##store initial results
	results<-c(geno.freqs,f_e=(sum(geno.freqs1[grep("e",names(geno.freqs1))])+sum(geno.freqs2[grep("e",names(geno.freqs2))]))/2,recomb=ri,fitness=prod(1+selec))

##START GENERATIONS
for (i in (1:num.gen)) {
	
	##apply a shift in environment if required
	
	if ((fluct)&(i%in%switch.points[[1]])) {
	opt.1[3:4]<-!opt.1[3:4]}
	if ((fluct)&(i%in%switch.points[[2]])) {
	opt.2[3:4]<-!opt.2[3:4]}
	
	##apply mutation
	geno.freqs<-geno.freqs%*%mut.matrix		
	
	##Calculate probabilities of joining mating pools under assortative mating
	##Only those that are ab and E join mating group 1, and AB and E join mating group 2
	group1<-which((rowSums(geno.combos[,1:2])==0)&(geno.combos[,5]==1))
	group2<-which((rowSums(geno.combos[,1:2])==2)&(geno.combos[,5]==1))
	p.group<-matrix(0.5,ncol=length(geno.ids),nrow=2)
	p.group[1,group1]<-1
	p.group[1,group2]<-0
	p.group[2,group1]<-0
	p.group[2,group2]<-1
	colnames(p.group)<-geno.ids
	rownames(p.group)<-c("mating1","mating2")

	##mating pool freqs
	mating.pools<-sweep(p.group,2,geno.freqs,"*")
	
	##random mating within mating pool 1
	mating.freqs1<-(mating.pools[1,]%o%mating.pools[1,])
	mating.freqs1<-mating.freqs1/sum(mating.freqs1)
	geno.freqs1<-colSums(RM.table*as.vector(mating.freqs1))

	##random mating in mating pool 2
	mating.freqs2<-(mating.pools[2,]%o%mating.pools[2,])
	mating.freqs2<-mating.freqs2/sum(mating.freqs2)
	geno.freqs2<-colSums(RM.table*as.vector(mating.freqs2))

	##calculate the frequency of assortative mating, i.e. abXab and ABxAB
	recomb.tmp<-(sum(mating.freqs1[grep("ab",rownames(mating.freqs1)),grep("ab",rownames(mating.freqs1))])+sum(mating.freqs2[grep("AB",rownames(mating.freqs2)),grep("AB",rownames(mating.freqs2))]))/2

	##Random settlement or stay in their own niche
	if (random.settlement) {
	geno.freqs.tmp<-(geno.freqs1+geno.freqs2)/2
	geno.freqs1<-geno.freqs.tmp
	geno.freqs2<-geno.freqs.tmp }

	##calculate fitness, selection operates separately in each group
		optim1<-sweep(geno.combos[,1:4],2,opt.1,"==")
		fitness1<-matrix(1,nrow=nrow(optim1),ncol=ncol(optim1))
		fitness1[optim1]<-fitness1[optim1]+sweep(optim1,2,selec,"*")[optim1]
		fitness1<-apply(fitness1,1,prod)

		optim2<-sweep(geno.combos[,1:4],2,opt.2,"==")
		fitness2<-matrix(1,nrow=nrow(optim2),ncol=ncol(optim2))
		fitness2[optim2]<-fitness2[optim2]+sweep(optim2,2,selec,"*")[optim2]
		fitness2<-apply(fitness2,1,prod)

	mean.fitness<-mean(c(weighted.mean(fitness1,geno.freqs1),weighted.mean(fitness2,geno.freqs2)))
		
	##store gen results
	results<-rbind(results,c(geno.freqs,f_e=(sum(geno.freqs1[grep("e",names(geno.freqs1))])+sum(geno.freqs2[grep("e",names(geno.freqs2))]))/2,recomb=recomb.tmp,fitness=mean(c(weighted.mean(fitness1,geno.freqs1),weighted.mean(fitness2,geno.freqs2)))))

	##apply selection
	geno.freqs1<-geno.freqs1*fitness1
	geno.freqs1<-geno.freqs1/sum(geno.freqs1)

	geno.freqs2<-geno.freqs2*fitness2
	geno.freqs2<-geno.freqs2/sum(geno.freqs2)
	
	##average to get overall frequences
	geno.freqs<-(geno.freqs1+geno.freqs2)/2
		
}	##end of gen loop

}  ##end of reps loop

##set up plotting layout

layout.matrix=matrix(c(1,2,1,2),nrow=2)

layout(mat = layout.matrix,
       heights = c(10, 4), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns

par(xpd = NA, mar = c(1, 5, 8, 2) + 0.1)

##plot the frequency of assortative mating

plot(results[,ncol(results)-1],type="l",lty=1,xlab="",ylab="Frequency of Assortative Mating",ylim=c(0.96,1),las=1,xaxt="n")
#plot(1-results[,ncol(results)-1],type="l",lty=1,xlab="Time (generations)",ylab="P(Assortative Mating)",ylim=c(0,0.01),las=1)

##add a dashed line showing the average
lines(c(0,1000),rep(mean(results[,ncol(results)-1]),2),lty=2)

# This plots the changes in physical environment in rectangles at the top
plot.rect.1<-unique(c(0,switch.points[[1]],1000))
plot.rect.2<-unique(c(0,switch.points[[2]],1000))
rect(plot.rect.1[-length(plot.rect.1)], 1.005, plot.rect.1[-1], 1.01, col=c("white","grey"))
rect(plot.rect.2[-length(plot.rect.2)], 1.01, plot.rect.2[-1], 1.015, col=c("grey","white"))
text(x=-70,y=1.0075,"patch 2")
text(x=-70,y=1.0125,"patch 1")

##plot the frequency of e alleles, i.e. those leading to random joining of mating pools for ab and AB
par(xpd = NA, mar = c(4,5,0,2) + 0.1)
plot(results[,ncol(results)-2]*1000,type="l",ylab="F(e) x 10-3",xlab="Time (generations)",ylim=c(0,max(results[,ncol(results)-2]*1000)),las=1)
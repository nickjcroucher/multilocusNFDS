# load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(rdist)
library(ape)
library(ggtree)
library(phytools)
library(cowplot)

# read matrices
read_matrices<-function(prefix, r, real = NULL, snp.calc.proc=FALSE) {
  # read data
  in.matrices<-list()
  snp.matrices<-list()
  for (i in 1:r) {
    fn<-paste(prefix,i,"finalPopGenotypes.tab",sep = '.')
    in.matrices[[i]]<-read.table(file=fn,header=T,comment.char = "")
    in.matrices[[i]]<-in.matrices[[i]][,-1]
    if (snp.calc.proc==TRUE) {
      snp_fn<-paste(prefix,i,"markers.finalPopGenotypes.tab",sep = '.')
      snp.matrices[[i]]<-read.table(file=snp_fn,header=T,comment.char = "")
      snp.matrices[[i]]<-snp.matrices[[i]][,-1]
    }
  }
  # process inputs - gene frequencies
  gene.freqs<-lapply(in.matrices,colMeans)
  gene.stats<-do.call(cbind, gene.freqs)
  # process inputs - genome content
  genome.sizes<-lapply(in.matrices,rowMeans)
  genome.stats<-do.call(cbind, genome.sizes)
  # process inputs - distances between isolates
  genome.dists<-lapply(in.matrices,rdist,metric="jaccard")
  dist.stats<-do.call(cbind, genome.dists)
  # process inputs - SNP frequencies
  snp.stats<-matrix()
  snp.dist.example<-NULL
  if (snp.calc.proc==TRUE) {
    snp.freqs<-lapply(snp.matrices,colMeans)
    snp.stats<-do.call(cbind, snp.freqs)
    snp.dists<-lapply(snp.matrices,rdist,metric="hamming")
    snp.dist.example<-snp.dists[[1]]
  }
  snp.dist.stats<-do.call(cbind, snp.dists)
  # compare distances to input data
  comp.dists<-lapply(in.matrices,get_comparative_dists,real[,6:ncol(real)])
  comp.dist.stats<-do.call(cbind, comp.dists)
  # output matrices
  output.list<-list(gene.stats,genome.stats,dist.stats,snp.stats,snp.dist.stats,snp.dist.example,comp.dist.stats)
  return(output.list)
}

get_comparative_dists<-function(i,r) {
  mod.r<-r[,match(colnames(i),colnames(r))]
  d<-cdist(mod.r,i,metric="jaccard") # rows are real data, columns are per query
  min.d<-apply(d,2,min)
  return(min.d)
}

# convert data to tidyverse format
convert_to_df<-function(outlist,el=NULL,snp.info=NULL,rec=NULL,select=NULL,summary=NULL,snp.calc.proc=FALSE) {
  # process gene stats
  message("Processing gene statistics\n")
  gene.data<-matrix(nrow=0,ncol=2)
  for (c in 1:ncol(outlist[[1]])) {
    gene.data<-rbind(gene.data,cbind(c,outlist[[1]][,c]))
  }
  ordered.el<-el[match(rownames(gene.data),el[,1]),2]
  gene.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Locus"=rownames(gene.data),"Replicate"=gene.data[,1],"Frequency"=gene.data[,2],"Equilibrium"=ordered.el)
  # parse genome sizes
  message("Parsing genome sizes\n")
  genome.data<-matrix(nrow=0,ncol=2)
  for (c in 1:ncol(outlist[[2]])) {
    genome.data<-rbind(genome.data,cbind(c,outlist[[2]][,c]))
  }
  genome.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(outlist[[2]]),"Replicate"=genome.data[,1],"Frequency"=genome.data[,2])
  # process pairwise distances
  message("Processing pairwise distances\n")
  dists.data<-matrix(nrow=0,ncol=2)
  for (c in 1:ncol(outlist[[3]])) {
    dists.data<-rbind(dists.data,cbind(c,outlist[[3]][,c]))
  }
  distance.sampling<-sample(nrow(dists.data),size=round(nrow(dists.data)/100),replace=FALSE)
  dists.data<-dists.data[distance.sampling,] # downsample for memory efficiency
  message("Combining pairwise distances\n")
  #dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(outlist[[3]]),"Replicate"=dists.data[,1],"Distances"=dists.data[,2])
  dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(dists.data),"Replicate"=dists.data[,1],"Distances"=dists.data[,2])
  # process comparisons with starting data
  message("Processing comparisons with starting data\n")
  comp.dists.data<-matrix(nrow=0,ncol=2)
  for (c in 1:ncol(outlist[[7]])) {
    comp.dists.data<-rbind(comp.dists.data,cbind(c,outlist[[7]][,c]))
  }
  comp.distance.sampling<-sample(nrow(comp.dists.data),size=round(nrow(comp.dists.data)/100),replace=FALSE)
  comp.dists.data<-comp.dists.data[comp.distance.sampling,] # downsample for memory efficiency
  comp.dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(comp.dists.data),"Replicate"=comp.dists.data[,1],"Distances"=comp.dists.data[,2])
  
  snp.df<-data.frame()
  snp.dists.df<-data.frame()
  
  if (snp.calc.proc==TRUE) {
    message("Processing SNPs\n")
    # process input SNPs
    snp.pos<-snp.info[,1]
    snp.freq<-snp.info[,2]
    # process gene stats
    snp.data<-matrix(nrow=0,ncol=2)
    for (c in 1:ncol(outlist[[4]])) {
      snp.data<-rbind(snp.data,cbind(c,outlist[[4]][,c]))
    }
    ordered.snp.freq<-snp.freq[match(rownames(snp.data),snp.pos)]
    #message(paste("SNP pos: ",head(snp.pos,n=1), " SNP data: ",head(snp.data,n=1)," SNP info: ",head(snp.info,n=1)," ordered SNPs: ",head(ordered.snp.freq,n=1),"\n"))
    snp.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Locus"=rownames(snp.data),"Replicate"=snp.data[,1],"Frequency"=snp.data[,2],"Equilibrium"=ordered.snp.freq)
    # SNP distance matrix
    snp.dists.data<-matrix(nrow=0,ncol=2)
    for (c in 1:ncol(outlist[[5]])) {
      snp.dists.data<-rbind(snp.dists.data,cbind(c,outlist[[5]][,c]))
    }
    snp.dists.data<-snp.dists.data[distance.sampling,] # downsample for memory efficiency - same rows as above
    message("Combining SNP pairwise distances\n")
    #dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(outlist[[3]]),"Replicate"=dists.data[,1],"Distances"=dists.data[,2])
    snp.dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(snp.dists.data),"Replicate"=snp.dists.data[,1],"Distances"=snp.dists.data[,2])
    
    # debug
    #snp.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Locus"=rownames(snp.data),"Replicate"=snp.data[,1],"Frequency"=snp.data[,2],"Equilibrium"=snp[match(rownames(snp.data),snp[,1]),2])
  }

  # complete
  return(list(gene.df,genome.df,dists.df,snp.df,snp.dists.df,comp.dists.df))
}

relevel_df<-function(df) {
  r<-match("Recombination",colnames(df))
  s<-match("Selection",colnames(df))
  df[,r]<-factor(df[,r],levels = c("No transformation","Symmetrical\ntransformation","Asymmetrical\ntransformation"))
  df[,s]<-factor(df[,s],levels = c("No selection","Multi-locus NFDS"))
  return(df)
}

process_simulation_data<-function(prefixes=NULL,recs=NULL,selects=NULL,summaries=NULL,eq.loci=NULL,real.data=NULL,snp.loci=NULL,snp.calc=FALSE) {
  # store plots
  snp.calc.val<-snp.calc
  out_plots<-list()
  locus.data<-list()
  genotype.data<-list()
  dist.data<-list()
  snp.data<-list()
  snp.dist.data<-list()
  tree.data<-list()
  tree.df<-data.frame()
  comp.dist.data<-list()
  
  # processing
  message("Initial processing of data frames\n")
  for (x in 1:length(prefixes)) {
    message(paste0("Processing ",prefixes[x],"\n"))
    processed.summaries<-read_matrices(prefixes[x],nrep,real=real.data,snp.calc.proc=snp.calc.val)
    message(paste0("Matrices read\n"))
    out.list<-convert_to_df(processed.summaries,el=eq.loci,snp.info=snp.loci,rec=recs[x],select=selects[x],summary=summaries[x],snp.calc.proc=snp.calc.val)
    if (snp.calc) {
      message(paste0("Calculating trees\n"))
      tree.data[[x]]<-make_nj_tree(processed.summaries[[6]],col="black")
    }
    message(paste0("Converted to data frame\n"))
    processed.summaries<-NULL # reduce memory
    locus.data[[x]]<-out.list[[1]]
    genotype.data[[x]]<-out.list[[2]]
    dist.data[[x]]<-out.list[[3]]
    if (snp.calc.val==TRUE) {
      snp.data[[x]]<-out.list[[4]]
      snp.dist.data[[x]]<-out.list[[5]]
    }
    comp.dist.data[[x]]<-out.list[[6]]
    out.list<-NULL # reduce memory
  }
  
  # merge data
  message("Merging data frames\n")
  raw.locus.df<-bind_rows(locus.data)
  mean.locus.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.locus.df,mean)
  lb.locus.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.locus.df,min)
  ub.locus.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.locus.df,max)
  locus.df<-cbind(mean.locus.df,lb.locus.df$Frequency,ub.locus.df$Frequency)
  mean.locus.df<-NULL # reduce memory
  ub.locus.df<-NULL # reduce memory
  lb.locus.df<-NULL # reduce memory
  locus.df<-relevel_df(locus.df)
  colnames(locus.df)[colnames(locus.df)=="lb.locus.df$Frequency"]<-"LowerBound"
  colnames(locus.df)[colnames(locus.df)=="ub.locus.df$Frequency"]<-"UpperBound"
  # genotype data
  genotype.df<-bind_rows(genotype.data)
  genotype.df<-relevel_df(genotype.df)
  # distance data
  dist.df<-bind_rows(dist.data)
  dist.df<-relevel_df(dist.df)
  # comparative distance data
  comp.dist.df<-bind_rows(comp.dist.data)
  comp.dist.df<-relevel_df(comp.dist.df)
  
  # SNP data
  if (snp.calc==TRUE) {
    # SNP frequencies
    message("Merging SNP data frames\n")
    raw.snp.df<-bind_rows(snp.data)
    mean.snp.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,mean)
    lb.snp.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,min)
    ub.snp.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,max)
    snp.df<-cbind(mean.snp.df,lb.snp.df$Frequency,ub.snp.df$Frequency)
    mean.snp.df<-NULL # reduce memory
    ub.snp.df<-NULL # reduce memory
    lb.snp.df<-NULL # reduce memory
    snp.df<-relevel_df(snp.df)
    colnames(snp.df)[colnames(snp.df)=="lb.snp.df$Frequency"]<-"LowerBound"
    colnames(snp.df)[colnames(snp.df)=="ub.snp.df$Frequency"]<-"UpperBound"
    # SNP distances
    snp.dist.df<-bind_rows(snp.dist.data)
    snp.dist.df<-relevel_df(snp.dist.df)
    # Tree plots
    tree.df<-data.frame("Recombination" = recs, "Selection" = selects, "Tree" = 1:length(tree.data))
  }
  
  # plot gene frequencies
  message("Plotting data\n")
  freq_plot<-ggplot(data=locus.df,aes(x=Equilibrium,y=Frequency))+
    geom_point()+
    geom_errorbar(aes(ymin=LowerBound,ymax=UpperBound),alpha=0.5)+
    xlab("Equilibrium frequency")+
    ylab("Simulated frequency at final timepoint")+
    facet_grid(Recombination~Selection)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red")
  out_plots[[1]]<-freq_plot

  # plot genome sizes
  genome_plot<-ggplot(data=genotype.df,aes(x=Frequency))+
    geom_histogram(aes(y=..density..),bins=50)+
    facet_grid(Recombination~Selection)+
    xlab("Proportion of loci in genome")+
    geom_density(data=genome_sizes, aes(x=Sizes), col = "red", trim=TRUE)
  out_plots[[2]]<-genome_plot

  # plot pairwise distances
  distance_plot<-ggplot(data=dist.df,aes(x=Distances))+
    geom_histogram(aes(y=..ndensity..),bins=50)+
    facet_grid(Recombination~Selection)+
    xlab("Pairwise binary Jaccard distances between genomes")+
    ylab("Density (log(density+1) scale)")+
    geom_density(data=genome_distances, aes(x=Distances, y=..scaled..), col = "red", trim=TRUE)+
    scale_y_continuous(trans="log1p")
  out_plots[[3]]<-distance_plot

  # plot distances to real data
  comp_dist_plot<-ggplot(data=comp.dist.df,aes(x=Distances))+
    geom_histogram(aes(y=..ndensity..),bins=50)+
    facet_grid(Recombination~Selection)+
    xlab("Minimum pairwise binary Jaccard distances between simulated and initial data")+
    ylab("Density (log(density+1) scale)")+
#    geom_density(data=genome_distances, aes(x=Distances, y=..scaled..), col = "red", trim=TRUE)+
    scale_y_continuous(trans="log1p")
  out_plots[[4]]<-comp_dist_plot
  
  # plot SNP data
  if (snp.calc==TRUE) {
    # plot SNP frequencies
    snp_plot<-ggplot(data=snp.df,aes(x=Equilibrium,y=Frequency))+
      geom_point()+
      geom_errorbar(aes(ymin=LowerBound,ymax=UpperBound),alpha=0.5)+
      xlab("Frequency in genomic data")+
      ylab("Simulated frequency at final timepoint")+
      facet_grid(Recombination~Selection)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red")
    out_plots[[5]]<-snp_plot
    # plot SNP distances
    snp_distance_plot<-ggplot(data=snp.dist.df,aes(x=Distances))+
      geom_histogram(aes(y=..ndensity..),bins=50)+
      facet_grid(Recombination~Selection)+
      xlab("Pairwise binary Jaccard SNP distances between genomes")+
      ylab("Density (log(density+1) scale)")+
      geom_density(data=genome_distances, aes(x=SNPDistances, y=..scaled..), col = "red", trim=TRUE)+
      scale_y_continuous(trans="log1p")
    out_plots[[6]]<-snp_distance_plot
    
    # plot SNP v accessory distances
    snp.dist.df$Accessory<-dist.df$Distances
    comparison_distance_plot<-ggplot(data=snp.dist.df,aes(x=Distances,y=Accessory))+
      geom_point(alpha = 0.05)+
      facet_grid(Recombination~Selection)+
      xlab("Pairwise binary Jaccard SNP distances between genomes")+
      ylab("Pairwise binary Jaccard distances between genomes")+
      geom_density_2d(data=genome_distances, aes(x=SNPDistances, y=Distances), colour = "red", alpha = 0.5)
    out_plots[[7]]<-comparison_distance_plot
    
    # plot trees
    tree.plot.list<-list()
    for (r in levels(snp.dist.df$Recombination)) {
      for (s in levels(snp.dist.df$Selection)) {
        tree.plot.list[[length(tree.plot.list)+1]]<-tree.data[[tree.df$Tree[tree.df$Recombination==r & tree.df$Selection==s]]]
      }
    }
    tree_comparison_plot<-cowplot::plot_grid(plotlist = tree.plot.list,ncol = 2, nrow = 3, align = "none")
    out_plots[[8]]<-tree_comparison_plot
  }
  
  # return final plots
  return(out_plots)
}

make_nj_tree<-function(s,col="black") {
  t<-nj(rdist(s,metric="hamming"))
  t<-midpoint.root(t)
  ggtree(t, layout="circular")
}

#########
# BEGIN #
#########

# data structures
recs<-c("Asymmetrical\ntransformation","Asymmetrical\ntransformation","No transformation","No transformation","Symmetrical\ntransformation","Symmetrical\ntransformation")
selects<-c("Multi-locus NFDS","No selection","Multi-locus NFDS","No selection","Multi-locus NFDS","No selection")
summaries<-c("Asymmetrical transformation, NFDS","Asymmetrical transformation, no selection","No transformation, NFDS","No transformation, no selection","Symmetrical transformation, NFDS","Symmetrical transformation, no selection")
nrep<-100

# read input data
eq.loci<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/mass.locusFreq")
colnames(eq.loci)<-c("Locus","Frequency")
snp.input<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/filtered.mass.snp.input",header=TRUE,comment.char = "")
snp.loci<-data.frame("Locus"=colnames(snp.input[6:ncol(snp.input)]),"Frequency"=as.numeric(colMeans(snp.input[,6:ncol(snp.input)])))

# read input file
ma.data<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/filtered.mass.input",header=TRUE,comment.char = "")
genome_sizes<-as.data.frame(rowMeans(ma.data[,6:ncol(ma.data)]))
colnames(genome_sizes)[1]<-"Sizes"
genome_distances<-as.data.frame(as.numeric(rdist(ma.data[,6:ncol(ma.data)],metric="jaccard")))
colnames(genome_distances)[1]<-"Distances"

# calculate SNP input data
snp.ma.data<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/filtered.mass.snp.input",header=TRUE,comment.char = "")
snp.ma.data<-snp.ma.data[match(ma.data$Taxon,snp.ma.data$Taxon),]
snp_distances<-as.data.frame(as.numeric(rdist(snp.ma.data[,6:ncol(snp.ma.data)],metric="hamming")))
genome_distances$SNPDistances<-snp_distances[,1]

# calculate core tree
ma.tree<-make_nj_tree(snp.ma.data[,6:ncol(snp.ma.data)])

# analysis with markers
marker_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/rec_v_select_with_markers_output/rec_v_select_with_markers",
                   "~/Documents/evoMLNFDS/rec_v_select/rec_v_no_select_with_markers_output/rec_v_no_select_with_markers",
                   "~/Documents/evoMLNFDS/rec_v_select/no_rec_v_select_with_markers_output/no_rec_v_select_with_markers",
                   "~/Documents/evoMLNFDS/rec_v_select/no_rec_v_no_select_with_markers_output/no_rec_v_no_select_with_markers",
                   "~/Documents/evoMLNFDS/rec_v_select/sym_rec_v_select_with_markers_output/sym_rec_v_select_with_markers",
                   "~/Documents/evoMLNFDS/rec_v_select/sym_rec_v_no_select_with_markers_output/sym_rec_v_no_select_with_markers"
)
marker_plots<-process_simulation_data(prefixes=marker_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE)

# save figures
for (i in 1:length(marker_plots)) {
  fn<-paste0("~/Documents/evoMLNFDS/rec_v_select/figures/Figure_A",i,".png")
  ggsave(marker_plots[[i]],file=fn)
}

# run processing for migration samples
no_migration_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/no_mig_rec_v_select_output/no_mig_rec_v_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_rec_v_no_select_output/no_mig_rec_v_no_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_no_rec_v_select_output/no_mig_no_rec_v_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_no_rec_v_no_select_output/no_mig_no_rec_v_no_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_sym_rec_v_select_output/no_mig_sym_rec_v_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_sym_rec_v_no_select_output/no_mig_sym_rec_v_no_select"
)
no_mig_plots<-process_simulation_data(prefixes=no_migration_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE)

# save figures
for (i in 1:length(no_mig_plots)) {
  fn<-paste0("~/Documents/evoMLNFDS/rec_v_select/figures/Figure_",i,".png")
  ggsave(no_mig_plots[[i]],file=fn)
}

# weak selection analysis
weak_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/no_mig_rec_v_weak_select_output/no_mig_rec_v_weak_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_rec_v_no_select_output/no_mig_rec_v_no_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_no_rec_v_weak_select_output/no_mig_no_rec_v_weak_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_no_rec_v_no_select_output/no_mig_no_rec_v_no_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_sym_rec_v_weak_select_output/no_mig_sym_rec_v_weak_select",
                         "~/Documents/evoMLNFDS/rec_v_select/no_mig_sym_rec_v_no_select_output/no_mig_sym_rec_v_no_select"
)

weak_plots<-process_simulation_data(prefixes=weak_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE)

# save figures
for (i in 1:length(weak_plots)) {
  fn<-paste0("~/Documents/evoMLNFDS/rec_v_select/figures/Figure_B",i,".png")
  ggsave(weak_plots[[i]],file=fn)
}

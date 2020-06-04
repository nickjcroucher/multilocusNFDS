# load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(rdist)
library(ape)
library(ggtree)
library(phytools)
library(cowplot)
library(stringr)
library(igraph)
library(reshape2)
library(ggthemes)
library(ggsci)
library(ggpubr)
library(propagate)
#library(pracma)
#library(pdist)

# calculate strain frequencies
get_strain_frequencies<-function(dists,slope=strain_divide_slope,intercept=strain_divide_intercept) {
  node_list<-unique(c(as.character(dists$Taxon1),as.character(dists$Taxon2)))
  edge_data<-dists[dists$Distances<(intercept+slope*dists$SNPDistances),1:2]
  strain_network<-graph_from_data_frame(d = edge_data, vertices = node_list, directed = FALSE)
  strains<-components(strain_network)
  strain_counts<-strains$csize/sum(strains$csize)
  strain_counts<-strain_counts[rev(order(strain_counts))]
  return(data.frame("Rank"=1:length(strain_counts),"Frequencies"=strain_counts))
}

# subsample tree - not used currently
# subsampling badly affects the Gamma statistic
subsample_tree<-function(r,t=ma.tree,m=616,n=100) {
  tip.sample<-sample(1:m,n,replace=FALSE)
  sub.tree<-keep.tip(t,tip.sample)
  return(sub.tree)
}

# generate tree
make_nj_tree<-function(i) {
  # Subsample tree down to 100 isolates
  # for a presentable tree
  # Use fan layout to lose the ultrametric
  # presentation of circle
  # Hamming distance more appropriate for
  # SNP comparisons
  # Process distance matrix
  gtree<-NULL
  if (class(i)=="data.frame") {
    j<-make.names(i[,1],unique = TRUE)
    s<-i[,2:ncol(i)]
    rownames(s)<-j
    colnames(s)<-rownames(s)
    # Make tree
    gtree<-nj(rdist(s,metric="hamming"))
  } else if (class(i)=="dist") {
    gtree<-nj(i)
    gtree$tip.label<-make.names(gtree$tip.label,unique = TRUE)
  } else {
    message("Need to provide data frame or distance matrix for tree generation\n")
    return(NULL)
  }
  gtree<-midpoint.root(gtree)
  return(gtree)
}

plot_nj_tree<-function(gtree,anno=ma.tip.labels,point.size=0.5,text.size=3,expansion=1.35,offset.size=10,palette=tree_palette) {
  # Update annotation for matching
  anno$taxa<-make.names(anno$taxa)
  # Process annotation
  k<-NULL
  k<-data.frame("taxa"=gtree$tip.label,stringsAsFactors = FALSE)
  tmp.k<-as.data.frame(str_split_fixed(k$taxa,"\\.[0-9]{1,3}$",2),stringsAsFactors = FALSE)
  k<-cbind(k,as.data.frame(str_split_fixed(tmp.k[,1],"_",2),stringsAsFactors = FALSE))
  colnames(k)[2:3]<-c("Original","Suffix")
  # fix anomaly in naming scheme
  k %<>%
    mutate(Original = case_when(
        taxa == "DH8X5.1" ~ "DH8X5.1",
        taxa == "DH8X5.2" ~ "DH8X5.2",
        TRUE ~ Original
      )
    )
  k$SC<-anno$SC[match(k$Original,anno$taxa)]
  # plot tree
  t.plot<-ggtree(gtree, layout="fan") +
    geom_treescale(x=expansion*max(branching.times(gtree)),y=5*616/8,width=0.025,offset=offset.size/point.size,linesize = 2,fontsize = text.size)
  
  t.plot <- t.plot %<+% k + geom_tippoint(aes(colour=SC),size=point.size) + palette + theme(plot.title = element_text(vjust = -35, hjust = 0.5)) #scale_colour_manual(values=c("#000000",rainbow(15)))
  return(t.plot)
}

return_distances<-function(input,snp=FALSE) {
  gene_dists<-NULL
  if (snp) {
    gene_dists<-rdist(input,metric="hamming")
  } else {
    gene_dists<-rdist(input,metric="jaccard")
  }
  return(gene_dists)
}

reformat_distances<-function(input,snp=FALSE) {
  genome_distances<-reshape2::melt(as.matrix(input), varnames = c("Taxon1", "Taxon2"))
  if (snp) {
    colnames(genome_distances)[3]<-"SNPDistances"
  } else {
    colnames(genome_distances)[3]<-"Distances"
  }
  return(genome_distances)
}

getGamma<-function(t) {
  u<-force.ultrametric(t)
  u$edge.length[u$edge.length<0]<-0
  return(gammaStat(u))
}

tree_analysis<-function(t) {
  gamma<-getGamma(t)
  pd_n<-mean(t$edge.length)
  t_df<-data.frame("PD"=pd_n,"Gamma"=gamma)
  return(t_df)
}

# read matrices
read_matrices<-function(prefix, r, real = NULL, snp.calc.proc=FALSE) {
  # read data
  in.matrices<-list()
  snp.matrices<-list()
  for (i in 1:r) {
    fn<-paste(prefix,i,"finalPopGenotypes.tab",sep = '.')
    in.matrices[[i]]<-read.table(file=fn,header=T,comment.char = "")
    rownames(in.matrices[[i]])<-make.unique(as.character(in.matrices[[i]][,1]))
    in.matrices[[i]]<-in.matrices[[i]][,-1]
    if (snp.calc.proc==TRUE) {
      snp_fn<-paste(prefix,i,"markers.finalPopGenotypes.tab",sep = '.')
      snp.matrices[[i]]<-read.table(file=snp_fn,header=T,comment.char = "")
      snp.matrices[[i]]<-snp.matrices[[i]][,-1]
      rownames(snp.matrices[[i]])<-rownames(in.matrices[[i]])
    }
  }
  # process inputs - gene frequencies
  gene.freqs<-lapply(in.matrices,colMeans)
  gene.stats<-do.call(cbind, gene.freqs)
  # process inputs - genome content
  genome.sizes<-lapply(in.matrices,rowMeans)
  genome.stats<-do.call(cbind, genome.sizes)
  # process inputs - distances between isolates
  dist.stats<-lapply(in.matrices,return_distances)
  dist.stats<-lapply(dist.stats,reformat_distances)
  #genome.dists<-lapply(in.matrices,rdist,metric="jaccard")
  #dist.stats<-do.call(cbind, genome.dists)
  # process inputs - SNP frequencies
  snp.dists<-matrix()
  snp.trees<-list()
  tree.stats<-list()
  #snp.dist.example<-NULL
  if (snp.calc.proc==TRUE) {
    snp.freqs<-lapply(snp.matrices,colMeans)
    snp.stats<-do.call(cbind, snp.freqs)
    #snp.dists<-lapply(snp.matrices,rdist,metric="hamming")
    snp.dists<-lapply(snp.matrices,return_distances,snp=TRUE)
    snp.trees<-lapply(snp.dists,make_nj_tree)
    tree.stats<-lapply(snp.trees,tree_analysis) # need to write tree_analysis; then make plots internally not return
    snp.dists<-lapply(snp.dists,reformat_distances,snp=TRUE)
    #snp.dist.example<-snp.dists[[1]]
  }
  #snp.dist.stats<-do.call(cbind, snp.dists)
  # compare distances to input data
  comp.dists<-lapply(in.matrices,get_comparative_dists,as.matrix(real[,6:ncol(real)]))
  comp.dist.stats<-do.call(cbind, comp.dists)
  # output matrices
  output.list<-list(gene.stats,genome.stats,dist.stats,snp.stats,snp.dists,snp.trees,comp.dist.stats,tree.stats)
  return(output.list)
}

get_comparative_dists<-function(i,r) {
  mod.r<-r[,match(colnames(i),colnames(r))]
  d<-cdist(mod.r,i,metric="jaccard") # rows are real data, columns are per query
  min.d<-apply(d,2,min)
  return(min.d)
}

# convert data to tidyverse format
convert_to_df<-function(outlist,el=NULL,snp.info=NULL,rec=NULL,select=NULL,summary=NULL,snp.calc.proc=FALSE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept) {
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
  #save(outlist,file="~/Desktop/outlist.Rdata")
  #dists.data<-suppressWarnings(bind_rows(outlist[[3]],.id="Replicate")) # Eliminated because sometimes could not deal with Taxon name columns and char/num values
  for (i in 1:length(outlist[[3]])) {outlist[[3]][[i]]$Replicate<-i}
  dists.data<-plyr::rbind.fill(outlist[[3]])
  
  # Process SNP information if provided
  snp.df<-data.frame()
  snp.dists.df<-data.frame()
  tree.stats.df<-data.frame()
  strain.frequencies<-NULL
  
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
    snp.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Locus"=rownames(snp.data),"Replicate"=snp.data[,1],"Frequency"=snp.data[,2],"Equilibrium"=ordered.snp.freq)
    # SNP distance matrix
    snp.dists.data<-matrix(nrow=0,ncol=2)
    combined.dist.list<-as.list(rep(0,times=length(outlist[[5]])))
    for (c in 1:length(outlist[[5]])) {
      #combined.dist.list[[c]]<-merge(outlist[[5]][[c]],outlist[[3]][[c]],by=c("Taxon1","Taxon2"))
      combined.dist.list[[c]]<-dplyr::inner_join(outlist[[5]][[c]],outlist[[3]][[c]],by=c("Taxon1","Taxon2"))
    }
    # Create joint data frame
    #snp.dists.data<-suppressWarnings(bind_rows(combined.dist.list,.id="Replicate"))
    for (i in 1:length(combined.dist.list)) {combined.dist.list[[i]]$Replicate<-i} # same avoidance of bind_rows as above
    snp.dists.data<-plyr::rbind.fill(combined.dist.list)
    snp.dists.data$Recombination<-rec
    snp.dists.data$Selection<-select
    # Calculate number of strains in population
    strain.frequency.list<-lapply(combined.dist.list,get_strain_frequencies,slope=strain_slope,intercept=strain_intercept)
    for (i in 1:length(strain.frequency.list)) {strain.frequency.list[[i]]$Replicate<-i} # same avoidance of bind_rows as above
    strain.frequencies<-plyr::rbind.fill(strain.frequency.list)
    #strain.frequencies<-suppressWarnings(bind_rows(strain.frequency.list,.id="Replicate"))
    strain.frequencies$Recombination<-rec
    strain.frequencies$Selection<-select
    # process tree statistics
    tree.stats.df<-plyr::rbind.fill(outlist[[8]])
    tree.stats.df$Recombination<-rec
    tree.stats.df$Selection<-select
    tree.stats.df<-relevel_df(tree.stats.df)
#    tree.stats.df$Simulation<-paste(tree.stats.df$Selection,tree.stats.df$Recombination,sep="\n")
  }
  
  # Downsample distance data for memory efficiency
  distance.sampling<-sample(nrow(dists.data),size=round(nrow(dists.data)/100),replace=FALSE)
  dists.data<-dists.data[distance.sampling,] # downsample for memory efficiency
  message("Combining pairwise distances\n")
  dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(dists.data),"Replicate"=dists.data$Replicate,"Distances"=dists.data$Distances)
  if (snp.calc.proc==TRUE) {
    snp.dists.data<-snp.dists.data[distance.sampling,] # downsample for memory efficiency - same rows as above
    message("Combining SNP pairwise distances\n")
   # snp.dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(snp.dists.data),"Replicate"=snp.dists.data[,1],"Distances"=snp.dists.data[,2])
  }
  
  # process comparisons with starting data
  message("Processing comparisons with starting data\n")
  comp.dists.data<-matrix(nrow=0,ncol=2)
  for (c in 1:ncol(outlist[[7]])) {
    comp.dists.data<-rbind(comp.dists.data,cbind(c,outlist[[7]][,c]))
  }
  comp.distance.sampling<-sample(nrow(comp.dists.data),size=round(nrow(comp.dists.data)/100),replace=FALSE)
  comp.dists.data<-comp.dists.data[comp.distance.sampling,] # downsample for memory efficiency
  comp.dists.df<-data.frame("Type"=summary,"Recombination"=rec,"Selection"=select,"Index"=1:nrow(comp.dists.data),"Replicate"=comp.dists.data[,1],"Distances"=comp.dists.data[,2])
  
  # complete
  return(list(gene.df,genome.df,dists.df,snp.df,snp.dists.data,comp.dists.df,strain.frequencies,tree.stats.df))
}

relevel_df<-function(df) {
  # rec_levels<-c("No transformation","Symmetrical transformation","Asymmetrical transformation")
  # select_levels<-c("Neutral","Multi-locus NFDS")
  r<-match("Recombination",colnames(df))
  s<-match("Selection",colnames(df))
  extra_rec_levels<-as.character(unique(df[,r])[!(unique(df[,r]) %in% rec_levels)])
  extra_select_levels<-as.character(unique(df[,s])[!(unique(df[,s]) %in% select_levels)])
  if (length(extra_rec_levels) > 0) {rec_levels<-c(rec_levels,extra_rec_levels)}
  if (length(extra_select_levels) > 0) {select_levels<-c(select_levels,extra_select_levels)}
  df[,r]<-factor(df[,r],levels = rec_levels)
  df[,s]<-factor(df[,s],levels = select_levels)
  return(df)
}

getDistStats<-function(in.dist,varname) {
  summary_formula<-as.formula(paste0(varname,"~Recombination+Selection"))
  kurtosis.df<-aggregate(summary_formula,data=in.dist,propagate::kurtosis)
  colnames(kurtosis.df)[ncol(kurtosis.df)]<-"Kurtosis"
  skewness.df<-aggregate(summary_formula,data=in.dist,propagate::skewness)
  colnames(skewness.df)[ncol(skewness.df)]<-"Skewness"
  mean.df<-aggregate(summary_formula,data=in.dist,mean)
  colnames(mean.df)[ncol(mean.df)]<-"Mean"
  var.df<-aggregate(summary_formula,data=in.dist,var)
  colnames(var.df)[ncol(var.df)]<-"Variance"
  summary.df<-dplyr::inner_join(kurtosis.df,skewness.df,by=c("Recombination","Selection"))
  kurtosis.df<-NULL
  skewness.df<-NULL
  summary.df<-dplyr::inner_join(summary.df,mean.df,by=c("Recombination","Selection"))
  mean.df<-NULL
  summary.df<-dplyr::inner_join(summary.df,var.df,by=c("Recombination","Selection"))
  mean.df<-NULL
  summary.df<-relevel_df(summary.df)
  return(summary.df)
}

process_simulation_data<-function(prefixes=NULL,recs=NULL,selects=NULL,summaries=NULL,eq.loci=NULL,real.data=NULL,snp.loci=NULL,snp.calc=FALSE,strain_slope=0,strain_intercept=0,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs,tree_data=ma.tree.stats) {
  # store plots
  snp.calc.val<-snp.calc
  out_plots<-list()
  locus.data<-list()
  genotype.data<-list()
  dist.data<-list()
  snp.data<-list()
  snp.dist.data<-list()
  strain.freq.data<-list()
  strain.freq.df<-data.frame()
  tree.data<-list()
  tree.df<-data.frame()
  comp.dist.data<-list()
  tree.stats<-list()
  
  # processing
  message("Initial processing of data frames\n")
  for (x in 1:length(prefixes)) {
    message(paste0("Processing ",prefixes[x],"\n"))
    processed.summaries<-read_matrices(prefixes[x],nrep,real=real.data,snp.calc.proc=snp.calc.val)
    message(paste0("Matrices read\n"))
    out.list<-convert_to_df(processed.summaries,el=eq.loci,snp.info=snp.loci,rec=recs[x],select=selects[x],summary=summaries[x],snp.calc.proc=snp.calc.val,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept)
    if (snp.calc) {
      message(paste0("Plotting trees\n"))
      # in.fn<-paste(prefixes[x],1,"markers.finalPopGenotypes.tab",sep = '.')
      # tree.input<-read.table(in.fn,header=TRUE,comment.char = "")
      # nj_tree<-make_nj_tree(tree.input)
      tree.data[[x]]<-plot_nj_tree(processed.summaries[[6]][[1]],anno=tip.anno)
    }
    message(paste0("Converted to data frame\n"))
    processed.summaries<-NULL # reduce memory
    locus.data[[x]]<-out.list[[1]]
    genotype.data[[x]]<-out.list[[2]]
    dist.data[[x]]<-out.list[[3]]
    if (snp.calc.val==TRUE) {
      snp.data[[x]]<- out.list[[4]]
      snp.dist.data[[x]]<-out.list[[5]]
      strain.freq.data[[x]]<-out.list[[7]]
      tree.stats[[x]]<-out.list[[8]]
    }
    comp.dist.data[[x]]<-out.list[[6]]
    out.list<-NULL # reduce memory
  }
  
  # process accessory locus data
  message("Merging data frames\n")
  #raw.locus.df<-bind_rows(locus.data)
  raw.locus.df<-do.call(rbind,locus.data)
  mean.locus.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.locus.df,mean)
  lb.locus.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.locus.df,min)
  ub.locus.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.locus.df,max)
  fixed.locus.df<-aggregate(cbind(Fixed=(Frequency==1|Frequency==0))~Type+Recombination+Selection+Locus+Equilibrium, data=raw.locus.df,  sum)
  locus.df<-cbind(mean.locus.df,lb.locus.df$Frequency,ub.locus.df$Frequency,fixed.locus.df$Fixed/nrep)
  mean.locus.df<-NULL # reduce memory
  ub.locus.df<-NULL # reduce memory
  lb.locus.df<-NULL # reduce memory
  locus.df<-relevel_df(locus.df)
  colnames(locus.df)[colnames(locus.df)=="lb.locus.df$Frequency"]<-"LowerBound"
  colnames(locus.df)[colnames(locus.df)=="ub.locus.df$Frequency"]<-"UpperBound"
  colnames(locus.df)[colnames(locus.df)=="fixed.locus.df$Fixed/nrep"]<-"Fixed"
  
  # process genotype data
  genotype.df<-do.call(rbind,genotype.data)
  genotype.df<-relevel_df(genotype.df)
  genotype.summary.df<-getDistStats(genotype.df,"Frequency")
  
  # distance data
  dist.df<-do.call(rbind,dist.data)
  dist.df<-relevel_df(dist.df)
  dist.summary.df<-getDistStats(dist.df,"Distances")
  
  # comparative distance data
  comp.dist.df<-do.call(rbind,comp.dist.data)
  comp.dist.df<-relevel_df(comp.dist.df)
  
  # SNP data
  if (snp.calc==TRUE) {
    # SNP frequencies
    message("Merging SNP data frames\n")
    #raw.snp.df<-bind_rows(snp.data)
    raw.snp.df<-do.call(rbind,snp.data)
    mean.snp.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,mean)
    lb.snp.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,min)
    ub.snp.df<-aggregate(Frequency~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,max)
    fixed.snp.df<-aggregate(cbind(Fixed=(Frequency==1|Frequency==0))~Type+Recombination+Selection+Locus+Equilibrium,data=raw.snp.df,sum)
    snp.df<-cbind(mean.snp.df,lb.snp.df$Frequency,ub.snp.df$Frequency,fixed.snp.df$Fixed/nrep)
    mean.snp.df<-NULL # reduce memory
    ub.snp.df<-NULL # reduce memory
    lb.snp.df<-NULL # reduce memory
    snp.df<-relevel_df(snp.df)
    colnames(snp.df)[colnames(snp.df)=="lb.snp.df$Frequency"]<-"LowerBound"
    colnames(snp.df)[colnames(snp.df)=="ub.snp.df$Frequency"]<-"UpperBound"
    colnames(snp.df)[colnames(snp.df)=="fixed.snp.df$Fixed/nrep"]<-"Fixed"
    # SNP distances
    #snp.dist.df<-bind_rows(snp.dist.data)
    snp.dist.df<-do.call(rbind,snp.dist.data)
    snp.dist.df<-relevel_df(snp.dist.df)
    snp.dist.summary.df<-getDistStats(snp.dist.df,"SNPDistances")
    # Tree plots
    tree.df<-data.frame("Recombination" = recs, "Selection" = selects, "Tree" = 1:length(tree.data))
    tree.df<-relevel_df(tree.df)
    # Strain frequencies
    raw.strain.freq.df<-do.call(rbind,strain.freq.data)
    save(raw.strain.freq.df,file="~/Desktop/strain_freq.Rdata")
    strain.freq.df<-raw.strain.freq.df %>%
      group_by(Recombination,Selection,Rank) %>%
      dplyr::summarize(Frequency = mean(Frequencies, na.rm = TRUE), MinFrequency = min(Frequencies, na.rm = TRUE), MaxFrequency = max(Frequencies, na.rm = TRUE))
    strain.freq.df<-relevel_df(as.data.frame(strain.freq.df))
    # tree statistics
    tree.stats[[length(tree.stats)+1]]<-tree_data
    tree.stats.df<-do.call(rbind,tree.stats)
    tree.stats.df<-relevel_df(tree.stats.df)
    tree.stats.df$Simulation<-paste(tree.stats.df$Selection,sub("\n"," ",tolower(tree.stats.df$Recombination)),sep=",\n")
    tree.stats.df$Simulation[grep("enomic",tree.stats.df$Simulation)]<-"Genomic data"
  }
  
  ############
  # Plotting #
  ############
  
  # plot gene frequencies
  message("Plotting data\n")
  freq_plot<-ggplot(data=locus.df,aes(x=Equilibrium,y=Frequency,colour=Fixed,fill=Fixed))+
    geom_point(alpha=my_scatter_alpha,size=0.5,shape=19,position=position_jitter(width=my_jitter_width,seed=my_jitter_seed))+
    geom_linerange(aes(ymin=LowerBound,ymax=UpperBound),size=0.5,alpha=my_line_alpha,position=position_jitter(width=my_jitter_width,seed=my_jitter_seed))+
    xlab(expression(bold("Initial frequency in genomic data"~(italic(f)[italic(l)*",0"]))))+
    ylab(expression(bold("Simulated frequency at final timepoint"~(italic(f)[italic(l)*",60000"]))))+
    facet_grid(Recombination~Selection)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = my_line_col)+
    scale_colour_gradient2(low = my_fill_low, mid = my_fill_mid, high = my_fill_high, midpoint = 0.5,name="Frequency of\nfixation")+
    scale_fill_gradient2(low = my_fill_low, mid = my_fill_mid, high = my_fill_high, midpoint = 0.5,name="Frequency of\nfixation")+
    stat_cor(aes(label = paste("rho ==",bquote(.(..r..)))),colour="black",label.x=0.05,label.y=0.9,hjust = 0,method="spearman",cor.coef.name="rho")
  out_plots[[1]]<-freq_plot

  # plot genome sizes
  genome_plot<-ggplot(data=genotype.df,aes(x=Frequency))+
    geom_density(aes(y=..density..),fill=my_fill_col,alpha=0.5,linetype=0)+
    facet_grid(Recombination~Selection)+
    xlab(expression("Proportion of loci in genome ("*Sigma^italic(L)*italic(g)[italic("i,l")]*")"))+
    ylab("Density")+
    geom_density(data=genome_sizes, aes(x=Sizes), col = my_line_col, trim=TRUE)+
    geom_text(data=genotype.summary.df,x=0.125,y=25,hjust=0,size=my_geom_text_size,aes(label=sprintf("Mean: %.2f\nVariance: %.3g\nSkewness: %.2f\n",Mean,Variance,Skewness)))
  out_plots[[2]]<-genome_plot

  # plot pairwise distances
  distance_plot<-ggplot(data=dist.df,aes(x=Distances))+
    geom_density(aes(y=..ndensity..),fill=my_fill_col,alpha=0.5,linetype=0)+
    facet_grid(Recombination~Selection)+
    xlab("Pairwise binary Jaccard distances between genomes")+
    ylab("Density (log(density+1) scale)")+
    geom_density(data=genome_distances, aes(x=Distances, y=..scaled..), col = my_line_col, trim=TRUE)+
    scale_y_continuous(trans="log1p")+
    geom_text(data=dist.summary.df,x=0.0,y=0.6,hjust=0,size=my_geom_text_size,aes(label=sprintf("Mean: %.2f\nVariance: %.3g\nSkewness: %.2f\n",Mean,Variance,Skewness)))
  out_plots[[3]]<-distance_plot

  # plot distances to real data
  comp_dist_plot<-ggplot(data=comp.dist.df,aes(x=Distances))+
    geom_density(aes(y=..ndensity..),fill=my_fill_col,alpha=0.5,linetype=0,adjust=1e-2)+
    facet_grid(Recombination~Selection)+
    xlab("Minimum pairwise binary Jaccard distances between simulated and initial data")+
    ylab("Density (log(density+1) scale)")+
    scale_y_continuous(trans="log1p")
  out_plots[[4]]<-comp_dist_plot
  
  # plot SNP data
  if (snp.calc==TRUE) {
    # plot SNP frequencies
    snp_plot<-ggplot(data=snp.df,aes(x=Equilibrium,y=Frequency,fill=Fixed,colour=Fixed))+
      geom_point(alpha=my_scatter_alpha,shape=19,size=0.5,position=position_jitter(width=my_jitter_width,seed=my_jitter_seed))+
      geom_linerange(aes(ymin=LowerBound,ymax=UpperBound),size=0.5,alpha=my_line_alpha,position=position_jitter(width=my_jitter_width,seed=my_jitter_seed))+
      xlab(expression(bold("Initial frequency in genomic data"~(italic(f)[italic(s)*",0"]))))+
      ylab(expression(bold("Simulated frequency at final timepoint"~(italic(f)[italic(s)*",60000"]))))+
      facet_grid(Recombination~Selection)+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red")+
      scale_colour_gradient2(low = my_fill_low,  high = my_fill_high, mid = my_fill_mid, midpoint = 0.5,name="Frequency of\nfixation")+
      scale_fill_gradient2(low = my_fill_low, high = my_fill_high, mid = my_fill_mid, midpoint = 0.5,name="Frequency of\nfixation")+
      stat_cor(aes(label = paste("rho ==",bquote(.(..r..)))),colour="black",label.x=0.05,label.y=0.9,hjust = 0,method="spearman",cor.coef.name="rho")
    out_plots[[5]]<-snp_plot
    
    # plot SNP distances
    snp_distance_plot<-ggplot(data=snp.dist.df,aes(x=SNPDistances))+
      geom_density(aes(y=..ndensity..),fill=my_fill_col,alpha=0.5,linetype=0)+
      facet_grid(Recombination~Selection)+
      xlab("Pairwise Hamming SNP distances between genomes")+
      ylab("Density (log(density+1) scale)")+
      geom_density(data=genome_distances, aes(x=SNPDistances, y=..scaled..), col = my_line_col, trim=TRUE)+
      scale_y_continuous(trans="log1p")+
      geom_text(data=snp.dist.summary.df,x=0.0,y=0.6,hjust=0,size=my_geom_text_size,aes(label=sprintf("Mean: %.2f\nVariance: %.3g\nSkewness: %.2f\n",Mean,Variance,Skewness)))
    out_plots[[6]]<-snp_distance_plot
    
    # plot SNP v accessory distances
    #snp.dist.df$Accessory<-dist.df$Distances
    comparison_distance_plot<-ggplot(data=snp.dist.df,aes(x=SNPDistances,y=Distances))+
      geom_point(alpha = 5000/nrow(snp.dist.df),size=0.5)+
      facet_grid(Recombination~Selection)+
      xlab("Pairwise Hamming SNP distances between genomes")+
      ylab("Pairwise binary Jaccard distances between genomes")+
      stat_cor(aes(label = paste("rho ==",bquote(.(..r..)))),colour="black",label.x=0.3,label.y=0.2,hjust = 0,method="spearman",cor.coef.name="rho")+
      #stat_cor(aes(label = paste(..rr.label..)),colour="black",label.x=0.3,label.y=0.2,hjust = 0)+
      #geom_density_2d(data=genome_distances, aes(x=SNPDistances, y=Distances), colour = "red", alpha = 0.5)
      geom_abline(colour=my_line_col, slope = strain_slope, intercept = strain_intercept, alpha=0.5, linetype = 2)+
      geom_density_2d(data=genome_distances[genome_distances$Distances>=(strain_intercept+strain_slope*genome_distances$SNPDistances),], aes(x=SNPDistances, y=Distances), colour = pal_npg()(2)[2], alpha = 0.5,size=0.5)+
      geom_density_2d(data=genome_distances[genome_distances$Distances<(strain_intercept+strain_slope*genome_distances$SNPDistances),], aes(x=SNPDistances, y=Distances), colour = my_fill_high, alpha = 0.5,size=0.5)
    
    out_plots[[7]]<-comparison_distance_plot
    
    # plot trees
    tree.plot.list<-list()
    new.tree.plot.list<-list()
    remove_vector<-c("tag","caption","title","subtitle","ylab-r","ylab-l","xlab-b","xlab-t","spacer","background","axis-t","axis-r","axis-l","axis-b")
    for (r in levels(tree.df$Recombination)) {
      for (s in levels(tree.df$Selection)) {
        tree.plot.list[[length(tree.plot.list)+1]]<-tree.data[[tree.df$Tree[tree.df$Recombination==r & tree.df$Selection==s]]]+theme(panel.background=element_rect(fill = "transparent",colour = NA))
      }
    }
    for (i in 1:length(tree.plot.list)) {
      test.gg<-ggplot_build(tree.plot.list[[i]])
      for (l in 1:4) {
        test.gg$data[[l]]$colour<-tree_pal[i]
      }
      new.tree.plot.list[[i]]<-gridExtra::arrangeGrob(gtable::gtable_filter((ggplot_gtable(test.gg)),"panel",trim=TRUE))
    }
    #out_plots[[8]]<-plot_grid(plotlist = new.tree.plot.list,nrow=3,ncol=2,scale=c(1.4,1.4,1.4,1.4,1.4,1.4))
    out_plots[[8]]<-cowplot::plot_grid(plotlist = new.tree.plot.list,
                                       nrow=3,
                                       ncol=2,
                                       scale=c(1.3,1.3,1.3,1.3,1.3,1.3),
                                       labels=paste0(rep(c("Neutral","Multi-locus NFDS"),times=3),",\n ",rep(tolower(c("No transformation","Symmetrical transformation","Asymmetrical transformation")),each=2)),
                                       label_fontface="plain",
                                       label_size = 9,
                                       hjust = 0)
    
    # plot strain frequencies
    max_rank<-20
    strain_frequency_plot<-ggplot(data=strain.freq.df[strain.freq.df$Rank<=max_rank,],aes(x=Rank,y=Frequency))+
      geom_point(color=my_fill_low)+
      geom_linerange(aes(ymin=MinFrequency,ymax=MaxFrequency),alpha=0.5,color=my_fill_low)+
      facet_grid(Recombination~Selection)+
      geom_point(data=real_strains[real_strains$Rank<=max_rank,],aes(x=Rank,y=Frequencies),colour=my_line_col,shape=4)+
      xlim(0,max_rank)+
      xlab("Ranked frequency of strain")+
      ylab("Frequency in population")+
      scale_y_continuous(trans="log1p")
    out_plots[[9]]<-strain_frequency_plot
  }
  
  # plot tree statistics
  tree.stats.df$Simulation<-factor(tree.stats.df$Simulation,levels=unique(tree.stats.df$Simulation[order(tree.stats.df$Selection,tree.stats.df$Recombination,tree.stats.df$Simulation)], ordered=TRUE))
  #new_tree_pal<-c(tree_pal[c(2,4,6)],tree_pal[c(1,3,5)],"#000000")
  tree_stats_plot<-ggplot(data=tree.stats.df,aes(x=PD,y=Gamma,colour=Simulation,shape=Simulation))+
    geom_point(alpha=0.75)+
    xlab("Phylogenetic diversity per tip")+
    ylab(expression(bold("Pybus & Harvey ")~bold(gamma)))+
    scale_fill_manual(values = new_tree_pal)+
    scale_colour_manual(values = new_tree_pal)+
    scale_shape_manual(values=new_tree_shapes)+
    stat_ellipse(data=tree.stats.df[tree.stats.df$Simulation!="Genomic data",],aes(group=Simulation),type="norm")+
    theme(legend.position = "bottom")
  out_plots[[10]]<-tree_stats_plot
  
  # compare genome sizes to genomic data
  genome_size_comparison_plot<-plot_simulation_comparison(genotype.df,genome_sizes$Sizes,vals_colname="Frequency")   
  out_plots[[11]]<-genome_size_comparison_plot
  
  # compare pairwise accessory distances
  accessory_dist_comparison_plot<-plot_simulation_comparison(dist.df,genome_distances$Distances,vals_colname="Distances")   
  out_plots[[12]]<-accessory_dist_comparison_plot
  
  # compare pairwise SNP distances
  snp_dist_comparison_plot<-plot_simulation_comparison(snp.dist.df,genome_distances$SNPDistances,vals_colname="SNPDistances")   
  out_plots[[13]]<-snp_dist_comparison_plot
  
  # return outputs
  output_list<-list(out_plots,processed.summaries[[6]])
  return(output_list)
}

# format graphs for publication
format_plots<-function(raw_plot_list,prefix="test") {
  
  # scale function
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # frequency plots
  freq_plots<-egg::ggarrange(raw_plot_list[[1]][[1]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                       panel.spacing = unit(0, "cm"),
                                                                       legend.position="none",
                                                                       panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                       axis.text.y = element_text(size = 8),
                                                                       axis.text.x = element_text(angle = 90, hjust = 1),
                                                                       axis.title.x = element_text(size = 10, face = "bold"),
                                                                       axis.title.y = element_text(size = 10, face = "bold"),
                                                                       strip.text.x = element_text(size = my_facet_label_size),
                                                                       strip.text.y = element_text(size = my_facet_label_size)
  ),
  raw_plot_list[[1]][[5]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                            panel.spacing = unit(0, "cm"),
                                            panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                            axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                            axis.text.y = element_text(size = 8),
                                            axis.title.x = element_text(size = 10, face = "bold"),
                                            axis.title.y = element_text(size = 10, face = "bold"),
                                            strip.text.x = element_text(size = my_facet_label_size),
                                            strip.text.y = element_text(size = my_facet_label_size),
                                            legend.position="bottom"),
  nrow=1,ncol=2,labels=c("a","b"),label.args = list(gp=grid::gpar(fontface = "bold")))
  ggsave(freq_plots, file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"_freq_plot.png"), width = 21, height = 21, units="cm")
  # genome sizes
  size_plots<-ggarrange(raw_plot_list[[1]][[2]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                  panel.spacing = unit(0, "cm"),
                                                                  panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                  axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                                                  axis.text.y = element_text(size = 8),
                                                                  axis.title.x = element_text(size = 10, face = "bold"),
                                                                  axis.title.y = element_text(size = 10, face = "bold"),
                                                                  strip.text.x = element_text(size = my_facet_label_size),
                                                                  strip.text.y = element_text(size = my_facet_label_size)))
  ggsave(size_plots,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"size_plot.png"),width = 15,height = 17.5, units="cm")
  # distance plots
  dist_plots<-egg::ggarrange(raw_plot_list[[1]][[3]]+ scale_x_continuous(labels=scaleFUN)+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                                                            panel.spacing = unit(0, "cm"),
                                                                                                            panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                                                            axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                                                                                            axis.text.y = element_text(size = 8),
                                                                                                            axis.title.x = element_text(size = 10, face = "bold"),
                                                                                                            axis.title.y = element_text(size = 10, face = "bold"),
                                                                                                            strip.text.x = element_text(size = my_facet_label_size),
                                                                                                            strip.text.y = element_text(size = my_facet_label_size),
                                                                                                            legend.position="none"),
                             raw_plot_list[[1]][[6]]+ scale_x_continuous(labels=scaleFUN)+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                                                            panel.spacing = unit(0, "cm"),
                                                                                                            panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                                                            axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                                                                                            axis.text.y = element_text(size = 8),
                                                                                                            axis.title.x = element_text(size = 10, face = "bold"),
                                                                                                            axis.title.y = element_text(size = 10, face = "bold"),
                                                                                                            strip.text.x = element_text(size = my_facet_label_size),
                                                                                                            strip.text.y = element_text(size = my_facet_label_size),
                                                                                                            legend.position="none"),
                             nrow=1,ncol=2,labels=c("a","b"),label.args = list(gp=grid::gpar(fontface = "bold")))
  ggsave(dist_plots,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"dist_plot.png"),width = 21,height = 21, units="cm")
  # tree plots
  tree_plots<-ggpubr::ggarrange(raw_plot_list[[1]][[10]]+theme_minimal()+theme(legend.position = "bottom",aspect.ratio=2,
                                                                               axis.line.x = element_line(color="black", size = 0.25, linetype='solid'),
                                                                               axis.line.y = element_line(color="black", size = 0.25, linetype='solid'),
                                                                               axis.title.x = element_text(size = 10, face = "bold"),
                                                                               axis.title.y = element_text(size = 10, face = "bold")
                                                                              )+panel_border(remove=FALSE)+guides(colour=guide_legend(title.position="top",title.hjust =0.5)),
                        raw_plot_list[[1]][[8]]+theme(aspect.ratio=2)+panel_border(remove=TRUE),
                        labels = "auto",common.legend=TRUE,legend="bottom")
  ggsave(tree_plots,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"tree_plot.png"), width = 21, height = 21, units="cm")
  # pairwise distance plots
  pairwise_plot<-ggarrange(raw_plot_list[[1]][[7]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                     panel.spacing = unit(0, "cm"),
                                                                     panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                     axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                                                     axis.text.y = element_text(size = 8),
                                                                     axis.title.x = element_text(size = 10, face = "bold"),
                                                                     axis.title.y = element_text(size = 10, face = "bold"),
                                                                     strip.text.x = element_text(size = my_facet_label_size),
                                                                     strip.text.y = element_text(size = my_facet_label_size)))
  ggsave(pairwise_plot,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"pairwise_plot.png"),width = 15,height = 17.5, units="cm")
  # strain frequency plots
  freq_plot<-ggarrange(raw_plot_list[[1]][[9]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                 panel.spacing = unit(0, "cm"),
                                                                 panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                 axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                                                 axis.text.y = element_text(size = 8),
                                                                 axis.title.x = element_text(size = 10, face = "bold"),
                                                                 axis.title.y = element_text(size = 10, face = "bold"),
                                                                 strip.text.x = element_text(size = my_facet_label_size),
                                                                 strip.text.y = element_text(size = my_facet_label_size)))
  ggsave(freq_plot,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"freq_plot.png"),width = 15,height = 17.5, units="cm")
  
  # divergence plots
  div_plot<-ggarrange(raw_plot_list[[1]][[4]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                panel.spacing = unit(0, "cm"),
                                                                panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                                                                axis.text.y = element_text(size = 8),
                                                                axis.title.x = element_text(size = 10, face = "bold"),
                                                                axis.title.y = element_text(size = 10, face = "bold"),
                                                                strip.text.x = element_text(size = my_facet_label_size),
                                                                strip.text.y = element_text(size = my_facet_label_size)))
  ggsave(div_plot,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"divergence_plot_test.png"),width = 15,height = 17.5, units="cm")
  
  # genome size comparison
  gsize_comp_plot<-raw_plot_list[[1]][[11]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                 panel.spacing = unit(0, "cm"),
                                                 panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                 axis.text.x = element_text(size = 8, angle = 90, hjust=1, vjust=0.5),
                                                 axis.text.y = element_text(size = 8),
                                                 axis.title.x = element_text(size = 10, face = "bold"),
                                                 axis.title.y = element_text(size = 10, face = "bold"),
                                                 strip.text.x = element_text(size = my_facet_label_size),
                                                 strip.text.y = element_text(size = my_facet_label_size),
                                                 legend.text = element_text(size = 7),
                                                 legend.position = "bottom",
                                                 legend.title = element_blank())
  ggsave(gsize_comp_plot,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"gsize_comp_plot.png"),width = 15,height = 17.5, units="cm")
  
  # compare pairwise accessory distances
  accdist_comp_plot<-raw_plot_list[[1]][[12]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                  panel.spacing = unit(0, "cm"),
                                                                  panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                  axis.text.x = element_text(size = 8, angle = 90, hjust=0, vjust=0.5),
                                                                  axis.text.y = element_text(size = 8),
                                                                  axis.title.x = element_text(size = 10, face = "bold"),
                                                                  axis.title.y = element_text(size = 10, face = "bold"),
                                                                  strip.text.x = element_text(size = my_facet_label_size),
                                                                  strip.text.y = element_text(size = my_facet_label_size),
                                                                  legend.text = element_text(size = 7),
                                                                  legend.position = "bottom",
                                                                  legend.title = element_blank())
  ggsave(accdist_comp_plot,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"accdist_comp_plot.png"),width = 15,height = 17.5, units="cm")
  
  # compare pairwise SNP distances
  snpdist_comp_plot<-raw_plot_list[[1]][[13]]+theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                                                                    panel.spacing = unit(0, "cm"),
                                                                    panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                                                                    axis.text.x = element_text(size = 8, angle = 90, hjust=1, vjust=0.5),
                                                                    axis.text.y = element_text(size = 8),
                                                                    axis.title.x = element_text(size = 10, face = "bold"),
                                                                    axis.title.y = element_text(size = 10, face = "bold"),
                                                                    strip.text.x = element_text(size = my_facet_label_size),
                                                                    strip.text.y = element_text(size = my_facet_label_size),
                                                                    legend.text = element_text(size = 7),
                                                                    legend.position = "bottom",
                                                                    legend.title = element_blank())
  ggsave(snpdist_comp_plot,file=paste0("~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/",prefix,"snpdist_comp_plot.png"),width = 15,height = 17.5, units="cm")
}

plot_simulation_comparison<-function(sims,actual,vals_colname="Frequency") {

  rec_vals<-NULL
  select_vals<-NULL
  divergence_vals<-NULL
  
  for (r in rec_levels) {
    for (s in select_levels) {
      sim.set<-sims[sims$Recombination==r & sims$Selection==s,colnames(sims) %in% c(vals_colname,"Replicate")]
      for (i in unique(sim.set$Replicate)) {
        ktest<-suppressWarnings(ks.test(actual,sim.set[sim.set$Replicate==i,colnames(sim.set)==vals_colname]))
        divergence_vals<-c(divergence_vals,as.numeric(ktest$statistic))
        rec_vals<-c(rec_vals,r)
        select_vals<-c(select_vals,s)
      }
    }
  }
  
  divergence_df<-data.frame("Recombination"=rec_vals,"Selection"=select_vals,"D"=divergence_vals)
  divergence_df<-relevel_df(divergence_df)
  divergence_df<-divergence_df[order(divergence_df$Selection,divergence_df$Recombination),]
  divergence_df$Simulation<-paste(divergence_df$Selection,sub("\n"," ",tolower(divergence_df$Recombination)),sep=",\n")
  divergence_df$Simulation<-factor(divergence_df$Simulation,levels=unique(divergence_df$Simulation))
  divergence_df$Colour<-new_tree_pal[as.numeric(divergence_df$Simulation)]
  
  comparison_list<-list(c(unique(divergence_df$Simulation)[1],unique(divergence_df$Simulation)[2]),
                        c(unique(divergence_df$Simulation)[1],unique(divergence_df$Simulation)[3]),
                        c(unique(divergence_df$Simulation)[2],unique(divergence_df$Simulation)[3]),
                        c(unique(divergence_df$Simulation)[4],unique(divergence_df$Simulation)[5]),
                        c(unique(divergence_df$Simulation)[4],unique(divergence_df$Simulation)[6]),
                        c(unique(divergence_df$Simulation)[5],unique(divergence_df$Simulation)[6]),
                        c(unique(divergence_df$Simulation)[1],unique(divergence_df$Simulation)[4]),
                        c(unique(divergence_df$Simulation)[2],unique(divergence_df$Simulation)[5]),
                        c(unique(divergence_df$Simulation)[3],unique(divergence_df$Simulation)[6])
                        )
  
  vio_comp_plot<-ggplot(divergence_df,aes(x=Simulation,y=D,colour=Simulation,shape=Simulation,fill=Simulation))+
    geom_violin(alpha=0.1)+
    geom_point(position=position_jitter(width=0.075),size=0.5)+
    scale_y_continuous(trans="log")+
    scale_color_manual(values = new_tree_pal)+
    scale_fill_manual(values = new_tree_pal)+
    scale_shape_manual(values = new_tree_shapes)+
    stat_compare_means(
                       label = "p.format",
                       #label = format.pval("..p..",eps=0.001), # tough to alter the format with comparisons - see https://stackoverflow.com/questions/59494698/r-reformat-p-value-in-ggplot-using-stat-compare-means
                       method="wilcox",
                       comparisons=comparison_list,
                       symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                       bracket.size = 0.25, step.increase = 0.06,tip.length=0.01,
                       label.y.npc = -1,
                       size=2.5, hjust = -1,
                       hide.ns=TRUE,show.legend = FALSE)+
    ylab(expression("Kolmogorov-Smirnov"~italic(D)))+
    xlab("Simulation type")+
    theme(
      legend.position="bottom",
      axis.text.x = element_text(angle = 270,hjust = 0, vjust=0.5)
    )
  
  return(vio_comp_plot)
}

#########
# BEGIN #
#########

# data structures
rec_levels<-c("No transformation","Symmetrical transformation","Asymmetrical transformation")
select_levels<-c("Neutral","Multi-locus NFDS")
recs<-c("Asymmetrical transformation","Asymmetrical transformation","No transformation","No transformation","Symmetrical transformation","Symmetrical transformation")
selects<-c("Multi-locus NFDS","Neutral","Multi-locus NFDS","Neutral","Multi-locus NFDS","Neutral")
summaries<-c("Asymmetrical transformation, NFDS","Asymmetrical transformation, neutral","No transformation, NFDS","No transformation, neutral","Symmetrical transformation, NFDS","Symmetrical transformation, neutral")
nrep<-100

# colour palette
npg_pal<-pal_npg("nrc")(9)
my_fill_col<-npg_pal[6]
my_line_col<-npg_pal[1]
my_fill_low<-npg_pal[4]
my_fill_mid<-npg_pal[7]
my_fill_high<-npg_pal[5]
my_jitter_width=0.25/113
my_jitter_seed=runif(1)
my_scatter_alpha<-0.9
my_line_alpha<-0.1
my_geom_text_size<-3
my_facet_label_size<-8

# palettes
tree_pal<-pal_npg()(8)[c(2,5,6,1,4,8)]
tree_pal[4]<-"#E37DA3"
tree_pal[3]<-"#6A0DAD"
new_tree_pal<-c(tree_pal[c(1,3,5)],tree_pal[c(2,4,6)],"#000000")
new_tree_shapes<-c(1,2,4,1,2,4,19)


# read input data
eq.loci<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/mass.locusFreq")
colnames(eq.loci)<-c("Locus","Frequency")
snp.input<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/filtered.mass.snp.input",header=TRUE,comment.char = "")
snp.loci<-data.frame("Locus"=colnames(snp.input[6:ncol(snp.input)]),"Frequency"=as.numeric(colMeans(snp.input[snp.input$Time==0,6:ncol(snp.input)])))

# read input file
ma.data<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/filtered.mass.input",header=TRUE,comment.char = "")

# calculate genome sizes
genome_sizes<-as.data.frame(rowMeans(ma.data[,6:ncol(ma.data)]))
colnames(genome_sizes)[1]<-"Sizes"
genome_sizes.stats<-genome_sizes
genome_sizes.stats$Recombination<-"Genome data"
genome_sizes.stats$Selection<-"Genome data"
getDistStats(genome_sizes.stats,"Sizes")

# calculate SNP input data
snp.ma.data<-read.table(file="~/Documents/evoMLNFDS/rec_v_select/filtered.mass.snp.input",header=TRUE,comment.char = "")
snp.ma.data<-snp.ma.data[match(ma.data$Taxon,snp.ma.data$Taxon),]

# generate tree annotation
ma.tip.labels<-
  snp.ma.data %>%
    dplyr::select(Taxon,Time,Serotype,VT,SC) %>%
    dplyr::rename(taxa = Taxon) %>%
    mutate(new_SC = SC) %>%
    group_by(SC) %>%
    mutate(new_SC = ifelse(n() <= 10,0,SC)) %>%
    ungroup() %>%
    dplyr::select(-SC) %>%
    dplyr::rename(SC = new_SC)
ma.tip.labels$SC<-as.factor(ma.tip.labels$SC)
tree_palette<-scale_discrete_manual("colour",values=c("#000000",pal_d3("category20")(length(levels(ma.tip.labels$SC))-1)),labels = levels(ma.tip.labels$SC), drop = FALSE)

# calculate core tree
ma.tree<-make_nj_tree(as.data.frame(snp.ma.data[,c(1,6:ncol(snp.ma.data))]))
ma.tree.plot<-plot_nj_tree(ma.tree,anno = ma.tip.labels,point.size=2,text.size=5,expansion=1.5,offset.size=20,palette=tree_palette)
ggsave(ma.tree.plot,file="~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/ma_tree.png",width=20,height=20,units = "cm")

# calculate subsampled trees from data
ma.tree.stats<-tree_analysis(ma.tree)
ma.tree.stats$Recombination<-"Genomic data"
ma.tree.stats$Selection<-"Genomic data"

# calculate pairwise gene content distances
actual_gene_data<-as.matrix(ma.data[,6:ncol(ma.data)])
rownames(actual_gene_data)<-ma.data[,1]
actual_gene_dists<-rdist(actual_gene_data,metric="jaccard")
genome_distances<-melt(as.matrix(actual_gene_dists), varnames = c("row", "col"))
colnames(genome_distances)<-c("Taxon1","Taxon2","Distances")

# calculate pairwise SNP distances
actual_snp_data<-as.matrix(snp.ma.data[,6:ncol(snp.ma.data)])
rownames(actual_snp_data)<-snp.ma.data[,1]
actual_snp_dists<-rdist(actual_snp_data,metric="hamming")
snp_distances<-melt(as.matrix(actual_snp_dists), varnames = c("row", "col"))
colnames(snp_distances)<-c("Taxon1","Taxon2","SNPDistances")

# join distance matrices
genome_distances<-merge(genome_distances,snp_distances,by=c("Taxon1","Taxon2"))
genome_distances<-genome_distances[as.numeric(genome_distances$Taxon1) > as.numeric(genome_distances$Taxon2), ]

# calculate variation distribution for real data
strain_divide_intercept = 0.9
strain_divide_slope = -3

ma.dists<-ggplot(genome_distances,aes(x=SNPDistances,y=Distances))+geom_point(alpha = 5000/nrow(genome_distances),size=0.5)+
  geom_abline(colour=my_line_col, slope = strain_divide_slope, intercept = strain_divide_intercept, alpha=0.5, linetype = 2)+
  xlab("Pairwise Hamming SNP distances between genomes")+
  ylab("Pairwise binary Jaccard distances between genomes")+
  stat_cor(aes(label = paste("rho ==",bquote(.(..r..)))),colour="black",label.x=0.3,label.y=0.2,hjust = 0,method="spearman",cor.coef.name="rho")+
  geom_density_2d(data=genome_distances[genome_distances$Distances>=(strain_divide_intercept+strain_divide_slope*genome_distances$SNPDistances),], aes(x=SNPDistances, y=Distances), colour = pal_npg()(2)[2], alpha = 0.5, size = 0.5)+
  geom_density_2d(data=genome_distances[genome_distances$Distances<(strain_divide_intercept+strain_divide_slope*genome_distances$SNPDistances),], aes(x=SNPDistances, y=Distances), colour = my_fill_high, alpha = 0.5, size = 0.5)+
  theme_minimal()+theme(plot.margin = unit(c(0,0,0,0), "cm"),
                         panel.spacing = unit(0, "cm"),
                         panel.border = element_rect(color = "black", fill = NA, size = 0.25),
                         axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
                         axis.text.y = element_text(size = 8),
                         axis.title.x = element_text(size = 8, face = "bold"),
                         axis.title.y = element_text(size = 8, face = "bold"),
                         strip.text.x = element_text(size = my_facet_label_size),
                         strip.text.y = element_text(size = my_facet_label_size))
ggsave(ma.dists,file="~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/ma_dists.png",height=8.5,width=8.5,units = "cm")

# calculate distance distribution statistics
genome_distances.stats<-genome_distances
genome_distances.stats$Recombination<-"Genome data"
genome_distances.stats$Selection<-"Genome data"
getDistStats(genome_distances.stats,"Distances")
getDistStats(genome_distances.stats,"SNPDistances")

# calculate strain frequencies from genomic data
ma.strain.freqs<-get_strain_frequencies(genome_distances,slope=strain_divide_slope,intercept=strain_divide_intercept)

# now recalculate distances, excluding the MGE COGs
cog_annotation<-read.table(file="/Users/nicholascroucher/Documents/frequencyDependence/FOR_PUB/population_comparison/annotated_intermediate_cogs.tab",header=FALSE)
cog_annotation<-cog_annotation[cog_annotation[,1] %in% colnames(ma.data),]
nomge.ma.data<-ma.data[,-1*match(cog_annotation[cog_annotation[,2]=="pink",1],colnames(ma.data))]
nomge_gene_data<-as.matrix(nomge.ma.data[,6:ncol(nomge.ma.data)])
rownames(nomge_gene_data)<-ma.data[,1]
nomge_gene_dists<-rdist(nomge_gene_data,metric="jaccard")
nomge_genome_distances<-melt(as.matrix(nomge_gene_dists), varnames = c("row", "col"))
colnames(nomge_genome_distances)<-c("Taxon1","Taxon2","Distances")
nomge_genome_distances<-merge(nomge_genome_distances,snp_distances,by=c("Taxon1","Taxon2"))
nomge_genome_distances<-nomge_genome_distances[as.numeric(nomge_genome_distances$Taxon1) > as.numeric(nomge_genome_distances$Taxon2), ]

nomge.ma.dists<-ggplot(nomge_genome_distances,aes(x=SNPDistances,y=Distances))+geom_point(alpha=5000/nrow(nomge_genome_distances))+
  geom_abline(colour=my_line_col, slope = strain_divide_slope, intercept = strain_divide_intercept, alpha=0.5)+
  xlab("Pairwise Hamming SNP distances between genomes")+
  ylab("Pairwise binary Jaccard distances between genomes")+
  geom_density_2d(data=nomge_genome_distances[nomge_genome_distances$Distances>=(strain_divide_intercept+strain_divide_slope*nomge_genome_distances$SNPDistances),], aes(x=SNPDistances, y=Distances), colour = my_fill_high, alpha = 0.5)+
  geom_density_2d(data=nomge_genome_distances[nomge_genome_distances$Distances<(strain_divide_intercept+strain_divide_slope*nomge_genome_distances$SNPDistances),], aes(x=SNPDistances, y=Distances), colour = my_fill_low, alpha = 0.5)+
  theme_few()+
  theme(axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"))

ggsave(nomge.ma.dists,file="~/Documents/evoMLNFDS/rec_v_select/zero_time/figures/nomge.ma_dists.png",height=20,width=20,units = "cm")

# analysis with markers
marker_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_select_mig/rec_v_select_mig",
                   "~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_no_select_mig/rec_v_no_select_mig",
                   "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_select_mig/no_rec_v_select_mig",
                   "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_no_select_mig/no_rec_v_no_select_mig",
                   "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_select_mig/sym_rec_v_select_mig",
                   "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_no_select_mig/sym_rec_v_no_select_mig"
)

marker_plots<-process_simulation_data(prefixes=marker_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs,tree_data=ma.tree.stats)

# save figures
format_plots(marker_plots,"multistrain_with_migration_")

# run processing for no migration samples
no_migration_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_select_no_mig/rec_v_select_no_mig",
                         "~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_no_select_no_mig/rec_v_no_select_no_mig",
                         "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_select_no_mig/no_rec_v_select_no_mig",
                         "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_no_select_no_mig/no_rec_v_no_select_no_mig",
                         "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_select_no_mig/sym_rec_v_select_no_mig",
                         "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_no_select_no_mig/sym_rec_v_no_select_no_mig"
)

no_mig_plots<-process_simulation_data(prefixes=no_migration_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs,tree_data=ma.tree.stats)

# save figures
format_plots(no_mig_plots,"multistrain_without_migration_")

# weak selection analysis
# weak_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_weak_select_mig/rec_v_weak_select_mig",
#                  "~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_no_select_mig/rec_v_no_select_mig",
#                  "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_weak_select_mig/no_rec_v_weak_select_mig",
#                  "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_no_select_mig/no_rec_v_no_select_mig",
#                  "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_weak_select_mig/sym_rec_v_weak_select_mig",
#                  "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_no_select_mig/sym_rec_v_no_select_mig"
# )
# 
# weak_plots<-process_simulation_data(prefixes=weak_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs)
# 
# # save figures
# format_plots(weak_plots,"multistrain_with_weak_selection_with_migration_")

# weak selection no migration analysis
weak_no_mig_prefixes<-c("~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_weak_select_no_mig/rec_v_weak_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/rec_v_no_select_no_mig/rec_v_no_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_weak_select_no_mig/no_rec_v_weak_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_no_select_no_mig/no_rec_v_no_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_weak_select_no_mig/sym_rec_v_weak_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/sym_rec_v_no_select_no_mig/sym_rec_v_no_select_no_mig"
)

weak_no_mig_plots<-process_simulation_data(prefixes=weak_no_mig_prefixes,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs)

# save figures
format_plots(weak_no_mig_plots,"multistrain_with_weak_selection_without_migration_")

# saltational recombination - no migration
saltational_prefixes_no_mig<-c("~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_rec_v_select_no_mig/saltational_rec_v_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_rec_v_no_select_no_mig/saltational_rec_v_no_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_select_no_mig/no_rec_v_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_no_select_no_mig/no_rec_v_no_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_sym_rec_v_select_no_mig/saltational_sym_rec_v_select_no_mig",
                 "~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_sym_rec_v_no_select_no_mig/saltational_sym_rec_v_no_select_no_mig"
)

saltational_no_mig_plots<-process_simulation_data(prefixes=saltational_prefixes_no_mig,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs)

# save figures
format_plots(saltational_no_mig_plots,"multistrain_with_saltational_trans_without_migration_")
# 
# # saltational recombination - with migration
# saltational_prefixes_with_mig<-c("~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_rec_v_select_mig/saltational_rec_v_select_mig",
#                                "~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_rec_v_no_select_mig/saltational_rec_v_no_select_mig",
#                                "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_select_mig/no_rec_v_select_mig",
#                                "~/Documents/evoMLNFDS/rec_v_select/zero_time/no_rec_v_no_select_mig/no_rec_v_no_select_mig",
#                                "~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_sym_rec_v_select_mig/saltational_sym_rec_v_select_mig",
#                                "~/Documents/evoMLNFDS/rec_v_select/zero_time/saltational_sym_rec_v_no_select_mig/saltational_sym_rec_v_no_select_mig"
# )
# 
# saltational_with_mig_plots<-process_simulation_data(prefixes=saltational_prefixes_with_mig,recs=recs,selects=selects,summaries=summaries,eq.loci=eq.loci,real.data=ma.data,snp.loci=snp.loci,snp.calc=TRUE,strain_slope=strain_divide_slope,strain_intercept=strain_divide_intercept,tip.anno=ma.tip.labels,real_strains=ma.strain.freqs)
# 
# # save figures
# format_plots(saltational_with_mig_plots,"multistrain_with_saltational_trans_with_migration_")

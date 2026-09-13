args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Usage: Rscript 09j_es01_mapping_multiverse.R <vp_root>')
vp <- normalizePath(args[[1]], winslash='/', mustWork=TRUE)
root <- normalizePath(file.path(vp,'../..'), winslash='/', mustWork=TRUE)
out <- file.path(vp,'10_exploratory_hypoxia_sensitivity','attempt_20260910_02')
.libPaths(c(file.path(vp,'02_env','Rlib'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(fgsea);library(jsonlite);library(digest)})
targets <- c('06_hypoxia_mapping_qc.tsv','06_analysis_multiverse.tsv','06_specification_concordance.tsv')
if(any(file.exists(file.path(out,targets)))) stop('Refusing to overwrite ES01-06 outputs.')
scope <- fromJSON(file.path(out,'01_scope_freeze.json'),simplifyVector=FALSE)
if(!identical(scope$status,'FROZEN_EXPLORATORY_SPECIFICATION')) stop('Scope is not frozen.')
universe <- fread(file.path(vp,'04_prepared','gene_universe.tsv'))$gene
if(length(universe)!=1388L||anyDuplicated(universe)) stop('Fixed gene universe mismatch.')
gmt_path <- file.path(root,'analysis','data','processed','gene_sets','Hallmark.gmt')
freeze <- fromJSON(file.path(vp,'01_protocol','VP_G06_enrichment_rule_freeze_v1.json'),simplifyVector=FALSE)
if(!identical(digest(file=gmt_path,algo='sha256'),freeze$input_sha256[['analysis/data/processed/gene_sets/Hallmark.gmt']])) stop('HOLD_INPUT_DRIFT: Hallmark GMT.')
parts <- strsplit(readLines(gmt_path,warn=FALSE,encoding='UTF-8'),'\t',fixed=TRUE)
sets <- lapply(parts,function(x) unique(x[-c(1L,2L)])); names(sets) <- vapply(parts,`[[`,character(1),1L)
hypoxia <- sets[['HALLMARK_HYPOXIA']]
mapping <- data.table(gene_set='HALLMARK_HYPOXIA',original_member_count=length(hypoxia),network_universe_member_count=length(intersect(hypoxia,universe)),mapping_fraction=length(intersect(hypoxia,universe))/length(hypoxia),missing_members=paste(sort(setdiff(hypoxia,universe)),collapse=';'),duplicate_members=length(hypoxia)-length(unique(hypoxia)),mapping_issue_status='PASS')
eligible <- lapply(sets,intersect,y=universe); eligible <- eligible[lengths(eligible)>=15L & lengths(eligible)<=500L]
formal_eligible <- fread(file.path(vp,'07_enrichment','v1','gene_set_eligibility.tsv'))[database=='Hallmark'&status=='TESTED',pathway]
if(!identical(sort(names(eligible)),sort(formal_eligible))) stop('HOLD_POTENTIAL_IMPLEMENTATION_DEFECT: eligibility mismatch.')
rank_scores <- function(x){o<-order(-x,names(x));setNames((length(x)-seq_along(x)+1)/length(x),names(x)[o])}
fg <- function(r,s){set.seed(s);z<-suppressWarnings(fgseaMultilevel(eligible,r,minSize=15L,maxSize=500L,eps=0,scoreType='std',nproc=1L));z[,pathway:=as.character(pathway)];z[,le:=vapply(leadingEdge,paste,collapse=';',character(1))];z[pathway=='HALLMARK_HYPOXIA']}
registry <- fread(file.path(vp,'05_runs','run_registry.tsv'));setorder(registry,gene,condition,donor,seed)
formal <- fread(file.path(vp,'07_enrichment','v1','hypoxia_primary_summary.tsv'))
rows <- list(); counter <- as.integer(freeze$method$random_seed)
for(gname in sort(unique(registry$gene))){
 cname <- unique(registry[gene==gname]$condition)
 base <- formal[target_gene==gname&condition==cname]
 rows[[length(rows)+1]] <- data.table(specification='S0',evidence_label='FORMAL_BASELINE',target_gene=gname,condition=cname,ES=base$ES,NES=base$NES,p_value=base$pval,FDR_within_Hallmark=base$padj,leading_edge=base$leadingEdge,result_status='FORMAL_BASELINE_READ_ONLY')
 for(spec in c('S1','S2')){
  donor_scores <- list()
  for(dname in sort(unique(registry[gene==gname&condition==cname]$donor))){
   rr <- registry[gene==gname&condition==cname&donor==dname]
   metric <- if(spec=='S1') 'distance' else 'Z'
   ranks <- sapply(seq_len(nrow(rr)),function(i){d<-fread(file.path(rr$output_attempt[i],'differential_regulation.tsv'),select=c('gene',metric));v<-setNames(d[[metric]],d$gene);rank_scores(v)})
   donor_scores[[dname]] <- if(spec=='S1') rowMeans(ranks) else apply(ranks,1,median)
  }
  score <- if(spec=='S1') rowMeans(do.call(cbind,donor_scores)) else apply(do.call(cbind,donor_scores),1,median)
  counter <- counter+1L; z <- fg(rank_scores(score),counter)
  rows[[length(rows)+1]] <- data.table(specification=spec,evidence_label='EXPLORATORY',target_gene=gname,condition=cname,ES=z$ES,NES=z$NES,p_value=z$pval,FDR_within_Hallmark=z$padj,leading_edge=z$le,result_status='EXPLORATORY')
 }
}
donor <- fread(file.path(vp,'07_enrichment','v1','hypoxia_donor_enrichment.tsv'))[pathway=='HALLMARK_HYPOXIA']
for(gname in sort(unique(donor$target_gene))){z<-donor[target_gene==gname];rows[[length(rows)+1]]<-data.table(specification='S3',evidence_label='EXPLORATORY',target_gene=gname,condition=z$condition[1],ES=NA_real_,NES=median(z$NES),p_value=NA_real_,FDR_within_Hallmark=NA_real_,leading_edge='NOT_APPLICABLE_DONOR_DESCRIPTIVE',result_status='EXPLORATORY_DONOR_DESCRIPTIVE')}
mv <- rbindlist(rows); setorder(mv,target_gene,specification)
con <- rbindlist(lapply(unique(mv$target_gene), function(g) {
  z <- mv[target_gene == g & specification %in% c('S0','S1','S2')]
  pairs <- combn(seq_len(nrow(z)), 2L)
  rbindlist(lapply(seq_len(ncol(pairs)), function(i) data.table(
    target_gene=g, specification_a=z$specification[pairs[1,i]], specification_b=z$specification[pairs[2,i]],
    NES_difference=abs(z$NES[pairs[1,i]]-z$NES[pairs[2,i]]),
    direction_concordant=sign(z$NES[pairs[1,i]])==sign(z$NES[pairs[2,i]])
  )))
}))
fwrite(mapping,file.path(out,targets[1]),sep='\t');fwrite(mv,file.path(out,targets[2]),sep='\t');fwrite(con,file.path(out,targets[3]),sep='\t')
cat(toJSON(list(gate='VP-ES01-06',status='PASS_COMPUTED_PENDING_AUDIT',mapping_rows=nrow(mapping),multiverse_rows=nrow(mv)),auto_unbox=TRUE),'\n')

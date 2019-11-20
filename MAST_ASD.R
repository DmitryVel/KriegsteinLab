fc<-function(data.sel,meta){
library("lme4")
library("lmtest");
out=matrix(,nrow=0,ncol=2);
exp=cbind(data.sel,meta);
colnames(exp)[1]<-"expression";
exp1=2^exp[exp$diagnosis=="Control",1]-1;
exp2=2^exp[exp$diagnosis=="ASD",1]-1;
exp.control=mean(exp1);
exp.asd=mean(exp2);
non.zero.cntr=length(exp1[exp1>0])/length(exp1);
non.zero.asd=length(exp2[exp2>0])/length(exp2);
samples.cntr=meta[meta$diagnosis=="Control",]$sample;
samples.asd=meta[meta$diagnosis=="ASD",]$sample;
sample.exp.cntr=mean(aggregate(exp1,by=list(as.character(samples.cntr)),FUN=mean)[,2]);
sample.exp.asd=mean(aggregate(exp2,by=list(as.character(samples.asd)),FUN=mean)[,2]);
FC=log(exp.asd/exp.control,2);
FC.sample=log(sample.exp.asd/sample.exp.cntr,2);
return(c(exp.control,exp.asd,non.zero.cntr,non.zero.asd,FC,FC.sample));
}
TEST.MAST.ASD<-function(data,meta,clusters="ALL",percentage=0.1,column=3,column2=12,group1="Control",group2="ASD",FDRTHRESHOLD=.05,FCTHRESHOLD=0,ENS_ref="GRCh38",cores=8){
  if(clusters!="ALL"){
  meta=meta[meta[,column]%in%clusters,];
  data=data[,rownames(meta)];
  }
  for(test.cluster in unique(meta[,column])){
    test.tsne=meta[meta[,column]==test.cluster,];
    test.data=data[,rownames(test.tsne)];
    rows=rowSums(as.matrix(test.data)> 0) > ncol(test.data)*percentage;
    test.data = test.data[rows,];
    print(paste("Testing for differential expression between",group1,"and",group2,"in cluster",test.cluster,sep=" "));
    print(paste("Testing",nrow(test.data),"genes"));
    sca=FromMatrix(as.matrix(test.data),as.data.frame(colnames(test.data)),as.data.frame(rownames(test.data)));
    groups=as.character(test.tsne[,column2]);
    cond<-factor(groups);
    cond<-relevel(cond,group1);
    cdr = colSums(assay(sca)>0);
    age=test.tsne$age;
    sex=as.character(test.tsne$sex);
    RIN=test.tsne$RIN;
    PMI=test.tsne$PMI;
    region=as.character(test.tsne$region);
    batch=as.character(test.tsne$batch);
    Capbatch=as.character(test.tsne$Capbatch);
    Seqbatch=as.character(test.tsne$Seqbatch);
    ind=as.character(test.tsne$individual);
    ribo_perc=test.tsne$ribo_perc;
    mito_perc=test.tsne$mito_perc;
    mito_nucl_perc=test.tsne$mito_nucl_perc;
    colData(sca)$diagnosis<-cond;
    colData(sca)$cngeneson <- scale(cdr);
    colData(sca)$age<-scale(age);
    colData(sca)$sex<-sex;
    colData(sca)$RIN<-scale(RIN);
    colData(sca)$PMI<-scale(PMI);
    colData(sca)$region<-region;
    colData(sca)$batch<-batch;
    colData(sca)$Capbatch<-Capbatch;
    colData(sca)$Seqbatch<-Seqbatch;
    colData(sca)$batch<-batch;
    colData(sca)$ind<-ind;
    colData(sca)$ribo_perc<-scale(ribo_perc);
    colData(sca)$mito_perc<-scale(mito_perc);
    colData(sca)$mito_nucl_perc<-scale(mito_nucl_perc);
    colData(sca)$PMI<-scale(PMI);
    form=as.formula("~diagnosis + (1|ind) + Capbatch + Seqbatch + region + cngeneson + age + sex + RIN + PMI + ribo_perc");
    zlmCond = zlm(form, sca, method = "glmer", ebayes = F, silent=T);
    summaryCond <- summary(zlmCond, doLRT='diagnosisASD');
    summaryDt <- summaryCond$datatable;
    save(summaryDt,file=paste("LMM_",test.cluster,".R",sep=""));
    ENS_ref_filename=paste("/home/dvelmeshev/ref/","GRCh38","_genes_ENS.txt",sep="");
    ENS_ref=read.table(ENS_ref_filename,header=TRUE,sep="\t");
    fcHurdle <- merge(summaryDt[contrast=='diagnosisASD' & component=='H',.(primerid,`Pr(>Chisq)`)],summaryDt[contrast=='diagnosisASD' & component=='logFC', .(primerid, coef,ci.hi, ci.lo)], by='primerid');
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')];
    fcHurdleSig <- merge(fcHurdle[fdr<FDRTHRESHOLD & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)),by='primerid');
    setorder(fcHurdleSig, fdr);
    gene.IDs<-as.character(fcHurdleSig$primerid);
    FCs<-as.character(fcHurdleSig$coef);
    FDRs<-as.character(fcHurdleSig$fdr);
    gene.names<-c();
    gene.biotypes<-c();
    for(gene.ID in gene.IDs){
      if(gene.ID%in%as.character(as.character(ENS_ref[,1]))){
        gene.name=as.character(ENS_ref[ENS_ref[,1]==gene.ID,2]);
        gene.biotype=as.character(ENS_ref[ENS_ref[,1]==gene.ID,3]);
      }
    else{
    gene.name="NA";
    gene.biotype="NA"
    }
    gene.names=append(gene.names,gene.name);
    gene.biotypes=append(gene.biotypes,gene.biotype);
    }
    out=cbind(gene.IDs,gene.names,gene.biotypes,FCs,FDRs);
    write.table(out,paste("LMM_",test.cluster,".txt",sep=""),sep="\t")
    zlmCond = zlm(~diagnosis, sca);
    summaryCond <- summary(zlmCond, doLRT='diagnosisASD');
    summaryDt <- summaryCond$datatable;
    ENS_ref_filename=paste("/home/dvelmeshev/ref/","GRCh38","_genes_ENS.txt",sep="");
    ENS_ref=read.table(ENS_ref_filename,header=TRUE,sep="\t");
    fcHurdle <- merge(summaryDt[contrast=='diagnosisASD' & component=='H',.(primerid,`Pr(>Chisq)`)],summaryDt[contrast=='diagnosisASD' & component=='logFC', .(primerid, coef,ci.hi, ci.lo)], by='primerid');
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')];
    fcHurdleSig <- merge(fcHurdle[fdr<=1 & abs(coef)>=0], as.data.table(mcols(sca)),by='primerid');
    setorder(fcHurdleSig, fdr);
    gene.IDs<-as.character(fcHurdleSig$primerid);
    FCs<-as.character(fcHurdleSig$coef);
    gene.names<-c();
    gene.biotypes<-c();
    for(gene.ID in gene.IDs){
      if(gene.ID%in%as.character(as.character(ENS_ref[,1]))){
        gene.name=as.character(ENS_ref[ENS_ref[,1]==gene.ID,2]);
        gene.biotype=as.character(ENS_ref[ENS_ref[,1]==gene.ID,3]);
      }
    else{
    gene.name="NA";
    gene.biotype="NA"
    }
    gene.names=append(gene.names,gene.name);
    gene.biotypes=append(gene.biotypes,gene.biotype);
    }
    out=cbind(gene.IDs,gene.names,gene.biotypes,FCs);
    write.table(out,paste(test.cluster,"MAST_FC.txt",sep="_"),sep="\t")
    cl <- makeCluster(cores);
    out=pbapply(cl=cl,test.data,1,fc,meta=test.tsne);
    res=cbind(rownames(test.data),t(out));
    colnames(res)<-c("gene_ID","exp_control","exp_asd","nz_control","nz_asd","FC","FC_sample");
    write.table(res,paste(test.cluster,"_FC.txt",sep=""),sep="\t",quote=F,row.names=F);
  }
}

#' Univariate analysis
#'
#' @param dataset dataset of bins or quantifications.
#' @param metadata Metadata to be used during the univariate analysis.
#'
#' @return p values from dataset of bins or quantificaitons
#' @export p_values
#'
#' @examples
#' setwd(paste(system.file(package = "Dolphin"),"extdata",sep='/'))
#' imported_data=import_data("Parameters_MTBLS242_15spectra_5groups.csv")
#' pval=p_values(imported_data$dataset,imported_data$Metadata)


p_values=function(dataset,metadata) {


types=unique(unlist(metadata[,-1]))[which(unique(unlist(metadata[,-1]))>=0)]
types2=abs(unique(unlist(metadata[,-1]))[which(unique(unlist(metadata[,-1]))<0)])
datasetlist=list()
paireddata=F
if (identical(types,types2)) {
  print ('Analysis of differences')


for (i in seq_along(types)) {
ind1=which(metadata[,-1] ==types[i])%%nrow(metadata)
ind2=which(metadata[,-1] ==-types[i])%%nrow(metadata)
ind1[ind1==0]=ind2[ind2==0]=nrow(metadata)
datasetlist[[i]]=dataset[ind2,,drop=F]-dataset[ind1,,drop=F]
# if (all(metadata[ind1,1]==metadata[ind2,1])==F)   paireddata=F
paireddata=F

}
} else {
  print ('Analysis of groups')
  paireddata=T
  lal=metadata[which(metadata[,-1] ==types[1])%%nrow(metadata),1]
  for (i in 1:length(types)) {
  ind1=which(metadata[,-1] ==types[i])%%nrow(metadata)
  ind1[ind1==0]=nrow(metadata)
  if (!identical(metadata[ind1,1],lal)) paireddata=F
  datasetlist[[i]]=dataset[ind1,,drop=F]
  lal=metadata[ind1,1]
  }
}
if (paireddata==F) {
  print('Unpaired data')
} else {
  print('Paired data')
}

if (length(datasetlist)==2) {
tt=matrix(NA,length(datasetlist),dim(dataset)[2])
for (ind in 1:length(datasetlist)) {
  for (k in 1:dim(dataset)[2]) {
    tt[ind,k]=tryCatch(shapiro.test(datasetlist[[ind]][,k])$p.value,error=function(e) NA)
  }

}
p_value=rep(NA,dim(dataset)[2])
for (k in 1:dim(dataset)[2]) {
  # if (!any(is.na(dataset[,k]))) {
  if (!any(tt[,k]<0.05,na.rm=T)) {

    p_value[k]=tryCatch(wilcox.test(datasetlist[[1]][,k],datasetlist[[2]][,k],paired=paireddata)$p.value,error=function(e) NA)
  } else {
    p_value[k]=tryCatch(t.test(datasetlist[[1]][,k],datasetlist[[2]][,k],paired=paireddata,var.equal=F)$p.value,error=function(e) NA)
  }

  # }
}

} else {

  # datasetlist=array(unlist(datasetlist), dim = c(nrow(datasetlist[[1]]), ncol(datasetlist[[1]]), length(datasetlist)))
  # tt=rep(NA,dim(dataset)[2])
  #   for (k in 1:dim(dataset)[2]) {
  #     fa=tryCatch(bartlett.test(as.data.frame(datasetlist[,k,]))$p.value,error=function(e) NA)
  #     fa2=tryCatch(fligner.test(as.data.frame(datasetlist[,k,]))$p.value,error=function(e) NA)
  #     tt[k]=suppressWarnings(min(c(fa,fa2),na.rm=T))
  #   }
  tt=rep(NA,dim(dataset)[2])
  for (k in 1:dim(dataset)[2]) {
    # fa=tryCatch(bartlett.test(as.data.frame(datasetlist[,k,]))$p.value,error=function(e) NA)

    tt[k]=tryCatch(fligner.test(dataset[,k],metadata[,2])$p.value,error=function(e) 0)
    #     tt[k]=suppressWarnings(min(c(fa,fa2),na.rm=T))
  }
  tt[is.na(tt)]=0
  p_value=rep(NA,dim(dataset)[2])
  for (k in 1:dim(dataset)[2]) {
    # if (!any(is.na(dataset[,k]))) {
    if (tt[k]<0.05) {

      # if (paireddata==T) {
      #   p_value[k]=withCallingHandlers({ friedman.test(datasetlist[,k,])$p.value}, warning = function(w) {})
      # } else {
      # p_value[k]=withCallingHandlers({ kruskal.test(as.data.frame(datasetlist[,k,]))$p.value}, warning = function(w) {})
      p_value[k]=tryCatch(kruskal.test(dataset[,k] ~ metadata[,2])$p.value,error=function(e) NA)
      # }
    } else {
      if (paireddata==T) {
        # no=stack(as.data.frame(datasetlist[,k,]))
        # no$subject = rep(rownames(as.data.frame(datasetlist[,k,])), length(types))
        # p_value[k]=tryCatch(summary(aov(values ~ ind + Error(subject/ind), data=no))[[2]][[1]][1,5],error=function(e) NA)
        p_value[k]=tryCatch(summary(aov(dataset[,k] ~ metadata[,2]+ Error(metadata[,1]/metadata[,2])))[[3]][[1]]$`Pr(>F)`[1],error=function(e) NA)


      } else {
        p_value[k]=tryCatch(summary(aov(dataset[,k] ~ metadata[,2]))[[1]]$`Pr(>F)`[1],error=function(e) NA)
      # no=data.frame(y=as.vector(as.matrix(datasetlist[,k,])),group=sort(rep(1:dim(datasetlist)[3],dim(datasetlist)[1])))
      # p_value[k]=tryCatch(anova(lm(y ~ group,no))[1,5],error=function(e) NA)

        }
    }

    # }
  }

}
p_value_final=rep(1,length(p_value))
p_value_final[which(!is.na(p_value))]=round(t(as.matrix(p.adjust(p_value[which(!is.na(p_value))],method="none"))),3)
names(p_value_final)=colnames(dataset)
return(p_value_final)
}

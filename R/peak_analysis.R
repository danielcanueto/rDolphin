peak_analysis=function(dataset,ppm,freq,export_path,metadata,repository,originaldataset) {
  print('Processing. Talk some gossip, meanwhile.')
  
  
  
  if (ppm[1]<ppm[2]) {
    setupRSPA(ppm)
    
  } else {
    setupRSPA(-ppm)
    
  }
  peakParam$ampThr=quantile(dataset,0.65,na.rm=T)
  
# refdataset<-dataset[suppressWarnings(selectRefSp(dataset,recursion$step)),]
refdataset<-apply(dataset,2,function(x)median(x,na.rm=T))
refSegments<- segmentateSp(refdataset, peakParam)



peak_ppm=c()
for (i in 1:length(refSegments$Peaks)) peak_ppm=c(peak_ppm,ppm[unlist(refSegments$Peaks[[i]][1])])

peak_info=matrix(NA,0,11)
for (i in 1:length(refSegments$Peaks)) peak_info=rbind(peak_info,refSegments$Peaks[[i]])
peak_ppm=peak_ppm[which(peak_info[,11]-peak_info[,10]>1)]
peak_info=peak_info[which(peak_info[,11]-peak_info[,10]>1),]



peak_shape_corr=matrix(NA,dim(dataset)[1],dim(peak_info)[1])
for (j in 1:dim(peak_info)[1]) {
  med=try(apply(dataset[,peak_info[j,10]:peak_info[j,11],drop=F],2,median),silent=T)
  med2=try(dataset[,peak_info[j,10]:peak_info[j,11]],silent=T)
  peak_shape_corr[,j]=try(cor(t(med2),med),silent=T)
}
peak_shape_corr[is.na(peak_shape_corr)]=0
colnames(peak_shape_corr)=peak_ppm

valid_peak_ind=apply(peak_shape_corr,2,function(x)quantile(x,0.5))
peak_info=peak_info[valid_peak_ind>0.8,]
spectra_lag = matrix(NA, nrow(dataset),nrow(peak_info))
spectra_position = matrix(NA, nrow(dataset),nrow(peak_info))

for (i in 1:nrow(dataset)) {
  for (j in 1:nrow(peak_info)) {
    d <-
      ccf(
        originaldataset[i, (-20:20)+peak_info$maxPos[j]],
        dataset[i, (-20:20)+peak_info$maxPos[j]],
        type = 'covariance',
        plot = FALSE)
    spectra_lag[i,j] = d$lag[which.max(d$acf)]
    spectra_position[i,j]=ppm[peak_info$maxPos[j]-d$lag[which.max(d$acf)]]

    
  }}
colnames(spectra_lag)=ppm[peak_info$maxPos]

peak_intensity=peak_quantification=matrix(NA,dim(dataset)[1],dim(peak_info)[1])
colnames(peak_quantification)=colnames(peak_intensity)=peak_ppm[valid_peak_ind>0.8]
for (i in 1:dim(dataset)[1]) {
  for (j in 1:dim(peak_info)[1]) {
    # peak_quantification[i,j]=sum(dataset[i,peak_info[j,10]:peak_info[j,11]])-peak_info[j,9]*length(dataset[i,peak_info[j,10]:peak_info[j,11]])
    peak_quantification[i,j]=sum(dataset[i,peak_info[j,10]:peak_info[j,11]]-seq(dataset[i,peak_info[j,10]],dataset[i,peak_info[j,11]],length.out = (peak_info[j,11]-peak_info[j,10])+1))
    # peak_intensity[i,j]=max(dataset[i,peak_info[j,10]:peak_info[j,11]]-seq(dataset[i,peak_info[j,10]],dataset[i,peak_info[j,11]],length.out = (peak_info[j,11]-peak_info[j,10])+1))
    peak_intensity[i,j]=max(dataset[i,peak_info[j,10]:peak_info[j,11]])
  }}

half_band_width=matrix(NA,nrow(dataset),nrow(peak_info))
for (i in 1:nrow(peak_info)) {
  for (j in 1:nrow(dataset)){
    
    fa=dataset[j,peak_info[i,10]:peak_info[i,11]]-rep(min(dataset[j,peak_info[i,10]:peak_info[i,11]]),peak_info[i,11]-peak_info[i,10]+1)
    nn=approx(ppm[peak_info[i,10]:peak_info[i,11]],fa/max(max(fa,1e-10)),seq(ppm[peak_info[i,10]],ppm[peak_info[i,11]],length.out = 10000))
    ss=diff(nn$x[which(diff(sign(nn$y-0.5))!=0)])*-freq
    if (length(ss)==1) half_band_width[j,i]=ss/2
  }}

aa=apply(half_band_width,1,function(x)median(x,na.rm=T))/median(apply(half_band_width,1,function(x)median(x,na.rm=T)))
bb=apply(half_band_width,2,function(x)median(x,na.rm=T))
expected_width=matrix(NA,nrow(half_band_width),ncol(half_band_width))
for (i in 1:nrow(expected_width)) {
  for (j in 1:ncol(expected_width)) {
    expected_width[i,j]=bb[j]*aa[i]
  }}

peak_quantification[peak_shape_corr[,valid_peak_ind>0.8]<0.8]=NA
# peak_intensity[peak_shape_corr[,valid_peak_ind>0.8]<0.8]=NA


corr_matrix=round(cor(peak_quantification,use='pairwise.complete.obs',method='spearman'),2)

threshold_corr_matrix=corr_matrix
threshold_corr_matrix[threshold_corr_matrix<0.8]=0

correlated_signals_indexes=threshold_corr_matrix[which(colSums(threshold_corr_matrix,na.rm=T)>1),which(colSums(threshold_corr_matrix)>1)]
secure_multiplets=list()
k=1
for (i in 1:dim(correlated_signals_indexes)[1]) {
  
  corr_ind=which(correlated_signals_indexes[i,]>0.75)
  cc=c()
  
  for ( j in 1:length(corr_ind)) {
    corr_ind2=which(correlated_signals_indexes[corr_ind[j],]>=0.75)
    if (identical(corr_ind2,corr_ind)) {
      cc=c(cc,corr_ind[j])   
    } else {
      cc=0
      break
    }
  }
  
  if (cc!=0) {
    secure_multiplets[[k]]=corr_ind
    k=k+1
  }
  
}
secure_multiplets=unique(secure_multiplets)
ab2=unlist(unique(secure_multiplets))
if (length(secure_multiplets)>0) {
for (i in 1:length(secure_multiplets)) {
  secure_multiplets[[i]]=which(colnames(threshold_corr_matrix) %in% names(secure_multiplets[[i]])==T)
  names(secure_multiplets[[i]])=colnames(threshold_corr_matrix)[secure_multiplets[[i]]]
}}

clustering_corr_matrix=corr_matrix[which(colSums(threshold_corr_matrix)>1),which(colSums(threshold_corr_matrix)>1)]
if (length(ab2)>0) clustering_corr_matrix=clustering_corr_matrix[-ab2,-ab2]
clustering_quality_indicator=matrix(NA,3,9)
for (i in seq(0.1,0.9,0.1)) {
  for (j in 1:3) {
    clustering_quality_indicator[j,i*10]=apcluster(negDistMat(r=j), clustering_corr_matrix,q=i)@netsim
  }}
dummy=apply(clustering_quality_indicator,1,function(x) diff(x)/mean(diff(x)))
optimal_clustering_values=(which(dummy==max(dummy), arr.ind = TRUE) +c(1,0))/c(10,1)
q_clustering_value=which.max(diff(clustering_quality_indicator))+1
signals_clusters <- apcluster(negDistMat(r=optimal_clustering_values[2]), clustering_corr_matrix,q=optimal_clustering_values[1])@clusters

for (i in 1:length(signals_clusters)) {
  signals_clusters[[i]]=which(colnames(threshold_corr_matrix) %in% names(signals_clusters[[i]])==T)
  names(signals_clusters[[i]])=colnames(threshold_corr_matrix)[signals_clusters[[i]]]
}

signals_list=unique(append(secure_multiplets,signals_clusters))
for (i in 1:length(which(colSums(threshold_corr_matrix)<=1)))  signals_list[[length(signals_list)+1]]=which(colSums(threshold_corr_matrix)<=1)[i]

peak_intensity_median=apply(peak_intensity,2,function(x)median(x,na.rm=T))
expected_width_median=apply(expected_width,2,function(x)median(x,na.rm=T))


CV <- function(x){
(sd(x)/mean(x))
}

p_value_final=p_values(peak_quantification,metadata)
# p_value_final=p_value_final[sort(as.numeric(names(p_value_final)),index.return=T)$ix]
signals_position=matrix(NA,nrow(dataset),0)
signals_width=matrix(NA,nrow(dataset),0)
signals_intensity=matrix(NA,nrow(dataset),0)


ROI_profile_suggestion=matrix(NA,0,10)
for (i in 1:length(signals_list)) {
  so=as.numeric(names(signals_list[[i]]))
  tt=which(abs(diff(so))>0.02)
  
  if(length(tt)>0) {
    dr=cbind(c(1,tt+1),c(tt,length(so)))
    ds=dr[which(dr[,2]-dr[,1]==0),]
    
    cvb=which(dr[,2]-dr[,1]>0)
    if (length(cvb)>0) {
      
      for (j in 1:length(cvb)) {
        comp=peak_intensity_median[signals_list[[i]][dr[cvb[j],1]:dr[cvb[j],2]]]
        vb=ifelse(is.na(CV(comp)),0,CV(comp))
        if (vb>0.1){
        fff=length(which(t(diff(t(spectra_lag[,signals_list[[i]][dr[cvb[j],1]:dr[cvb[j],2]]])))!=0  ))
        if (length(fff)<0.1*nrow(dataset)) {
          ds=rbind(ds,dr[cvb[j],])
          # dt=
          
        } else {
  ds=rbind(ds,cbind(dr[cvb[j],1]:dr[cvb[j],2],dr[cvb[j],1]:dr[cvb[j],2]))
        }
          # dr[cvb[j],which.max(dr[cvb[j],])]
        } else {
          ds=rbind(ds,dr[cvb[j],])
        }
      }
    } 
    # else {
    #   ds=dr
    # }
  }else {
    ds=t(as.matrix(c(1,length(so))))
  }
  signals_patterns=ds
  # signals_patterns=cbind(1:length(signals_list[[i]]),1:length(signals_list[[i]]))
  signals_intensities=c()
  # for (k in 1:dim(signals_patterns)[1]) signals_intensities=c(signals_intensities,mean(peak_intensity_median[signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]]))
  for (k in 1:dim(signals_patterns)[1]) signals_intensities=c(signals_intensities,peak_intensity_median[signals_list[[i]][signals_patterns[k,1]]])
  
  signals_intensities2=c()
  # for (k in 1:dim(signals_patterns)[1]) signals_intensities=c(signals_intensities,mean(peak_intensity_median[signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]]))
  for (k in 1:dim(signals_patterns)[1]) signals_intensities2=c(signals_intensities2,peak_intensity_median[signals_list[[i]][signals_patterns[k,2]]])
  roof_effect = (signals_intensities2/signals_intensities)-1
  for (k in 1:dim(signals_patterns)[1]) signals_position=cbind(signals_position,apply(spectra_position[,signals_list[[i]][signals_patterns[k,]]],1,mean))
  for (k in 1:dim(signals_patterns)[1]) signals_width=cbind(signals_width,apply(expected_width[,signals_list[[i]][signals_patterns[k,]]],1,mean))
  for (k in 1:dim(signals_patterns)[1]) signals_intensity=cbind(signals_intensity,apply(peak_intensity[,signals_list[[i]][signals_patterns[k,]]],1,mean))
  
  for (k in 1:dim(signals_patterns)[1]) {
    p_v=mean(p_value_final[signals_list[[i]]][signals_patterns[k,1]:signals_patterns[k,2]])
    # ROI_profile_suggestion=rbind(ROI_profile_suggestion,c(paste(i,k,sep='_'),mean(as.numeric(names(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]))),1,1,length(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]),abs((as.numeric(names(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]))[2]-as.numeric(names(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]))[1])*500.3),0,0.001,signals_intensities[k]/max(signals_intensities)))
    ROI_profile_suggestion=rbind(ROI_profile_suggestion,c(i,mean(as.numeric(names(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]))),0.7,k,length(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]),abs((as.numeric(names(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]))[2]-as.numeric(names(signals_list[[i]][signals_patterns[k,1]:signals_patterns[k,2]]))[1])*freq),roof_effect[k],0.001,signals_intensities[k]/max(signals_intensities),p_v))
    
  }
}



ROI_profile_suggestion[is.na(ROI_profile_suggestion)]=0 
ROI_profile_suggestion=as.data.frame(ROI_profile_suggestion[sort(ROI_profile_suggestion[,2],decreasing=T,index.return=T)$ix,],stringsAsFactors = F)
ROI_profile_suggestion[,-1]=lapply(ROI_profile_suggestion[,-1],function(x) as.numeric(x))

signals_diff=c(0,which(abs(diff(ROI_profile_suggestion[,2]))>0.03),dim(ROI_profile_suggestion)[1])
ROI_limits_suggestion=matrix(NA,dim(ROI_profile_suggestion)[1],2)
for (i in 1:dim(ROI_profile_suggestion)[1]) {
  relevant_signals_diff=which(signals_diff %in% i ==T)
  if (length(relevant_signals_diff)>0) {
    ROI_limits_suggestion[(signals_diff[relevant_signals_diff-1]+1):signals_diff[relevant_signals_diff],]=matrix(rep(c(ROI_profile_suggestion[(signals_diff[relevant_signals_diff-1]+1),2]+0.02,ROI_profile_suggestion[signals_diff[relevant_signals_diff],2]-0.02),2),length((signals_diff[relevant_signals_diff-1]+1):signals_diff[relevant_signals_diff]),2,byrow=T)
  }
}

ROI_profile_suggestion=cbind(round(ROI_limits_suggestion,3),rep('Baseline Fitting',dim(ROI_profile_suggestion)[1]),ROI_profile_suggestion)
colnames(ROI_profile_suggestion)=c('ROI_left','ROI_right','Q.Mode','Signal','Position..ppm.','Width','Q.Signal','Multiplicity','J.coupling..Hz.','Roof.effect','Shift.tolerance','Intensity','p_value')

suggested_metabolites=matrix(NA,dim(ROI_profile_suggestion)[1],5)
colnames(suggested_metabolites)=paste('Suggested Metabolite',1:5)
for (i in 1:dim(suggested_metabolites)[1]) {
  a=which(abs(repository[,3]-ROI_profile_suggestion[i,5])<0.05)
  b=which(as.numeric(repository[a,7])==ROI_profile_suggestion[i,8])
  d=a[sort(abs(repository[a,3]-ROI_profile_suggestion[i,5]),index.return=T)$ix[1:5]]
  e=c(setdiff(a[b],d),d)
  e=unique(e)[1:5]
  suggested_metabolites[i,]=paste(repository[e,1],' (', round(abs(repository[e,3]-ROI_profile_suggestion[i,5])[1:5],3),')',sep='')
}

# suggested_metabolites[i,]=paste(repository[!is.na(repository[,5]),][sort(abs(repository[,5]-ROI_profile_suggestion[i,5]),index.return=T)$ix[1:5],1],' (', round(sort(abs(repository[,5]-ROI_profile_suggestion[i,5]),index.return=T)$x[1:5],3),')',sep='')
ROI_profile_suggestion=cbind(ROI_profile_suggestion,suggested_metabolites)

try(write.csv(ROI_profile_suggestion,paste(export_path,'ROI_profile_suggestion.csv',sep='/'),row.names = F),silent=T)
try(write.csv(peak_quantification,paste(export_path,'peak_quantification.csv',sep='/'),row.names = F),silent=T)
print('Done! Look on your export folder, you should have a new csv with ROI profiles suggestions and a new CSV with quantifications of integrated peaks.')
# dumy=list(signals_intensity=signals_intensity,signals_width=signals_width,signals_position=signals_position,spectra_lag=spectra_lag)
# return(dumy)
}

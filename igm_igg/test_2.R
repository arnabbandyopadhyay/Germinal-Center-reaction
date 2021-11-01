library(gridExtra)
library(ggplot2)
library(stringr)
library(reshape2)
library("readxl")




files <-list.files(pattern = "output_mem_cls_swt_test_")
dat<-read_excel(files[1], col_types = "numeric")
col<-2
dd<-dat[1:9900,col]
for (i in 2:9){
  #  dat2<-read_excel(paste0('antg_tfh/',files[i]), col_types = "numeric")
  dat2<-read_excel(files[i], col_types = "numeric")
  dd2<-dat2[1:9900,col]
  dd<-cbind(dd,dd2)
}

mean<-rowMeans(dd)
sd<-apply(dd, 1, sd)
nd<-data.frame(time=0.01*(1:9900), mean=mean, sd=sd)


p<- ggplot(nd, aes(x=time, y=mean)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))












################################

dat<-read_excel("output_mem_cls_swt_40.xlsx", col_types = "numeric")

#library("ggplot2")

# basic plots, DZ, LZ, DZ/LZ ratio, avg div and affinity maturation
colors <- c("DZ+LZ" = "red", "DZ" = "blue", "LZ" = "green")

my_gc<-ggplot(dat, aes(x = time/24)) +
  geom_line(aes(y = DZ+LZ, color = "DZ+LZ"), size = 1) +
  geom_line(aes(y = DZ, color = "DZ"), size = 1) +
  geom_line(aes(y = LZ, color = "LZ"), size = 1) +
  labs(x = "Days",
       y = "B cells",
       color = "Legend") +
  scale_color_manual(values = colors)
my_gc

ggplot(dat, aes(x = time/24,y=dat$`mean affinity DZ`))+geom_line()

ggplot(dat, aes(x = time/24,y=dat$aff_igm_in_igm))+geom_line(color='red')+
  geom_line(aes(y=dat$aff_igg_in_igm),color='blue')+geom_line(aes(y=dat$aff_igm_in_igg),color='red')+
  geom_line(aes(y=dat$aff_igg_in_igg),color='blue')

ggplot(dat, aes(x = time/24,y=dat$aff_igm_in_ps))+geom_line(color='red')+
  geom_line(aes(y=dat$aff_igg_in_ps),color='blue')


nd<-dat[,c('time','avg_div')]
nd<-nd[!is.na(nd$avg_div),]

k1<-seq.int(1,500,1)
result <- data.frame(matrix(nrow = length(k1), ncol = 2))
for (i in 1:length(k1)){
  print(c(k1[i],k1[i+1]))
  dd<-nd[nd$time>k1[i] & nd$time<k1[i+1],]
  result[i,1]<-colMeans(dd)[1]
  result[i,2]<-colMeans(dd)[2]
  
}
#pdf('Selected_no_of_div.pdf')
gc_div <- ggplot(data = result, mapping = aes(x = result$X1/24, y = result$X2))+geom_line(color='red')+labs(x = "days",y='No of Division')
gc_div
#dev.off()

gc_ratio <- ggplot(data = dat, mapping = aes(x = dat$time/24, y = dat$DZ/dat$LZ))+
  ylim(0,10)+ geom_line(color='red')+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  geom_hline(yintercept=2, linetype="dashed", color = "red")+
  labs(x = "Days",y = "DZ/LZ ratio")

gc_ratio

gc_aff <- ggplot(data = dat, mapping = aes(x = dat$time/24, y = dat$`mean affinity DZ`))+geom_line()+
  labs(x = "Days",y = "Affinity of B cells")

gc_aff

#pdf('basic_gc.pdf')
grid.arrange(my_gc, gc_ratio,gc_div, gc_aff,nrow = 2)
dev.off()

######## igm and igg

my_plot <- ggplot(data = dat, mapping = aes(x = dat$time, y = dat$igm_in_igm+dat$igm_in_igg))+ geom_line(color='red')#+
geom_line(data = dat, mapping = aes(x = dat$time, y = dat$igg_in_igm+dat$igg_in_igg),color='green')

my_plot

my_plot1 <- ggplot(data = dat, mapping = aes(x = dat$time, y = dat$igm_in_ps))+ geom_line(color='red')#+
geom_line(data = dat, mapping = aes(x = dat$time, y = dat$igg_in_igm),color='green')

my_plot1

### switching 
colors <- c("igm_in_LZ" = "blue", "igg_in_LZ" = "red", "LZ" = "green")
my_plot2 <- ggplot(data = dat, aes(x = time/24))+ geom_line(aes(y = igm_in_LZ+igm_in_DZ,color="igm_in_LZ"))+
  geom_line(data = dat, mapping = aes(y = igg_in_LZ+igg_in_DZ,color="igg_in_LZ"))+
  geom_line(data = dat, mapping = aes(y = DZ+LZ,color="LZ"))+
  labs(x = "Days",y = "B cells",color = "Legend") +
  scale_color_manual(values = colors)

#pdf('class_switch.pdf')
my_plot2
#dev.off()


ggplot(data=dat, aes(x=time/24, y=100*(igm_in_LZ+igm_in_DZ)/(DZ+LZ))) + geom_line(color='red',size = 1)+
  geom_line(aes(y=100*(igg_in_LZ+igg_in_DZ)/(DZ+LZ)),color='blue')




ggplot(data=dat, aes(x=time/24, y=100*(igm_in_LZ)/(DZ+LZ))) + geom_line(color='red',size = 1)+
  geom_line(aes(y=100*(igg_in_LZ)/(DZ+LZ)),color='blue')





aff_ig <- ggplot(data = dat, mapping = aes(x = dat$time, y = dat$aff_igm_in_igm))+ geom_line(color='red')+
  geom_line(data = dat, mapping = aes(x = dat$time, y = dat$aff_igg_in_igm),color='green')+
  geom_line(data = dat, mapping = aes(x = dat$time, y = dat$aff_igm_in_igg),color='blue')+
  geom_line(data = dat, mapping = aes(x = dat$time, y = dat$aff_igg_in_igg),color='black')+
  geom_line(data = dat, mapping = aes(x = dat$time, y = dat$aff_igm_in_ps),color='orange')+
  geom_line(data = dat, mapping = aes(x = dat$time, y = dat$aff_igg_in_ps),color='magenta')
aff_ig

ggplot(data = dat, mapping = aes(x = dat$time, y = dat$`free Tfh`))+ geom_line(color='red')+ ylim(150,200)

pdf('ig_dyn.pdf')
grid.arrange(my_plot2, aff_ig,nrow = 1)
dev.off()

nd<-dat[,c('time','igm_in_igm','igm_in_igg','igg_in_igm','igg_in_igg','igm_in_ps','igg_in_ps')]
nd$total_igm<-nd$igm_in_igm+nd$igm_in_igg
nd$total_igg<-nd$igg_in_igm+nd$igg_in_igg
nd$total_ps<-nd$igm_in_ps+nd$igg_in_ps


k1<-seq.int(1,900,1)
result <- data.frame(matrix(nrow = length(k1), ncol = 4))
for (i in 1:length(k1)){
  print(c(k1[i],k1[i+1]))
  dd<-nd[nd$time>k1[i] & nd$time<k1[i+1],]
  result[i,1]<-colMeans(dd)[1]
  result[i,2]<-colSums(dd)[8]
  result[i,3]<-colSums(dd)[9]
  result[i,4]<-colSums(dd)[10]
  
  
}
colors<-c("igm"="blue","igg"="red","ps"="black")
p<-ggplot(data=result,mapping = aes(x=X1/24))+geom_line(aes(y=X2,color='igm'))+
  geom_line(data=result,mapping = aes(y=X3,color='igg'))+
  geom_line(data=result,mapping = aes(y=X4,color='ps'))+labs(x = "Days",y = "B cells",color = "Legend") +
  scale_color_manual(values = colors)
#pdf('per_hour_memory_long.pdf')
p
#dev.off()

cd<-result[1,2:4]
for (i in 2:900){
  cd<-rbind(cd,result[i,2:4]+cd[nrow(cd),])
  
}
cd$X1<-result$X1
colors<-c("igm"="blue","igg"="red","ps"="black")
#pdf('cumulative_memory_long.pdf')
ggplot(data = cd, mapping = aes(x = X1/24))+ geom_line(aes(y=X2,color='igm'),size=1)+
  geom_line(aes(y=X3,color='igg'),size=1)+geom_line(aes(y=X4,color='ps'),size=1)+labs(x = "Days",y = "B cells",color = "Legend") +
  scale_color_manual(values = colors)
#dev.off()

###################### no-normalize
######## wt data
nd<-dat[,c('time','igm_in_igm','igm_in_igg','igg_in_igm','igg_in_igg','igm_in_ps','igg_in_ps')]
nd$total_igm<-nd$igm_in_igm+nd$igm_in_igg
nd$total_igg<-nd$igg_in_igm+nd$igg_in_igg
nd$total_ps<-nd$igm_in_ps+nd$igg_in_ps
nd$time<-nd$time
#nd$time<-nd$time/24

#k1<-seq.int(1,40,1)
k1<-seq.int(1,900,1)
result <- data.frame(matrix(nrow = length(k1), ncol = 4))
for (i in 1:length(k1)){
  print(c(k1[i],k1[i+1]))
  dd<-nd[nd$time>k1[i] & nd$time<k1[i+1],]
  result[i,1]<-colMeans(dd)[1]
  result[i,2]<-colSums(dd)[8]
  result[i,3]<-colSums(dd)[9]
  result[i,4]<-colSums(dd)[10]
  
}

igd<-read.table('/home/abp19/Projects/GC/ODE_GC/weisel_datasets.csv',sep=',')

fun<-function(x){as.numeric(as.character(x))}
igd<-sapply(igd,fun)
igd2<-data.frame(time=igd[2:10,1],igm=igd[2:10,2],igg=igd[2:10,4],ps=igd[2:10,6])

cd<-result
cd$X1<-cd$X1/24

colors<-c("igm"="blue","igg"="red","ps"="black")
#pdf('no-normalize_mem.pdf')
ggplot()+geom_point(aes(x=igd2$time,y=igd2$igm,color='igm',size=1))+
  geom_point(aes(x=igd2$time,y=igd2$igg,color='igg',size=1))+geom_point(aes(x=igd2$time,y=igd2$ps,color='ps',size=1))+
  geom_line(aes(x=cd$X1,y=cd$X2*100/600,color='igm'))+
  geom_line(aes(x=cd$X1,y=cd$X3*100/600,color='igg'))+
  geom_line(aes(x=cd$X1,y=cd$X4*100/600,color='ps'))+labs(x = "Days",y = "B cells",color = "Legend") +
  scale_color_manual(values = colors)

#dev.off()

###### adoptive transfer data
nd<-dat[,c('time','igm_in_igm','igm_in_igg','igg_in_igm','igg_in_igg','igm_in_ps','igg_in_ps','DZ+LZ')]
nd$total_igm<-nd$igm_in_igm+nd$igm_in_igg
nd$total_igg<-nd$igg_in_igm+nd$igg_in_igg
nd$total_ps<-nd$igm_in_ps+nd$igg_in_ps
nd$time<-nd$time
#nd$time<-nd$time/24

ggplot()+geom_line(aes(x=nd$time,y=nd$total_igg_per))

#k1<-seq.int(1,40,1)
k1<-seq.int(1,900,1)
result <- data.frame(matrix(nrow = length(k1), ncol = 5))
for (i in 1:length(k1)){
  print(c(k1[i],k1[i+1]))
  dd<-nd[nd$time>k1[i] & nd$time<k1[i+1],]
  result[i,1]<-colMeans(dd)[1]
  result[i,2]<-colSums(dd)[9]
  result[i,3]<-colSums(dd)[10]
  result[i,4]<-colSums(dd)[11]
  result[i,5]<-colMeans(dd)[8]
  
  
}

fig1d<-read.table('/home/abp19/Projects/GC/ODE_GC/fig_1D.csv',sep=',')
fig1f<-read.table('/home/abp19/Projects/GC/ODE_GC/fig_1F.csv',sep=',')

fig1f$igg<-fig1d$V2*100/fig1f$V2

result$X1<-result$X1/24
result$X1<-result$X1+1

ggplot()+geom_point(aes(x=fig1d$V1,y=fig1d$V2,size=0.1))+
  geom_line(aes(x=result$X1,y=100*result$X3/result$X5))

result$per<-100*result$X3/result$X5
result<-result[-900,]
nk<-seq.int(2,37,3)
nk2<-data.frame(matrix(nrow = length(nk), ncol = 2))


for (i in 1:length(nk)){
  print(c(nk[i],nk[i+1]))
  dd<-result[result$X1>nk[i] & result$X1<nk[i+1],]
  
  nk2[i,1]<-nk[i+1]
  nk2[i,2]<-colSums(dd[sample(nrow(dd),10),])[6]
  
}

for (i in 1:length(nk)){
  print(c(nk[i],nk[i+1]))
  dd<-result[result$X1>nk[i] & result$X1<nk[i+1],]
  
  nk2[i,1]<-nk[i+1]
  nk2[i,2]<-colSums(dd)[6]
  
}

ggplot()+geom_point(aes(x=fig1f$V1,y=fig1f$igg,size=0.1))+
  geom_line(aes(x=nk2$X1,y=nk2$X2*0.13))

ggplot()+geom_point(aes(x=fig1f$V1,y=fig1f$igg,size=0.1))+
  geom_line(aes(x=nk2$X1,y=nk2$X2))


#igd<-read.table('/home/abp19/Projects/GC/ODE_GC/weisel_datasets.csv',sep=',')
igd<-read.table('/home/abp19/Projects/GC/ODE_GC/weisel_datasets_adoptive_transfer.csv',sep=',')

fun<-function(x){as.numeric(as.character(x))}
igd<-sapply(igd,fun)
igd2<-data.frame(time=igd[3:242,1],igm=igd[3:242,2],igg=igd[3:242,4],ps=igd[3:242,6])

cd<-result
cd$X1<-cd$X1/24
# no<-igd[2:10,1]
# pos<-NULL
# for (i in 1:length(no)){
#   pos[i]<-which.min(abs(cd$X1-no[i]))
# }
# cd2<-cd[pos,]
# 
# cd2[,2:4]<-cd2[,2:4]*100/600

#nd<-cbind(cd2,igd2)

# cd3<-cd2[1,2:4]
# for (i in 2:13){
#   cd3<-rbind(cd3,cd2[i,2:4]+cd3[nrow(cd3),])
#   
# }
# 
# 
# nd<-cbind(cd2$X1,cd3,igd2)

colors<-c("igm"="blue","igg"="red","ps"="black")
#pdf('no-normalize_mem.pdf')
# ggplot(data=nd,mapping = aes(x=nd$time))+geom_point(aes(y=nd$igm,color='igm'))+
#   geom_point(aes(x=nd$time,y=nd$igg,color='igg'))+geom_point(aes(x=nd$time,y=nd$ps,color='ps'))+
#   geom_line(aes(x=nd$time,y=nd$X2,color='igm'))+
#   geom_line(aes(x=nd$time,y=nd$X3,color='igg'))+
#   geom_line(aes(x=nd$time,y=nd$X4,color='ps'))+labs(x = "Days",y = "B cells",color = "Legend") +
#   scale_color_manual(values = colors)
#dev.off()
ggplot()+geom_point(aes(x=igd2$time,y=igd2$igm,color='igm',size=0.1))+
  geom_point(aes(x=igd2$time,y=igd2$igg,color='igg',size=0.1))+geom_point(aes(x=igd2$time,y=igd2$ps,color='ps',size=0.1))+
  geom_line(aes(x=cd$X1,y=cd$X2*100/600,color='igm'))+
  geom_line(aes(x=cd$X1,y=cd$X3*100/600,color='igg'))+
  geom_line(aes(x=cd$X1,y=cd$X4*100/600,color='ps'))+labs(x = "Days",y = "B cells",color = "Legend") +
  scale_color_manual(values = colors)


result<-result[-900,]
da<-seq.int(0,38,3)

md<-data.frame(matrix(nrow = length(k1), ncol = 4))
for (i in 1:length(da)){
  dd<-result[result$X1>da[i] & result$X1<da[i+1],]
  md[i,1]<-da[i]+2
  md[i,2]<-colSums(dd)[2]
  md[i,3]<-colSums(dd)[3]
  md[i,4]<-colSums(dd)[4]
}

ggplot()+geom_line(aes(x=md$X1,y=md$X2),color='red')+
  geom_line(aes(x=md$X1,y=md$X3),color='blue')+geom_line(aes(x=md$X1,y=md$X4),color='green')

############ normalize wrt igg


igd<-read.table('/home/abp19/Projects/GC/ODE_GC/wpd_datasets.csv',sep=',')

fun<-function(x){as.numeric(as.character(x))}
igd<-sapply(igd,fun)
igd2<-data.frame(time=igd[3:15,1])
igd2$igm=igd[3:15,2]/igd[3:15,4]

igd2$ps=igd[3:15,6]/igd[3:15,4]
igd2$igg=igd[3:15,4]/igd[3:15,4]

cd$X1<-cd$X1/24
cd$X2<-cd$X2/cd$X3
cd$X4<-cd$X4/cd$X3
cd$X3<-cd$X3/cd$X3




nd<-data.frame(time=rep('na',900-13),igm=rep('na',900-13),igg=rep('na',900-13),ps=rep('na',900-13))
nd<-rbind(igd2,nd)
nd<-cbind(nd,cd)

nd<-sapply(nd,fun)
nd<-as.data.frame(nd)
colors<-c("igm"="blue","igg"="red","ps"="black")
#pdf('mem_data_sim.pdf')
ggplot(data=nd,mapping = aes(x=time))+geom_point(aes(y=igm,color='igm'))+
  geom_point(aes(x=time,y=igg,color='igg'))+geom_point(aes(x=time,y=ps,color='ps'))+
  geom_line(aes(x=X1,y=X2,color='igm'))+
  geom_line(aes(x=X1,y=X3,color='igg'))+
  geom_line(aes(x=X1,y=X4,color='ps'))+labs(x = "Days",y = "B cells",color = "Legend") +
  scale_color_manual(values = colors)

#dev.off()

###################### normalize wrt igm

cd<-result[1,2:4]
for (i in 2:500){
  cd<-rbind(cd,result[i,2:4]+cd[nrow(cd),])
  
}
cd$X1<-result$X1
ggplot(data = cd, mapping = aes(x = X1, y = X2))+ geom_line(color='red')+
  geom_line(aes(y=X3),color='green')+geom_line(aes(y=X4),color='blue')

igd<-read.table('/home/abp19/Projects/GC/ODE_GC/wpd_datasets.csv',sep=',')

fun<-function(x){as.numeric(as.character(x))}
igd<-sapply(igd,fun)
igd2<-data.frame(time=igd[3:15,1])


igd2$ps=igd[3:15,6]/igd[3:15,2]
igd2$igg=igd[3:15,4]/igd[3:15,2]

igd2$igm=igd[3:15,2]/igd[3:15,2]

cd$X1<-cd$X1/24

cd$X4<-cd$X4/cd$X2
cd$X3<-cd$X3/cd$X2
cd$X2<-cd$X2/cd$X2



nd<-data.frame(time=rep('na',500-13),igm=rep('na',500-13),igg=rep('na',500-13),ps=rep('na',500-13))
nd<-rbind(igd2,nd)
nd<-cbind(nd,cd)

nd<-sapply(nd,fun)
nd<-as.data.frame(nd)
ggplot(data=nd,mapping = aes(x=nd$time,y=nd$igm))+geom_point(color='red')+
  geom_point(aes(x=nd$time,y=nd$igg),color='green')+geom_point(aes(x=nd$time,y=nd$ps),color='blue')+
  geom_line(aes(x=nd$X1,y=nd$X2),color='red')+
  geom_line(aes(x=nd$X1,y=nd$X3),color='green')+
  geom_line(aes(x=nd$X1,y=nd$X4),color='blue')




################# normalize wrt ps
cd<-result[1,2:4]
for (i in 2:500){
  cd<-rbind(cd,result[i,2:4]+cd[nrow(cd),])
  
}
cd$X1<-result$X1
ggplot(data = cd, mapping = aes(x = X1, y = X2))+ geom_line(color='red')+
  geom_line(aes(y=X3),color='green')+geom_line(aes(y=X4),color='blue')

igd<-read.table('/home/abp19/Projects/GC/ODE_GC/wpd_datasets.csv',sep=',')

fun<-function(x){as.numeric(as.character(x))}
igd<-sapply(igd,fun)
igd2<-data.frame(time=igd[3:15,1])



igd2$igg=igd[3:15,4]/igd[3:15,6]

igd2$igm=igd[3:15,2]/igd[3:15,6]
igd2$ps=igd[3:15,6]/igd[3:15,6]

cd$X1<-cd$X1/24


cd$X3<-cd$X3/cd$X4
cd$X2<-cd$X2/cd$X4
cd$X4<-cd$X4/cd$X4


nd<-data.frame(time=rep('na',500-13),igm=rep('na',500-13),igg=rep('na',500-13),ps=rep('na',500-13))
nd<-rbind(igd2,nd)
nd<-cbind(nd,cd)

nd<-sapply(nd,fun)
nd<-as.data.frame(nd)
ggplot(data=nd,mapping = aes(x=nd$time,y=nd$igm))+geom_point(color='red')+
  geom_point(aes(x=nd$time,y=nd$igg),color='green')+geom_point(aes(x=nd$time,y=nd$ps),color='blue')+
  geom_line(aes(x=nd$X1,y=nd$X2),color='red')+
  geom_line(aes(x=nd$X1,y=nd$X3),color='green')+
  geom_line(aes(x=nd$X1,y=nd$X4),color='blue')+ylim(0,10)

























my_plot <- ggplot(data = dat, 
                  mapping = aes(x = dat$`0`, y = dat$`3`)) +
  geom_line(color='red') + geom_line(data = dat, 
                                     mapping = aes(x = dat$`0`, y = dat$`4`),color='blue')+ geom_line(data = dat, 
                                                                                                      mapping = aes(x = dat$`0`, y = dat$`5`),color='green')
my_plot

igm<-dat[,c(1,4,5,6)]

k1<-seq.int(1,500,1)
result <- data.frame(matrix(nrow = length(k1), ncol = 4))
for (i in 1:length(k1)){
  print(c(k1[i],k1[i+1]))
  dd<-igm[igm$`0`>k1[i] & igm$`0`<k1[i+1],]
  result[i,1]<-colMeans(dd[,1])
  result[i,2]<-colSums(dd[,2])
  result[i,3]<-colSums(dd[,3])
  result[i,4]<-colSums(dd[,4])
  
}

my_plot <- ggplot(data = result, 
                  mapping = aes(x = result[,1], y = result[,2])) +
  geom_line(color='red') + geom_line(data = result, 
                                     mapping = aes(x = result[,1], y = result[,3]),color='blue')+ geom_line(data = result, 
                                                                                                            mapping = aes(x = result[,1], y = result[,4]),color='green')
my_plot


d <- density(igm) # returns the density data
plot(d) # plots the results 

my_plot <- ggplot(data = dat, mapping = aes(x = dat$`0`, y = dat$`12`))+ geom_line(color='red')

my_plot

p <- ggplot(data = dat, mapping = aes(x = dat$`0`, y = dat$`9`))+ geom_line(color='red')+
  geom_line(color='green',mapping = aes(x = dat$`0`, y = dat$`10`))
p

my_plot <- ggplot(data = dat, mapping = aes(x = dat$`0`, y = dat$`7`+dat$`9`))+ geom_line(color='red')+
  geom_line(color='green',mapping = aes(x = dat$`0`, y = dat$`8`+dat$`10`))

my_plot

my_plot <- ggplot(data = dat, mapping = aes(x = dat$`0`, y = dat$`10`))+ geom_line(color='red')+
  geom_line(color='green',mapping = aes(x = dat$`0`, y = dat$`11`))

my_plot



my_plot <- ggplot(data = dat, 
                  mapping = aes(x = dat$`0`, y = dat$`11`))+geom_line()
my_plot
#########################

dat1<-read_excel('affinity.xlsx', col_types = "numeric")
dat2<-read_excel('tfh_help.xlsx', col_types = "numeric")

tfh_help<-1.0
affinity <-0.7

k1<-seq.int(1,500,1)
result <- data.frame(matrix(nrow = length(k1), ncol = 4))
kp1<-NULL
kp2<-NULL
kp3<-NULL
kp4<-NULL

for (i in 1:(length(k1)-1)){
  print(c(k1[i],k1[i+1]))
  dd1<-dat1[dat1$`0`>k1[i] & dat1$`0`<k1[i+1],]
  dd2<-dat2[dat2$`0`>k1[i] & dat2$`0`<k1[i+1],]
  s1<-dim(dd2)[1]
  s2<-dim(dd2)[2]
  bb1<-NULL
  bb2<-NULL
  bb3<-NULL
  bb4<-NULL
  for (j in 1:s1){
    pos1<-which(dd2[j,c(2:s2)]<tfh_help)
    pos2<-which(dd2[j,c(2:s2)]>tfh_help)
    
    bb1[j]<-length(which(dd1[j,c(2:s2)][pos1]<affinity))
    # bb1[j]<-length(dd1[j,c(2:s2)][pos1])
    bb2[j]<-length(which(dd1[j,c(2:s2)][pos1]>affinity))
    
    bb3[j]<-length(which(dd1[j,c(2:s2)][pos2]<affinity))
    bb4[j]<-length(which(dd1[j,c(2:s2)][pos2]>affinity))
    
  }
  kp1[i]<-sum(bb1)
  kp2[i]<-sum(bb2)*0.3
  kp3[i]<-sum(bb3)*0.3
  kp4[i]<-sum(bb4)
  # result[i,1]<-colMeans(dd[,1])
  # result[i,2]<-colSums(dd[,2])
  # result[i,3]<-colSums(dd[,3])
  # result[i,4]<-colSums(dd[,4])
  
}

nd<-data.frame(time=seq.int(1.5,500,1),igm=kp1,igg=kp2,pls=kp3,pls2=kp4)

ggplot(data=nd, aes(x=time, y=igm)) + geom_line(color='red',size = 1)+
  geom_line(aes(x=time, y=igg),color='blue')+
  geom_line(aes(x=time, y=pls),color='green')+
  geom_line(aes(x=time, y=pls2),color='magenta')
labs(x = "days",y='production/hr')




###############

dd<-read.table('/home/abp19/Projects/hyphasma_ar/hyphasma_ar/gc_numbers_arA.out')
dd<-read.table('/home/abp19/Projects/hyphasma_ar/hyphasma_ar/apolog.out')

igms<-dd[dd$V2==2 & dd$V8==0,]
iggs<-dd[dd$V2==2 & dd$V8==1,]

k1<-seq.int(1,500,1)
result1 <- data.frame(matrix(nrow = length(k1), ncol = 2))
result2 <- data.frame(matrix(nrow = length(k1), ncol = 2))
for (i in 1:length(k1)){
  print(c(k1[i],k1[i+1]))
  ddm<-igms[igms$V1>k1[i] & igms$V1<k1[i+1],]
  ddg<-igms[iggs$V1>k1[i] & iggs$V1<k1[i+1],]
  if (nrow(ddm)==0){
    result1[i,1]<-mean(c(k1[i],k1[i+1]))
    result1[i,2]<-0
  }
  else{
    result1[i,1]<-mean(c(k1[i],k1[i+1]))
    result1[i,2]<-nrow(ddm)
    
  }
  
  if (nrow(ddg)==0){
    result2[i,1]<-mean(c(k1[i],k1[i+1]))
    result2[i,2]<-0
  }
  else{
    result2[i,1]<-mean(c(k1[i],k1[i+1]))
    result2[i,2]<-nrow(ddg)
    
  }
  
  
}

ggplot(data=dd, aes(x=V1/24, y=V2)) + geom_line(color='red',size = 1)+
  geom_line(aes(y=V3),color='blue')+
  geom_line(aes(y=V4),color='green')#+geom_line(aes(y=V5+V6+V7),color='magenta')
dd$n<-dd$V5+dd$V6+dd$V7

ggplot(data=dd, aes(x=V1/24, y=100*V5/n)) + geom_line(color='red',size = 1)+
  geom_line(aes(y=100*V6/n),color='blue')+
  geom_line(aes(y=100*V7/n),color='green')

ggplot() + geom_line(aes(x=result1$X1,y=result1$X2),color='red')+
  geom_line(aes(x=result2$X1,y=result2$X2),color='blue')



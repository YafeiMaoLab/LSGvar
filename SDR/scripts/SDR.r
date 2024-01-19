
library(dbscan)
library(dplyr)
library(data.table)
library(IRanges)
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
print(args[2])
pos <- read.delim(args[2], header = FALSE, col.names = c("ref_chr", "ref_start", "ref_end", "ref_pos", "query_chr", "query_start", "query_end", "query_pos", "orient"))
chrpc<-fread(args[3])
chrpc<-distinct(chrpc)
colnames(chrpc)<-c("ref_chr","ref_len","query_chr","query_len")

chrnames <- unique(pos$ref_chr)
numeric_part <- as.numeric(gsub("\\D", "", chrnames))
sorted_chrnames <- chrnames[order(numeric_part)] 



# pos$ref_start<-pos$ref_start+1
# pos$query_start<-pos$query_start+1

# store$len1<-store$ref_end-store$ref_start
# store$len2<-store$query_end-store$query_start



cluster0paras<-700000 
cluster1paras<-700000
#sorted_chrnames<- sorted_chrnames[14:23]
for(chrid in sorted_chrnames){

  store<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss')
  storesmall<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,anno='ss',orient="+")
  duplication<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster=0)
  inversion<-data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0)
  smalltrans<- data.frame(ref_chr=0,ref_start = 0, ref_end = 0,query_chr=0,query_start = 0, query_end = 0,orient=0,cluster=0)
  pos.chr<-pos[pos$ref_chr==chrid,]
  pos.chr<-pos.chr[order(pos.chr$ref_start),]
  
  

  endcluster0<-split_region(pos.chr,0,cluster0paras)
  # Dupex<-clusterall(endcluster0)
  # Dupex<-Dupex[order(Dupex$ref_start),]
  # if(nrow(Dupex)!=1){
  #   for(i in 2:dim(Dupex)[1]){
  #     if(Dupex[i,]$ref_start<Dupex[i-1,]$ref_end & Dupex[i,]$ref_start>Dupex[i-1,]$ref_start){
  #       a<-Dupex[i,]$ref_start
  #       b<-min(Dupex[c(i,i-1),]$ref_end)
  #       if(nrow(endcluster0[endcluster0$ref_start>=a & endcluster0$ref_end<=b,])!=0)
  #       {
  #         c<-clusterall(endcluster0[endcluster0$ref_start>=a & endcluster0$ref_end<=b,])
  #         duplication<-rbind(duplication,c)
  #       }
  #       else{
  #         if((b-a)>10000){
  #           c<-Dupex[i,]
  #           c$ref_start<-a
  #           c$ref_end<-b
  #           c$query_end<-Dupex[i,]$query_start
  #           cc<-Dupex[which((c$query_end-Dupex$query_end)>0),]
  #           c$query_start<-cc$query_end[which.min(c$query_end - cc$query_end)]
  #           duplication<-rbind(duplication,c)
  #         }
  #       }
  #       
  #     }
  #   }
  # }
  # 
  aftertans<-smalltransf(endcluster0)
  # if(length(aftertans$tran)!=0){
  #   data1<-aftertans$tran
  #   for(j in 2:dim(duplication)[1]){
  #     data2<-data1[data1$ref_start>=duplication[j,]$ref_start & data1$ref_end<=duplication[j,]$ref_end & data1$query_start>=duplication[j,]$query_start &data1$query_end<=duplication[j,]$query_end,]
  #     if(nrow(data2)!=0){
  #       duplication<-duplication[-j,]
  #     }
  #   }
  # }
  # 
  
  
  
  
  #dotplot_cluster(endcluster0)
  for(i in unique((aftertans$tranbefore)$cluster)){
    a<-aftertans$tranbefore[(aftertans$tranbefore)$cluster==i,]
    a$cluster<-as.character(1:dim(a)[1])
    if(sum(a$orient == "-") > 0.5* nrow(a)){
      reverse_end<-reverse.region(a,chrid,3)
    }
    if(sum(a$orient == "+") > 0.5 * nrow(a) |sum(a$orient == "-") == 0.5 * nrow(a)){
      reverse_end<-reverse.region(a,chrid,2)
    }
    if(exists("reverse_end")){
      middle<-reverse_end$reverse
      if(nrow(middle)!=0){
        middle<-intersect.unit(middle,unique(middle$query_chr)) 
      }
      if(sum(a$orient == "+") > 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"+"
      }
      if(sum(a$orient == "-") > 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"-"
      }
      if(sum(a$orient == "-") == 0.5 * nrow(a) & nrow(middle)!=0){
        middle$orient<-"no"
      }
      
      storesmall<-rbind(storesmall,middle)
      minusmiddle<-reverse_end$reverse[reverse_end$reverse$ref_start<=reverse_end$reverse$ref_end & reverse_end$reverse$query_start<=reverse_end$reverse$query_end,]
      if(nrow(minusmiddle)!=0){
        minusmiddle$orient<-"no"
      }
      
      storesmall<-rbind(storesmall,minusmiddle)
      storesmall<-distinct(storesmall)
    }
  }
  if(length(aftertans$tran)!=0){
    smalltrans<-rbind(smalltrans,aftertans$tran)
  }
  if(!isEmpty(aftertans$tran)){
    inversion<-rbind(inversion,aftertans$tranbefore[aftertans$tranbefore$orient=="-",colnames(inversion)])
    transdup<-transduplication_extract(endcluster0,aftertans$tran)
    if(nrow(transdup)!=0){
      duplication<-rbind(duplication,transdup[,colnames(duplication)])
    }
    
 }
  
  
  
  endcluster0<-aftertans$endcluster0
  endcluster1<-split_region(endcluster0,1,cluster1paras)

  endcluster1<-minus.next(endcluster1)
  inver<-inversion.extract(endcluster1,chrid)
  list<-which(rle(rle(endcluster1$cluster)$lengths)$values==1 &rle(rle(endcluster1$cluster)$lengths)$lengths >5)
  cross.region<-cross.calcu(endcluster1)
  initclsuer<-1000
  for(k in list){
    chaoval<-rle(rle(endcluster1$cluster)$lengths)$values
    chaolen<-rle(rle(endcluster1$cluster)$lengths)$lengths
    start<-sum(chaoval[1:k-1]*chaolen[1:k-1])+1
    
    endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),]$cluster<-as.character(initclsuer)
    print(endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),])
    chaos<-endcluster1[intersect(which(!endcluster1$cluster %in% cross.region),start:(start+chaolen[k]-1)),]
    chaos$query_start<-abs(chaos$query_start)
    chaos$query_end<-abs(chaos$query_end)
    if(nrow(chaos)!=0){
      chaos<-clusterall(chaos)
    }
    chaos$anno<-"COMPLEX"
    store<-rbind(store,chaos[,colnames(store)])
    initclsuer<-initclsuer+1
  }
  
  reverse_end<-reverse.region(endcluster1,chrid,0)
  duplic<-reverse_end$dup
  if(!is.character(duplic)){
    if(nrow(duplic)!=0){
      duplication<-rbind(duplication,duplic[,colnames(duplication)])
    }
    
  }
  minimap<-reverse_end$minimaploc 
  
  
  if ("delsytenic" %in% names(reverse_end)){
    for(cluster in unique(reverse_end$delsytenic$cluster)){
      datalist<-reverse_end$delsytenic[reverse_end$delsytenic$cluster==cluster,]
      datalist$cluster<-1:dim(datalist)[1]
      reverse_end1<-reverse.region(datalist,unique(datalist$ref_chr),2)
      middle<-reverse_end1$reverse
      if(nrow(middle)!=0){
        middle<-intersect.unit(middle,unique(middle$query_chr))
        middle$orient<-"no"
        storesmall<-rbind(storesmall,middle)
      }
    }
  }
  
  
  for (value in unique(reverse_end$reverse$query_chr)) {
    datarev<-reverse_end$reverse[reverse_end$reverse$query_chr==value,]
    if(nrow(datarev)!=0){
      rows_to_remove <- which(datarev$query_chr == value)
      if (length(rows_to_remove) > 1) {
        startend<-datarev[c(rows_to_remove[1], tail(rows_to_remove, 1)), ]
        store<-rbind(store,datarev[-c(rows_to_remove[1], tail(rows_to_remove, 1)), ]) 
        startend$orient<-"+"
        storesmall<-rbind(storesmall,startend)
      }
    }
  }

  vectore<-which(store$ref_start-store$ref_end==2)
  if(length(vectore)!=0){
    for(l in vectore){
      store[l,]$ref_start<-store[l,]$ref_start-1
      store[l,]$ref_end<-store[l,]$ref_end+1
    }
  }
  vectore<-which((store$query_start-store$query_end)==2)
  if(length(vectore)!=0){
    for(l in vectore){
      store[l,]$query_start<-store[l,]$query_start-1
      store[l,]$query_end<-store[l,]$query_end+1
    }
  }
  
  if(nrow(store[store$ref_start>store$ref_end,])!=0){
    store[store$ref_start>store$ref_end,]$ref_end=store[store$ref_start>store$ref_end,]$ref_start 
  }
  if(length(unique(store$query_chr)[-1])!=0){
    for(chr_child in unique(store$query_chr)[-1]){
      new_row<-intersect.unit(store,chr_child) 
      if(nrow(new_row)!=1){
        new_row<-complex(new_row)
      }
      percen<-max(new_row[new_row$anno=="COMPLEX",]$ref_end-new_row[new_row$anno=="COMPLEX",]$ref_start)
      if(percen>100000000){ 
        new_row<-intersect.unit(store,chr_child) 
      }
      assign(chr_child,new_row)
    }
    store<-docall(store)
    rm(list=unique(store$query_chr)) 
  }
  
  if(dim(inver)[1]!=0){
    inver<-inver[,-1]
    inver<-inver[,-7]
    inversion<-rbind(inversion,inver)
  }

  
  endcluster1<-reverse_end$endcluster1


  
  m=0
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster%in% as.numeric(names(count_minus))])
  ##cluster
  for(clusid in as.numeric(names(count_minus[(count_minus / count_total) >0.6])) ){
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,100000)
    minusduplic<-duplication_extract(endcluster1,miuscluster)  
    if(!isEmpty(minusduplic$dupli)){
      duplication<-rbind(duplication,minusduplic$dupli[,colnames(duplication)])
    }
    
   
    reverse_end<-reverse.region(miuscluster,chrid,3)
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){
  
      local<-intersect.unit(middle,unique(middle$query_chr))
      changed_rows <- anti_join(local,middle)
      store<-rbind(store,inner_join(local,middle))
      if(nrow(changed_rows)!=0){
        changed_rows$anno<-'COMPLEX'
        store<-rbind(store,changed_rows[,colnames(store)])
      }}
    
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster$cluster<-paste(miuscluster$cluster, "8",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    
    m=m+1
  }
  
  
  count_minus <- table(endcluster1$cluster[endcluster1$orient == '-'])
  count_total <- table(endcluster1$cluster[endcluster1$cluster%in% as.numeric(names(count_minus))])
  
  
  
  SDRminudlist<-as.numeric(names(count_minus[(count_minus / count_total) >0.6])) ##cluster

  for(k in SDRminudlist){
    region<-endcluster1[endcluster1$cluster==k,]
    if(dim(region)[1]==1){
      next
    }
    #endcluster2<-split_region(region,2)
    region$cluster<-1:dim(region)[1]
   
    if(dim(region)[1]==1){
      next
    }
    
    endcluster2<-duplication_extract(region,region)
    if(nrow(endcluster2$dupli)!=0){
      duplication<-rbind(duplication,endcluster2$dupli[,colnames(duplication)])
    }

    endcluster2before<-endcluster2$pos_end
    orientid="+"
    endcluster2<-smallcluster(endcluster2$pos_end,orientid) 
    reverse_end<-reverse.region(endcluster2,chrid,3)
    middle<-reverse_end$reverse
    if(nrow(middle)!=0){
      #middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start 
      local<-intersect.unit(middle,unique(middle$query_chr)) 
      changed_rows <- anti_join(local,middle)
      if(nrow(changed_rows)!=0){
        changed_rows$anno<-'COMPLEX'
        store<-rbind(store,changed_rows[,colnames(store)])
      }
      
      if(length(which(middle$query_start>middle$query_end))!=0){
        middle<-middle[-(which(middle$query_start>middle$query_end)),]
        
      }
      if(nrow(middle)!=0){
        middle$orient<-"-"
        storesmall<-rbind(storesmall,middle[,colnames(storesmall)]) 
      }
      
    }
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
  }
  
  
  
  
  for(k in setdiff(unique(endcluster1$cluster),SDRminudlist)){
    region<-endcluster1[endcluster1$cluster==k,]
    region$cluster<-1:dim(region)[1]
    if(dim(region)[1]==1){
      next
    }
    
    endcluster2<-duplication_extract(region,region)
    if(nrow(endcluster2$dupli)!=0){
      duplication<-rbind(duplication,endcluster2$dupli[,colnames(duplication)])
    }

    #claster2<-(repeat.integrate(region,1))
    #endcluster2<-claster2$afterdup
    # if(!is.null(claster2$repeat.region)){
    #   dup<-distinct(claster2$repeat.region)
    #   duplication<-rbind(duplication,dup)
    # }
    endcluster2before<-endcluster2$pos_end
    orientid="-"
    endcluster2<-smallcluster(endcluster2$pos_end,orientid) 
    reverse_end<-reverse.region(endcluster2,chrid,2)
    middle<-reverse_end$reverse
    local<-intersect.unit(middle,unique(middle$query_chr)) 
    changed_rows <- anti_join(local,middle)
    if(nrow(changed_rows)!=0){
      print(k)
      changed_rows$anno<-'COMPLEX'
      store<-rbind(store,changed_rows[,colnames(store)])
    }
    # middle[middle$ref_start>middle$ref_end,]$ref_end=middle[middle$ref_start>middle$ref_end,]$ref_start 
    # middle<-intersect.unit(middle,unique(middle$query_chr)) 
    if(length(which(middle$query_start>middle$query_end))!=0){
      middle<-middle[-(which(middle$query_start>middle$query_end)),]
    }
    if(nrow(middle)!=0){
      middle$orient<-"+"
      storesmall<-rbind(storesmall,middle[,colnames(storesmall)])
    }
    
    storesmall<-insertsmall(endcluster2before,storesmall,orientid)
    
  }
  inversion$anno<-"INV"
  inversion$orient<-"-"
  duplication<-distinct(duplication)
  duplication$anno<-"DUP"
  duplication$orient<-"no"
  smalltrans$anno<-"TRANS"
  smalltrans$orient<-"no"
  storesmall$anno<-"SDR_NM"
  store$orient<-"no"
  duplication<-duplication[,colnames(store)]
  smalltrans<-smalltrans[,colnames(store)]
  storesmall<-storesmall[-1,colnames(store)]
  store<-store[store$ref_chr!=0,]
  all<-rbind(store,inversion[-1,],distinct(duplication[-1,]),smalltrans[-1,],storesmall)
  all<-all[order(all$ref_start),]
  for (chr_child in unique(all$query_chr[all$query_chr!=0])){
    assign(chr_child,endfilter(all[all$query_chr==chr_child,],chrid,chr_child))
    
  }
  data<-docall(all)
  if(length(which(data$ref_start>data$ref_end))!=0){
    data<-data[-which(data$ref_start>data$ref_end),]
  }
  data<-distinct(data)
  if(length(which(data$reflen==0 & data$querylen==0))!=0){
    data<-data[-which(data$reflen==0 & data$querylen==0),]
  }
  COMPLEX<-data[data$anno=="COMPLEX",]
  if(nrow(data[data$anno=="COMPLEX",])!=0){
    data<-data[-which(data$anno=="COMPLEX"),]
    newcomplex<-complexinte(COMPLEX)
    newcomplex$anno<-"COMPLEX"
    newcomplex$orient<-"no"
    newcomplex$reflen<-newcomplex$ref_end-newcomplex$ref_start
    newcomplex$querylen<-newcomplex$query_end-newcomplex$query_start
    data<-rbind(data,newcomplex[,colnames(data)])
  }
  
  
  write.table(data, paste(args[4],chrid,"SDRend.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  #write.table(minimap, paste("D:/MS/saffire测试/R/minimap/",chrid,"minimap.tsv",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
}

#write.csv(result_data[result_data$SV=="DUP",], paste("D:/MS/saffire测试/R/result/SDRend.csv",sep = ""), quote = FALSE, row.names = FALSE)


exchange<-function(data){
  data <- data[data$ref_start != 0, ]
  data[data$ref_start > data$ref_end, c('ref_start', 'ref_end')] <-  data[data$ref_start > data$ref_end, c('ref_end','ref_start')]
  data[data$query_start > data$query_end, c('query_start', 'query_end')] <-  data[data$query_start > data$query_end, c('query_end','query_start')] 
  return(data)
}

intersect.unit<-function(data,chr_child){
  rm(inte)
  data<-data[data$query_chr==chr_child,]
  data<-exchange(data)
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref) 
  overlapsque<-findOverlaps(ir1que,ir1que) 
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
  }
  if(length(which(overlapsque@from!=overlapsque@to))!=0){
    if(exists("inte")){
      inte<-rbind(inte,as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)]))
    }
    else{
      inte<-as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)])
    }
  }
  if(exists("inte")){
    df_sorted <- t(apply(inte, 1, function(x) sort(x)))
    unique_df <- unique(df_sorted)
    inte <- as.data.frame(unique_df)
    for(j in 1:dim(inte)[1]){
      a<-c()
      a<-append(a,inte[j,1])
      a<-append(a,inte[j,2])
      for(k in 1:dim(inte)[1]){
        if(inte[k,1] %in% ((min(a)-1):(max(a)+1))){
          a<-unique(append(a,inte[k,1]))
        }
        if(inte[k,2] %in% ((min(a)-1):(max(a)+1))){
          a<-unique(append(a,inte[k,2]))
        }
      }
      if(j==1){
        x<-a
        data[min(x),]$ref_start<-min(data[x,]$ref_start)
        data[min(x),]$ref_end<-max(data[x,]$ref_end)
        data[min(x),]$query_start<-min(data[x,]$query_start)
        data[min(x),]$query_end<-max(data[x,]$query_end)
        data[x,]<-data[min(x),]
        data[x,]$anno<-"COMPLEX"
        print(x)
      }
      else{
        if(length(intersect(a,x))==0){
          x<-a
          data[min(x),]$ref_start<-min(data[x,]$ref_start)
          data[min(x),]$ref_end<-max(data[x,]$ref_end)
          data[min(x),]$query_start<-min(data[x,]$query_start)
          data[min(x),]$query_end<-max(data[x,]$query_end)
          data[x,]<-data[min(x),]
          data[x,]$anno<-"COMPLEX"
        }
      }
    }
    data<-distinct(data)
  }
  data<-data[order(data$ref_start),]
  return(data)
}


alignintersect<-function(align1,align2,align3){
  del.list<-c()
  ir1ref <- IRanges(start = align1$ref_start, end = align1$ref_end)
  ir1que  <- IRanges(start = align1$query_start, end = align1$query_end)
  ir2 <- IRanges(start =align2$ref_start , end =align2$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir2)
  if(length(overlapsref@from)!=0){  
    xx<-cbind(overlapsref@from,
              align1[overlapsref@from,]$ref_chr,
              align1[overlapsref@from,]$ref_start,
              align1[overlapsref@from,]$ref_end,
              overlapsref@to,
              align2[overlapsref@to,]$ref_chr,
              align2[overlapsref@to,]$ref_start,
              align2[overlapsref@to,]$ref_end)
    xx<-as.data.frame(xx)
    for(i in 1:dim(xx)[1]){
      if(xx$V3[i]>=xx$V7[i] & xx$V4[i]<=xx$V8[i]){
        del.list<-append(del.list,xx$V1[i])
      }
    }
    
  }
  
  ir3 <- IRanges(start = align3$query_start, end = align3$query_end)
  overlapsque<-findOverlaps(ir1que,ir3)
  if(length(overlapsque@from)!=0){  
    xx<-cbind(overlapsque@from,
              align1[overlapsque@from,]$query_chr,
              align1[overlapsque@from,]$query_start,
              align1[overlapsque@from,]$query_end,
              overlapsque@to,
              align3[overlapsque@to,]$query_chr,
              align3[overlapsque@to,]$query_start,
              align3[overlapsque@to,]$query_end)
    xx<-as.data.frame(xx)
    for(i in 1:dim(xx)[1]){
      if(xx$V3[i]>=xx$V7[i] & xx$V4[i]<=xx$V8[i]){ 
        del.list<-append(del.list,xx$V1[i])
      }
    }
    
  }
  if(length(del.list)!=0)
  {align1<-align1[-as.numeric(del.list),]}
  return(align1)
}

reverse_xy<-function(data){
  selected_rows <- data$orient == '-'
  new_query_start <- data[selected_rows, ]$query_end
  new_query_end <- data[selected_rows, ]$query_start
  data[selected_rows, ]$query_start <- new_query_start
  data[selected_rows, ]$query_end <- new_query_end
  return(data)
}
split_region<-function(pos.chr.region,cluster.id,clusterparas){
  pos.chr.region<-reverse_xy(pos.chr.region)  #负链反向
  df <- pos.chr.region %>%
    rowwise() %>%
    mutate(length = sqrt((ref_end - ref_start)^2 + (query_end- query_start)^2),
           segments = ceiling(length / 5000))
  df_segments <- df[rep(seq_len(nrow(df)), df$segments), ]
  df_segments <- df_segments %>%
    group_by(ref_start, ref_end, query_start, query_end) %>%
    mutate(segment_id = row_number(),
           segment_mid_x = ref_start + (ref_end - ref_start) * (segment_id - 0.5) / segments,
           segment_mid_y = query_start + (query_end - query_start) * (segment_id - 0.5) / segments)
  
  if(cluster.id==2){
    dbs.para<-10000
  }
  if(cluster.id==0){
    dbs.para<-clusterparas
    
  }
  if(cluster.id==1){
    df_segments[df_segments$orient=='-',]$segment_mid_y<--df_segments[df_segments$orient=='-',]$segment_mid_y
    df_segments[df_segments$orient=='-',]$query_start<--df_segments[df_segments$orient=='-',]$query_start
    df_segments[df_segments$orient=='-',]$query_end<--df_segments[df_segments$orient=='-',]$query_end
    pos.chr.region[pos.chr.region$orient=='-',]$query_start<--pos.chr.region[pos.chr.region$orient=='-',]$query_start
    pos.chr.region[pos.chr.region$orient=='-',]$query_end<--pos.chr.region[pos.chr.region$orient=='-',]$query_end
    dbs.para<-clusterparas}
  dbscan_result <- dbscan(df_segments[,c("segment_mid_x", "segment_mid_y")], eps = dbs.para, minPts =1)
  df_segments$cluster <- dbscan_result$cluster
  df_segments <- df_segments %>%
    group_by(ref_chr,ref_start, ref_end, ref_pos,query_chr,query_start, query_end,query_pos,orient) %>%
    summarise(cluster = names(which.max(table(cluster))))
  df_segments$query_start<-abs(df_segments$query_start)
  df_segments$query_end<-abs(df_segments$query_end)
  df_segments<-reverse_xy(df_segments)  #负链反向
  return(df_segments)
}


inte.minud<-function(endcluster1,id){
  if(id==1){
    if(length( which(endcluster1$orient=="-"))!=0){
      endcluster1<-endcluster1[-(which(endcluster1$orient=="-")),]
    }
  }
  if(id==2){
    if(length( which(endcluster1$orient=="+"))!=0){
      endcluster1<-endcluster1[-(which(endcluster1$orient=="+")),]
    }
  }
  
  
  # len<-rle(endcluster1[['orient']])$lengths
  # mask<-rle(endcluster1[['orient']])$values
  # num<-which(mask=='+' & len<3)[which(mask=='+' & len<3)!=1 & which(mask=='+' & len<3)!=length(len)]
  # for(i in num){
  #   len[i-1]<-len[i-1]+len[i]
  # }
  # minus_rows <- which(mask=='-' & len>=1)
  # minus_rows <-minus_rows[minus_rows!=1 & minus_rows!=length(len)]
  # for(i in minus_rows){
  #   i = minus_rows
  #   from=sum(len[1:i-1])+1
  #   loc<- seq(from = from, by = 1, length.out = len[i])
  #   endcluster1[loc,]$ref_start<-min(endcluster1[loc,]$ref_start)
  #   endcluster1[loc,]$ref_end<-max(endcluster1[loc,]$ref_end)
  #   endcluster1[loc,]$query_start<-min(endcluster1[loc,]$query_start)
  #   endcluster1[loc,]$query_end<-max(endcluster1[loc,]$query_end)
  #   endcluster1[loc,]$cluster<-0
  #   endcluster1[loc,]$orient<-'-'
  # }
  # endcluster1<-distinct(endcluster1)
  return(endcluster1)
}



docall<-function(data){
  if(length(unique(data$query_chr)[unique(data$query_chr)!='0'])!=1){
    data<-do.call(rbind,mget(unique(data$query_chr)[unique(data$query_chr)!='0'], envir = .GlobalEnv))
  }
  else{
    data<-get( unique(data$query_chr)[unique(data$query_chr)!='0'])
  }
  return(data)
}




inversion.extract<-function(endcluster1,chrid){
  endcluster1$query_start<-abs(endcluster1$query_start)
  endcluster1$query_end<-abs(endcluster1$query_end)
  inversion<-endcluster1[endcluster1$orient=="-",] %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),
              ref_end=max(ref_end),
              query_chr=query_chr,
              query_starttem=min(pmin(query_end,query_start)),
              query_endtem=max(pmax(query_start,query_end)),
              (names(which.max(table(orient)))))
  colnames(inversion)[6]<-"query_start"
  colnames(inversion)[7]<-"query_end"
  inversion<-distinct(inversion) 
  return(inversion)
  
}


clusterbigminus<-function(endcluster1,cluster.id){
  rm("xx")
  pos_end<-endcluster1 %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))) ###正链负链？？重新算一下
  # pos_end<-pos_end[,-1]
  # pos_end<-pos_end[,-7]
  pos_end<-distinct(pos_end)
  pos_end<-pos_end[order(pos_end$ref_start),]
  if(cluster.id==0){
    newpos_end<-pos_end
    for(chr_child in unique(pos_end$query_chr)){
      minus<-which(pos_end$`(names(which.max(table(orient))))`=="-" & pos_end$query_chr==chr_child)
      if(length(minus)!=1 &length(minus)!=0){
        for(i in 2:length(minus)){
          if((minus[i]==minus[i-1]+1 |(minus[i]==minus[i-1]+2)) & pos_end[minus[i],]$query_end<=pos_end[minus[i-1],]$query_start ){ ## 说明这两个-相邻且为反转
            if((minus[i]==minus[i-1]+2)){
              endcluster1[endcluster1$cluster==newpos_end[minus[i-1]+1,]$cluster,]$cluster<-newpos_end[minus[i-1],]$cluster
              newpos_end[minus[i-1]+1,]$cluster<-newpos_end[minus[i-1],]$cluster
            }
            endcluster1[endcluster1$cluster==newpos_end[minus[i],]$cluster,]$cluster<-newpos_end[minus[i-1],]$cluster
            newpos_end[minus[i],]$cluster<-newpos_end[minus[i-1],]$cluster
            colnames(newpos_end)[8]<-"orient"
            xx<-newpos_end %>% group_by(cluster)%>%  #聚cluster
              summarise(ref_chr=ref_chr,ref_start=min(ref_start),ref_end=max(ref_end),query_chr=query_chr,query_start=min(query_start),query_end=max(query_end),(names(which.max(table(orient))))) ###
            xx<-distinct(xx)
            
          }
        }
        
      }
    }
  }
  if(exists("xx")){
    return(list(endcluster1=endcluster1,pos_end=xx))
  }else{
    return(list(endcluster1=endcluster1,pos_end=pos_end))
  }
}
reverse.region<-function(endcluster1,chrid,cluster.id){
  endcluster1$query_start<-abs(endcluster1$query_start)
  endcluster1$query_end<-abs(endcluster1$query_end)
  ## cluster.id筛选的为大聚类还是小聚类
  if(cluster.id==0){
    ## 如果大聚类的cluster中包含了一些小聚类我们把它划为大的
    ## 如果交叉就不管
    for(chr_child in unique(endcluster1$query_chr)){
      x=table(endcluster1[endcluster1$query_chr==chr_child,]$cluster)
      thereshod<-quantile(sort(as.vector(x)),0.8)/5
      if((thereshod)<5){ ## eg，1  10  11   2   3   4   5   6   7   8   9 
        # 390   1   2   1   2 838   1   1   2   1   4 
        thereshod=5
      }
      library(stringr)
      cross.region<-cross.calcu(endcluster1)
      # pattern <- "(\\d)\\1{2,}((\\d)\\2{0,1}\\1{1,2})*\\2{3,}"
      # ## 1112211222
      # cross.region<-c()
      # matches <- str_match_all(paste(endcluster1$cluster, collapse = ''), pattern)
      # if (length(matches[[1]]) > 0) {
      #   for (row_index in 1:nrow(matches[[1]])) {
      #     if(matches[[1]][row_index,][2]!=matches[[1]][row_index,][3]){
      #       cross.region<-append(cross.region,matches[[1]][row_index,][2])
      #       cross.region<-append(cross.region,matches[[1]][row_index,][3])
      #     }
      #   }
      # }
      print(cross.region)
      for(k in setdiff(names(x)[x<=thereshod],cross.region)){
        #print(k)
        if(k==1){
          next
        }
        for(j in names(x)[x>thereshod]){
          a<-which(endcluster1$cluster==k & endcluster1$query_chr==chr_child )
          b<-which(endcluster1$cluster==j & endcluster1$query_chr==chr_child)
          if(min(b)<min(a) & max(b)>max(a)){
            cat(k,j,"\n")
            endcluster1[(which(endcluster1$cluster==k & endcluster1$query_chr==chr_child)),]$cluster<-j
            endcluster1<-distinct(endcluster1)
            break
          }
          
        }
        
      }
      
      
      
      
      x=table(endcluster1[endcluster1$query_chr==chr_child,]$cluster)
      thereshod<-quantile(sort(as.vector(x)),0.8)/5
      if((thereshod)<5){ 
        thereshod=5
      }
      rm("delsytenic")
      
      for(k in names(x)[x<=thereshod]){
        if(k==1){
          next
        }
        if(max(which(endcluster1$cluster==k))==dim(endcluster1)[1]){
          next
        }
        
        idd<-(which(endcluster1$cluster==k & endcluster1$query_chr==chr_child))
        if(min(idd)==1 |max(idd)==dim(endcluster1)[1]){
          next
        }
        if(endcluster1[min(idd),]$orient=="-"){next}
        if((endcluster1[min(idd)-1,]$orient=="-" & endcluster1[max(idd)+1,]$orient=="-") |(endcluster1[max(idd)+1,]$orient=="-" & endcluster1[min(idd),]$ref_end>endcluster1[max(idd)+1,]$ref_start)){
          if((endcluster1[min(idd),]$query_end<=endcluster1[min(idd)-1,]$query_start) &(endcluster1[max(idd),]$query_start>=endcluster1[max(idd)+1,]$query_end)){
            if(abs(endcluster1[min(idd),]$query_end-endcluster1[min(idd)-1,]$query_start)<abs(endcluster1[max(idd),]$query_start-endcluster1[max(idd)+1,]$query_end)){
              endcluster1[idd,]$cluster<-endcluster1[min(idd)-1,]$cluster
            }
            else{
              endcluster1[idd,]$cluster<-endcluster1[max(idd)+1,]$cluster
            }
          }
          
          else if((endcluster1[min(idd),]$query_start>=endcluster1[min(idd)-1,]$query_end) &(endcluster1[max(idd),]$query_end<=endcluster1[max(idd)+1,]$query_start)){
            if(abs(endcluster1[min(idd),]$query_start-endcluster1[min(idd)-1,]$query_end)<abs(endcluster1[max(idd),]$query_end-endcluster1[max(idd)+1,]$query_start)){
              endcluster1[idd,]$cluster<-endcluster1[min(idd)-1,]$cluster
            }
            else{
              endcluster1[idd,]$cluster<-endcluster1[max(idd)+1,]$cluster
            }
          }
          
          else{
            if(exists("delsytenic")){
              delsytenic<-rbind(delsytenic,endcluster1[idd,])
            }
            else{
              delsytenic<-endcluster1[idd,]
            }
            print(delsytenic)
            endcluster1<-endcluster1[-idd,]
          }
        }
        
        
      }
    }
    
    

    rm("dupli")
    for (iteration in 1:3) {
      del.list<-c()
      for(chr_child in unique(endcluster1$query_chr)){
        for(k in unique(endcluster1[endcluster1$query_chr==chr_child,]$cluster)[-1]){
          datalist<-which(endcluster1$cluster==k  & endcluster1$query_chr==chr_child)
          if(endcluster1[min(datalist),]$ref_start<endcluster1[min(datalist)-1,]$ref_end){
            intervalue<-abs(endcluster1[min(datalist),]$ref_start-endcluster1[min(datalist)-1,]$ref_end)
            allval<-endcluster1[min(datalist),]$ref_end-endcluster1[min(datalist),]$ref_start
            allval2<-endcluster1[min(datalist)-1,]$ref_end-endcluster1[min(datalist)-1,]$ref_start
            if(intervalue/allval>0.5){
              del.list<-append(del.list,min(datalist))
            }
            else if(intervalue/allval2>0.5){
              del.list<-append(del.list,min(datalist)-1)
            }
            else{
              
            }
          }
        }
        
      }
      if(length(del.list)!=0){
        if(exists("dupli")){
          dupli<-rbind(dupli,endcluster1[del.list,])
        }else{dupli<-endcluster1[del.list,]}
        #endcluster1<-endcluster1[-del.list,]
      }
    }
    if(!exists("dupli")){dupli<-"no"}
    
    # for(k in unique(endcluster1[endcluster1$orient=='-',]$cluster)){
    #   if(k==1){
    #     next
    #   }
    #   a<-endcluster1[min(which(endcluster1$cluster==k))-1,]$cluster
    #   print('start')
    #   print(k)
    #   print(a)
    #   b<-endcluster1[max(which(endcluster1$cluster==k))+1,]$cluster
    #   print(b)
    #   if((!is.na(b)& length(b)!=0)  &  (!is.na(a) & length(a)!=0 ) ){
    #     if(a==b){
    #       
    #     }
    #   }
    #   
    #   #if(nrow(endcluster1[endcluster1$cluster==k,])!=0){endcluster1[min(which(endcluster1$cluster==k,)):max(which(endcluster1$cluster==k,)),"cluster"]=k}
    # }
    # for(k in unique(endcluster1[endcluster1$orient=='+',]$cluster)){
    #   if(k==1){
    #     next
    #   }
    #   a<-endcluster1[min(which(endcluster1$cluster==k))-1,]$cluster
    #   print('start')
    #   print(k)
    #   print(a)
    #   b<-endcluster1[max(which(endcluster1$cluster==k))+1,]$cluster
    #   print(b)
    #   if((!is.na(b)& length(b)!=0)  &  (!is.na(a) & length(a)!=0 ) ){
    #     if(a==b){
    #       endcluster1[(which(endcluster1$cluster==k)),]$cluster<-a
    #     }
    #   }
    #   #if(nrow(endcluster1[endcluster1$cluster==k,])!=0){endcluster1[min(which(endcluster1$cluster==k,)):max(which(endcluster1$cluster==k,)),"cluster"]=k}
    # }
    if(cluster.id==1){
      invdata<-endcluster1[endcluster1$orient=="-",]
      for(k in unique(invdata$cluster)){
        if(k==1){
          next
        }
        temdata<-invdata[invdata$cluster==k,]
        mark<-which(endcluster1$cluster==k & endcluster1$orient=="-")
        if(endcluster1[min(mark),]$ref_start<endcluster1[min(mark)-1,]$ref_end | endcluster1[min(mark),]$query_start<endcluster1[min(mark)-1,]$query_end){
          endcluster1[mark,]$cluster<-endcluster1[min(mark)-1,]$cluster
        }
        if(endcluster1[max(mark),]$ref_end>endcluster1[max(mark)+1,]$ref_start | endcluster1[min(mark),]$query_end>endcluster1[min(mark)+1,]$query_start){
          endcluster1[mark,]$cluster<-endcluster1[max(mark)+1,]$cluster
        }
      }
    }
    
    #dotplot_cluster(endcluster1)
  }
  
  testminus<-clusterbigminus(endcluster1,cluster.id)
  endcluster1<-testminus$endcluster1
  pos_end<-testminus$pos_end
  pos_end<-pos_end[order(pos_end$ref_start),]
 
  if(nrow(endcluster1)!=1){
    middledata<-duplication_extract(endcluster1,pos_end) ##对大片段之间的重复进行处理
    pos_end<-middledata$pos_end
    
    if(exists("dupli")){
      if(is.character(dupli)){
        dupli<-middledata$dupli
      }else{
        dupli<-rbind(dupli,middledata$dupli)
      }
    }else{
      dupli<-middledata$dupli
    }
  }

  
 

  
  
  reall<-merge(pos_end,chrpc,by=c("ref_chr","query_chr"))
  reall<-reall[order(reall$ref_start),]
  for(chr_child in unique(endcluster1$query_chr)){

    chrchch<-reall[reall$query_chr==chr_child,]
    new_row <- data.frame(ref_chr=chrid,ref_start = 0, ref_end = as.numeric(chrchch$ref_start[1])-1,query_chr=chr_child,query_start = 0, query_end = as.numeric(min(as.numeric(chrchch$query_start))-1)) ##改为初始queryend值为querystart最小值减1
    
    
    if(cluster.id==0){
      #note<-c()
      # if(dim(chrchch)[1]!=1){
      #   for(i in 1:(dim(chrchch)[1]-1)){
      #     if(chrchch[i+1,]$ref_start<chrchch[i,]$ref_end & chrchch[i+1,]$ref_start>chrchch[i,]$ref_start){
      #       intersectrefstart<-chrchch[i+1,]$ref_start
      #       intersectrefend<-min(chrchch[c(i,i+1),]$ref_end)
      #       if(nrow(endcluster1[endcluster1$query_chr==chr_child & endcluster1$ref_start>=intersectrefstart &endcluster1$ref_end<= intersectrefend,])!=0)
      #       {
      #         if(is.character(dupli)){
      #           dupli<-clusterall(endcluster1[endcluster1$query_chr==chr_child & endcluster1$ref_start>=intersectrefstart &endcluster1$ref_end<= intersectrefend,])
      #         }else{
      #           dupli<-rbind(dupli,clusterall(endcluster1[endcluster1$query_chr==chr_child & endcluster1$ref_start>=intersectrefstart &endcluster1$ref_end<= intersectrefend,]))
      #         }
      #         note<-append(note,i)}
      #     }
      #   }
      # }
      
      if(dim(chrchch)[1]!=1){
        for (i in 1:dim(chrchch)[1]){
          # if(length(note)!=0){
          #   if(i %in% note){
          #     next
          #   }
          # }
         
          nowvalue=chrchch$query_end[i]
        
          if(nowvalue==max(chrchch$query_end)){
            if(i==dim(chrchch)[1]){ 
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]),  query_end = chrchch$query_len[i])
            }
            else{
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]),  query_end = chrchch$query_len[i])
            }
          }else{
            query_endlist<-chrchch$query_start[which(chrchch$query_start>=nowvalue)]
            if(length(query_endlist)==0){
              next
            }
            queendval<-query_endlist[which.min(abs(query_endlist -nowvalue))]
            
            if(i==dim(chrchch)[1]){
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = queendval)
            }else{
              xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = queendval)
            }
          }
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
        }
        
      }
      else{
        for (i in 1:dim(chrchch)[1]){
          xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = chrchch$ref_len[i],query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]),  query_end = chrchch$query_len[i])
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
        }
      }
      
    }
    if(cluster.id==2){
      if(dim(chrchch)[1]!=1){
        for (i in 2:dim(chrchch)[1]-1){
          
          if(as.numeric(chrchch$ref_end[i])==chrchch$ref_start[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end =as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = as.numeric(chrchch$query_start[i+1]))
          }else if(as.numeric(chrchch$query_end[i])==chrchch$query_start[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end = as.numeric(chrchch$query_start[i+1]))
          }else{
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i]), query_end =as.numeric(chrchch$query_start[i+1]))
          }
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
        }
      }
    }
    
    
    
    
    if(cluster.id==3){
     
      if(dim(chrchch)[1]!=1){
        for (i in 2:dim(chrchch)[1]-1){
          
          if(as.numeric(chrchch$ref_end[i])==chrchch$ref_start[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end =as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
          }else if(as.numeric(chrchch$query_start[i])==chrchch$query_end[i+1]){
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
          }else{
            xx=data.frame(ref_chr=chrchch$ref_chr[i],ref_start = as.numeric(chrchch$ref_end[i]), ref_end = as.numeric(chrchch$ref_start[i+1]),query_chr=chrchch$query_chr[i],query_start = as.numeric(chrchch$query_end[i+1]), query_end = as.numeric(chrchch$query_start[i]))
          }
          colnames(xx)<-colnames(new_row)
          new_row <-rbind(new_row,xx)
          
        }
        
      }
      
    }
    assign(chr_child,new_row)
  }
  allchr.reverse<-do.call(rbind,mget(unique(endcluster1$query_chr)))
  pointer<-which(allchr.reverse$ref_start==1)[which(allchr.reverse$ref_start==1)!=1]
  allchr.reverse[pointer,]$ref_start<-allchr.reverse[pointer-1,]$ref_start
  allchr.reverse[pointer-1,]$ref_end<-allchr.reverse[pointer,]$ref_end
  rm(list=unique(endcluster1$query_chr)) 
  
  
  if(cluster.id==2 |cluster.id==3){
    allchr.reverse<-allchr.reverse[-1,]
  }
  if(nrow(allchr.reverse)!=0){
    allchr.reverse$anno<-"SDR_NM"
  }
  if(cluster.id==0){
    if(exists("delsytenic")){
      return(list(reverse=allchr.reverse,endcluster1=endcluster1,minimaploc=reall,dup=dupli,delsytenic=delsytenic))
    }else{
      return(list(reverse=allchr.reverse,endcluster1=endcluster1,minimaploc=reall,dup=dupli))
    }
    
  }else{
    return(list(reverse=allchr.reverse,endcluster1=endcluster1))
  }
  
}



repeat.integrate<-function(data,repeatid){
  repeat_region<-data.frame()
  data$query_start<-abs(data$query_start)
  data$query_end<-abs(data$query_end)
  data<-exchange(data)
  for(kk in 1:2){
    ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
    overlapsref<-findOverlaps(ir1ref,ir1ref) 
    if(length(which(overlapsref@from!=overlapsref@to))!=0){
      inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
      if(exists("inte")){
        df_sorted <- t(apply(inte, 1, function(x) sort(x)))
        unique_df <- unique(df_sorted)
        inte <- as.data.frame(unique_df)
        deletelist<-c()
        if(nrow(inte)!=0){
          for(j in 1:dim(inte)[1]){
            if(data[inte[j,2],]$ref_start==data[inte[j,1],]$ref_end){
              deletelist<-append(deletelist,j)
            }
          }
          if(length(deletelist)!=0){
            inte<-inte[-deletelist,]
          }
        }
        
        if(nrow(inte)!=0){
          for(j in 1:dim(inte)[1]){
            a<-c()
            a<-append(a,inte[j,1])
            a<-append(a,inte[j,2])
            for(k in 1:dim(inte)[1]){
              if(inte[k,1] %in% ((min(a)):(max(a)))){
                a<-unique(append(a,inte[k,1]))
              }
              if(inte[k,2] %in% ((min(a)):(max(a)))){
                a<-unique(append(a,inte[k,2]))
              }
            }
            x<-a
            x<-sort(x,decreasing = TRUE)
            for(mm in 1:(length(x)-1)){
              data[x[mm],]  ##下一列
              data[x[mm+1],] ##上一列
              if(data[x[mm],]$ref_start<data[x[mm+1],]$ref_end & data[x[mm],]$ref_start>=data[x[mm+1],]$ref_start){
                v<-min(data[x[c(mm,mm+1)],]$ref_end)-data[x[mm],]$ref_start
                o<-data[x[mm],]$ref_end-data[x[mm],]$ref_start
                o2<-data[x[mm+1],]$ref_end-data[x[mm+1],]$ref_start
                if(o2==0 | o==0){
                  next
                }
                if(v/o>0.5){
                  repeat_region<-rbind(repeat_region,data[x[mm],])
                }
                if(v/o2>0.5){
                  repeat_region<-rbind(repeat_region,data[x[mm+1],])
                }
                if(dim(data[data$ref_start>=min(data[x[c(mm,mm+1)],]$ref_end) & data$ref_end<=-data[x[mm],]$ref_start,])[1]!=0){
                  repeat_region<-rbind(data[data$ref_start>=min(data[x[c(mm,mm+1)],]$ref_end) & data$ref_end<=-data[x[mm],]$ref_start,],)
                }
                data[x[mm+1],]$ref_end<-data[x[mm+1],]$ref_start
               
              }
            }
            
            
            
            # if(abs(max(data[x,]$ref_start)-min(data[x,]$ref_end))>10000 & length(unique(data[x,]$orient))!=2){
            #   repeat_region<-append(repeat_region,min(x))
            #   if(data[min(x),]$orient == data[max(x),]$orient){
            #     data[min(x),]$ref_start<-min(data[x,]$ref_start)
            #     data[min(x),]$ref_end<-max(data[x,]$ref_end)
            #     data[min(x),]$query_start<-min(data[x,]$query_start)
            #     data[min(x),]$query_end<-max(data[x,]$query_end)
            #     data[x,]<-data[min(x),]
            #   }
            #   
            # }
            # if(abs(max(data[x,]$ref_start)-min(data[x,]$ref_end))<10000 & length(unique(data[x,]$orient))!=2){
            #   if(max(data[x,]$ref_start)!=min(data[x,]$ref_end)){
            #     data[min(x)+1,]$ref_start<- data[min(x),]$ref_end
            #   }
            #   
            # }
          }
        }
        
 
      }
    }
    if(kk==1){
      dup<-repeat_region
    }
    data<-distinct(data)
  }
  
  # ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  # overlapsque<-findOverlaps(ir1que,ir1que)
  # if(length(which(overlapsque@from!=overlapsque@to))!=0){
  #   inte<-as.data.frame(overlapsque[which(overlapsque@from!=overlapsque@to)])
  #   if(exists("inte")){
  #     df_sorted <- t(apply(inte, 1, function(x) sort(x)))
  #     unique_df <- unique(df_sorted)
  #     inte <- as.data.frame(unique_df)
  #     for(j in 1:dim(inte)[1]){
  #       a<-c()
  #       a<-append(a,inte[j,1])
  #       a<-append(a,inte[j,2])
  #       for(k in 1:dim(inte)[1]){
  #         if(inte[k,1] %in% ((min(a)):(max(a)))){
  #           a<-unique(append(a,inte[k,1]))
  #         }
  #         if(inte[k,2] %in% ((min(a)):(max(a)))){
  #           a<-unique(append(a,inte[k,2]))
  #         }
  #       }
  #       x<-a
  #       if(abs(max(data[x,]$query_start)-min(data[x,]$query_end))>10000){
  #         #repeat_region<-append(repeat_region,min(x))
  #         data[min(x),]$ref_start<-min(data[x,]$ref_start)
  #         data[min(x),]$ref_end<-max(data[x,]$ref_end)
  #         data[min(x),]$query_start<-min(data[x,]$query_start)
  #         data[min(x),]$query_end<-max(data[x,]$query_end)
  #         data[x,]<-data[min(x),]
  #       }
  #       if(abs(max(data[x,]$query_start)-min(data[x,]$query_end))<10000){
  #         if(repeatid==2){
  #           data[min(x)+1,]$query_end<- data[min(x),]$query_start
  #         }
  #         if(repeatid==1){
  #           data[min(x)+1,]$query_start<- data[min(x),]$query_end
  #         }
  #        
  #       }
  #     }
  #   }
  # }
  
  
  if(repeatid==1){
    i=2
    while (i<=dim(data)[1]){
      list<-which(data$query_start[i]<data[1:(i-1),]$query_end)
      
      if(length(list)!=0){
        if(data$orient[i]=="-" & all(data$orient[list] == '-')){
          i<-i+1
          next
        }
        if(length(list)>0.5*dim(data)[1]){
          data<-data[-i,]
          next
        }
        print(i)
        data[i,]$ref_start<-min(data[c(list,i),]$ref_start)
        data[i,]$ref_end<-max(data[c(list,i),]$ref_end)
        data[i,]$query_start<-min(data[c(list,i),]$query_start)
        data[i,]$query_end<-max(data[c(list,i),]$query_end)
        data<-data[-list,]
 
      }
      else{
        i<-i+1
      }
      
    }
  }
  
  
  data<-data[order(data$ref_start),]
  data<-distinct(data)
  if(nrow(dup)!=0){
    dup<-dup[c("ref_chr","ref_start","ref_end","query_chr","query_start","query_end","orient","cluster")]
  }
  
  return(list(afterdup=data,repeat.region=dup))
  
}


endfilter<-function(all,chrid,chr_child){
  data<-all
  data$reflen<-data$ref_end-data$ref_start
  data$querylen<-data$query_end-data$query_start
  if(nrow(data[(data$reflen<10000 & data$querylen<10000) & data$anno=="SDR_NM",])!=0){
    data[(data$reflen<10000 & data$querylen<10000) & data$anno=="SDR_NM",]$anno<-"SV_NM"
  }
  
  if(length(which(data$reflen==0 & data$querylen>=10000))!=0){
    data[data$reflen==0 & data$querylen>=10000,]$anno<-"SDR_INS"
  }
  if(length(which(data$reflen==0 & data$querylen<10000))!=0){
    data[data$reflen==0 & data$querylen<10000,]$anno<-"SV_INS"
  }
  if(length(which(data$querylen==0 & data$reflen>=10000))!=0){
    data[data$querylen==0 & data$reflen>=10000,]$anno<-"SDR_DEL"
  }
  if(length(which(data$querylen==0 & data$reflen<10000))!=0){
    data[data$querylen==0 & data$reflen<10000,]$anno<-"SV_DEL"
  }
  
  if(length(which(data$anno=="INV" & (data$reflen<10000 | data$querylen<10000)))!=0){
    data[data$anno=="INV" & (data$reflen<10000 | data$querylen<10000),]$anno<-"SV_INV"
  }
  
  if(length(which(data$anno=="INV" & (data$reflen>=10000 | data$querylen>=10000)))!=0){
    data[data$anno=="INV" & (data$reflen>=10000 | data$querylen>=10000),]$anno<-"SDR_INV"
  }
  if(length(which(data$anno=="TRANS" & (data$reflen>=10000 | data$querylen>=10000)))!=0){
    data[data$anno=="TRANS" & (data$reflen>=10000 | data$querylen>=10000),]$anno<-"SDR_TRANS"
  }
  if(length(which(data$anno=="TRANS" & (data$reflen<10000 | data$querylen<10000)))!=0){
    data[data$anno=="TRANS" & (data$reflen<10000 | data$querylen<10000),]$anno<-"SV_TRANS"
  }
  if(length(which(data$anno=="DUP" & (data$reflen>=10000 | data$querylen>=10000)))!=0){
    data[data$anno=="DUP" & (data$reflen>=10000 | data$querylen>=10000),]$anno<-"SDR_DUP"
  }
  if(length(which(data$anno=="DUP" & (data$reflen<10000 | data$querylen<10000)))!=0){
    data[data$anno=="DUP" & (data$reflen<10000 | data$querylen<10000),]$anno<-"SV_DUP"
  }
  
  
  # teloquery<-result_chantelo[result_chantelo$chr==chr_child,]
  # m<-teloquery[teloquery$id=="end",]
  # n<-teloquery[teloquery$id=="first",]
  # if(nrow(m)!=0){ ##说明是端粒末端
  #   data[nrow(data),6]<-m$query_start
  # }
  # if(nrow(n)!=0){
  #   data[1,5]<-n$query_end
  # }
  # teloref<-result_hmtelo[result_hmtelo$chr==chrid,]
  # m<-teloref[teloref$id=="end",]
  # n<-teloref[teloref$id=="first",]
  # if(nrow(m)!=0){ 
  #   data[nrow(data),3]<-m$ref_start
  # }
  # if(nrow(n)!=0){
  #   data[1,2]<-n$ref_end
  # }
  # centroref<-result_hmcentr[result_hmcentr$chr==chrid,]
  # centroquery<-result_chancentr[result_chancentr$chr==chr_child,]
  # #data<-data[!((data$ref_start %in% (centroref$ref_start:centroref$ref_end)) & (data$ref_end %in% (centroref$ref_start:centroref$ref_end))),]
  # #data<-data[!((data$query_start %in% (centroquery$query_start:centroquery$query_end)) & (data$query_end %in% (centroquery$query_start:centroquery$query_end))),]
  # ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  # ir2 <- IRanges(start =centroref$ref_start , end =centroref$ref_end)
  # overlapsref<-findOverlaps(ir1ref,ir2) 
  # ir1que  <- IRanges(start = data$query_start, end = data$query_end)
  # ir3 <- IRanges(start =centroquery$query_start , end =centroquery$query_end)
  # overlapsquery<-findOverlaps(ir1que,ir3) 
  # ref<-overlapsref@from
  # que<-overlapsquery@from
  # print(overlapsquery)
  # data$note<-"a"
  # ## 
  # if(length(intersect(ref,que))!=0){
  #   for(l in intersect(ref,que)){
  #     
  #     if(data[l,]$ref_start> centroref$ref_start & data[l,]$ref_end<centroref$ref_end)
  #     {data$note[l]<-'del'}
  # 
  #     if(data[l,]$ref_start< centroref$ref_start & data[l,]$ref_end>centroref$ref_end)
  #     {data$note[l]<-'del'
  #     data<-rbind(data,data[l,])
  #     data<-rbind(data,data[l,])
  #     data[dim(data)[1]-1,]$ref_end<-centroref$ref_start-1
  #     data[dim(data)[1],]$ref_start<-centroref$ref_end+1
  #     if(data[l,]$query_start< centroquery$query_start & data[l,]$query_end>centroquery$query_end){
  #       data[dim(data)[1]-1,]$query_end<-centroquery$query_start-1
  #       data[dim(data)[1]-1,]$note<-"a"
  #       data[dim(data)[1],]$query_start<-centroquery$query_end+1
  #       data[dim(data)[1],]$note<-"a"
  #     }
  #     }
  #   }
  # }
  # ## 
  # setdiff(ref,que)
  # ## 
  # setdiff(que,ref)
  # 
  # if(length(which(data$note=="del"))!=0){
  #   data<-data[-which(data$note=="del"),]
  # }
  # 
  # data<-data[order(data$ref_start),]
  # data<-data[,-10]
  # 
  
  return(data)
}

insertsmall<-function(endcluster2before,storesmall,orientid){
  orientlist<-rle(endcluster2before[["orient"]] == orientid)
  if(orientid=="+"){ 
    reverseid=2
  }
  else{
    reverseid=3
  }
  for (i in which(orientlist$values)){
    if(i!=1){
      start<-(sum(orientlist$lengths[1:(i-1)])+1)
      end<-(sum(orientlist$lengths[1:(i-1)]))+orientlist$lengths[i]
      insertpos=endcluster2before[start:end,]
      reverse_end<-reverse.region(insertpos,chrid,reverseid)
      if(nrow(reverse_end$reverse)!=0){
        xx<-reverse_end$reverse
        if(reverseid==2){
          xx$orient<-"+"
        }else{
          xx$orient<-"-"
        }
        
        storesmall<-rbind(storesmall,xx[,colnames(storesmall)])
        print(storesmall)
      }
      
    }
    # else{
    #   insertpos=data[1:orientlist$lengths[i],]
    #   reverse_end<-reverse.region(insertpos,chrid,reverseid)
    #   if(nrow(reverse_end$reverse)!=0){
    #     storesmall<-rbind(storesmall,reverse_end$reverse[,colnames(storesmall)])
    #     print(storesmall)
    #   }
    #   
    # }
  }
  return(storesmall)
}


smalltransf<-function(endcluster0){
  x=table(endcluster0$cluster)
  thereshod<-quantile(sort(as.vector(x)),0.8)/5
  dellist<-c()
  for(k in names(x)[x<=thereshod]){
    if(k==1){
      next
    }
    numberslist<-rle(endcluster0$cluster)$values
    biglist<-names(x)[x>thereshod]
    tranloc<-which(numberslist==k)
    bigloc<-which(numberslist %in% biglist)
    positive_positions <- which(min(tranloc)-bigloc > 0)
    aclus<-bigloc[positive_positions][which.min(min(tranloc)-bigloc[positive_positions])]
    positive_positions <- which(bigloc-max(tranloc) > 0)
    bclus<-bigloc[positive_positions][which.min(min(tranloc)-bigloc[positive_positions])]
    s<-numberslist[aclus]
    e<-numberslist[bclus]
    a<-which(endcluster0$cluster==k)
    loc1<-max(intersect(1:min(a), which(endcluster0$cluster==s)))
    loc2<-min(intersect(max(a):dim(endcluster0)[1], which(endcluster0$cluster==e)))
    x1<-abs(endcluster0[min(a),]$query_start-endcluster0[loc1,]$query_start)
    x2<-abs(endcluster0[max(a),]$query_end-endcluster0[loc2,]$query_end)
    if(unique(is.na(x2)) |length(x2)!=1){
      if(x1>2000000){
        dellist<-append(dellist,k)
      }
    }
    if( unique(is.na(x1)) |length(x1)!=1){
      if(x2>2000000){
        dellist<-append(dellist,k)
      }
    }
    if(unique(!is.na(x1)) & unique(!is.na(x2)) &length(x2)==1  &length(x1)==1){
      if(min(x1,x2)>2000000){
        dellist<-append(dellist,k)
      }
    }
       
  }
  if(length(dellist)!=0){
    tran<-endcluster0[endcluster0$cluster %in% dellist,]
    tranbefore<-tran
    tran<-tran %>% group_by(cluster)%>%  #聚cluster
      summarise(ref_chr=ref_chr,ref_start=min(ref_start),
                ref_end=max(ref_end),
                query_chr=query_chr,
                query_starttem=min(pmin(query_end,query_start)),
                query_endtem=max(pmax(query_start,query_end)),
                (names(which.max(table(orient))))) ###正链负链？？重新算一下
    colnames(tran)[6]<-"query_start"
    colnames(tran)[7]<-"query_end"
    endcluster0<-endcluster0[!endcluster0$cluster %in% dellist,]
    colnames(tran)[8]<-"orient"
    return(list(tran=distinct(tran),endcluster0=endcluster0,tranbefore=tranbefore))
  }
  else{
    return(list(tran=NULL,endcluster0=endcluster0))
  }
  
}
clusterall<-function(data){
  data<-data %>% group_by(cluster)%>%  #聚cluster
    summarise(ref_chr=ref_chr,ref_start=min(ref_start),
              ref_end=max(ref_end),
              query_chr=query_chr,
              query_starttem=min(pmin(query_end,query_start)),
              query_endtem=max(pmax(query_start,query_end)),
              (names(which.max(table(orient))))) ###正链负链？？重新算一下
  colnames(data)[6]<-"query_start"
  colnames(data)[7]<-"query_end"
  colnames(data)[8]<-"orient"
  data<-distinct(data)
  return(data)
}

smallcluster<-function(data,orientid){
  orientlist<-rle(endcluster2[["orient"]] == orientid)
  for (i in which(orientlist$values)){
    if(i!=1){
      start<-(sum(orientlist$lengths[1:(i-1)])+1)
      end<-(sum(orientlist$lengths[1:(i-1)]))+orientlist$lengths[i]
    }
    else{
      start<-1
      end<-orientlist$lengths[i]
    }
    data[start,]$ref_start<-min(data[start:end,]$ref_start)
    data[start,]$ref_end<-max(data[start:end,]$ref_end)
    data[start,]$query_start<-min(data[start:end,]$query_start)
    data[start,]$query_end<-max(data[start:end,]$query_end)
    data[start,]$ref_pos<-(data[start,]$ref_start+data[start,]$ref_end)/2
    data[start,]$query_pos<-(data[start,]$query_start+data[start,]$query_end)/2
    data[start:end,]<-data[start,]
  }
  return(distinct(data))
}

complex<-function(data){
  for(i in 2:dim(data)[1]){
    if(data[i,]$ref_start<data[i-1,]$ref_end){
      data[i,]$ref_start<-min(data[i,]$ref_start,data[i-1,]$ref_start)
      data[i,]$ref_end<-max(data[i,]$ref_end,data[i-1,]$ref_end)
      data[i,]$query_start<-min(data[i,]$query_start,data[i-1,]$query_start)
      data[i,]$query_end<-max(data[i,]$query_end,data[i-1,]$query_end)
      data[c(i,i-1),]<-data[i,]
      data[c(i,i-1),]$anno<-'COMPLEX'
    }
  }
  data<-distinct(data)
  for(i in 2:dim(data)[1]){
    if(data[i,]$query_start<data[i-1,]$query_end){
      data[i,]$ref_start<-min(data[i,]$ref_start,data[i-1,]$ref_start)
      data[i,]$ref_end<-max(data[i,]$ref_end,data[i-1,]$ref_end)
      data[i,]$query_start<-min(data[i,]$query_start,data[i-1,]$query_start)
      data[i,]$query_end<-max(data[i,]$query_end,data[i-1,]$query_end)
      data[c(i,i-1),]<-data[i,]
      data[c(i,i-1),]$anno<-'COMPLEX'
    }
  }
  data<-distinct(data)
  list<-rle(data$anno)
  for (j in which(list$values=="COMPLEX")){
    if(j!=1){
      start<-sum(list$lengths[1:(j-1)])+1
      end<-start+list$lengths[j]-1
      
    }
    else{
      start<-1
      end<-list$lengths[j]
    }
    data[start,]$ref_start<-min(data[start:end,]$ref_start)
    data[start,]$ref_end<-max(data[start:end,]$ref_end)
    data[start,]$query_start<-min(data[start:end,]$query_start)
    data[start,]$query_end<-max(data[start:end,]$query_end)
    data[start:end,]<-data[start,]
  }
  return(distinct(data))
}
minus.next<-function(endcluster1){
  m=0
  for(clusid in unique(endcluster1[endcluster1$orient=="-",]$cluster)){
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,200000)
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster$cluster<-paste(miuscluster$cluster, "1",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    m=m+1
  }
  for(clusid in unique(endcluster1[endcluster1$orient=="+",]$cluster)){
    miuscluster<-split_region(endcluster1[endcluster1$cluster==clusid,],1,200000)
    if(length(unique(miuscluster$cluster))!=1){
      miuscluster$cluster<-paste(miuscluster$cluster, "1",as.character(m), sep = "")
      endcluster1<-endcluster1[endcluster1$cluster!=clusid,]
      endcluster1<-rbind(endcluster1,miuscluster)
    }
    m=m+1
  }
  endcluster1<-endcluster1[order(endcluster1$ref_start),]
  return(endcluster1)
}
cross.calcu<-function(endcluster1){
  del<-c()
  list<-unique(endcluster1$cluster)
  for(i in 2:(length(unique(endcluster1$cluster))-1)){
    a<-which(endcluster1$cluster==list[i])
    b<-which(endcluster1$cluster==list[i-1])
    if(max(b)>max(a) & min(b)>min(a) &  max(a)>min(b)){
      del<-append(del,list[i])
      del<-append(del,list[i-1])
    }
    c<-which(endcluster1$cluster==list[i+1])
    if(max(c)>max(a)& min(c)>min(a) & max(a)>min(c)){
      del<-append(del,list[i])
      del<-append(del,list[i+1])
    }
    
  }
  return(unique(del))
}
duplication_extract<-function(endcluster1,pos_end){
  data<-pos_end
  rm("inte")
  ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
  overlapsref<-findOverlaps(ir1ref,ir1ref) 
  if(length(which(overlapsref@from!=overlapsref@to))!=0){
    inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    }
  duplilist<-data.frame()
  if(exists("inte")){
    df_sorted <- t(apply(inte, 1, function(x) sort(x)))
    unique_df <- unique(df_sorted)
    inte <- as.data.frame(unique_df)
    if(nrow(inte)!=0){
      for(j in 1:dim(inte)[1]){
        if(data[inte[j,2],]$ref_start==data[inte[j,1],]$ref_end){
          next
        }
        if(data[inte[j,2],]$ref_start<data[inte[j,1],]$ref_end & data[inte[j,2],]$ref_start>=data[inte[j,1],]$ref_start){
          startregion<-data[inte[j,2],]$ref_start
          endregion<-min(data[inte[j,2],]$ref_end,data[inte[j,1],]$ref_end)
          smalldup<-endcluster1[endcluster1$ref_start>=startregion && endcluster1$ref_end<=endregion,]
          if(nrow(smalldup)!=0){
            duplilist<-rbind(duplilist,smalldup) 
          }else{
            data[inte[j,1],]$ref_end=data[inte[j,2],]$ref_start
          }
          ir1ref <- IRanges(start = endcluster1$ref_start, end = endcluster1$ref_end)
          ir1que <- IRanges(start=startregion,end=endregion)
          overlapsref<-findOverlaps(ir1ref,ir1que) 
          for (m in overlapsref@from) {
            if(endcluster1[m,]$ref_end-endcluster1[m,]$ref_start==0){
              next
            }
            if((endregion-startregion)/(endcluster1[m,]$ref_end-endcluster1[m,]$ref_start)>0.5){
              duplilist<-rbind(duplilist,endcluster1[m,]) 
            }
          }
        }
        
      }
    }
    
    duplilist<-distinct(duplilist)
  
  }
 
  return(list(pos_end=data,dupli=duplilist))}

  
  complexinte<-function(data){
    data$cluster<-1:dim(data)[1]
    for(kk in 1:5){
      rm("inte")
      ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
      overlapsref<-findOverlaps(ir1ref,ir1ref) 
      if(length(which(overlapsref@from!=overlapsref@to))!=0){
        inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
      }
      if(exists("inte")){
        df_sorted <- t(apply(inte, 1, function(x) sort(x)))
        unique_df <- unique(df_sorted)
        inte <- as.data.frame(unique_df)
        if(nrow(inte)!=0){
          for(j in 1:dim(inte)[1]){
            data[inte[j,2],]$cluster=data[inte[j,1],]$cluster
          }
          data<-clusterall(data)
        }
      }
      
    
    }
    return(data)
    }
  
  transduplication_extract<-function(endcluster1,pos_end){
    data<-pos_end
    duplilist<-data.frame()
    for(j in 1:dim(data)[1]){
      x<-endcluster1[which(endcluster1$ref_start>=data[j,]$ref_start & endcluster1$ref_end<=data[j,]$ref_end),]
      duplilist<-rbind(duplilist,duplication_extract(x,x)$dupli)
      endcluster1<-endcluster1[-which(endcluster1$ref_start>=data[j,]$ref_start & endcluster1$ref_end<=data[j,]$ref_end),]
      }
    
    rm("inte")
    ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
    ir1req <- IRanges(start = endcluster1$ref_start, end = endcluster1$ref_end)
    overlapsref<-findOverlaps(ir1ref,ir1req) 
    if(length(which(overlapsref@from!=overlapsref@to))!=0){
      inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
    }
   
    if(exists("inte")){
      df_sorted <- t(apply(inte, 1, function(x) sort(x)))
      unique_df <- unique(df_sorted)
      inte <- as.data.frame(unique_df)
      if(nrow(inte)!=0){
        for(j in 1:dim(inte)[1]){
          if(endcluster1[inte[j,2],]$ref_start==data[inte[j,1],]$ref_end){
            next
          }
          if((endcluster1[inte[j,2],]$ref_start<data[inte[j,1],]$ref_end & endcluster1[inte[j,2],]$ref_start>data[inte[j,1],]$ref_start)|(data[inte[j,1],]$ref_start<endcluster1[inte[j,2],]$ref_end & data[inte[j,1],]$ref_start>endcluster1[inte[j,2],]$ref_start)){
            startregion<-max(endcluster1[inte[j,2],]$ref_start,data[inte[j,1],]$ref_start)
            endregion<-min(endcluster1[inte[j,2],]$ref_end,data[inte[j,1],]$ref_end)
            smalldup<-endcluster1[endcluster1$ref_start>=startregion && endcluster1$ref_end<=endregion,]
            if(nrow(smalldup)!=0){
              duplilist<-rbind(duplilist,smalldup) ##结合完全在其中的
            }
  
            if((endregion-startregion)/(data[inte[j,1],]$ref_end-data[inte[j,1],]$ref_start)>0.5){
                duplilist<-rbind(duplilist,data[inte[j,1],]) ##如果有交集的也放进去
            }
            if((endregion-startregion)/(endcluster1[inte[j,2],]$ref_end-endcluster1[inte[j,2],]$ref_start)>0.5){
              duplilist<-rbind(duplilist,endcluster1[inte[j,2],]) ##如果有交集的也放进去
            }
          }
          
        }
      }
      
      duplilist<-distinct(duplilist)
      
    }
    
    return(duplilist)}
  
  
  complexinte<-function(data){
    data$cluster<-1:dim(data)[1]
    for(kk in 1:5){
      rm("inte")
      ir1ref <- IRanges(start = data$ref_start, end = data$ref_end)
      overlapsref<-findOverlaps(ir1ref,ir1ref) 
      if(length(which(overlapsref@from!=overlapsref@to))!=0){
        inte<-as.data.frame(overlapsref[which(overlapsref@from!=overlapsref@to)])
      }
      if(exists("inte")){
        df_sorted <- t(apply(inte, 1, function(x) sort(x)))
        unique_df <- unique(df_sorted)
        inte <- as.data.frame(unique_df)
        if(nrow(inte)!=0){
          for(j in 1:dim(inte)[1]){
            data[inte[j,2],]$cluster=data[inte[j,1],]$cluster
          }
          data<-clusterall(data)
        }
      }
      
      
    }
    return(data)
  }
  
  
  
  


rm(list=ls())
setwd("C:/Users/HP/Desktop/Site_extraction/")
library(dplyr)

# 读csv表
table_all <- read.csv(file = "HQN9.csv",header = T,quote="",sep=",",encoding = 'utf8')
table_PB2 <- table_all[,c(2,3,6)]


PB2_position <- c(74,191,402,511,535,559,570,627,647,701)


site_extraction <- function(table_name,segname){
  write.table(paste("===============segment===============:",segname),file="results.txt",append = T,col.names = F,row.names = F)
  #print(paste("===============segment===============:",segname))
  for(i in get(paste(segname,"_position",sep = ""))){
    position = i
    #print(paste("===============position===============:",position))
    write.table(paste("===============position===============:",position),file="results.txt",append = T,col.names = F,row.names = F)
    # 筛选所有Host=="Human"的序列，重新生成数据框
    table_human <- table_name %>% select(Host,segname) %>% filter(Host=="Human")

    # 提取待分析位点列
    x = ''
    for(i in 1:length(table_human[,2])){
      x[i] <- substr(table_human[,2][i],position,position)
    }
    # 输出结果
    table_human_count <- as.data.frame(x)
    table_human_count <- table(table_human_count$x)
    table_human_count <- as.data.frame(table_human_count)
    write.table("table_human_summary:",file="results.txt",append=T,col.names = F,row.names = F)
    write.table(table_human_count,file="results.txt",append=T,col.names = F,row.names = F)
    #print("table_human_summary:")
    #print(table_human_count)
    
    
    table_bird <- table_name %>% select(Host,segname) %>% filter(Host %in% c('Chicken','Duck','Goose','Pigeon'))
    x_2 = ''
    for(i in 1:length(table_bird[,2])){
      x_2[i] <- substr(table_bird[,2][i],position,position)
    }
    table_bird_count <- as.data.frame(x_2)
    table_bird_count <- table(table_bird_count$x_2)
    table_bird_count <- as.data.frame(table_bird_count)
    write.table("table_bird_summary:",file="results.txt",append=T,col.names = F,row.names = F)
    write.table(table_bird_count,file="results.txt",append=T,col.names = F,row.names = F)
    #print("table_bird_summary:")
    #print(table_bird_count)
    
    
    table_en <- table_name %>% select(Host,segname) %>% filter(Host == 'Environment')
    x_3 = ''
    for(i in 1:length(table_en[,2])){
      x_3[i] <- substr(table_en[,2][i],position,position)
    }
    table_en_count <- as.data.frame(x_3)
    table_en_count <- table(table_en_count$x_3)
    table_en_count <- as.data.frame(table_en_count)
    write.table("table_en_summary:",file="results.txt",append=T,col.names = F,row.names = F)
    write.table(table_en_count,file="results.txt",append=T,col.names = F,row.names = F)
    #print("table_en_summary:")
    #print(table_en_count)
  }
}

site_extraction(table_PB2,"PB2")



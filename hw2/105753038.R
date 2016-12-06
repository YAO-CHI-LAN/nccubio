######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript 105753038.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<-args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-args[i+1]
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

title=paste("PARAMETERS",(paste("\ninput file         :", i_f)),(paste("\noutput file        :", o_f)))
title=paste(title,"\nscore file         :", s_f,"\nalignments         :",aln_mode)
title=paste(title,"\ngap open penalty   :", g_o,paste("\ngap extend penalty :", g_e))
######################################
# main
######################################
#part1讀取相關參數
content=as.matrix(read.table(file=i_f))#source
scorematrix=as.matrix(read.table(s_f))#評分table
#從此可任意讀入PAM250 或者 PAM100
seq1=toString(content[2])#序列1
seq2=toString(content[4])#序列2
gop=as.numeric(g_o)#扣分機制gop
gep=as.numeric(g_e)#扣分機制gep
alnm=toString(aln_mode)#對應模式
colnames(scorematrix)[24]="*" #修正讀取錯誤
sl1=nchar(seq1,type="chars")#序列1長度
sl2=nchar(seq2,type="chars")#序列2長度
R_seq1=sapply(lapply(strsplit(as.character(seq1), NULL), rev), paste, collapse="")#反轉序列1
R_seq2=sapply(lapply(strsplit(as.character(seq2), NULL), rev), paste, collapse="")#反轉序列2
local_max=0#找到的區域最大值
local=c(0,0)#區域最大值位於matrix的值1=列2=行

#part2建立Alignments matrix
sc_table=matrix(data = 0,nrow=sl2+1,ncol=sl1+1,dimnames=list(c(1:(sl2+1)),c(1:(sl1+1))))#matrix大小設定
#matrix行列名稱設定
colnames(sc_table)[1]="-"
rownames(sc_table)[1]="-"
for(i in 2:(sl1+1)){colnames(sc_table)[i]=substr(R_seq1,i-1,i-1)}
for(j in 2:(sl2+1)){rownames(sc_table)[j]=substr(R_seq2,j-1,j-1)}


#part3填入Matrix值
for(i in 1:(sl1+1)){
  for(j in 1:(sl2+1)){
    if(i==1){
      sc_table[j,i]=gep*(j-1)#第一行的設定
    }else if(j==1){
      sc_table[j,i]=gep*(i-1)#第一列設定
    }else{
      #其他j表示列,i表示行
      sc_1=sc_table[j-1,i]+gep #從右上來的Insertion
      sc_2=sc_table[j-1,i-1]+(scorematrix[rownames(sc_table)[j],colnames(sc_table)[i]])
      #從左上來的Alignment，左上原本的分數加上對應該兩個字元的分數
      sc_3=sc_table[j,i-1]+gep#從左下來的Deletion
      
      if(alnm %in% "global"){#如果alnm的模式為全域則直接比較該三個
        sc_table[j,i]=max(sc_1,sc_2,sc_3)#選擇best score
      }else{#比較模式為local與0比較
        sc_table[j,i]=max(0,max(sc_1,sc_2,sc_3))
        #找到現在最好的score位置與值
        if(sc_table[j,i]>=local_max){
          local_max=sc_table[j,i]
          local[1]=j
          local[2]=i
        }
      }
    }
  }
}

#part4開始trace-back

#開始
newseq1=""#對應後的seq1
newseq2=""#對應後的seq2
if(alnm %in% "global"){#先判斷模式為何
  #global
  #表示我們現在有sl1+1乘上sl2+1那麼大的matrix現在要從[j,i]回到[1,1]
  i=sl1+1#行數
  j=sl2+1#列數
  while(i!=1 || j!=1){#只要還沒回到[1,1]就繼續
    now=sc_table[j,i]#現在的值
    if((i!=1 && j!=1) && (now==sc_table[j-1,i-1]+scorematrix[rownames(sc_table)[j],colnames(sc_table)[i]])){
      #如果現在的值是由[i-1,j-1]來的，則XX，
      #但若當兩者(i=1 或者j=1)其中一個滿足的話會找不到從左上[i-1,j-1]的值故要排除
      #也就是該情控只發生在i,j均大於1
      newseq1=paste(newseq1,colnames(sc_table)[i],sep="")
      newseq2=paste(newseq2,rownames(sc_table)[j],sep="")
      i=i-1
      j=j-1
    }else if((i>1) && (now==sc_table[j,i-1]+gep)){
      #如果現在的值是由[j,i-1]來的則seq1 - seq2 X  
      newseq1=paste(newseq1,colnames(sc_table)[i],sep="")
      newseq2=paste(newseq2,"-",sep="")
      i=i-1
    }else if((j>1) && (now==sc_table[j-1,i]+gep)){
      #如果現在的值是由[j-1,i]來的則seq1 X seq2 -
      newseq1=paste(newseq1,"-",sep="")
      newseq2=paste(newseq2,rownames(sc_table)[j],sep="")
      j=j-1
    }
  }
}else{
  #local現在要從最高分的score trace-back到某個[j,i]=0的起點此配對為最好的local配對
  j=local[1]
  i=local[2]
  
  while(!(sc_table[j,i]==0)){
    now=sc_table[j,i]
    if((i>1 && j>1) && (now==sc_table[j-1,i-1]+scorematrix[rownames(sc_table)[j],colnames(sc_table)[i]])){
      newseq1=paste(newseq1,colnames(sc_table)[i],sep="")
      newseq2=paste(newseq2,rownames(sc_table)[j],sep="")
      i=i-1
      j=j-1
    }else if((i>1) && (now==sc_table[j,i-1]+gep)){
      #如果現在的值是由[j,i-1]來的則seq1 - seq2 X  
      newseq1=paste(newseq1,colnames(sc_table)[i],sep="")
      newseq2=paste(newseq2,"-",sep="")
      i=i-1
    }else if((j>1) && (now==sc_table[j-1,i]+gep)){
      #如果現在的值是由[j-1,i]來的則seq1 X seq2 -
      newseq1=paste(newseq1,"-",sep="")
      newseq2=paste(newseq2,rownames(sc_table)[j],sep="")
      j=j-1
    }
  }
}
#part5不管是global還是local得到新的對應序列將其寫入檔案
data=paste(">1aboA\n",newseq1,"\n>1ycsB\n",newseq2,sep = "")
write(paste(data),file=o_f)

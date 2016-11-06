pam1=read.table("pam1.txt",sep=",",header=T)
pam1=as.matrix(pam1)/10000
rownames(pam1)=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
pam250=pam1
for(x in 1:250){
  pam250=pam250%*%pam1
}
pam250=round(pam250*100,3)
write.table(pam250,"pam250.txt",sep=",")

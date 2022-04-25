
#### This code requires PLINK bed,bim,fam files and takes fam file (sample information in PLINK) as input for shuffling sample labels to break genotype-phenotype relation and then perform fisher's exact test between permuted labels. test.fam file contains 36 samples (complete population)
##### Fix the number of permutations to be performed
##### 1st iteration is performed separately and output is taken in "first_trimmed.txt", then subsequent iterations output is appended in this file
### This code performs 10,000 iterations 

shuffled_fam<-read.table("test.fam",header=F,sep="",stringsAsFactors=F)
set.seed(1) ## for reproducibility
keep<-sample(shuffled_fam$V1,size=36,replace=F) ## randomly sample 36 samples (that need to be tested as case:control)
make_pheno<-keep[1:18] ### create file with 18 samples to be latter assigned as case
make_pheno<-cbind(make_pheno,make_pheno) ### create phenotype dataframe
write.table(make_pheno,"./make_pheno.txt",row.names=F,quote=F,col.names=F)
write.table(keep,"./keep.fam",row.names=F,quote=F) 

random_out<-paste("./random",sep="")  ### Name of random permuted datasets
 system(paste("plink --keep-fam keep.fam --bfile ./test --make-bed --make-pheno make_pheno.txt '*' --out ",random_out,sep="")) ### PLINK command to assign case-control
system(paste("plink --bfile ",random_out," --allow-no-sex --assoc fisher --out ",random_out,sep="")) ### PLINK command to do fisher test between case-control
cmd<-'awk -F" "  \'{print $2 \"\t"\  $8}\' ./random.assoc.fisher > first_trimmed.txt'   ### extracting SNP and P value from fisher output file
system(cmd)
#### Now performing remaining iterations
for(i in 2:10000){
set.seed(i) ## for reproducibility
keep<-sample(shuffled_fam$V1,size=36,replace=F)
make_pheno<-keep[1:18]
make_pheno<-cbind(make_pheno,make_pheno)
write.table(make_pheno,"./make_pheno.txt",row.names=F,quote=F,col.names=F)
write.table(keep,"./keep.fam",row.names=F,quote=F) 
index<-which(shuffled_fam$V7 %in% keep)
random_out<-paste("./random",sep="")
 system(paste("plink --keep-fam keep.fam --bfile ./test --make-bed --make-pheno make_pheno.txt '*' --out ",random_out,sep=""))
system(paste("plink --bfile ",random_out," --allow-no-sex --assoc fisher --out ",random_out,sep=""))
cmd<-'awk -F" " \'{print $8}\' ./random.assoc.fisher|paste -d, first_trimmed.txt - > second_trimmed.txt'
system(cmd)
system('mv second_trimmed.txt first_trimmed.txt')
print(i)
}

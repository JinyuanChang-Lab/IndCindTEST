x <- read.csv('D:/0-xnc-du-dr/1-paper/0-code-new/new_FNN/res_red681.csv')
y <- read.csv('D:/0-xnc-du-dr/1-paper/0-code-new/new_FNN/res_red021.csv')





name <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
          'Ind', 'IT', 'Mat', 'RE', 'Uti')
alpha <- 0.02
A=matrix(0,11,11)

for (i in 1:55) {
  if (x$pv[i] <= alpha*i/55){
    index1 <- which(x$name1[i]==name)
    index2 <- which(x$name2[i]==name)
    A[index1, index2] <-1
    A[index2, index1] <-1
  }
}




rownames(A) <-  c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
                  'Ind', 'IT', 'Mat', 'RE', 'Uti')
colnames(A) <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
                 'Ind', 'IT', 'Mat', 'RE', 'Uti')
all_1 <- apply(A, 2, sum)
A1 =cbind(A,all_1)


B=matrix(0,11,11)

for (i in 1:55) {
  if (y$pv[i] <= alpha*i/55){
    index1 <- which(y$name1[i]==name)
    index2 <- which(y$name2[i]==name)
    B[index1, index2] <-1
    B[index2, index1] <-1
  }
}


rownames(B) <-  c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
                  'Ind', 'IT', 'Mat', 'RE', 'Uti')
colnames(B) <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
                 'Ind', 'IT', 'Mat', 'RE', 'Uti')

all_2 <- apply(B, 2, sum)
B1 =cbind(B,all_2)


A1
B1

write.csv(cbind(A1,B1),  paste0("D:/0-xnc-du-dr/1-paper/0-code-new/new_FNN/real_red1",".csv"))


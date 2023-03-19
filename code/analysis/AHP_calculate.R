data <- read.csv('/Users/duhyeonglee/Downloads/survey.csv')
data <- as.matrix(data)
data <- data[,4:17]
data <- as.numeric(data)
## data 가중치 재표현
data[which(data == 1)] <- 1/9
data[which(data == 2)] <- 1/7
data[which(data == 3)] <- 1/5
data[which(data == 4)] <- 1/3
data[which(data == 5)] <- 1
data[which(data == 6)] <- 3
data[which(data == 7)] <- 5
data[which(data == 8)] <- 7
data <- 1/data ## 처음부터 1:9, 2:7 ... 이렇게 안넣은 이유는 5:1 부터 데이터가 망가짐
data <- matrix(as.numeric(data), ncol = 14) 
coln <- c("AB", "AC", "AD", "BC", "BD", "CD",
          "A12", "A13", "A23",
          "B12",
          "C12",
          "D12", "D13", "D14")
colnames(data) <- coln

w <- matrix(0, ncol = 35, nrow = 10) #가중치 저장
CR <- matrix(0, ncol = 35, nrow = 5) #일치성 비율 저장
CR.test <- matrix(0, ncol = 35, nrow = 5) # 일치성 비율 검정 통과 유무 (통과 : 1)
alpha <- 0.15 #일치성비율 검정 기준값
w11 <- matrix(0, nrow = 35, ncol = 4)
w22 <- matrix(0, nrow = 35, ncol = 3)
w33 <- matrix(0, nrow = 35, ncol = 2)
w44 <- matrix(0, nrow = 35, ncol = 2)
w55 <- matrix(0, nrow = 35, ncol = 3)
for (i in 1:35){
  temp <- data[i,]
  temp1 <- diag(x=1, nrow = 4, ncol = 4)
  temp1[1,2] <- temp[1]; temp1[2,1] <- 1/temp1[1,2]
  temp1[1,3] <- temp[2]; temp1[3,1] <- 1/temp1[1,3]
  temp1[1,4] <- temp[3]; temp1[4,1] <- 1/temp1[1,4]
  temp1[2,3] <- temp[4]; temp1[3,2] <- 1/temp1[2,3]
  temp1[2,4] <- temp[5]; temp1[4,2] <- 1/temp1[2,4]
  temp1[3,4] <- temp[6]; temp1[4,3] <- 1/temp1[3,4]
  temp.ei <- max(Mod(eigen(temp1)$values))
  CR[1,i] <- (temp.ei - 4)/(4-1)/0.9
  CR.test[1,i] <- CR[1,i] < alpha
  #대문항 가중치
  w1 <- apply(temp1/matrix(rep(apply(temp1, 2, sum), 4), nrow =4, byrow = T), 1, mean)
  
  w11[i,] <- w1
  
  #소문항 가중치 및 소문항 별 CR계산
  temp.a <- diag(x = 1, nrow =3, ncol =3)
  temp.a[1,2] <- temp[7]; temp.a[2,1] <- 1/temp.a[1,2]
  temp.a[1,3] <- temp[8]; temp.a[3,1] <- 1/temp.a[1,3]
  temp.a[2,3] <- temp[9]; temp.a[3,2] <- 1/temp.a[2,3]
  temp.ei <- max(Mod(eigen(temp.a)$values))
  CR[2,i] <- (temp.ei - 3)/(3-1)/0.58
  CR.test[2,i] <- CR[2,i] < alpha
  #A소문항 가중치
  w.a <- apply(temp.a/matrix(rep(apply(temp.a, 2, sum), 3), nrow =3, byrow = T), 1, mean)
  w22[i,] <- w.a
  
  
  temp.b <- diag(x = 1, nrow = 2, ncol = 2)
  temp.b[1,2] <- temp[10]; temp.b[2,1] <- 1/temp.b[1,2]
  temp.ei <- max(Mod(eigen(temp.b)$values))
  CR[3,i] <- 0
  CR.test[3,i] <- CR[3,i] < alpha
  #B소문항 가중치
  w.b <- apply(temp.b/matrix(rep(apply(temp.b, 2, sum), 2), nrow =2, byrow = T), 1, mean)
  w33[i,] <- w.b
  
  temp.c <- diag(x = 1, nrow = 2, ncol = 2)
  temp.c[1,2] <- temp[11]; temp.c[2,1] <- 1/temp.c[1,2]
  temp.ei <- max(Mod(eigen(temp.c)$values))
  CR[4,i] <- 0
  CR.test[4,i] <- CR[4,i] < alpha
  #C소문항 가중치
  w.c <- apply(temp.c/matrix(rep(apply(temp.c, 2, sum), 2), nrow =2, byrow = T), 1, mean)
  w44[i,] <- w.c
  
  temp.d <- diag(x = 1, nrow =3, ncol =3)
  temp.d[1,2] <- temp[12]; temp.d[2,1] <- 1/temp.d[1,2]
  temp.d[1,3] <- temp[13]; temp.d[3,1] <- 1/temp.d[1,3]
  temp.d[2,3] <- temp[14]; temp.d[3,2] <- 1/temp.d[2,3]
  temp.ei <- max(Mod(eigen(temp.d)$values))
  CR[5,i] <- (temp.ei - 3)/(3-1)/0.58
  CR.test[5,i] <- CR[5,i] < alpha
  #D소문항 가중치
  w.d <- apply(temp.d/matrix(rep(apply(temp.d, 2, sum), 3), nrow =3, byrow = T), 1, mean)
  w55[i,] <- w.d
  temp.w <- c(w.a * w1[1], w.b * w1[2], w.c * w1[3], w.d * w1[4])
  
  
  ## 논문바탕 수정가중치 적용
  temp.adjw <- c(rep(w1[1], 3), rep(w1[2],2), rep(w1[3],2), rep(w1[4],3)) * temp.w * c(rep(3,3), rep(2,2), rep(2,2), rep(3,3))
  w[,i] <- temp.adjw/sum(temp.adjw)
}
# 대문항 CR.test 를 통과한 표본 중 0인 소문항은 제외하고 1인 소문항만 사용하는 코드
## CR.test를 통과한 가중치 계산에 필요한 행렬 가져오기
w.calc <- w[,which(CR.test[1,]==1)]
CR.test.calc <- CR.test[,which(CR.test[1,]==1)]

Indicator.calc <- matrix(1, nrow = 10, ncol = 8)
Indicator.calc[1:3,which(CR.test.calc[2,] == 0)] = 0
Indicator.calc[4:5,which(CR.test.calc[3,] == 0)] = 0
Indicator.calc[6:7,which(CR.test.calc[4,] == 0)] = 0
Indicator.calc[8:10,which(CR.test.calc[5,] == 0)] = 0

sum(CR.test.calc[2,]==1)


# 대문항과 소문항 모두 CR.test 통과한 표본만 쓰는 코드

matrix(apply(w[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)/sum(matrix(apply(w[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)
)

#matrix(apply(t(w11)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)/sum(matrix(apply(t(w11)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1))
#matrix(apply(t(w22)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)/sum(matrix(apply(t(w22)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)
#)
#matrix(apply(t(w33)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)/sum(matrix(apply(t(w33)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)
#)
#matrix(apply(t(w44)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)/sum(matrix(apply(t(w44)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)
#)
#matrix(apply(t(w55)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)/sum(matrix(apply(t(w55)[,which(apply(CR.test, 2, function(x){prod(x)})==1)],1,function(x){prod(x)^{1/length(x)}}),ncol = 1)
#)







library(umap)
#----data read ----
RNA1 <- read.csv('RNA_c.csv',header = TRUE)
rownames(RNA1) <- RNA1[,1]
RNA1 <- RNA1[,-1]
Z <- as.matrix(RNA1)
G <- nrow(Z)
M <- ncol(Z)
#----hyper-paramters----
eta_h <- 0
tau_h <- 5
alpha_1 <- 3
beta_1 <- 0.1

#----iteration----
num_iter = 5000
set.seed(1234)
#---- K=5 ----
K = 5
gam_1 <- rep(2, K)
C_t <- sample(1:K, M, replace = TRUE)
h_t <- matrix(NA, G, K)
for(k in 1:K){
  h_t[, k] <- rowMeans(Z[ , C_t == k])
}
sgm_sq_t_1 <- apply(Z - h_t[, C_t], 1, var)

H_T1 <- array(0, dim = c(G, K,num_iter))
Sgm_sq_t1 <- matrix(0, G,num_iter)
C_T1 <- matrix(0, M, num_iter)
for(t in 1:num_iter){
  h_t <- h_update(Z,C_t, sgm_sq_t_1, eta_h, tau_h, G,K)
  sgm_sq_t_1 <- sgm_sq_update(Z,h_t,C_t,alpha_1,beta_1,G,M)
  C_t <- C_update(Z, h_t, C_t, sgm_sq_t_1, gam_1, M, K)
  H_T1[,,t] <- h_t
  Sgm_sq_t1[,t] <- sgm_sq_t_1
  C_T1[,t] <- C_t
}

h_sim1 <- H_T1[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
sgm_sq_sim1 <- rowMeans(Sgm_sq_t1[,(4*num_iter/5):num_iter])
C_sim1 <- rowMeans(C_T1[,(4*num_iter/5):num_iter]) %>% floor()
BIC1 <- BIC_k(Z,h_sim1, C_sim1,sgm_sq_sim1)
BIC1
table(C_sim1)

#---- K=6 ----
K = 6
gam_1 <- rep(2, K)
C_t <- sample(1:K, M, replace = TRUE)
h_t <- matrix(NA, G, K)
for(k in 1:K){
  h_t[, k] <- rowMeans(Z[ , C_t == k])
}
sgm_sq_t_1 <- apply(Z - h_t[, C_t], 1, var)

H_T2 <- array(0, dim = c(G, K,num_iter))
Sgm_sq_t2 <- matrix(0, G,num_iter)
C_T2 <- matrix(0, M, num_iter)

for(t in 1:num_iter){
  h_t <- h_update(Z,C_t, sgm_sq_t_1, eta_h, tau_h, G,K)
  sgm_sq_t_1 <- sgm_sq_update(Z,h_t,C_t,alpha_1,beta_1,G,M)
  C_t <- C_update(Z, h_t, C_t, sgm_sq_t_1, gam_1, M, K)
  H_T2[,,t] <- h_t
  Sgm_sq_t2[,t] <- sgm_sq_t_1
  C_T2[,t] <- C_t
}

h_sim2 <- H_T2[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
sgm_sq_sim2 <- rowMeans(Sgm_sq_t2[,(4*num_iter/5):num_iter])
C_sim2 <- rowMeans(C_T2[,(4*num_iter/5):num_iter]) %>% floor()
BIC2 <- BIC_k(Z,h_sim2, C_sim2,sgm_sq_sim2)
BIC2
table(C_sim2)
#---- K=7 ----
K = 7
gam_1 <- rep(2, K)
C_t <- sample(1:K, M, replace = TRUE)
h_t <- matrix(NA, G, K)
for(k in 1:K){
  h_t[, k] <- rowMeans(Z[ , C_t == k])
}
sgm_sq_t_1 <- apply(Z - h_t[, C_t], 1, var)

H_T3 <- array(0, dim = c(G, K,num_iter))
Sgm_sq_t3 <- matrix(0, G,num_iter)
C_T3 <- matrix(0, M, num_iter)

for(t in 1:num_iter){
  h_t <- h_update(Z,C_t, sgm_sq_t_1, eta_h, tau_h, G,K)
  sgm_sq_t_1 <- sgm_sq_update(Z,h_t,C_t,alpha_1,beta_1,G,M)
  C_t <- C_update(Z, h_t, C_t, sgm_sq_t_1, gam_1, M, K)
  H_T3[,,t] <- h_t
  Sgm_sq_t3[,t] <- sgm_sq_t_1
  C_T3[,t] <- C_t
}

h_sim3 <- H_T3[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
sgm_sq_sim3 <- rowMeans(Sgm_sq_t3[,(4*num_iter/5):num_iter])
C_sim3 <- rowMeans(C_T3[,(4*num_iter/5):num_iter]) %>% floor()
BIC3 <- BIC_k(Z,h_sim3, C_sim3,sgm_sq_sim3)
BIC3
table(C_sim3)
#---- K=8 ----
K = 8
gam_1 <- rep(2, K)
C_t <- sample(1:K, M, replace = TRUE)
h_t <- matrix(NA, G, K)
for(k in 1:K){
  h_t[, k] <- rowMeans(Z[ , C_t == k])
}
sgm_sq_t_1 <- apply(Z - h_t[, C_t], 1, var)

H_T4 <- array(0, dim = c(G, K,num_iter))
Sgm_sq_t4 <- matrix(0, G,num_iter)
C_T4 <- matrix(0, M, num_iter)

for(t in 1:num_iter){
  h_t <- h_update(Z,C_t, sgm_sq_t_1, eta_h, tau_h, G,K)
  sgm_sq_t_1 <- sgm_sq_update(Z,h_t,C_t,alpha_1,beta_1,G,M)
  C_t <- C_update(Z, h_t, C_t, sgm_sq_t_1, gam_1, M, K)
  H_T4[,,t] <- h_t
  Sgm_sq_t4[,t] <- sgm_sq_t_1
  C_T4[,t] <- C_t
}

h_sim4 <- H_T4[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
sgm_sq_sim4 <- rowMeans(Sgm_sq_t4[,(4*num_iter/5):num_iter])
C_sim4 <- rowMeans(C_T4[,(4*num_iter/5):num_iter]) %>% floor()
BIC4 <- BIC_k(Z,h_sim4, C_sim4,sgm_sq_sim4)
BIC4
table(C_sim4)

#---- K = 4----
K = 4
gam_1 <- rep(2, K)
C_t <- sample(1:K, M, replace = TRUE)
h_t <- matrix(NA, G, K)
for(k in 1:K){
  h_t[, k] <- rowMeans(Z[ , C_t == k])
}
sgm_sq_t_1 <- apply(Z - h_t[, C_t], 1, var)

H_T5 <- array(0, dim = c(G, K,num_iter))
Sgm_sq_t5 <- matrix(0, G,num_iter)
C_T5 <- matrix(0, M, num_iter)

for(t in 1:num_iter){
  h_t <- h_update(Z,C_t, sgm_sq_t_1, eta_h, tau_h, G,K)
  sgm_sq_t_1 <- sgm_sq_update(Z,h_t,C_t,alpha_1,beta_1,G,M)
  C_t <- C_update(Z, h_t, C_t, sgm_sq_t_1, gam_1, M, K)
  H_T5[,,t] <- h_t
  Sgm_sq_t5[,t] <- sgm_sq_t_1
  C_T5[,t] <- C_t
}

h_sim5 <- H_T5[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
sgm_sq_sim5 <- rowMeans(Sgm_sq_t5[,(4*num_iter/5):num_iter])
C_sim5 <- rowMeans(C_T5[,(4*num_iter/5):num_iter]) %>% floor()
BIC5 <- BIC_k(Z,h_sim5, C_sim5,sgm_sq_sim5)
BIC5
table(C_sim5)

#------choose best K----
BIC_s <- as.data.frame(cbind(c(4:8),c(BIC5, BIC1,BIC2,BIC3,BIC4)))
ggplot(BIC_s,aes(x = BIC_s[,1], y = BIC_s[,2])) + 
  geom_line(color = 'blue', alpha = 0.6) + geom_point(color = 'blue', alpha = 0.6) + 
  labs(x = 'cluster numbers', y = 'BIC') +
  geom_vline(xintercept = 7,linetype = 'dotted') 
BIC_s
K_opt = BIC_s[which.min(BIC_s[,2]),1]
bx <- as.numeric(H_T3[20, 3 ,][2000:5000])
by <- cbind(c(2000:5000),bx)
tmp_0 <- as.data.frame(cbind(c(2000:5000),H_T3[20, 3 ,][2000:5000]))
names(tmp_0) <- c('iter', 'h') 
ggplot(tmp_0,aes(x = iter, y = h)) + 
  geom_line(color = 'blue') + 
  labs(x = 'iteration', y = 'h') +
  #scale_y_continuous(limits = c(0, 2)) +
  geom_hline(yintercept = h_sim3[20,3], color = 'red')

tmp_2 <- as.data.frame(cbind(c(2000:5000), Sgm_sq_t3[100,][2000:5000]))
ggplot(tmp_2,aes(x = tmp_2[,1], y = tmp_2[,2])) + 
  geom_line(color = 'blue') + 
  labs(x = 'iteration', y = expression(sigma^2)) +
  scale_y_continuous(limits = c(0,0.2)) +
  geom_hline(yintercept = sgm_sq_sim3[100], color = 'red')
table(C_sim3)

rna_umap <- umap(t(Z))
rna_7 <- rna_umap$layout
rna_7 <- as.data.frame(rna_7)
rna_7$CT <- C_sim3
ggplot(rna_7,aes(x = rna_7[,1], y = rna_7[,2], color = letters[CT])) + 
  geom_point() + xlab('umap1') + ylab('umap2')


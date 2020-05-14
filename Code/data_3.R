rm(list = ls(all = TRUE))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
#----data read ----
ST1 <- read.csv('ST_complete.csv',header = TRUE)
rownames(ST1) <- ST1[,1]
ST1 <- ST1[,-1]
X <- as.matrix(ST1)
L = 31
W = 33
G = dim(ST1)[1]
N = L * W

ind_com <- matrix(0,L*W,2)
for(j in 1:W){
  for(i in 1:L){
    ind_com[i+L*(j-1),] = c(i,j)
  }
}
ind_com <- as.data.frame(ind_com)
names(ind_com) <- c('row_ind','col_ind')
null_na <- read.csv('null_index.csv')[,-1]
index_null <- c()
for(i in 1:16){
  index_null[i] <- L * (null_na[i,2] -1) + null_na[i,1]
}

#----hyper-paramters----
eta_theta <- 0
tau_theta <- 10
eta_mu <- 0
tau_mu <- 10
par_alpha <- 5
par_beta <- 0.1
#tuning parameter in the proposal distribution
tau_0 <- 0.1

#----choose different region number S
set.seed(1996)
#----S = 5----
S = 5
km1 <- kmeans(t(X), S)
R_t_vec <- km1$cluster
# R_t_vec <- sample(2:S, N, replace = TRUE)
# R_t_vec[index_null] <- 1
R_t <- matrix(R_t_vec, L, W)

mu_t <- matrix(NA, G, S)
for(s in 1:S){
  mu_t[ ,s] <- rowMeans(X[ , R_t_vec == s])
}

sgm_sq_t <- as.numeric(apply(X - mu_t[ ,R_t_vec], 1, var))
theta_t <- matrix(rnorm(S*S, eta_theta, 0.01),S,S)
theta_t <- (theta_t + t(theta_t))/2
diag(theta_t) <- 0
#iteration
num_iter <- 5000
Mu <- array(0, dim = c(G,S, num_iter))
Sgm_sq <- matrix(0,G,num_iter)
R_T <- array(0, dim = c(L,W,num_iter))
Theta <- array(0, dim = c(S,S,num_iter))
ptm <- proc.time()
  for(t in 1:num_iter){
    mu_t = mu_update(X, sgm_sq_t, R_t, tau_mu, eta_mu,S, G)
    sgm_sq_t = sgm_sq_star_update(X, R_t, mu_t, S, G, N, par_alpha, par_beta)
    R_t = R_update(X, R_t, mu_t, theta_t, sgm_sq_t,S,G,L,W)
    theta_t = theta_update(X, R_t, mu_t, theta_t, sgm_sq_t,
                           S, G, L, W,tau_0, eta_theta,tau_theta)
    Mu[,,t] <- mu_t
    Sgm_sq[,t] <- sgm_sq_t
    R_T[,,t] <- R_t
    Theta[,,t] <- theta_t
  }
print(proc.time()-ptm)

mu_sim <- Mu[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
R_sim <- R_T[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean) %>% floor()
sgm_sq_sim <- rowMeans(Sgm_sq[,(4*num_iter/5):num_iter])
theta_sim <- Theta[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
BIC_5 <- BIC_k(X, mu_sim, c(R_sim), sgm_sq_sim)
BIC_5
tmp1 <- ind_com
tmp1$Re <- as.numeric(c(R_sim))
ggplot(tmp1, aes(col_ind, row_ind,color = letters[Re])) + geom_point(alpha = 0.6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
table(R_sim)

#---- S = 6----
S = 6
km2 <- kmeans(t(X), S)
R_t_vec <- km2$cluster
# R_t_vec <- sample(2:S, N, replace = TRUE)
# R_t_vec[index_null] <- 1
R_t <- matrix(R_t_vec, L, W)
mu_t <- matrix(NA, G, S)
for(s in 1:S){
  mu_t[ ,s] <- rowMeans(X[ , R_t_vec == s])
}
sgm_sq_t <- as.numeric(apply(X - mu_t[ ,R_t_vec], 1, var))
theta_t <- matrix(rnorm(S*S, eta_theta, 0.01),S,S)
theta_t <- (theta_t + t(theta_t))/2
diag(theta_t) <- 0

num_iter <- 5000
Mu2 <- array(0, dim = c(G,S, num_iter))
Sgm_sq2 <- matrix(0,G,num_iter)
R_T2 <- array(0, dim = c(L,W,num_iter))
Theta2 <- array(0, dim = c(S,S,num_iter))
#iteration
for(t in 1:num_iter){
  mu_t = mu_update(X, sgm_sq_t, R_t, tau_mu, eta_mu,S, G)
  sgm_sq_t = sgm_sq_star_update(X, R_t, mu_t, S, G, N, par_alpha, par_beta)
  R_t = R_update(X, R_t, mu_t, theta_t, sgm_sq_t,S,G,L,W)
  theta_t = theta_update(X, R_t, mu_t, theta_t, sgm_sq_t,
                         S, G, L, W,tau_0, eta_theta,tau_theta)
  Mu2[,,t] <- mu_t
  Sgm_sq2[,t] <- sgm_sq_t
  R_T2[,,t] <- R_t
  Theta2[,,t] <- theta_t
}

mu_sim2 <- Mu2[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
R_sim2 <- R_T2[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean) %>% floor()
sgm_sq_sim2 <- rowMeans(Sgm_sq2[,(4*num_iter/5):num_iter])
theta_sim2 <- Theta2[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
BIC_6 <- BIC_k(X, mu_sim2, c(R_sim2), sgm_sq_sim2)
BIC_6

tmp2 <- ind_com
tmp2$Re <- as.numeric(c(R_sim2))
ggplot(tmp2, aes(col_ind, row_ind,color = letters[Re])) + geom_point(alpha = 0.6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
plot(Mu2[10,5,],type = 'l')
plot(Theta2[1,3,], type = 'l')
plot(Sgm_sq2[100,], type = 'l')

#----S = 7----
S = 7
km3 <- kmeans(t(X), S)
R_t_vec <- km3$cluster
# R_t_vec <- sample(2:S, N, replace = TRUE)
# R_t_vec[index_null] <- 1
R_t <- matrix(R_t_vec, L, W)
mu_t <- matrix(NA, G, S)
for(s in 1:S){
  mu_t[ ,s] <- rowMeans(X[ , R_t_vec == s])
}
sgm_sq_t <- as.numeric(apply(X - mu_t[ ,R_t_vec], 1, var))
theta_t <- matrix(rnorm(S*S, eta_theta, 0.01),S,S)
theta_t <- (theta_t + t(theta_t))/2
diag(theta_t) <- 0

num_iter <- 5000
Mu3 <- array(0, dim = c(G,S, num_iter))
Sgm_sq3 <- matrix(0,G,num_iter)
R_T3 <- array(0, dim = c(L,W,num_iter))
Theta3 <- array(0, dim = c(S,S,num_iter))

#iteration
ptm <- proc.time()
for(t in 3403:num_iter){
  mu_t = mu_update(X, sgm_sq_t, R_t, tau_mu, eta_mu,S, G)
  sgm_sq_t = sgm_sq_star_update(X, R_t, mu_t, S, G, N, par_alpha, par_beta)
  R_t = R_update(X, R_t, mu_t, theta_t, sgm_sq_t,S,G,L,W)
  theta_t = theta_update(X, R_t, mu_t, theta_t, sgm_sq_t,
                         S, G, L, W,tau_0, eta_theta,tau_theta)
  Mu3[,,t] <- mu_t
  Sgm_sq3[,t] <- sgm_sq_t
  R_T3[,,t] <- R_t
  Theta3[,,t] <- theta_t
}
print(proc.time()-ptm)
mu_sim3 <- Mu3[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
R_sim3 <- R_T3[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean) %>% floor()
sgm_sq_sim3 <- rowMeans(Sgm_sq3[,(4*num_iter/5):num_iter])
theta_sim3 <- Theta3[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
BIC_7 <- BIC_k(X, mu_sim3, c(R_sim3), sgm_sq_sim3)
BIC_7

tmp3 <- ind_com
tmp3$Re <- as.numeric(c(R_sim3))
ggplot(tmp3, aes(col_ind, row_ind,color = letters[Re])) + geom_point(alpha = 0.6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

#---- S=3 ----
S = 3
km4 <- kmeans(t(X), S)
R_t_vec <- km4$cluster
# R_t_vec <- sample(1:S, N, replace = TRUE)
# R_t_vec[index_null] <- 1
R_t <- matrix(R_t_vec, L, W)
mu_t <- matrix(NA, G, S)
for(s in 1:S){
  mu_t[ ,s] <- rowMeans(X[ , R_t_vec == s])
}
sgm_sq_t <- as.numeric(apply(X - mu_t[ ,R_t_vec], 1, var))
theta_t <- matrix(rnorm(S*S, eta_theta, 0.01),S,S)
theta_t <- (theta_t + t(theta_t))/2
diag(theta_t) <- 0

num_iter <- 5000
Mu4 <- array(0, dim = c(G,S, num_iter))
Sgm_sq4 <- matrix(0,G,num_iter)
R_T4 <- array(0, dim = c(L,W,num_iter))
Theta4 <- array(0, dim = c(S,S,num_iter))

#iteration
ptm <- proc.time()
for(t in 1:num_iter){
  mu_t = mu_update(X, sgm_sq_t, R_t, tau_mu, eta_mu,S, G)
  sgm_sq_t = sgm_sq_star_update(X, R_t, mu_t, S, G, N, par_alpha, par_beta)
  R_t = R_update(X, R_t, mu_t, theta_t, sgm_sq_t,S,G,L,W)
  theta_t = theta_update(X, R_t, mu_t, theta_t, sgm_sq_t,
                         S, G, L, W,tau_0, eta_theta,tau_theta)
  Mu4[,,t] <- mu_t
  Sgm_sq4[,t] <- sgm_sq_t
  R_T4[,,t] <- R_t
  Theta4[,,t] <- theta_t
}
print(proc.time()-ptm)
mu_sim4 <- Mu4[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
R_sim4 <- R_T4[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean) %>% floor()
sgm_sq_sim4 <- rowMeans(Sgm_sq4[,(4*num_iter/5):num_iter])
theta_sim4 <- Theta4[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
BIC_3 <- BIC_k(X, mu_sim4, c(R_sim4), sgm_sq_sim4)
BIC_3

tmp4 <- ind_com
tmp4$Re <- as.numeric(c(R_sim4))
ggplot(tmp4, aes(col_ind, row_ind,color = letters[Re])) + geom_point(alpha = 0.6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

#----- S = 4----
S = 4
km5 <- kmeans(t(X), S)
R_t_vec <- km5$cluster
# R_t_vec <- sample(2:S, N, replace = TRUE)
# R_t_vec[index_null] <- 1
R_t <- matrix(R_t_vec, L, W)
mu_t <- matrix(NA, G, S)
for(s in 1:S){
  mu_t[ ,s] <- rowMeans(X[ , R_t_vec == s])
}
sgm_sq_t <- as.numeric(apply(X - mu_t[ ,R_t_vec], 1, var))
theta_t <- matrix(rnorm(S*S, eta_theta, 0.01),S,S)
theta_t <- (theta_t + t(theta_t))/2
diag(theta_t) <- 0

num_iter <- 5000
Mu5 <- array(0, dim = c(G,S, num_iter))
Sgm_sq5 <- matrix(0,G,num_iter)
R_T5 <- array(0, dim = c(L,W,num_iter))
Theta5 <- array(0, dim = c(S,S,num_iter))

#iteration
ptm <- proc.time()
for(t in 1:num_iter){
  mu_t = mu_update(X, sgm_sq_t, R_t, tau_mu, eta_mu,S, G)
  sgm_sq_t = sgm_sq_star_update(X, R_t, mu_t, S, G, N, par_alpha, par_beta)
  R_t = R_update(X, R_t, mu_t, theta_t, sgm_sq_t,S,G,L,W)
  theta_t = theta_update(X, R_t, mu_t, theta_t, sgm_sq_t,
                         S, G, L, W,tau_0, eta_theta,tau_theta)
  Mu5[,,t] <- mu_t
  Sgm_sq5[,t] <- sgm_sq_t
  R_T5[,,t] <- R_t
  Theta5[,,t] <- theta_t
}
print(proc.time()-ptm)
mu_sim5 <- Mu5[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
R_sim5 <- R_T5[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean) %>% floor()
sgm_sq_sim5 <- rowMeans(Sgm_sq5[,(4*num_iter/5):num_iter])
theta_sim5 <- Theta5[,,(4*num_iter/5):num_iter] %>% apply(c(1,2),mean)
BIC_4 <- BIC_k(X, mu_sim5, c(R_sim5), sgm_sq_sim5)
BIC_4

tmp5 <- ind_com
tmp5$Re <- as.numeric(c(R_sim5))
ggplot(tmp5, aes(col_ind, row_ind,color = letters[Re])) + geom_point(alpha = 0.6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


#------intergrate--------
#BIC plot
BIC_total <- cbind(c(3,4,5,6,7),c(BIC_3, BIC_4, BIC_5, BIC_6, BIC_7))
BIC_total <- as.data.frame(BIC_total)
ggplot(BIC_total,aes(x = BIC_total[,1], y = BIC_total[,2])) + 
  geom_line(color = 'blue', alpha = 0.6) + geom_point(color = 'blue', alpha = 0.6) + 
  labs(x = 'Region numbers', y = 'BIC') +
  geom_vline(xintercept = 6,linetype = 'dotted') 

#K = 6
#parameters iteration
tmp1 <- as.data.frame(cbind(c(2000:4000), Mu2[10, 3 ,][2000:4000]))
ggplot(tmp1,aes(x = tmp1[,1], y = tmp1[,2])) +
  geom_line(color = 'blue') + 
  labs(x = 'iteration', y = expression(mu)) + 
  scale_y_continuous(limits = c(1.25,2.0)) +
  geom_hline(yintercept = mu_sim2[10,3], color = 'red')

tmp2 <- as.data.frame(cbind(c(2000:4000), Sgm_sq2[100,][2000:4000]))
ggplot(tmp2,aes(x = tmp2[,1], y = tmp2[,2])) + 
  geom_line(color = 'blue') + 
  labs(x = 'iteration', y = expression(sigma^2)) + 
  scale_y_continuous(limits = c(0.02,0.18)) +
  geom_hline(yintercept = sgm_sq_sim2[100], color = 'red')

tmp3 <- as.data.frame(cbind(c(2000:5000), Theta2[2, 3 ,][2000:5000]))
ggplot(tmp3,aes(x = tmp3[,1], y = tmp3[,2])) + 
  geom_line(color = 'blue') + 
  labs(x = 'iteration', y = expression(theta)) +
  scale_y_continuous(limits = c(-4,8)) +
  geom_hline(yintercept = theta_sim2[2,3], color = 'red')


# 
# write.csv(Mu2[10,,],'mu_10.csv')
# write.csv(Sgm_sq2, 'Sgm_sq.csv')
# write.csv(Theta2[2,,],'Theta_2.csv')


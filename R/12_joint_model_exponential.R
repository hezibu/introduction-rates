library(tidyverse)
library(rstan)
library(loo)
library(shinystan)
library(brms)
library(bayesplot)
library(mobsim)

model <- stan_model("Stan/stan_exp_model.stan")

new.sps=function(sps.pool, sps.not.entered, n.enter){
  sps.enter=sample(sps.pool, n.enter, replace = F)  #names of entering sps
  sps.new=sps.not.entered[which(sps.not.entered%in%sps.enter)] #names of NEW sps 
  return(sps.new) 
}

t = 1:150 #time seres
sps.pool = 1:2000 #assumed species pool of invasives
b0 = 1
b1 = 0.015
sim.data.list=list()
all.sps.entered.list=list()

sps.entered=c()
all.sps.entered=list()  

u=exp(b0+b1*t)
sim.data=cbind.data.frame(t,u)
plot(u~t, sim.data)

#number of species entering the system in year t
sim.data$y.t=apply(sim.data, 1, function(x) rpois(1,x["u"]))

sps.not.entered=sps.pool
sim.data$n.new.sps.enter=0
sim.data$n.tot.inv=0

for (i in 1:length(t)){
  # names & number of new species entering in t
  new.sps.enter=new.sps(sps.pool=sps.pool, sps.not.entered = sps.not.entered,
                        n.enter = sim.data$y.t[i])
  sim.data$n.new.sps.enter[i]=length(new.sps.enter)
  
  #update list of unknown species (sps.not.entered)
  if(sim.data$n.new.sps.enter[i]==0){
    sps.not.entered= sps.not.entered
  } else{
    sps.not.entered=sps.not.entered[-c(which(sps.not.entered%in%new.sps.enter))]
  }
  
  sps.entered=c(sps.entered, new.sps.enter)
  all.sps.entered[[i]]=sps.entered
}

sim.data$n.tot.inv=cumsum(sim.data$n.new.sps.enter)
plot(n.tot.inv~t, sim.data)

full_sim.data = sim.data
full_all.sps.entered = all.sps.entered

M=5000 #species pool, natives
nat.sps.pool=1:M

samp.time=sort(sample(t[-1], 60, replace = F))
sim.data=full_sim.data[c(1,samp.time),c("t", "n.tot.inv")]
all.sps.entered=full_all.sps.entered[c(1,samp.time)]

sim.data$n.Inv.discov=0
sim.data$n.Nat.discov=0
all.Inv.discov=c()
all.Nat.discov=c()

t.steps=dim(sim.data)[1]

for(j in 1:t.steps) {
  samp.effort=   rpois(1, 100) #10
  
  #given SAD and sampling-effort, how many species are sampled:
  # length(sad)
  sad <- sim_sad(s_pool = M+length(all.sps.entered[[j]]), 
                 n_sim = samp.effort, 
                 sad_type = "lnorm",
                 sad_coef = list("meanlog" = 5, "sdlog" = 0.5))
  
  # probability that sampled species id invasive
  pr_Inv= length(all.sps.entered[[j]]) / (M+length(all.sps.entered[[j]]))
  
  #number & id of invasive and native sps sampled
  n.sps.samp=rmultinom(1, size=length(sad), prob = c(pr_Inv, 1-pr_Inv)) 
  id.Inv=sample(all.sps.entered[[j]], n.sps.samp[1], replace = F)
  id.Nat=sample(nat.sps.pool, n.sps.samp[2], replace = F)
  
  #new discoveries = sample species not in current discovered_species list
  new.Inv=id.Inv[!id.Inv%in%all.Inv.discov]
  new.Nat=id.Nat[!id.Nat%in%all.Nat.discov]
  
  #number of new discoveries at t
  sim.data$n.Inv.discov[j]=length(new.Inv)   
  sim.data$n.Nat.discov[j]=length(new.Nat)
  
  # update discovered_species list
  all.Inv.discov=c(all.Inv.discov, new.Inv)
  all.Nat.discov=c(all.Nat.discov, new.Nat)
}

sim.data$tot.discov.Inv=cumsum(sim.data$n.Inv.discov)
sim.data$tot.discov.Nat=cumsum(sim.data$n.Nat.discov)

#shift up by 1 time-step to make total number at start of t+1
sim.data$tot.discov.Inv=c(0, sim.data$tot.discov.Inv[1:(t.steps-1)])
sim.data$tot.discov.Nat=c(0,sim.data$tot.discov.Nat[1:(t.steps-1)])
sim.data=sim.data[c("t","n.tot.inv","tot.discov.Inv","tot.discov.Nat", "n.Inv.discov","n.Nat.discov" )]
sim.data$dsps=sim.data$n.Inv.discov+sim.data$n.Nat.discov

g = ggplot(data = sim.data)  + 
  geom_point(aes(x = t, y = n.tot.inv), col='black') +
  xlim(0, 150) + ylim(0,M) + ylab("No. of Species") +
  geom_point(aes(x = t, y = tot.discov.Inv), color='blue') +
  geom_point(aes(x = t, y = tot.discov.Nat), color='red')
g


d=list(
  M = M,
  N=as.integer(dim(sim.data)[1]),
  dI=as.integer(sim.data$n.Inv.discov), #discovery.data$n.Inv_t
  d_Nativ=as.numeric(sim.data$n.Nat.discov),  #discovery.data$n.Nativ_t
  dsps=as.integer(sim.data$dsps), # 
  t=as.numeric(sim.data$t),  # ----
  n_Inv=as.numeric(sim.data$tot.discov.Inv), #$n.Inv
  n_Nativ=as.numeric(sim.data$tot.discov.Nat) #n.Nativ
)

files <- dir("~/invasionrate/Invasives/", pattern = ".stan", full.names = T)

# run model and store output

fit1.test <- stan(file = "Stan/stan_exp_model.stan",
                  data = d,
                  chains = 4,      # number of Markov chains
                  cores  = 2,       # number of cores   
                  warmup = 10000,   # number of warmup iterations per chain
                  iter = 20000,     # total number of iterations per chain
                  refresh = 0,     # show progress every 'refresh' iterations
                  thin = 20,
                  control = list(adapt_delta = 0.99)
)

fitted.data = rstan::extract(fit1.test, permuted=TRUE)
sim.data$av.Itot=apply(fitted.data$Itot, 2, median)
sim.data$av.dI=apply(fitted.data$dI_rep, 2, median)
sim.data$av.discov_t=apply(fitted.data$discov_t, 2, median)

sim.data$quant.01.Itot=apply(fitted.data$Itot, 2, 
                             function (x) quantile(x, probs=c(0.01)))
sim.data$quant.99.Itot=apply(fitted.data$Itot, 2, 
                             function (x) quantile(x, probs=c(0.99)))                   
sim.data$quant.01.dI=apply(fitted.data$dI_rep, 2, 
                           function (x) quantile(x, probs=c(0.01)))
sim.data$quant.99.dI=apply(fitted.data$dI_rep, 2, 
                           function (x) quantile(x, probs=c(0.99)))   
sim.data$quant.01.discov_t=apply(fitted.data$discov_t, 2, 
                                 function (x) quantile(x, probs=c(0.01)))
sim.data$quant.99.discov_t=apply(fitted.data$discov_t, 2, 
                                 function (x) quantile(x, probs=c(0.99))) 

f.Itot = ggplot(data = sim.data)  + 
  geom_point(aes(x = t, y = n.tot.inv), col='black') +
  geom_line(aes(x = t, y = av.Itot), col='red') +
  xlim(1,150) + ylim(0,length(sps.pool)+20) + ylab("Introduced Species") 

f.Itot = f.Itot + geom_ribbon(data = sim.data, 
                              aes(x = t, ymin=quant.01.Itot, ymax=quant.99.Itot), 
                              alpha=0.2, fill='red') 
f.Itot



f.discov_t = ggplot(data = sim.data)  + 
  geom_point(aes(x = t, y = tot.discov.Inv), col='black') +
  geom_line(aes(x = t, y = av.discov_t), col='red') +
  xlim(1,150) + ylim(0,length(sps.pool)+20) + ylab("Discovered Inv Species") 

f.discov_t = f.discov_t + geom_ribbon(data = sim.data, 
                                      aes(x = t, ymin=quant.01.discov_t, ymax=quant.99.discov_t), 
                                      alpha=0.2, fill='red') 
f.discov_t



f.dI = ggplot(data = sim.data)  + 
  geom_point(aes(x = t, y = n.Inv.discov), col='black') +
  geom_line(aes(x = t, y = av.dI), col='red') +
  ylab("New Inv Species") 

f.dI = f.dI + geom_ribbon(data = sim.data, 
                          aes(x = t, ymin=quant.01.dI, ymax=quant.99.dI), 
                          alpha=0.2, fill='red') 
f.dI


hist(fitted.data$b1)
abline(v = 0.015)


optim.fit=optimizing(model, 
                     data = d, hessian = TRUE)
optim.fit$par[1:2]
opt.Itot=optim.fit$par[str_detect(c(names(optim.fit$par)),"Itot")]
plot(opt.Itot~d$t)

y <- d$dI
yrep1 <- extract(fit1.test)[["dI_rep"]]
samp100 <- sample(nrow(yrep1), 100)
ppc_dens_overlay(y, yrep1[samp100, ])  
ppc_stat(y, yrep1, stat = 'mean')

fit1=brm(bf(#n.tot.inv ~ d+ (a-d) *(1 / (1 + exp(-(log(t)-c)*b))),
  n.tot.inv ~ exp(b0 + b1*t) ,
  b0 + b1 ~ 1 , 
  nl = TRUE),
  data = sim.data, #sim.data.list[[L]], 
  family = gaussian(),
  prior = c(prior(normal(1, 1), nlpar = "b0"), #b0
            prior(normal(0.1, 0.01), nlpar = "b1")),  #b1 
  chains = 3,
  iter = 10000,
  thin = 2,
  save_model= 'brms_5PL',
  control = list(adapt_delta = 0.9))

summary(fit1)
plot(fit1)
pp_check(fit1)

plot(conditional_effects(fit1), points = TRUE)
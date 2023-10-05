### Code to estimate cumulative incidence from cross-sectional PCR survey during
# BQ.1 wave in Salvador, Brazil
### Matt Hitchings, based on method of James Hay et al

library(devtools)
library(lazymcmc)
library(virosolver)
library(extraDistr)
require(dplyr)
require(ggplot2)
library(virosolver)
library(tidyverse)
library(patchwork)
library(lazymcmc)
library(foreach)
library(doParallel)
require(lubridate)

# Set parameters to define CT curve over time and probability of PCR positivity
rho = 0.06
prob_detect = 0.2
level_switch=38
viral_peak = 25
desired_mode=3.8
t_switch=7

file_suffix = paste0('omicronct_pd',
                     gsub('\\.','',round(prob_detect,digits = 2)),
                     '_viralpeak',viral_peak,
                     '_tswitch',t_switch,
                     '_levelswitch',level_switch,
                     '_rho',gsub('\\.','',rho))

# Number of days before the first PCR test that we will infer incidence
td = 28

pcr_data = read.csv('Data/PdL_CT_0110.csv')

colnames(pcr_data) = c("ID","HH_ID","PCRRes_D0","ct","Date_D0","Symp_D0","SympOnset_D0","PCRRes_D7","CT_D7","Date_D7","Symp_D7","SympOnset_D7")
pcr_data = pcr_data %>% mutate(Date_D0=as_date(Date_D0,format='%m/%d/%Y'),
                               SympOnset_D0=as_date(SympOnset_D0,format='%m/%d/%Y'))
pcr_data = pcr_data %>% mutate(w=
                                (ceiling(as.numeric((Date_D0-as.Date("2022-11-15"))/7))),
                               #epiweek(Date_D0),
                               epiweek = epiweek(Date_D0),
                               woy = week(Date_D0),
                               d = as.numeric(Date_D0-min(Date_D0)+td),
                               t_w = 7*(w-1+td/7),
                               ct = case_when(PCRRes_D0=="Positive" ~ ct,
                                              PCRRes_D0=="Negative" ~ 40),
                               symponset_to_testdate = as.numeric(Date_D0-SympOnset_D0))

pcr_data = pcr_data %>% mutate(t=d) 

ggplot(pcr_data %>% filter(ct<40)) + 
  geom_violin(aes(x=factor(w),group=factor(w),y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(aes(x=factor(w),y=ct),size=0.1,width=0.4,height=0) + 
  scale_y_continuous(trans="reverse") +
  theme_bw() +
  ylab("Ct value") +
  xlab("Epidemiological week of sample") +
  ggtitle("Observed Ct values < 40 over time")+
    theme(title = element_text(size=6),
          axis.text = element_text(size=6))

## Set up parameters for virosolver
times <- 0:max(pcr_data$t)
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
par_tab <- example_gp_partab
par_tab <- bind_rows(par_tab[par_tab$names != "prob",], par_tab[par_tab$names == "prob",][1:length(times),])
pars <- par_tab$values
names(pars) <- par_tab$names

par_tab$values[par_tab$names=="desired_mode"]=desired_mode
par_tab$values[par_tab$names=="viral_peak"]=viral_peak
par_tab$values[par_tab$names=="t_switch"]=t_switch
par_tab$values[par_tab$names=="level_switch"]=level_switch
par_tab$values[par_tab$names=="prob_detect"]=prob_detect
par_tab$values[par_tab$names=="rho"]=rho

# CT model
pars <- par_tab$values
names(pars) <- par_tab$names

## Solve the Ct model over a range of times since infection (referred to as "ages")
test_ages <- seq(1,50,by=1)

## This gives the modal Ct value
cts <- viral_load_func(pars, test_ages)

p_ct_model <- ggplot(data.frame(ct=c(40,cts),t=c(0,test_ages))) + 
  geom_line(aes(x=t,y=ct)) + 
  scale_y_continuous(trans="reverse",
                     limits=c(40,10)) +
  theme_bw() +
  ylab("Modal Ct value") +
  xlab("Days since infection")
p_ct_model
## Note that this model does not solve for t=0, 
## as it is always assumed that no one is detectable 0 days post infection
prop_detect <- prop_detectable(test_ages,pars, cts)
p_ct_model_detectable <- ggplot(data.frame(p=c(0,prop_detect),t=c(0,test_ages))) + 
  geom_line(aes(x=t,y=p)) + 
  theme_bw() +
  ylab("Proportion of infections\n still detectable") +
  xlab("Days since infection")
p_ct_model/p_ct_model_detectable

### Likelihood and priors for multiple cross-sections
## MCMC chain options
mcmc_pars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
               "thin"=2500,"adaptive_period"=200000,"save_block"=100)

# Set pointer to the Gaussian Process model as the incidence function
incidence_function <- gaussian_process_model
## Read in the GP model parameter control table
data(example_gp_partab)

## Pull out the current values for each parameter, and set these as the prior means
means <- par_tab$values
names(means) <- par_tab$names

## Set standard deviations of prior distribution
sds_gp <- c("obs_sd"=0.5,"viral_peak"=2,
            "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
            "prob_detect"=0.03,
            "incubation"=0.25, "infectious"=0.5)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
## Prior for GP version
prior_func_gp <- function(pars, ...){
  par_names <- names(pars)
  
  ## Viral kinetics parameters
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_gp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_gp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_gp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_gp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_gp["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_gp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ### VERY IMPORTANT
  ## Gaussian process prior, un-centered version
  k <- pars[which(par_names=="prob")]
  ## Leave this - correct for uncentered version as per Chapter 14 Statistical Rethinking
  prob_priors <- sum(dnorm(k, 0, 1, log=TRUE))
  #########
  
  nu_prior <- dexp(pars["nu"], 1/means[which(names(means) == "nu")],log=TRUE)
  rho_prior <- dexp(pars["rho"], 1/means[which(names(means) == "rho")],log=TRUE)
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + prob_priors +
    nu_prior + rho_prior
}

## Check that posterior function solves and returns a finite likelihood
posterior_function <- create_posterior_func(parTab=par_tab, 
                                            data=pcr_data, 
                                            PRIOR_FUNC=prior_func_gp, 
                                            INCIDENCE_FUNC=incidence_function,
                                            t_dist=t_dist)
posterior_function(par_tab$values)

####################################################################################
####################################################################################
####################################################################################
############# IF THE MODEL IS ALREADY RUN, DO NOT RUN THIS CODE ####################
####################################################################################
####################################################################################
####################################################################################

dir.create(paste0("mcmc_chains/multiple_cross_section",file_suffix),recursive=TRUE)

cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)
#########################
#
## RUN THE MCMC FRAMEWORK
## Run 3 MCMC chains. Note that it is possible to parallelize this loop with foreach and doPar
## Note the `use_pos` argument needs to be set here too
nchains <- 3
res <- foreach(chain_no=1:nchains,.packages = c("virosolver","lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
  ## Get random starting values
  start_tab <- generate_viable_start_pars(par_tab,pcr_data,
                                          create_posterior_func,
                                          incidence_function,
                                          prior_func_gp)
  
  output <- run_MCMC(parTab=start_tab,
                     data=pcr_data %>% select(t,ct),
                     INCIDENCE_FUNC=incidence_function,
                     PRIOR_FUNC=prior_func_gp,
                     mcmcPars=mcmc_pars,
                     filename=paste0("mcmc_chains/multiple_cross_section",file_suffix,"/gp_",chain_no),
                     CREATE_POSTERIOR_FUNC=create_posterior_func,
                     mvrPars=NULL,
                     OPT_TUNING=0.2,
                     use_pos=FALSE,
                     t_dist=t_dist)
}

####################################################################################
####################################################################################
####################################################################################


## Read in the MCMC chains
chains <- load_mcmc_chains(location=paste0("mcmc_chains/multiple_cross_section",file_suffix),
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=TRUE,
                           unfixed=TRUE,
                           multi=FALSE)
## Reshape for plotting
chains_melted <- chains$chain %>% as_tibble %>% group_by(chain) %>% mutate(sampno=1:n()) %>% pivot_longer(-c(sampno,chain))
## Look at trace plots
p_trace_gp <- chains_melted %>%
  filter(!(name %in% paste0("prob.",1:max(times)))) %>%
  ggplot() + 
  geom_line(aes(x=sampno,y=value,col=as.factor(chain))) + 
  facet_wrap(~name,scales="free_y") + 
  scale_color_viridis_d(name="Chain") + 
  theme_bw() +
  xlab("Iteration") +
  ylab("Value")
p_trace_gp
ggsave(paste0('./mcmc_plots/trace',file_suffix,'.png'),p_trace_gp,height=4,width=4,units='in',device='png')

## Load in MCMC chains again, but this time read in the fixed parameters too 
## to ensure that the posterior draws are compatible with the model functions
chains <- load_mcmc_chains(location=paste0("mcmc_chains/multiple_cross_section",file_suffix),
                           parTab=par_tab,
                           burnin=mcmc_pars["adaptive_period"],
                           chainNo=FALSE,
                           unfixed=FALSE,
                           multi=FALSE)
## Do some reshaping to allow correct subsampling (we need each sampno to correspond to one unique posterior draw)
chain_comb <- as.data.frame(chains$chain)
chain_comb$sampno <- 1:nrow(chain_comb)

## Load in true incidence curve to compare to our prediction
predictions <- plot_prob_infection(chain_comb,nsamps=603, INCIDENCE_FUNC=incidence_function,
                                   solve_times=0:max(pcr_data$t)+1,obs_dat=pcr_data,
                                   smooth=TRUE) ## Smooth the trajectories a bit

# x-axis scale
p_incidence_prediction <- predictions$plot + 
  scale_x_continuous(name='Epidemiological week',
                     limits=c(0,max(pcr_data$t)+2),
                     breaks = seq(min(pcr_data$t-23),max(pcr_data$t),7),
                     labels = (43:51),
                     )
p_incidence_prediction

##### DAY OF PEAK INCIDENCE
peak_incidence_day = sapply(1:max(predictions$predictions$sampno),function(x) {
    ps=predictions$predictions %>% filter(sampno==x & t>14 & t<56)
    return(ps$t[ps$prob_infection==max(ps$prob_infection)])}
)

peak_incidence_summ = quantile(peak_incidence_day,c(0.025,0.5,0.975)) - td + min(pcr_data$Date_D0)

cat('Predicted peak in incidence ',paste0(peak_incidence_summ[2], ' (95% CrI ',peak_incidence_summ[1],' to ',peak_incidence_summ[3],')'))
  
p = predictions$predictions %>% filter(sampno==1)
op = 1 - exp(sum(log(1-p$prob_infection)))
chain_comb$myoverall_prob = sapply(1:nrow(chain_comb),function(x) {
  p = predictions$predictions %>% filter(sampno==x)
  #return(1 - exp(sum(log(1-p$prob_infection))))
  return(sum(p$prob_infection))
})
  
lqs = sapply(1:max(predictions$predictions$t),
             function(x) quantile(predictions$predictions$prob_infection[predictions$predictions$t==x],0.025))
uqs = sapply(1:max(predictions$predictions$t),
             function(x) quantile(predictions$predictions$prob_infection[predictions$predictions$t==x],0.975))
meds = sapply(1:max(predictions$predictions$t),
             function(x) quantile(predictions$predictions$prob_infection[predictions$predictions$t==x],0.5))

d_forplot = data.frame('t' = rep(predictions$map_prediction$t),
                       'y' = meds,
                       'lq' = lqs,
                       'uq' = uqs
)

theme_set(theme_bw())
Panel3G = ggplot() + 
  geom_line(data=d_forplot,aes(x=t,y=y),col='blue') + 
  geom_ribbon(data=d_forplot,aes(x=t,ymin=lq,ymax=uq),alpha=0.2,col='white',fill='blue')+
  scale_y_continuous(name='Probability of infection',
                     #expand=c(0,0),
                     limits=c(0,0.04))+
  scale_x_continuous(name='Epidemiological week',
                     limits=c(0,max(pcr_data$t)+2),
                     #breaks = seq(1,95,1)
                     breaks = seq(min(pcr_data$t-23),max(pcr_data$t),7),
                     labels = (43:51),
  )+
  theme(axis.title.x = element_text(color = "grey20", size = 15, angle = 0, hjust = 0.5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),         
        axis.text.x = element_text(color = "grey20", size = 15, angle = 60, hjust = .5, vjust = .5, face = "plain"),
        axis.ticks.x= element_line(size = 1, linetype="longdash"),
        panel.grid = element_blank()
  )
for (d in unique(pcr_data$t)) {
  Panel3G = Panel3G + geom_vline(xintercept=d,col='red',linetype=2)
}
Panel3G
ggsave(paste0('./mcmc_plots/Panel3G',file_suffix,'.pdf'),Panel3G,height=4,width=6,units='in',device='pdf')

## Use create_posterior_func to return the predicted Ct distribution rather than the posterior probability
model_func_gp <- create_posterior_func(par_tab,pcr_data,NULL,incidence_function,"model")
## Pass model_func to a plotting function to observe predicted Ct distribution against data
p_distribution_fit_gp <- plot_distribution_fits(chain_comb, pcr_data, model_func_gp,100,pos_only=FALSE)
p_distribution_fit_gp$data

best_pars <- lazymcmc::get_best_pars(chain_comb)
samps <- sample(unique(chain_comb$sampno), 100)
all_res <- NULL
for (i in seq_along(samps)) {
  samp <- samps[i]
  tmp_pars <- lazymcmc::get_index_pars(chain_comb, samp)
  prob_infection_tmp <- incidence_function(tmp_pars, 0:max(pcr_data$t)+1)
  all_res[[i]] <- tibble(myoverall_prob = sum(prob_infection_tmp), overall_prob = tmp_pars["overall_prob"],
                         sampno = i)
}
posterior_overallprob <- do.call("bind_rows", all_res)
plot(density(chain_comb$overall_prob),type='l')

# Overall proportion who are PCR positive
posterior_dat <- predicted_distribution_fits(chain_comb, model_func_gp, 
                                             100)
if (grepl('day',file_suffix)) {
  summary_prop_detectable <- posterior_dat %>% mutate(
    w = ceiling((t-27)/7)
  ) %>%
    filter(ct == 
             best_pars["intercept"]) %>% group_by(w,sampno) %>% mutate(density = 1 - density)
  
  summary_prop_detectable$n = rep(unname(table(pcr_data$t)),max(summary_prop_detectable$sampno))
  
  summary_prop_detectable = summary_prop_detectable %>% group_by(w,sampno) %>% 
    mutate(wdensity = sum(density*n)/sum(n),
           meandensity = mean(density)) %>% 
    group_by(w,sampno) %>% summarize(wdensity=mean(wdensity),meandensity=mean(meandensity)) %>% 
    group_by(w) %>%
    summarize(lower = quantile(wdensity, 0.025), 
              median = quantile(wdensity, 0.5), upper = quantile(wdensity, 0.975),
              mlower = quantile(meandensity, 0.025), 
              mmedian = quantile(meandensity, 0.5), mupper = quantile(meandensity, 0.975))
} else {
  summary_prop_detectable <- posterior_dat %>% mutate(
    w = ceiling((t-27)/7)
  ) %>%
    filter(ct == 
             best_pars["intercept"]) %>% group_by(w) %>% mutate(density = 1 - density) %>% 
    summarize(lower = quantile(density, 0.025), 
              median = quantile(density, 0.5), upper = quantile(density, 0.975))
}

cat('Estimated overall PCR positivity: ',(sum(summary_prop_detectable$median * table(pcr_data$w))/nrow(pcr_data)),'\n')
cat('Observed overall PCR positivity: ',mean(pcr_data$PCRRes_D0=="Positive"),'\n')

### Plot observed and estimated PCR positivity by week of sampling
summary_prop_detectable$obs = sapply(summary_prop_detectable$w,function(x) mean(pcr_data$PCRRes_D0[pcr_data$w==x]=="Positive"))
summary_prop_detectable_long = pivot_longer(summary_prop_detectable,cols='lower':'obs')
Panel3D = ggplot() + 
  geom_point(data=summary_prop_detectable_long %>% filter(name %in% c('median','obs')),
               aes(x=as.numeric(w),y=value,col=relevel(factor(name),ref='obs')),shape=16,size=2,
             position = position_dodge(width=0.4))+
  geom_errorbar(data=summary_prop_detectable,aes(x=as.numeric(w)+0.1,ymin=lower,ymax=upper),width=0.1,col='blue')+
  scale_x_continuous(name='Epidemiological week',breaks=1:6,
  )+
  ylab('Prevalence of PCR positivity')+
  scale_color_manual(name='',values=c('grey20','blue'),labels=c('Observed','Fitted'))+
  theme_classic()+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 45, hjust = 0.5, vjust = 0.5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),         
        axis.title.x=element_blank(),
        axis.title.y=element_text(color = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"),
        axis.ticks.x= element_line(size = 1, linetype="longdash"))
Panel3D
ggsave(paste0('./mcmc_plots/Panel3D',file_suffix,'.pdf'),Panel3D,height=3,width=3,units='in',device='pdf')

  

# Overall probability of infection
cat("Cumulative incidence: ",paste0(round(median(chain_comb$overall_prob),2)," (",round(quantile(chain_comb$overall_prob,0.025),2),
       ",",round(quantile(chain_comb$overall_prob,0.975),2),")"),'\n')


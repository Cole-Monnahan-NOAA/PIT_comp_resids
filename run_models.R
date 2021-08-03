#### hacky code from a different experiment, repurposed to
#### explore residuals in composition data
## Cole | 2021-07

### Run and compare models
library(tidyverse)
library(r4ss)
library(snowfall)
library(ggplot2)
theme_set(theme_bw())
packageVersion('r4ss') #  '1.42.0'
source('functions.R')


### ------------------------------------------------------------
### Run bootstrapped data sets back through the model (matches
### perfectly).
## For each boostrap data set, run a retrospective analysis
Nreps <- 498 # messed up starter and put 500 total
reps <- 1:Nreps
Npeels <- 0
peels <- 0:-Npeels
## Setup to run parallel, saving a single core free.
cpus <- parallel::detectCores()-2
sfStop()
sfInit( parallel=TRUE, cpus=cpus)
sfExportAll()
## ## Run the models
## run_model(1:Nreps, model.name='BSAI_FHS', clean.files=FALSE)


### ------------------------------------------------------------
### Do the PIT calculations

## First get Pearson residuals, it's long here so pivot wider
## replist <- SS_output('models/BSAI_FHS', covar=FALSE,
##                      printstats=FALSE, verbose=FALSE)
## saveRDS(replist, 'results/replist.RDS')
replist <- readRDS('results/replist.RDS')
## This is split by sex so need to recombine to be in same format
## as the PIT residuals below.
lc.pearson.female <- replist$lendbase %>% filter(Sex==1) %>%
    select(Yr, FltSvy=Fleet, Gender=Sex,
         Nsamp=Nsamp_in, Pearson, Bin) %>%
  pivot_wider(names_from=Bin, names_prefix='f',
              values_from=Pearson)
lc.pearson.male <- replist$lendbase %>% filter(Sex==2) %>%
    select(Yr, FltSvy=Fleet, Gender=Sex,
         Nsamp=Nsamp_in, Pearson, Bin) %>%
  pivot_wider(names_from=Bin, names_prefix='m',
              values_from=Pearson)
lc.pearson <- cbind(lc.pearson.female, lc.pearson.male[,-(1:4)]) %>% select(-Gender)

## Grab the original data and process it a bit. The Nsamp won't
## necessary match the bootstrap data sets due to variance
## factors, the data may not sum to 1 or Nsamp depending on how
## user inputted them, and it may be in a different order b/c SS
## sorts the bootstrap data (fleet then year) but not the original.
## clean that up here before appending.
## dat0 <- SS_readdat('models/BSAI_FHS/data.ss_new', verbose=FALSE,
##                   section=1)
## datboot <- SS_readdat('models/BSAI_FHS/data.ss_new', verbose=FALSE,
##                   section=3)
## saveRDS(dat0, 'results/dat0.RDS')
## saveRDS(datboot, 'results/datboot.RDS')
dat0 <- readRDS('results/dat0.RDS')       # original data file
datboot <- readRDS('results/datboot.RDS') # single bootstrap data set

## The sample size from the bootstrapped data is correct b/c it's
## been weighted by the variance factors. These factors can be
## different by fleet, and it appears that SS rounds this and has
## a minimum of 1. So it's not going to match exactly, but we
## need the observed counts to reflect the simulated sample size
## as closely as possible. First have to get the order right.
if(any(dat0$lencomp$Gender!=3)) stop("Not setup to work with gender!=3")
lc0 <- dat0$lencomp %>% arrange(abs(FltSvy), Yr)
stopifnot(all.equal(lc0[,1:5], datboot$lencomp[,1:5]))
## these match variance adjustment factors when Nsamp large
## plot(datboot$lencomp$Nsamp, (1/(lc0$Nsamp/datboot$lencomp$Nsamp)))

## Since the rows are lined up can now scale the proportions up
## to the correct sample size so observed and simulated counts
## are equivalent
lc0[,-(1:6)] <- lc0[,-(1:6)]*datboot$lencomp$Nsamp/rowSums(lc0[,-(1:6)])
stopifnot(all.equal(datboot$lencomp$Nsamp, rowSums(lc0[,-(1:6)])))
## lc0 is now processed and comparable to the bootstrap data

## ## Read in each simulated data set and the original. Could read
## ## in from full ss_new but that'd be slower... this still takes a
## ## few minutes at least
## lc <- array(NA, dim=c(85, 48, 499))
## ## ac <- array(NA, dim=c(1085, 42,499))
## for(boot in 1:Nreps){
##   wd <- paste0('runs/BSAI_FHS/boot_', boot)
##   dat <- SS_readdat(file.path(wd, 'data.ss'), verbose=FALSE)
##   ## ac[,,boot+1] <- dat$agecomp[,-(1:9)] %>% as.matrix
##   lc[,,boot+1] <- dat$lencomp[,-(1:6)] %>% as.matrix
## }
## ## Now can set the first element to be the observed data for
## ## easier calculations
## lc[,,1] <- as.matrix(lc0[,-(1:6)])
## ## So first element is the observed data, rest are simulated.
## ## saveRDS(lc, file='results/lencomp_boots.RDS')
lc <- readRDS('results/lencomp_boots.RDS')
Nsamp <- rowSums(lc[,,1])
plot(ecdf(lc[5,5,-1])); abline(v=lc[5,5,1])
plot(ecdf(lc[43,15,-1])); abline(v=lc[43,15,1])

## Now calculate PIT resids
lc.pit <- array(NA, dim=c(dim(lc)[1:2]))
set.seed(1214124)
for(y in 1:nrow(lc.pit)){ # loop year
  for(b in 1:ncol(lc.pit)){ # loop length bins
    ## P(obs data>simulated data)
    lc.pit[y,b] <- qnorm(runif(n=1,
                          min=mean(lc[y,b,1]>lc[y,b,-1]),
                          max=mean(lc[y,b,1]>=lc[y,b,-1]) ))
  }
}
## Cleanup to match Pearson ones
lc.pit <- cbind(lc0[,1:6], lc.pit) %>%
  ## Drop the ghost fleets and negative years
  filter(Yr>0 & FltSvy>0) %>% select(-Seas, -Gender, -Part)
stopifnot(all(dim(lc.pearson) == dim(lc.pit)))
names(lc.pit ) <- names(lc.pearson)
## Make sure the two are lined up.
stopifnot(all.equal(lc.pit[,1:3], lc.pearson[,1:3]))

plot(as.vector(as.matrix(lc.pit[,-(1:3)])), as.numeric(as.matrix(lc.pearson[,-(1:3)])))


## Merge together for comparing the two types
lc.all <- bind_rows(cbind(type='PIT', lc.pit),
                    cbind(type='Pearson', lc.pearson)) %>%
  pivot_longer(-c(type, Yr, FltSvy, Nsamp), names_to='bin',
               values_to='residual') %>%
  mutate(FltSvy=paste0('Fleet',FltSvy),
         sex=ifelse(grepl('f', bin), 'female', 'male'),
         type=factor(type, levels=c('PIT', "Pearson"))) #%>%   filter(abs(residual)<10)
lc.all.wide <- lc.all %>%
  pivot_wider(id_cols=c(Yr, FltSvy, Nsamp, sex, bin),
              names_from='type',
              values_from='residual')
## Can carefully merge this back into the original Pearson resid
## database
x <- replist$lendbase %>%
  mutate(bin=paste0(ifelse(Sex==1, 'f', 'm'), Bin),
         sex=ifelse(Sex==1, 'female', 'male'),
         FltSvy=paste0('Fleet', Fleet))
test <- merge(x, lc.all.wide, by=c('Yr', 'FltSvy', 'bin'))
stopifnot(all.equal(test$Pearson.x, test$Pearson.y))
lc.final <- test %>% select(Yr, FltSvy, bin, Bin, sex=sex.x,
                            Pearson=Pearson.x, PIT,
                            Nsamp, Nsamp_adj, Obs, Exp) %>%
  filter(Pearson<10)

lc.final$tiny <- with(lc.final, paste0('tiny=',Obs<.99e-4 | Exp<.99e-4))

## Scatter plots
## g <- ggplot(lc.final, aes(Pearson, PIT, color=log10(Exp))) +
##   geom_abline(slope=1, intercept=0) + geom_point(alpha=.5)+ facet_grid(tiny~.)
g <- ggplot(lc.final, aes(Pearson, PIT, color=tiny)) +
  geom_abline(slope=1, intercept=0) + geom_point(alpha=.5)+ facet_grid(FltSvy~sex)
ggsave("plots/resids_scatter.png", g, width=7, height=5)
g <- ggplot(lc.final, aes(log10(Obs), log10(Exp), color=Pearson)) +
  geom_abline(slope=1, intercept=0) + geom_point(alpha=.5)+ facet_grid(FltSvy~sex)
ggsave("plots/obs_vs_exp.png", g, width=7, height=5)

## Histograms of each
g <- lc.final %>% pivot_longer(c(PIT,Pearson), names_to='type', values_to='residual') %>%
  ggplot(aes(residual, fill=type)) +
  geom_histogram(position='identity', bins=50, alpha=.5) + facet_grid(FltSvy~sex)
ggsave("plots/resids_by_type.png", g, width=7, height=3)


## ## Look at some bins closely
## g <- lc.final %>%  filter(bin %in% c('f6', 'f20', 'f52')) %>%
##   ggplot(aes(Pearson, PIT, size=Nsamp, color=FltSvy)) +
##   geom_hline(yintercept=0) + geom_vline(xintercept=0)+
##   facet_wrap('bin')+ geom_point(alpha=.7)
## ggsave("plots/compare_resids_bins.png", g, width=7, height=3)
## ## Look at some years closely
## g <- lc.all.wide%>%  filter(Yr %in% c(1973, 1983, 1999, 2020)) %>%
##   ggplot(aes(Pearson, PIT, size=Nsamp, color=FltSvy)) +
##   geom_hline(yintercept=0) + geom_vline(xintercept=0)+
##   facet_wrap('Yr')+ geom_point(alpha=.7)
## ggsave("plots/compare_resids_years.png", g, width=7, height=5)





### Part two is to look at the behavior of Pearson residuals with
### the correctly specified model

## ## Loop through and read in the comp residuals for all boottrap
## ## iterations.
## lc.list <- list()
## for(boot in 0:498){ # 0 is the original file
##   ff <- paste0('runs/BSAI_FHS/boot_',boot)
##   out <- SS_output(ff, covar=FALSE, printstats=FALSE, verbose=FALSE)
##   lc.list[[boot+1]] <- select(out$lendbase, Yr, Fleet, Sex, Bin,
##                               Pearson, Obs, Exp, Nsamp_adj) %>% cbind(boot=boot)
## }
## lcboot <- bind_rows(lc.list)
## saveRDS(lcboot, 'results/lcboot.RDS')

## Quick output prep for plotting and analysis
lcboot <- readRDS('results/lcboot.RDS')
lcboot <- mutate(lcboot, Fleet=paste0('fleet',Fleet)) %>%
  ## Drop half the years for easier plotting
  filter(Yr %% 2 == 0)
lcboot0 <- filter(lcboot, boot==0)
lcboot <- lcboot %>% filter(boot!=0)


## Make some big picture plots
g <- ggplot(lcboot, aes(Obs, Pearson, color=Nsamp_adj)) +
  geom_point(alpha=.25) + scale_x_log10() +
  facet_wrap("Fleet", scales='free_x',  ncol=1)
ggsave("plots/null_vs_obs.png", g, width=7, height=5)
g <- ggplot(lcboot, aes(Exp, Pearson, color=Nsamp_adj)) +
  geom_point(alpha=.25) + scale_x_log10() +
  facet_wrap("Fleet", scales='free_x',  ncol=1)
ggsave("plots/null_vs_exp.png", g, width=7, height=5)
g <- ggplot(lcboot, aes(Exp, Obs, color=abs(Pearson)>5))+
  geom_abline(slope=1, intercept=0)+
  geom_point(alpha=.25) + scale_x_log10() + scale_y_log10() +
  facet_wrap("Fleet", scales='free_x',  ncol=1)
ggsave("plots/obs_vs_exp.png", g, width=7, height=5)



g <- lcboot %>%
  ggplot(aes(x=factor(Yr), y=Pearson)) + geom_boxplot() +
  coord_cartesian(ylim=c(-2.5,5))+
  facet_wrap("Fleet", scales='free_x',  ncol=1) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data=lcboot0, mapping=aes(x=factor(Yr), y=Pearson),
             col=2, pch='-')
ggsave("plots/null_distribution_years.png", g, width=7, height=7)
g <- lcboot %>%
  ggplot(aes(x=factor(Bin), y=Pearson)) + geom_boxplot() +
  coord_cartesian(ylim=c(-2.5,5))+ facet_wrap("Fleet", ncol=1) +
  theme(axis.text.x = element_text(angle = 90))+
    geom_point(data=lcboot0, mapping=aes(x=factor(Bin), y=Pearson),
             col=2, pch='-')
ggsave("plots/null_distribution_bins.png", g, width=7, height=7)


## Zoom in a bit? These don't really help I don't think
lcboot2 <- group_by(lcboot, Yr, Fleet, Sex, Bin) %>%
  summarize(mean=mean(Pearson), sd=sd(Pearson), max=max(Pearson),
            range=range(Pearson),
            pvalue=ks.test(Pearson, 'pnorm', mean=0, sd=1)$p.value, .groups='drop')
## ggplot(lcboot2, aes(pvalue)) + geom_histogram()
##   scale_color_gradient(lim=c(0,1))
ggplot(lcboot2, aes(Yr, Bin, color=range)) + geom_point() +
  scale_color_gradient(low='red',  high='blue') +
  facet_wrap('Fleet')



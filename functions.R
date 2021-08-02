

#' Calculation the retrospective metrics for a single
#' bootstrapped data set. Designed to work with parallel
#' execution.
run_SS_boot_iteration <- function(boot, model.name,
                                  clean.files=TRUE, miller=FALSE){
  ## Some of these are global variables
  library(r4ss)
  if(!file.exists(file.path('models',model.name)))
    stop("model not found")
  ## The two types of bootstraps are implemented here
  if(!miller){
    dat <- SS_readdat(file.path('models', model.name,'data.ss_new'),
                      verbose=TRUE, section=2+boot)
    wd <- file.path('runs', model.name, paste0("boot_", boot))
  } else {
    ## Original data
    dat0 <- SS_readdat(file.path('models', model.name,'data.ss'),
                      verbose=TRUE, section=1)
    ## This resamples the observed data
    dat <- sample_miller_boot(boot=boot, datlist=dat0, test=FALSE)
    wd <- file.path('runs', model.name, paste0("millerboot_", boot))
  }

  dir.create(wd, showWarnings=TRUE, recursive=TRUE)
  blank.files <- list.files(file.path('models', model.name,'blank'), full.names=TRUE)
  test <- file.copy(from=blank.files, to=wd, overwrite=TRUE)
  if(!all(test)){
    message(paste(blank.files[!test], collapse= '\n'))
    stop("Some blank files failed to copy for iteration ", boot)
  }
  ## Write new data
  SS_writedat(dat, outfile=paste0(wd, '/data.ss'), verbose=FALSE,
              overwrite=TRUE)
  newfiles <- list.files(wd, full.names=TRUE) # mark for deletion
  setwd(wd)
  system('ss -nohess -iprint 1000 -nox')
  ## ## Run retro for this bootstrap one
  ## SS_doRetro(masterdir=getwd(), oldsubdir=wd,
  ##            newsubdir=paste0(wd, '/retros'),
  ##            years=peels, extras='-nohess -nox -iprint 1000')
  ## dirvec <- file.path(paste0(wd, '/retros'), paste("retro",peels,sep=""))
  ## if(length(dirvec)!=Npeels+1)
  ##   stop("Some retro runs missing in iteration ", boot)
  ## retroModels <- SSgetoutput(dirvec=dirvec, getcovar=FALSE)
  ## retroSummary <- SSsummarize(retroModels)
  ## if(model.name %in% c('BSAI_GT', 'BSAI_GT2')) {
  ##   ## What is going on here?? hack to fix error
  ##   retroSummary$startyrs <- rep(1961, length(peels))
  ## }
  ## saveRDS(retroSummary, file=paste0(file.path(wd,'retroSummary.RDS')))
  ## saveRDS(retroModels, file=paste0(file.path(wd,'retroModels.RDS')))
  ## endyrvec <- retroSummary$endyrs + peels
  ## SSplotComparisons(retroSummary, subplots=1:3, endyrvec=endyrvec, png=TRUE, plot=FALSE,
  ##                   plotdir=wd, legendlabels=paste("Data",peels,"years"))
  ## ## Calculate Mohn's rho. Not sure why I need to increase
  ## ## startyr here but errors out for some models if I don't. It
  ## ## looks like this is only relevant for WHOI rho calcs so seems
  ## ## fine to ignore.
  ## rho <- SSmohnsrho(retroSummary, endyrvec=endyrvec, startyr=retroSummary$startyr[1]+1)
  ## rhos <- data.frame(model=model.name, miller=miller, boot=boot, rho)
  ## tmp <- ifelse(miller, 'results_miller_rho.csv', 'results_rho.csv')
  ## write.csv(x=rhos, file=file.path(wd, tmp), row.names=FALSE)
  if(clean.files){
    unlink(file.path(wd, 'retros'), recursive=TRUE)
    file.remove(file.path(wd, 'retroSummary.RDS'))
    file.remove(file.path(wd, 'retroModels.RDS'))
    ## trash <-
    ##   file.remove(list.files(file.path(wd), pattern='.exe', full.names=TRUE))
    ## tmp <- file.path(wd, 'wtatage.ss')
    ## if(file.exists(tmp)) file.remove(tmp)
    ## ## this can be large so it'd be nice to keep but deleting for now
    ## file.remove(file.path(wd, 'data.ss'))
    file.remove(newfiles)
  }
    trash <-
      file.remove(list.files(file.path(getwd()), pattern='.exe',
                             full.names=TRUE))
  setwd("../../..")
}

#' Wrapper to run and save a single model
run_model <- function(reps, model.name, miller=FALSE, clean.files=TRUE){
  ## Run all bootstrap results. The clean.files argument is
  ## helpful b/c it's Nreplicates*Npeels SS runs which gets huge
  ## fast.
  trash <- sfLapply(reps, function(i)
    run_SS_boot_iteration(boot=i, model.name=model.name, clean.files=clean.files,
                          miller=miller))

  ## It fails on some. Not sure why this is happening. But hack is
  ## to loop through and figure out which failed and simply rerun
  ## them. Try 5 loops and break out if they all worked.
  tmp <- ifelse(miller, 'results_miller_rho', 'results_rho')

  for(i in 1:5){
    ff <- list.files(path=file.path('runs', model.name),
                     pattern=tmp, recursive=TRUE, full.names=TRUE)
    results <- lapply(ff, read.csv) %>% bind_rows()
    ind <- which(!reps %in% results$boot)
    if(length(ind)>0){
      warning("Rerunning failed ", model.name, " models= ", paste(ind, collapse=','))
      trash <- sfLapply(reps, function(i)
        run_SS_boot_iteration(boot=i, model.name=model.name, clean.files=clean.files,
                              miller=miller))
    }
  }
  ## Read in all final results, including those not necessarily
  ## in reps since they could have been run ea
  ff <- list.files(path=file.path('runs', model.name),
                   pattern=tmp, recursive=TRUE, full.names=TRUE)
  results <- lapply(ff, read.csv) %>% bind_rows()
  f <- paste0(model.name,ifelse(miller, '_millerboot_retros.RDS',
                                '_boot_retros.RDS'))
  saveRDS(results, file=file.path('results', f))
  message("All replicates finished for model=", model.name,
    " and miller=", miller)
}



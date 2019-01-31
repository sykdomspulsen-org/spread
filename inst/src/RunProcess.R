# 100 i oslo
# 1 i oslo
# 1 i bodø
# 1 i utsira

# R0=4/8

# 1 i oslo, 1 i bodø try w/ 3 days infectious vs 6 days infectious
# everything else 6 days infectious (r0=8)

# 12 uker



fhi::DashboardInitialiseOpinionated("commuter")

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(pbmcapply))


dirTemp <- SetupCPPAndStructure()

startVals <- CreateDataFiles()
startVals <- startVals$pop_wo_com
startVals[,value:=0]

startValsOslo1 <- copy(startVals)
startValsOslo1[location=="municip0301",value:=1]

startValsOslo100 <- copy(startVals)
startValsOslo100[location=="municip0301",value:=100]

startValsBodo1 <- copy(startVals)
startValsBodo1[location=="municip1804",value:=1]

startValsUtsira1 <- copy(startVals)
startValsUtsira1[location=="municip1151",value:=1]

m <- RunSim(
  startVals=startValsOslo100,
  R0=1.88,
  gammaTheoretical =3,
  a=1.9,
  asymptomaticProb=0.33,
  asymptomaticRelativeInfectiousness=0.5,
  verbose=T)

q <- ggplot(m,aes(x=day,y=INCIDENCE))
q <- q + geom_col()
q



stackUnique <- list(
  list(
    location="Oslo",
    scenario="1 person",
    startVals=startValsOslo1,
    gammaTheoretical=6,
    gammaEffective=3,
    a=12,
    asymptomaticProb=0,
    asymptomaticRelativeInfectiousness=1
  ),
  list(
    location="Oslo",
    scenario="1 person",
    startVals=startValsOslo1,
    gammaTheoretical=6,
    gammaEffective=6,
    a=12,
    asymptomaticProb=0,
    asymptomaticRelativeInfectiousness=1
  ),
  list(
    location="Oslo",
    scenario="100 personer",
    startVals=startValsOslo100,
    gammaTheoretical=6,
    gammaEffective=6,
    a=12,
    asymptomaticProb=0,
    asymptomaticRelativeInfectiousness=1
  ),
  list(
    location="Bodø",
    scenario="1 person",
    startVals=startValsBodo1,
    gammaTheoretical=6,
    gammaEffective=6,
    a=12,
    asymptomaticProb=0,
    asymptomaticRelativeInfectiousness=1
  ),
  list(
    location="Utsira",
    scenario="1 person",
    startVals=startValsUtsira1,
    gammaTheoretical=6,
    gammaEffective=6,
    a=12,
    asymptomaticProb=0,
    asymptomaticRelativeInfectiousness=1
  )
)

replicates <- 50
stack <- vector("list",length=length(stackUnique)*replicates)
index <- 1
for(i in seq_along(stackUnique)) for(j in 1:replicates) for(R0 in c(4,8)){
  stack[[index]] <- stackUnique[[i]]
  stack[[index]]$replicate <- j
  stack[[index]]$R0 <- R0
  stack[[index]]$id <- index
  index <- index + 1
}

res <- pbmclapply(
  stack,
  function(x){
    m <- RunSim(
      id=x$id,
      startVals=x$startVals,
      R0=x$R0,
      gammaTheoretical=x$gammaTheoretical,
      gammaEffective=x$gammaEffective,
      a=x$a,
      asymptomaticProb=x$asymptomaticProb,
      asymptomaticRelativeInfectiousness=x$asymptomaticRelativeInfectiousness,
      M=7*12,
      verbose = F)

    retval <- m[,.(
      S=sum(S),
      E=sum(E),
      I=sum(I),
      IA=sum(IA),
      R=sum(R),
      INCIDENCE=sum(INCIDENCE)
    ),by=day]

    retval[,location:=x$location]
    retval[,scenario:=x$scenario]
    retval[,R0:=x$R0]
    retval[,gammaTheoretical:=x$gammaTheoretical]
    retval[,gammaEffective:=x$gammaEffective]
    retval[,a:=x$a]
    retval[,replicate:=x$replicate]

    return(retval)
  },
  mc.cores = parallel::detectCores()
)

res <- rbindlist(res)
res[,week:=floor((day-1)/7+1)]
res[,dayOfWeek:=(day-1)%%7+1]

res[,R_initial:=min(R),
  by=.(
    location,
    scenario,
    R0,
    gammaTheoretical,
    gammaEffective,
    a,
    replicate
  )]

pd <- res[dayOfWeek==7 &
            location=="Oslo" &
            scenario=="1 person" &
            R0==4 &
            gammaEffective==3]

pd <- res[dayOfWeek==7,.(
        EI_p05=quantile(E+I,probs=c(0.05)),
        EI_p25=quantile(E+I,probs=c(0.25)),
        EI_p50=quantile(E+I,probs=c(0.5)),
        EI_p75=quantile(E+I,probs=c(0.75)),
        EI_p95=quantile(E+I,probs=c(0.95)),

        EIR_p05=quantile(E+I+R-R_initial,probs=c(0.05)),
        EIR_p25=quantile(E+I+R-R_initial,probs=c(0.25)),
        EIR_p50=quantile(E+I+R-R_initial,probs=c(0.5)),
        EIR_p75=quantile(E+I+R-R_initial,probs=c(0.75)),
        EIR_p95=quantile(E+I+R-R_initial,probs=c(0.95))
      ),keyby=.(
        day,
        location,
        scenario,
        R0,
        gammaEffective
        )]
pd[location=="Bodø"]

q <- ggplot(pd,aes(x=day,y=E+I,group=replicate))
q <- q + geom_line()
#q <- q + facet_grid(scenario~replicate,scales="free")
q



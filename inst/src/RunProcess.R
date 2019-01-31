fhi::DashboardInitialiseOpinionated("commuter")

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(doSNOW))
suppressMessages(library(iterators))


dirTemp <- SetupCPPAndStructure()
x <- CreateDataFiles()
x <- x$pop_wo_com
x[,pop:=0]
x[location=="municip0301",pop:=100]

fwrite(x[,"pop"],
       file=file.path(dirTemp,"start_infected.txt"),
       sep=" ",
       col.names=F)

# beta;  // infection parameter, 0.6
# gamma; // 1/infectious period, 1/3
# a;  // 1/latent period, 1/1.9
  res <- withr::with_dir(dirTemp,{
    processx::run(
      command=file.path(dirTemp,"infl_kommuner.exe"),
      args=c(
        as.character(0.75),
        as.character(3),
        as.character(1.9)
      ),
      echo_cmd = T,
      echo = T)
  })

  loc <- fread(file.path(dirTemp,"pop_wo_com.txt"))
  setnames(loc,c("location","pop"))
  loc[,kn:=1:.N-1]

  m <- fread(file.path(dirTemp,"cpp_res_series.txt"))
  setnames(m,c("kn","S","E","I","IA","R","INCIDENCE"))
  #m <- m[seq(1,nrow(d),2)]
  m[,day:=rep(1:(.N/2),each=2),by=kn]

  m <- merge(m,loc,by="kn")

  m <- m[,.(
    S=base::mean(S),
    E=base::mean(E),
    I=base::mean(I),
    IA=base::mean(IA),
    R=base::mean(R),
    INCIDENCE=sum(INCIDENCE)),
    by=.(
      location,day
    )]

q <- ggplot(m,aes(x=day,y=INCIDENCE))
q <- q + geom_col()
q


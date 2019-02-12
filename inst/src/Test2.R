library(data.table)
devtools::load_all()
dirTemp <- SetupCPPAndStructure()

startVals <- CreateDataFiles()
nrow(startVals$di_edge_list)
startVals$di_edge_list[from=="municip1151"]
startVals$di_edge_list[to=="municip1151"]

startValsUtsira1 <- copy(startVals$pop_wo_com)
startValsUtsira1[,pop:=NULL]
startValsUtsira1[,value:=0]
startValsUtsira1[municip=="municip1151",value:=1]

m <- RunSim(
  startVals=startValsUtsira1,
  R0=8,
  a=12,
  gammaTheoretical =6,
  gammaEffective =6,
  asymptomaticProb=0,
  asymptomaticRelativeInfectiousness=1,
  verbose=T)

setorder(m,day,location)
m[INCIDENCE>0]

m[INCIDENCE>0 & location!="municip1151"]

which(startVals$pop_wo_com$municip=="municip1151")
which(startVals$pop_wo_com$municip=="municip1223")

which(startVals$pop_wo_com$municip=="municip1424")
which(startVals$pop_wo_com$municip=="municip1412")

from          to   n
1: municip0101 municip0104  84
2: municip0101 municip0105 837
3: municip0101 municip0106 698
4: municip0101 municip0111  32
5: municip0101 municip0118 265
6: municip0101 municip0119  22
7: municip0101 municip0121   1
8: municip0101 municip0122  15
9: municip0101 municip0123   8
10: municip0101 municip0124  17

0 municip0101
1 municip0104
2 municip0105
3 municip0106
4 municip0111
5 municip0118
6 municip0119
7 municip0121
8 municip0122
9 municip0123
10 municip0124
11 municip0125
12 municip0127
13 municip0128
14 municip0135
15 municip0136
16 municip0137
17 municip0138
18 municip0211

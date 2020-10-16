# TODO
# - Bestem avgrensninger mellom MatrixCorrelation, ROSA, multiblock
# - Kompletter og funksjonaliser SO-PLS
# - Legg til mye brukte metoder
# - ASCA?
# - Documentation
#   - Kort beskrivelse av hver metode

# library(multiblock)
# library(pls) # Change to importFrom(pls, cvsegments, plsr, R2, RMSEP)
# library(SSBtools) # Change to importFrom(SSBtools, RowGroups)
# library(parallel) # Change to importFrom(parallel, makeCluser)
# cl <- makeCluster(c("localhost","localhost"))
# # cl <- makeForkCluster(nnodes=getOption("mc.cores",10L))
# pls.options(parallel=cl)
# 
# load("data/process.rda")
# B <- NULL
# set.seed(101)
# ADEm <- sopls_dirIndTotAddOver_multiple(process_data[c(1,4,5)], c("A","D","E"), c(4,3), sel.comp = c(4,3), sequential = TRUE, B=B, k=795)
# BDEm <- sopls_dirIndTotAddOver_multiple(process_data[c(2,4,5)], c("B","D","E"), c(12,4), sel.comp = c(12,4), sequential = TRUE, B=B, k=795)
# CDEm <- sopls_dirIndTotAddOver_multiple(process_data[c(3,4,5)], c("C","D","E"), c(9,2), sel.comp = c(9,2), sequential = TRUE, B=B, k=795)
# 
# pd1 <- process_data; pd1[[5]] <- pd1[[5]][,1]
# ADE1m <- sopls_dirIndTotAddOver_multiple(pd1[c(1,4,5)], c("A","D","E"), c(3,3), sel.comp = c(3,3), sequential = TRUE, B=B, k=795)
# BDE1m <- sopls_dirIndTotAddOver_multiple(pd1[c(2,4,5)], c("B","D","E"), c(5,3), sel.comp = c(5,3), sequential = TRUE, B=B, k=795)
# CDE1m <- sopls_dirIndTotAddOver_multiple(pd1[c(3,4,5)], c("C","D","E"), c(3,1), sel.comp = c(3,1), sequential = TRUE, B=B, k=795)
# 
# pd2 <- process_data; pd2[[5]] <- pd2[[5]][,2]
# ADE2m <- sopls_dirIndTotAddOver_multiple(pd2[c(1,4,5)], c("A","D","E"), c(2,2), sel.comp = c(2,2), sequential = TRUE, B=B, k=795)
# BDE2m <- sopls_dirIndTotAddOver_multiple(pd2[c(2,4,5)], c("B","D","E"), c(4,2), sel.comp = c(4,2), sequential = TRUE, B=B, k=795)
# CDE2m <- sopls_dirIndTotAddOver_multiple(pd2[c(3,4,5)], c("C","D","E"), c(7,1), sel.comp = c(7,1), sequential = TRUE, B=B, k=795)
# 
# pd3 <- process_data; pd3[[5]] <- pd3[[5]][,3]
# ADE3m <- sopls_dirIndTotAddOver_multiple(pd3[c(1,4,5)], c("A","D","E"), c(2,2), sel.comp = c(2,2), sequential = TRUE, B=B, k=795)
# BDE3m <- sopls_dirIndTotAddOver_multiple(pd3[c(2,4,5)], c("B","D","E"), c(3,1), sel.comp = c(3,1), sequential = TRUE, B=B, k=795)
# CDE3m <- sopls_dirIndTotAddOver_multiple(pd3[c(3,4,5)], c("C","D","E"), c(3,2), sel.comp = c(3,2), sequential = TRUE, B=B, k=795)
# 
# ADEm[]
# BDEm[]
# CDEm[]
# 
# ADE1m[]
# BDE1m[]
# CDE1m[]
# 
# ADE2m[]
# BDE2m[]
# CDE2m[]
# 
# ADE3m[]
# BDE3m[]
# CDE3m[]
# 
# 
# 
# # sds <- sopls_dirIndTotAddOver(process_data[1:4], process_data[[5]], c(4,2,2,3), max_comps = 11, sel.comp = c(4,2,2,3), sequential = TRUE, B = 50, k=795, type="consecutive")
# system.time(scvE <- sopls_cv(process_data[1:4], process_data[[5]][,2], c(5,8,8,4), max_comps = 10, cvsegments(795,795, type="consecutive")))
# maageSeq(scvE, expl_var=TRUE, compSeq=c(2,4,1,0), ylim=c(0,50))
# system.time(scvD <- sopls_cv(process_data[1:3], process_data[[4]], c(5,8,8), max_comps = 10, cvsegments(795,795, type="consecutive")))
# maageSeq(scvD, expl_var=TRUE, compSeq=c(2,4,2), ylim=c(0,20))
# 
# 
# pd2 <- process_data; pd2[[5]] <- pd2[[5]][,2]
# # ABCDE2 <- sopls_dirIndTotAddOver_multiple(pd2[1:5], c("A","B","C","D","E"), c(2,4,1,0), sel.comp = c(2,4,1,0), sequential = TRUE, B=B, k=795)
# ABCDE2 <- sopls_dirIndTotAddOver(pd2[1:4], pd2[[5]], c(2,4,1,0), sel.comp = c(2,4,1,0), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
# ABCDE2
# 
# ABCD2 <- sopls_dirIndTotAddOver(pd2[1:3], pd2[[4]], c(2,4,2), sel.comp = c(2,4,2), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
# ABCD2



# library(multiblock)
# library(pls) # Change to importFrom(pls, cvsegments, plsr, R2, RMSEP)
# library(SSBtools) # Change to importFrom(SSBtools, RowGroups)
# library(parallel) # Change to importFrom(parallel, makeCluser)
# cl <- makeCluster(c("localhost","localhost"))
# # cl <- makeForkCluster(nnodes=getOption("mc.cores",10L))
# pls.options(parallel=cl)
# 
# load("data/process.rda")
# B <- 100
# set.seed(101)
#ABCD2 <- multiblock:::sopls_dirIndTotAddOver(pd2[1:3], pd2[[4]], c(4,5,5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
#BCD2 <- multiblock:::sopls_dirIndTotAddOver(pd2[2:3], pd2[[4]], c(5,5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
#CD2 <- multiblock:::sopls_dirIndTotAddOver(pd2[3:3], pd2[[4]], c(5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
#ABCDE2 <- multiblock:::sopls_dirIndTotAddOver(pd2[1:4], pd2[[5]], c(4,5,5,5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
#BCDE2 <- multiblock:::sopls_dirIndTotAddOver(pd2[2:4], pd2[[5]], c(5,5,5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
#CDE2 <- multiblock:::sopls_dirIndTotAddOver(pd2[3:4], pd2[[5]], c(5,5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)
#DE2 <- multiblock:::sopls_dirIndTotAddOver(pd2[4:4], pd2[[5]], c(5), computeAdditional=TRUE, sequential = TRUE, B=B, k=795)

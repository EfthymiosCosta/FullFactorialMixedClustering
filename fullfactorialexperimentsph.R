require(MixSim)
require(cluster)
require(clustrd)
require(clustMixType)
require(mclust)
require(kmed)
require(FactoMineR)
require(fpc)
require(mclust)
require(kamila)

# categorized numerical variable function
intv <- function(vec, class) {
  nbase <- (1:(class-1))/class
  nq <- numeric(length(nbase))
  for (i in 1:length(nq)) {
    nq[i] <- quantile(vec, nbase[i])
  }
  res <- c(min(vec), nq, max(vec)) 
  res[1] <- res[1]-1
  for (i in 2:length(res)){
    if (res[i-1]==res[i]){
      res[i] <- res[i]+2e-15
    }
  }
  return(res)
}

# Full factorial simulation study - SPHERICAL CLUSTERS

# Empty data frame to store simulation results
fullfactorialsph <- data.frame(seed=numeric(),
                            nClust=numeric(),
                            overlap=numeric(),
                            nrows=numeric(),
                            ncols=numeric(),
                            pi=numeric(),
                            method=character(),
                            res=numeric(),
                            stringsAsFactors=FALSE)

cluster_sizes <- c(3:5)
overlap_levels <- c(0.01, seq(0.05, 0.20, by=0.05))
num_rows <- c(300, 600, 1200)
num_cols <- c(6, 10, 14)
pi_values <- c(0.01, 1.0)

nreps <- 50
tot <- nreps*length(cluster_sizes)*length(overlap_levels)*length(num_rows)*length(num_cols)*length(pi_values)

for (l in 0:(nreps-1)){
  # Random seed
  set.seed(1234+l)
  # Number of clusters
  for (clusters in cluster_sizes){
    # Overlap level
    for (overlap in overlap_levels){
      # Number of rows
      for (rows in num_rows){
        # Number of variables/columns
        for (columns in num_cols){
          # Balanced/unbalanced design
          for (pi_val in pi_values){
            # Construct artificial data set - sphericity here is set to TRUE
            ifelse((overlap==0.01 & clusters==5), maxom <- 7*overlap/4, maxom <- 3*overlap/2)
            mixsimaux <- MixSim(BarOmega = overlap, MaxOmega = maxom, PiLow=pi_val, K = clusters, p = columns, sph = TRUE, resN = 10000000)
            mixdtaux <- simdataset(n = rows, Pi = mixsimaux$Pi, Mu = mixsimaux$Mu, S = mixsimaux$S)
            # Discretise first half attributes using 4 levels
            for (k in 1:(columns/2)){
              mixdtaux$X[,k] <- 
                as.factor(cut(mixdtaux$X[,k], intv(mixdtaux$X[,k], 4), labels = (1:4)))
            }
            mixdt1df <- as.data.frame(mixdtaux$X)
            for (k in 1:(columns/2)){
              mixdt1df[,k] <- as.factor(mixdt1df[,k])
            }
            # this temporarily fixes a nasty bug of cluspcamix
            colnames(mixdt1df) <- sprintf("a%d", 1:columns)
            save(mixdt1df,file=paste("mixdt1dfsph",run,sep = "",".Rdata"))
            # Calculate Gower's dissimilarities and Ahmad distances
            gower_dist <- daisy(mixdt1df, metric = "gower")
            mix_sim <- distmix(mixdt1df, method = "ahmad", idcat = c(1:(columns/2)), idnum = c((columns/2+1):columns))
            # Specify continuous & categorical attributes
            conDf <- data.frame(scale(mixdt1df[,((columns/2+1):columns)]))
            catDf <- dummyCodeFactorDf(data.frame(mixdt1df[,1:(columns/2)]))
            # Extract principal components
            outpcamix <- FAMD(mixdt1df, ncp = clusters -1, graph=FALSE)
            # PAM with Gower's
            pam_fit <- pam(gower_dist, diss = TRUE, k = clusters, do.swap = FALSE, cluster.only = TRUE,
                           nstart=100)
            fullfactorialsph <- rbind(fullfactorialsph, data.frame(seed=l+1, nClust=clusters, overlap=overlap,
                                                             nrows=rows, ncols=columns, pi=pi_val,
                                                             method='PAM', res=adjustedRandIndex(pam_fit, mixdtaux$id)))
            save(fullfactorialsph, file='fullfactorialsph.RData')
            cat('PAM done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            # K-prototypes
            outk <- kproto(mixdt1df, clusters, nstart = 100,verbose = FALSE)
            fullfactorialsph <- rbind(fullfactorialsph, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                             ncols=columns, pi=pi_val, method='K-Prot',
                                                             res=adjustedRandIndex(outk$cluster, mixdtaux$id)))
            save(fullfactorialsph, file='fullfactorialsph.RData')
            cat('K-prototypes done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            # Mixed K-means
            kmedres <- fastkmed(mix_sim, clusters, iterate = 100, init = sample(1:nrow(mixdt1df), clusters))
            fullfactorialsph <- rbind(fullfactorialsph, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                             ncols=columns, pi=pi_val, method='Mixed',
                                                             res=adjustedRandIndex(kmedres$cluster, mixdtaux$id)))
            save(fullfactorialsph, file='fullfactorialsph.RData')
            cat('Mixed K-means done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            # Modha-Spangler
            msRes <- gmsClust(conDf, catDf, nclust = clusters, searchDensity = 5, nstart=100)
            fullfactorialsph <- rbind(fullfactorialsph, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                             ncols=columns, pi=pi_val, method='Modha-Spangler',
                                                             res=adjustedRandIndex(msRes$results$cluster, mixdtaux$id)))
            save(fullfactorialsph, file='fullfactorialsph.RData')
            cat('Modha-Spangler done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            # FAMD + K-means
            outkm <- kmeans(outpcamix$ind$coord, nstart=100, clusters, iter.max=1000, algorithm = "MacQueen")
            fullfactorialsph <- rbind(fullfactorialsph, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                             ncols=columns, pi=pi_val, method='FAMD',
                                                             res=adjustedRandIndex(outkm$cluster, mixdtaux$id)))
            save(fullfactorialsph, file='fullfactorialsph.RData')
            cat('FAMD done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            # Mixed RKM
            outmix = clustrd::cluspcamix(mixdt1df, nclus = clusters, ndim = clusters-1, nstart=100)
            fullfactorialsph <- rbind(fullfactorialsph, data.frame(seed=l+1, nClust=clusters, overlap=overlap, nrows=rows,
                                                             ncols=columns, pi=pi_val, method='RKM',
                                                             res=adjustedRandIndex(outmix$cluster, mixdtaux$id)))
            save(fullfactorialsph, file='fullfactorialsph.RData')
            cat('Mixed RKM done for dataset',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            cat('Simulations for data set',l+1,'with',clusters,'clusters,',rows,'rows',columns,
                'columns, an overlap of',overlap, 'and pi',pi_val,'\n')
            cat('*****  Run',run,'/',tot,' *****','\n')
            run <- run + 1
          }
        }
      }
    } 
  }
}
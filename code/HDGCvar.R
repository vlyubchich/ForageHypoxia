# devtools::install_github("Marga8/HDGCvar")

library(HDGCvar)
library(igraph)
x = HDGC_VAR_all(data=sample_dataset_I1, p=2, d=2, parallel=TRUE )

# matrix of test p-values
x$tests[,,2,1]

Plot_GC_all(x, Stat_type="FS_cor", #alpha=0.01,
            multip_corr=list(F),
            mode = "directed",
            layout=layout.circle,
            # main="Network",
            edge.arrow.size=.2,vertex.size=5,
            # vertex.color=c("lightblue"),
            # vertex.frame.color="blue",
            # vertex.label.size=2,
            # vertex.label.color="black",
            # vertex.label.cex=0.6,
            # vertex.label.dist=1,
            # edge.curved=0.
            cluster=list(T,5,"black",0.8,1,0)
            )


# X = HDGC_VAR_all(data=Dwide_mat[,1:30], p=30, d=2, parallel=TRUE )
#
# Error in checkForRemoteErrors(val) :
#     7 nodes produced errors; first error: the leading minor of order 332 is not positive definite

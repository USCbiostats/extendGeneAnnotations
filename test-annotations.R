library(aphylo)

set.seed(1231)
x <- raphylo(100, P = 2)
x <- rdrop_annotations(x, .6)
# GO:0000749 GO:0071444
ann <- x$tip.annotation
colnames(ann) <- c("GO:0000749", "GO:0071444")

write.csv(ann, file = "test-annotations.csv", row.names = FALSE)

library(snpStats)
data(testdata, package="snpStats")
x <- Autosomes[1:10,1:10]
y <- Autosomes[11:20,6:15]
z <- Autosomes[1:20,1:15]

test_that("snpmatrix.combine", {
    z2 <- snpmatrix.combine(x,y)
    expect_equal(ncol(z2),ncol(z))
    expect_equal(nrow(z2),nrow(z))
    expect_identical(z[1:10,1:10], z2[1:10,colnames(z)[1:10]])
    expect_identical(z[11:20,6:15], z2[11:20,colnames(z)[6:15]])
    expect_identical(z2[11:20,6:15],z[11:20,6:15])
    expect_false(identical(z2,z))
})

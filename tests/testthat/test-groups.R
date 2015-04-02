
gr <- new("groups",
          list(letters[c(1,5,9,15,21)],letters[2:4],c("u","i","o","p")),
          tags=c("a","b","u"))
tg <- new("tags",
          .Data=unlist(gr@.Data),
          tags=rep(gr@tags,times=sapply(gr@.Data,length)))

context("Conversion groups <-> tags")

test_that("conversion creates correct classes", {
  expect_is(as(gr,"tags"), "tags")
  expect_is(as(tg,"groups"), "groups")
})

test_that("conversion is lossless", {
  expect_identical(as(gr,"tags"),tg)
  ## check groups by slots to avoid conflict in names
  gr2 <- as(tg,"groups")
  expect_identical(gr@.Data,gr2@.Data)
  expect_identical(gr@tags,gr2@tags)
})


test_that("subsetting", {
  expect_is(gr[[2]],"character")
  expect_is(gr[["a"]],"character")
  expect_is(gr[1:2],"groups")
  expect_identical(gr[1:2],gr[c("a","b")])
  expect_equal(length(gr[["a"]]),5)
})

test_that("union", {
  ugr <- union(gr[1:2],gr[2:3])
  cgr <- c(gr[1:2],gr[2:3])
  expect_equal(length(ugr),2)
  expect_equal(length(cgr),4)
})

test_that("accessors", {
    expect_is(snps(gr),"list")
    expect_is(tags(gr),"character")
    expect_is(snps(tg),"character")
    expect_is(tags(tg),"character")
})


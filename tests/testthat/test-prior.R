library(ALDEx2)
set.seed(81)

## Setup
reads <- rbind(c(900, 12,  0, 18),
               c(  7, 10,  1, 42),
               c(  0,  0, 15,  0))

reads_2 <- rbind(c(900, 12,  0, 18),
                 c(  0,  0,  0,  0),
                 c(  7, 10,  1, 42),
                 c(  0,  0, 15,  0),
                 c(  0,  0,  0,  0))

row.names(reads) <- paste0("taxa_", 1:nrow(reads))
colnames(reads) <- paste0("sample_", 1:ncol(reads))
row.names(reads_2) <- paste0("taxa_", 1:nrow(reads_2))
colnames(reads_2) <- paste0("sample_", 1:ncol(reads_2))

## Test Errors
test_that("aldex2 gives error for prior vector", {
    expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=c(0.1, 0.2)),
                 "prior should be a single numeric")
})

test_that("aldex2 wrong input for prior", {
    expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior="A"),
                 "prior must be numeric")
})

test_that("aldex2 wrong input for prior NULL", {
    expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=NULL),
                 "prior must be numeric")
})

test_that("aldex2 gives error for prior matrix wrong rows", {
      prior <- matrix(0.5, nrow=3, ncol=5)
      expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                   "prior is a matrix")
})

test_that("aldex2 gives error for prior matrix wrong cols", {
      prior <- matrix(0.5, nrow=5, ncol=4)
      expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                   "prior is a matrix")
})

test_that("aldex2 gives error for prior array wrong dim 1", {
      prior <- array(0.5, c(5, 4, 150))
      expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                   "prior is an array")
})

test_that("aldex2 gives error for prior array wrong dim 2", {
      prior <- array(0.5, c(3, 5, 150))
      expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                   "prior is an array")
})

test_that("aldex2 gives error for prior array wrong dim 3", {
      prior <- array(0.5, c(3, 4, 151))
      expect_error(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                   "prior is an array")
})

## Test Warnings
test_that("aldex2 warning small prior", {
    prior <- matrix(0.5, nrow=3, ncol=4)
    prior[2, 3] <- 0.4
    expect_warning(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                  regexp="lead to high posterior sampling variance")
})

test_that("aldex2 warning extremly prior", {
    prior <- array(0.5, c(3,4,150))
    prior[1,2,12] <- 0.05
    expect_warning(aldex(reads, c(0,0,1,1), mc.samples=150, prior=prior),
                   regexp="lead to extremly high posterior sampling variance")
})

## Test Expected Functionality
test_that("aldex2 prior single numeric works", {
    aldex.clr.obj <- aldex.clr(reads, mc.samples=150, prior=1000000)
    dir.data <- getDirichletInstances(aldex.clr.obj)
    expect_true(all(round(unlist(dir.data),2)==0.33))
})

test_that("aldex2 prior numeric matrix works", {
    prior <- matrix(0.5, nrow=3, ncol=4)
    prior[1,2] <- 1000000
    aldex.clr.obj <- aldex.clr(reads, mc.samples=150, prior=prior)
    dir.data <- getDirichletInstances(aldex.clr.obj)
    expect_true(dir.data$sample_1[1,65]>0.9)
})

test_that("aldex2 prior numeric matrix works 2", {
    prior <- matrix(0.5, nrow=3, ncol=4)
    prior[,2] <- c(1000000, 1000000, 1000000)
    prior[1,4] <- 1000000
    aldex.clr.obj <- aldex.clr(reads, mc.samples=150, prior=prior)
    dir.data <- getDirichletInstances(aldex.clr.obj)
    expect_true(all(round(dir.data$sample_2, 2)==0.33))
    expect_true(all(dir.data$sample_4[1,]>0.98))
    expect_true(dir.data$sample_1[1,65]>0.9)
    expect_true(dir.data$sample_1[3,65]<0.1)
})

test_that("aldex2 prior numeric matrix works 3", {
    prior <- matrix(0.5, nrow=5, ncol=4)
    prior[,2] <- c(100000, 100000, 0.5, 0.5, 100000)
    prior[,3] <- c(1000000, 0.5, 1000000, 1000000, 0.5)
    aldex.clr.obj <- aldex.clr(reads_2, mc.samples=150, prior=prior)
    dir.data <- getDirichletInstances(aldex.clr.obj)
    expect_true(all(dir.data$sample_2[1,]>0.98))
    expect_true(all(round(dir.data$sample_3, 2)==0.33))
})

test_that("aldex2 prior numeric array works", {
    prior <- array(0.5, c(3, 4, 150))
    prior[3, 1, 122] <- 1000000
    aldex.clr.obj <- aldex.clr(reads, mc.samples=150, prior=prior)
    dir.data <- getDirichletInstances(aldex.clr.obj)
    expect_true(dir.data$sample_1[3,122]>0.99)
    expect_true(dir.data$sample_1[3,123]<0.1)
    expect_true(dir.data$sample_4[3,122]<0.1)
})

test_that("aldex2 prior numeric array works 2", {
    prior <- array(0.5, c(5, 4, 150))
    prior[3, 1, 122] <- 1000000
    prior[2, 1, 124] <- 1000000
    aldex.clr.obj <- aldex.clr(reads_2, mc.samples=150, prior=prior)
    dir.data <- getDirichletInstances(aldex.clr.obj)
    expect_true(dir.data$sample_1[2,122]>0.99)
    expect_true(dir.data$sample_1[2,123]<0.1)
    expect_true(dir.data$sample_4[3,122]<0.1)
    expect_true(dir.data$sample_4[3,124]<0.1)
})

## Test Getter
test_that("aldex2 getPrior works as expected", {
    aldex.clr.obj <- aldex.clr(reads, mc.samples=150, prior=0.92)
    expect_true(getPrior(aldex.clr.obj)==0.92)
})

test_that("aldex2 getPrior works as expected matrix", {
    prior <- matrix(1:20, nrow=5, ncol=4)
    aldex.clr.obj <- aldex.clr(reads_2, mc.samples=150, prior=prior)
    expected_prior <- rbind(c(1,6,11,16),
                            c(3,8,13,18),
                            c(4,9,14,19))
    expect_equal(getPrior(aldex.clr.obj), expected_prior)
})

test_that("aldex2 getPrior works as expected array", {
    prior <- array(1:3000, c(5,4,150))
    aldex.clr.obj <- aldex.clr(reads_2, mc.samples=150, prior=prior)
    expected_prior <- rbind(c(141,146,151,156),
                            c(143,148,153,158),
                            c(144,149,154,159))
    expect_equal(getPrior(aldex.clr.obj)[,,8], expected_prior)
})

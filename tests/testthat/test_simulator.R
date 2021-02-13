test_that("Simulator creates proper object",{
	output<-SimulateHisse(c(.3, .1), c(0, 0), matrix(c(NA, 0.2, 0.3, NA), nrow=2), max.taxa=15, x0=1)
	expect_equal(class(SimToPhylo(output$result)), "phylo")
})

test_that("Simulator creates proper object with multiple starting taxa",{
	output<-SimulateHisse(c(.3, .1), c(0, 0), matrix(c(NA, 0.2, 0.3, NA), nrow=2), max.taxa=15, x0=1, nstart=2)
	expect_equal(class(SimToPhylo(output$result)), "phylo")
})

test_that("Simulator generates correct number of taxa",{
	output<-SimulateHisse(c(.3, .1), c(0, 0), matrix(c(NA, 0.2, 0.3, NA), nrow=2), max.taxa=15, x0=1)
	expect_equal(output$n.surviving, 15)
})

test_that("Simulator stops at correct height",{
	output<-SimulateHisse(c(.3, .1), c(0, 0), matrix(c(NA, 0.2, 0.3, NA), nrow=2), max.t=15, x0=1)
	expect_equal(max(output$results$height), 15)
})

test_that("Simulator runs with mass extinctions",{
	output<-SimulateHisse(c(.3, .1), c(0, 0), matrix(c(NA, 0.2, 0.3, NA), nrow=2), max.t=15, x0=1,
	mass.extinction.heights=c(3,10), mass.extinction.magnitudes=c(0.5, 0.8))
	while(output$n.surviving==0) {
		output<-SimulateHisse(c(.3, .1), c(0, 0), matrix(c(NA, 0.2, 0.3, NA), nrow=2), max.t=15, x0=1,
	mass.extinction.heights=c(3,10), mass.extinction.magnitudes=c(0.5, 0.8))	
	}
	expect_equal(max(output$results$height), 15)
})

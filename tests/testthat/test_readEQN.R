
eqnfile <- system.file("MPTmodels/2htm.eqn", package = "TreeBUGS")
model <- "# 2HTM
Target    Hit    Do
Target    Hit    (1-Do)*g
Target    Miss   (1-Do)*(1-g)
Lure      FA     (1-Dn)*g
Lure      CR     (1-Dn)*(1-g)
Lure      CR     Dn
"


test_that("readEQN results in proper MPT model", {
  expect_is(readEQN(eqnfile), "data.frame")
  expect_is(readEQN(model), "data.frame")

  res <- readEQN(model)
  expect_named(res, c("Tree", "Category", "Equation", "EQN"))

  expect_is(
    readEQN(model, restrictions = list("g=.5", "Dn=Do")),
    "data.frame"
  )

  m1 <- readEQN(model)
  expect_identical(
    m1$EQN,
    c("Do", "(1-Do)*g", "(1-Do)*(1-g)", "(1-Dn)*g", "(1-Dn)*(1-g)", "Dn")
  )

  m2 <- readEQN(model, restrictions = list("g=.5", "Dn=Do"))
  expect_identical(
    m2$EQN,
    c("Do", "(1-Do)*.5", "(1-Do)*(1-.5)", "(1-Do)*.5", "(1-Do)*(1-.5)", "Do")
  )

  m3 <- readEQN(model, restrictions = list("Dn=Do=g"))
  expect_identical(
    m3$EQN,
    c("g", "(1-g)*g", "(1-g)*(1-g)", "(1-g)*g", "(1-g)*(1-g)", "g")
  )

  m4 <- readEQN(model, restrictions = list("Dn=0.273", "g=0.93"))
  expect_identical(
    m4$EQN,
    c("Do", "(1-Do)*0.93", "(1-Do)*(1-0.93)", "(1-0.273)*0.93", "(1-0.273)*(1-0.93)", "0.273")
  )

  # misspecified constraints
  expect_error(readEQN(model, restrictions = list("g=.5", "Dnasd=Do")))
  expect_error(readEQN(model, restrictions = list("g=.1=.5")))
  expect_warning(readEQN(model, restrictions = list("g=-4")))

  # parsed model matrices
  mod <- readEQN(model, restrictions = list("g=.5", "Dn=Do"), parse = TRUE)
  expect_named(mod$Table, c("Tree", "Category", "Equation", "EQN"))

})

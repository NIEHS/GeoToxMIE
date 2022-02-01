model <- c("Hill_three", "Hill_three", "Hill_three", "Hill_three")
param <- matrix(c(3.94e-5, 0.97, 1.9, 5.16e-4, 1.50, 1.9, 3.43e-6, 1.04, 1.9, 9.18e-6, 0.77, 1.9), nrow = 4, ncol = 3, byrow = TRUE)
rownames(param) <- c('Ni', 'Zn', 'Cu', 'Mn')
colnames(param) <- c('Alpha', 'Beta', 'Gamma')
## example 1
# using GCA to predict the mixtures designed by equal effect concentration ratio (eecr) at # the effect concentration of EC05 and EC50
# the eecr mixture design is based on four heavy metals (four factors).
mdl3 <- gcaHill(model, param, mixType = "eecr", effv = c(0.05, 0.5), rtype = 'continuous')


model <- c("Hill_two", "Hill_two", "Hill_two", "Hill_two")
param <- matrix(c(3.94e-5, 0.97, 0, 5.16e-4, 1.50, 0, 3.43e-6, 1.04, 0, 9.18e-6, 0.77, 0), nrow = 4, ncol = 3, byrow = TRUE)
rownames(param) <- c('Ni', 'Zn', 'Cu', 'Mn')
colnames(param) <- c('Alpha', 'Beta', 'Gamma')
## example 1
# using GCA to predict the mixtures designed by equal effect concentration ratio (eecr) at # the effect concentration of EC05 and EC50
# the eecr mixture design is based on four heavy metals (four factors).
mdl2 <- gcaHill(model, param, mixType = "eecr", effv = c(0.05, 0.5), rtype = 'continuous')

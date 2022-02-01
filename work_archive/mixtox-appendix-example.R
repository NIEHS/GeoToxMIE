par(mfrow = c(1, 2))
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
eeca <- caPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))
eeia <- iaPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))
# plot EE05 mixture
par(mar = c(5, 5, 1, 1))
x <- antibiotox$ee05$x
expr <- antibiotox$ee05$y
ee05fit <- curveFit(x, expr, eq = antibiotox$eecr.mix$model[1],
                    param = antibiotox$eecr.mix$param[1, ])
plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
     xlab = "log(c) mol/L", ylab = "Inhibition[]", cex = 1.8, cex.lab = 1.8,
     cex.axis = 1.8)
lines(log10(x), ee05fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
lines(log10(x), ee05fit$crcInfo[, 6] * 100, col = "green", lwd = 1.5, lty = 3)
lines(log10(x), ee05fit$crcInfo[, 7] * 100, col = "green", lwd = 1.5,  lty = 3)
lines(log10(eeia$ia[1, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
lines(log10(eeca$ca[1, ]), eeca$e * 100, col = "blue", lwd = 2.5, lty = 2)
legend("topleft", legend = "EE05", border = "white", cex = 1.4)
# plot EE50 mixture
par(mar = c(5, 5, 1, 1))
x <- antibiotox$ee50$x
expr <- antibiotox$ee50$y
ee50fit <- curveFit(x, expr, eq = antibiotox$eecr.mix$model[1],
                    param = antibiotox$eecr.mix$param[1, ])
plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
     xlab = "log(c) mol/L", ylab = "", cex = 1.8, cex.lab = 1.8, cex.axis = 1.8)
lines(log10(x), ee50fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
lines(log10(x), ee50fit$crcInfo[, 6] * 100, col = "green", lwd = 1.5, lty = 3)
lines(log10(x), ee50fit$crcInfo[, 7] * 100, col = "green", lwd = 1.5, lty = 3)
lines(log10(eeia$ia[2, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
lines(log10(eeca$ca[2, ]), eeca$e * 100, col = "blue", lwd = 2.5, lty = 2)
legend("topleft", legend = "EE50", border = "white", cex = 1.4)

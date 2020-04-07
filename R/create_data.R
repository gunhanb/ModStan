dose     = c(0, 300, 200, 300, 100, 300)
schedule = c(1, 1, 2, 2, 3, 3)
freq     = c(24 * 7, 24 * 7, 24 * 7 * 2,  24 * 7 * 2, 24 * 7 * 4, 24 * 7 * 4)
resp     = c(-18.1, -73.7, -65.4, -68.2, -44.8, -63.5)
sigma    = c(5.2, 5.2, 5.2, 5.1, 5.0, 4.9)
## For complete pooling
dose_rescaled = c(0, 600, 200, 300, 50, 150)

dat.Dupilumab = data.frame(dose = dose,
                           schedule = schedule,
                           resp = resp,
                           freq = freq,
                           sigma = sigma,
                           dose_rescaled = dose_rescaled)


save(dat.Dupilumab, file = "data/dat.Dupilumab.RData")

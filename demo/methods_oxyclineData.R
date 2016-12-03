
# Set directories where the Echopen's outputs are located
fileMode <- list(fish38_file   = system.file("extdata", "fish38.mat", package = "oXim"),
                 fluid120_file = system.file("extdata", "fluid120.mat", package = "oXim"),
                 blue38_file   = system.file("extdata", "blue38.mat", package = "oXim"))


# Read echograms (echoData object)
echoData <- readEchograms(fileMode = fileMode)

# Calculate oxycline limits (oxyclineData object)
oxyLimits <- getOxyrange(fluidMatrix = echoData)

# Print method
print(oxyLimits)

# Summary method
summaryOxyLimits <- summary(oxyLimits)

# Print.summary method
print(summaryOxyLimits)

# Plot method
plot(oxyLimits)

# echogramPlot method
echogramPlot(oxyLimits)

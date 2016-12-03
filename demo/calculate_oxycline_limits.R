
# Calculate oxycline limits from echoData object -------------------------------

# Set directories where the Echopen's outputs are located
fileMode <- list(fish38_file   = system.file("extdata", "fish38.mat", package = "oXim"),
                 fluid120_file = system.file("extdata", "fluid120.mat", package = "oXim"),
                 blue38_file   = system.file("extdata", "blue38.mat", package = "oXim"))


# Read echograms (echoData object)
echoData <- readEchograms(fileMode = fileMode)

# Calculate oxycline limits
oxyLimits <- getOxyrange(fluidMatrix = echoData)

# Show main information using print method
print(oxyLimits)

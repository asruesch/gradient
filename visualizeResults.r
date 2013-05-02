source("G:/gradient/code/gradientFunctions.r")

up = 71237
down = 482

new = elevProfile(up, down, shedEdges, "minElevFix2", "maxElevFix2")
old = elevProfile(up, down, shedEdges, "minElevFix1", "maxElevFix1")


plot(minelv ~ l
     , data = old
     , xlab = "distance downstream (m)"
     , ylab = "minimum segment elevation (m)"
     , type = "l"
     , col = "blue"
     , lwd = 2)
lines(minelv ~ l, data = new, col = "red", lwd = 2)

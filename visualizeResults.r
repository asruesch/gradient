source("G:/gradient/gradientFunctions.r")

up = 213452
down = 224856

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

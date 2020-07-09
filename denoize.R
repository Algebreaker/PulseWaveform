# Denoizing attempt:

# Example ear ppg data (1 wave):
data <- read.csv("~/Desktop/AK465-T0-__BH-R1-PHYS copy.csv", header = T)
data <- data[504500:507000, 3:5]

# Apply scaling, multiplier and offset:
test2 <- round(((data$Ear.PulseOX...PPG..OXY100E..Pulse.*800)-(-22.7093))/0.9765714)
test <- ((data$Ear.PulseOX...PPG..OXY100E..Pulse.*800)-(-22.7093))/0.9765714

# Define mode function:
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Apply rolling mode iteratively:
for(i in 1:13){
  for(j in i:length(test2)){
    test2[j] <- mode(test2[(j-i):(j+i)])
  }
}

# Plot:
plot(test)
lines(test2, col = "red")
plot(test2, type = "l")
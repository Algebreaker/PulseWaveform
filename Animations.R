# Animations of pulse, respiration and HR data using gganimate. 

# To run this, you will need to have run a time-series though the main R script so as to generate 
# both the undetrended_data dataframe as well as the pulse_stacked dataframe (and w$w_poly_peaks). 

library(stringr)
library(transformr)
library(gifski)
library(png)
library(gganimate)
library(magick)
library(dplyr)


##### Pulse Animation #####

# Convert pulse_stacked wave_ID column to numeric form, and remove NA values from all waves. 
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}
for(i in 1:nrow(pulse_stacked)){
  pulse_stacked[i, 2] <- as.numeric(numextract(pulse_stacked[i, 2]))
}
pulse_stacked$wave_ID <- as.numeric(pulse_stacked$wave_ID)
pulse_stacked_nonas <- pulse_stacked[!is.na(pulse_stacked$values), ]
# Create the plot to be animated:
p <- ggplot(pulse_stacked_nonas, aes(x = x, y = values), col = wave_ID) + geom_line() 
# Create the animation:
anim_p <- p + transition_states(wave_ID, transition_length = 2, state_length = 2) + ggtitle('Wave {closest_state}')
# Render as a gif:
a_gif <- animate(anim_p, renderer = gifski_renderer(), nframes = (length(w$w_poly_peaks)*4))



##### Respiration Animation #####

# Create a dataframe with RIP values:
RIP_z <- data_frame(1:length(w$w_poly_peaks))
RIP_z <- cbind(RIP_z, undetrended_data$RIP[w$w_poly_peaks])   # note respiration values are those that correspond in time to the PPGs w peaks
colnames(RIP_z)[1] <- "waveID"
colnames(RIP_z)[2] <- "RIP"
# Create plot:
r <- ggplot(RIP_z, aes(x = waveID, y = RIP)) + geom_point(size = 3) + ylim(-0.005, 0.005) 
# Create animation:
anim_r <- r + transition_states(waveID, transition_length = 2, state_length = 2) + view_follow(fixed_y = TRUE) + ggtitle('Respiration at wave {closest_state}')
# Render:
b_gif <- animate(anim_r, renderer = gifski_renderer(), nframes = (length(w$w_poly_peaks)*4))



##### Heart Rate Animation #####

# Create dataframe:
HR_z <- data_frame(1:length(w$w_poly_peaks))
HR_z <- cbind(HR_z, undetrended_data$Heart.Rate.PulseOx1[w$w_poly_peaks])
colnames(HR_z)[1] <- "waveID"
colnames(HR_z)[2] <- "HR"
# Create plot:
h <- ggplot(HR_z, aes(x = waveID, y = HR)) + geom_line() + geom_point()
# Create animation:
anim_h <- h + transition_reveal(waveID) + ggtitle('Heart Rate')
# Render:
c_gif <- animate(anim_h, renderer = gifski_renderer(), nframes = (length(w$w_poly_peaks)*4))



##### Combine animations #####

a_mgif <- image_read(a_gif)
b_mgif <- image_read(b_gif)
c_mgif <- image_read(c_gif)   
new_gif <- image_append(c(a_mgif[1], b_mgif[1], c_mgif))
for(i in 2:(length(w$w_poly_peaks)*4)){ 
  combined <- image_append(c(a_mgif[i], b_mgif[i], c_mgif[i]))
  new_gif <- c(new_gif, combined)
}
# View:
new_gif

# From here, you can use the anim_save() function to export the combined gif. I have found it simpler to 
# save from the webpage that automatically opens when clicking to expand the animation in R. 


# Useful links for exploring other options with gganimate:
https://gganimate.com/articles/gganimate.html#rendering-1
https://www.datanovia.com/en/blog/gganimate-how-to-create-plots-with-beautiful-animation-in-r/
https://ugoproto.github.io/ugo_r_doc/pdf/gganimate.pdf

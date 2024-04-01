#install.packages("chromoMap")
library(chromoMap)
#packageVersion("chromoMap")

#chromosome_length <- c(307041717, 244442276, 235667834, 246994605, 223902240, 174033170, 182381542, 181122637, 159769782, 150982314)

chromosome_data <- data.frame(
  chromosome = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  length = c(307041717, 244442276, 235667834, 246994605, 223902240, 174033170, 182381542, 181122637, 159769782, 150982314)
)

hotspot <- read.table("C:/Users/harel/Downloads/hotspot_data.txt", header = FALSE)

# Extracting the desired columns from the results data frame
hotspot_data <- results[, c("Chromosome","Hotspot.start.position..bp.", "Hotspot.end.position..bp.")]

# Renaming the columns to 'start' and 'end'
colnames(hotspot_data) <- c("chromosome", "start", "end")

# Display the first few rows of the hotspot_data to verify
head(hotspot_data)
write.table(hotspot_data,"C:/Users/harel/Downloads/hotspot_data.txt")

# Create the chromosome map- option 1
  #chromoMap("C:/Users/harel/Downloads/ChromosomeLength.txt", "C:/Users/harel/Downloads/hotspot_data.txt",
   #       labels=T, chr_color="grey",scale = 0.1,n_win.factor = 10,
   #      win.summary.display=T)
# Create the chromosome map
chromoMap("C:/Users/harel/Downloads/ChromosomeLength.txt", "C:/Users/harel/Downloads/hotspot_data.txt",
          anno_col= c("black"),
          chr_color = c("grey"),
          n_win.factor = 1,
          chr.2D.plot = T,
          labels=F,
          label_angle = -65,
          chr_length = 5,
          chr_width = 8,
          plot_height = 20,
          canvas_width = 800,          
          canvas_height = 1000,
          win.summary.display=T,
          export.options=TRUE,
          legend = T,
          ch2D.lg_y=100,
          plot.legend.labels = c("AC","AD","FA","HC","PA","WE"),
          #plot_y_domain = list(c(0, 2)),
          plot_filter = list(c("col","byCategory")),
          ch2D.colors = c('#82C9ED', '#CC6677', '#DCCA71', '#107632', '#413190', '#882255'))

# Create the chromosome map- option 2
#chromoMap("C:/Users/harel/Downloads/ChromosomeLength.txt", "C:/Users/harel/Downloads/hotspot_data.txt",
#       labels=T, chr_color="grey",scale = 0.1,n_win.factor = 10,
#      win.summary.display=T)

# Create the chromosome map- option 3
#chromoMap("C:/Users/harel/Downloads/ChromosomeLength.txt", "C:/Users/harel/Downloads/hotspot_data.txt",
 #         labels=F,
 #        label_angle = -65,
 #       chr_length = 6,
 #        chr_width = 25,
 #       canvas_width = 800,
 #        data_based_color_map = T,
 #       plot_filter = list(c("col","byCategory")),
 #      data_type = "categorical",
 #     data_colors = list(c('#82C9ED', '#CC6677', '#DCCA71', '#107632', '#413190', '#882255')))

library(ggplot2)
library(MASS)
library(patchwork)

###Load the functions
source("OCP_Procedures.R")
source("plot_generating_function.R")

#Generates Figure 4 of the paper "The Online Closure Principle"
ggsave("Figure4.pdf", plot=plot_generator(1000,4,0))

#Generates Figure 5 of the paper "The Online Closure Principle"
ggsave("Figure5.pdf", plot=plot_generator(1000,4,-2))

#Generates Figure 6 of the paper "The Online Closure Principle"
ggsave("Figure6.pdf", plot=plot_generator(100,4,0))

#Generates Figure 7 of the paper "The Online Closure Principle"
ggsave("Figure7.pdf", plot=plot_generator(1000,3,0))

#Generates Figure 8 of the paper "The Online Closure Principle"
ggsave("Figure8.pdf", plot=plot_generator(100,3,0))
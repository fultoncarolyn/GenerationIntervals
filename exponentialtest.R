# Generate exponentially distributed random variables and compare CDF 

library(ggplot2)
library(tidyverse)
library(dplyr)
library(hrbrthemes)

expvals = rexp(1000,0.8)
expvals <- sort(expvals)
timevals <- seq(from = 0,to = 10,length.out = 1000)
df <- data.frame(expvals = as.numeric(expvals),
                 timevals = as.numeric(timevals))

 df$expcdf <- 1 - exp(-(0.8)*df$expvals)
 df$expcdft <- 1 - exp(-(0.8)*df$timevals)
 df$id <- ave(df$expvals, FUN=seq_along)
 #df$csum <- ave(df$id, FUN=cumsum)
 #df$normcsum <- df$csum/max(df$csum, na.rm=TRUE)
 df$normexpvals <- df$id/max(df$id, na.rm=TRUE)

print(df)
 
 fig <-df %>%
   ggplot(aes(x = expvals)) +
   geom_line(aes(y = normexpvals, color = "red")) +
   geom_line(aes(y=expcdf)) +
   geom_line(aes(y=expcdft), linetype = "dashed") +
   xlab("Time") +
   ylab("Cumulative Secondary Infections")
 
 print(fig)
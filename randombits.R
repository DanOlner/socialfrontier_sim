#Random bits
library(stringr)

#Example string with number at end we want
string = c("thingwithnnumberatend1","thingwithnnumberatend2","thingwithnnumberatend3")

#Pull out range of characters
subz = str_sub(string,nchar(string),nchar(string))
#Or pull out a single character at the chosen position. Here it's the last position (given by nchar, number of characters)
subz = str_sub(string,nchar(string))

#Then convert to numeric so can be used in a plot
plot(as.numeric(subz))


#rep will repeat numbers
rep(1:5)

#To get repeated values, just use "each" it turns out!
rep(1:5, each = 3)

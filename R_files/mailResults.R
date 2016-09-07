library(mailR)

send.results <- function(sendTo, countries){
sender <- "jessicalynnfessler@gmail.com"  # Replace with a valid address
recipients <- c(sendTo)  # Replace with one or more valid addresses
email <- send.mail(from = sender,
                   to = recipients,
                   subject="The countries you can choose from are...",
                   body = paste0(countries,"!"),
                   smtp = list(host.name = "smtp.gmail.com", port = 465, 
                               user.name = "jessicalynnfessler@gmail.com", 
                               passwd = "featherF!$H8", ssl = TRUE),
                   authenticate = TRUE,
                   send = TRUE)
}
## Not run: email$send() # execute to send email

guestList <- c("Haileybrown@uchicago.edu", "stamper@uchicago.edu", "acjoslin@gmail.com", "tiffany.marchell@gmail.com")

countries <-read.csv("~/Desktop/new_countries.csv", header = FALSE, sep = ',')
countries <- as.list(countries$V1)
selected_countries <- sample(countries, (length(guestList)*3), replace = FALSE)

count <- 1
for (guest in guestList){
  country1 <- selected_countries[[count]]
  count = count + 1
  country2 <- selected_countries[[count]]
  count = count + 1
  country3 <- selected_countries[[count]]
  count = count + 1
  country_list <- paste(country1, country2, country3, sep = ", ")
  print(country_list)
  send.results(guest, country_list)
}

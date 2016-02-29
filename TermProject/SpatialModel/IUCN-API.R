require(XML)
require(maptools)
require(jpeg)
input = "panthera-uncia"
h <- htmlParse(paste("http://api.iucnredlist.org/go/",
               input, sep = ""))
distr1 <- xpathSApply(h, '//ul[@class="countries"]', xmlValue)
distr2 <- unlist(strsplit(distr1, "\n"))
distr2[distr2 == "Russian Federation"] <- "Russia"
pop <-xpathSApply(h, '//div[@id="population"]/text()[preceding-sibling::br]', xmlValue)
status <- xpathSApply(h, '//div[@id="red_list_category_code"]', xmlValue)
data(wrld_simpl)
pdf("IUCN_map.pdf", width = 10, height = 10, pointsize = 20)
par(mar = c(3, 3, 1, 1))
plot(wrld_simpl, col = "grey98", xlim=c(-170, 170), axes = T)
plot(wrld_simpl[wrld_simpl$NAME %in% distr2,], col = "grey75", add = T)
text(0, 150, gsub("-", " ", toupper(input)), font = 3)
text(0, 130, paste("--Status: ", status, "--", sep = ""))
text(0, -98, "--Population--", cex = 0.5, font = 2) 
text(0, -140, paste(strwrap(pop, width = 30), collapse = "\n"), cex = 0.4)
# download image:
graphics.off()
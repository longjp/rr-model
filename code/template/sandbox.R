load("make_template.RData")


tem$dust


SDSS    u 
SDSS    g 
SDSS    r 
SDSS    i 
SDSS    z

macri_dust <- c(4.799,3.737,2.587,1.922,1.430)
names(macri_dust) <- c("u","g","r","i","z")
tem$dust

macri_dust <- macri_dust[order(names(macri_dust))]


plot(macri_dust,tem$dust)

#STUB for code to summarize results for users. Also consider output from plot.hisse.states
#
# all.poly <- ls() #list all output hisse objects
# all.poly <- all.poly[grepl("poly", all.poly)]
# all.poly <- all.poly[-which(all.poly=="poly")]
# all.poly <- all.poly[-which(all.poly=="poly.dat")]
# all.poly.list <- list()
# for (i in sequence(length(all.poly))) {
#   all.poly.list[[i]] <- eval(parse(text=all.poly[i]))
#   names(all.poly.list)[i] <- all.poly[i]
# }
# all.AICc <- unlist(lapply(all.poly.list, "[[", "AICc"))
# all.AICc <- all.AICc-min(all.AICc)
#

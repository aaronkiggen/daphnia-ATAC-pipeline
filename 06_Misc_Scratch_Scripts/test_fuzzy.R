library(dplyr)
library(stringr)

tf_names <- c("FOXM1", "FOXA1", "FoxA-a", "FoxM", "FoxI-a", "FOXK1", "Foxf1", "FoxL2")

normalize_tf <- function(name) {
  name <- toupper(name)
  name <- str_replace_all(name, "[-_]", "")
  name <- str_replace_all(name, "[A-Z]$", "") # remove trailing letter if after number? no.
  # simplify FOX family
  name <- str_replace(name, "^FOX[A-Z0-9]+", str_extract(name, "^FOX[A-Z]"))
  return(name)
}

print(data.frame(orig=tf_names, norm=normalize_tf(tf_names)))

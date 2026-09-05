roundTbl <- function(R, k, seuil = 0.05, diag = "  ~~~") {
  # Affiche une matrice arrondie avec marquage des valeurs significatives.
  # Le triangle supérieur est marqué si la p-valeur (seuil > 0) ou la valeur
  # elle-même (seuil < 0) dépasse le seuil. Le triangle inférieur est marqué
  # "_" si la valeur miroir dans le triangle supérieur est significative.
  txt <- capture.output(write.table(round(R, k), sep = "\t", quote = FALSE))
  txt[1] <- paste0("\t", txt[1])
  txt <- gsub("NA", diag, txt)

  marque <- ifelse(seuil > 0, "*", "=")
  ss <- sign(seuil)

  txt <- sapply(seq_along(txt), function(li) {
    line  <- txt[li]
    parts <- strsplit(line, "\t")[[1]]

    dbl_pos <- if (li > 1 && li <= length(parts)) li else integer(0)
    if (length(dbl_pos) > 0) parts[dbl_pos] <- diag
    for (i in seq_along(parts)) {
      if (parts[i] == diag) next
      val <- suppressWarnings(as.numeric(parts[i]))
      if (!is.na(val)) {
        if (length(dbl_pos) > 0 && i < min(dbl_pos)) {
          row_idx <- li - 1
          col_idx <- i - 1
          p_val <- R[col_idx, row_idx]
          if (!is.na(p_val) && ss * p_val < seuil)
            parts[i] <- paste0(parts[i], "_")
        }
        if (val >= 0) parts[i] <- paste0("\u00A0", parts[i])
        if (length(dbl_pos) > 0 && i > max(dbl_pos) && ss * val < seuil)
          parts[i] <- paste0(parts[i], marque)
      }
    }
    paste(parts, collapse = "\t")
  })

  txt[1] <- gsub("\t", paste0("\t", strrep("\u00A0", ceiling((k + 2) / 2))), txt[1])
  cat(txt, sep = "\n")
}

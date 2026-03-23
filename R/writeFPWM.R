#' Writes the FPWM in transfac or FPWMtransfac format
#'
#' This function saves the FPWM in transfac format.
#' @param FPWM FPWM object.
#' @param format [character] the "transfac" option will output a standard transfac matrix per binding partner. The "FPWMtransfac" prints a single matrix with all the binding partners in FPWM format.
#' @param fileName [character] name of the file where the FPWM is going to be written.
#' @return A plain text file with the FPWM in transfac or FPWMtransfac format.
#' @keywords writeFPWM
#' @export
#' @examples
#' fpwm <- createFPWM(mainTF = "CEBPB", partners = c("ATF4","ATF7","ATF3"), cell = "K562", forkPosition = 5)
#' writeFPWM(FPWM = fpwm, format = "transfac", fileName = "FPWM.tf" )
#' writeFPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "FPWM.FPWMtf" )

writeFPWM <- function(FPWM = NULL,
                       format = "transfac",
                       fileName = NULL) {

  if (is.null(FPWM)) {
    stop("Please provide an FPWM object: write.TRANSFAC( FPWM = yourFPWM, ...)")
  }

  if (format != "transfac" && format != "FPWMtransfac") {
    stop("formats supported: transfac or FPWMtransfac ")
  }

  # Get format from forked
  diff_rowSums <- length(unique(rowSums(FPWM@forked[, 2:5])))
  first_rowsum <- rowSums(FPWM@forked[, 2:5])[1]
  matrix_format <- "Matrix"
  if (diff_rowSums != 1) {
    matrix_format <- "Count matrix"
  } else {
    if (first_rowsum > 1) {
      matrix_format <- "Scale Count matrix"
    }
    if (first_rowsum == 1) {
      matrix_format <- "Probability matrix"
    }
  }

  if (format == "FPWMtransfac") {
    ynames <- paste(unlist(FPWM@id), collapse = "_&_")

    if (is.null(fileName)) {
      fileName <- paste0(FPWM@xid, "_+_", ynames, ".FPWMtf")
    }

    fileConn <- file(fileName)

    transfac_vector <- character()  # Initialize an empty character vector
    transfac_vector <- c(
      transfac_vector,
      paste0("AC ", FPWM@xid, "_+_", ynames),
      "XX",
      paste0("parentLogo : ", FPWM@xid),
      paste0("leafLogos : ", paste(unlist(FPWM@id), collapse = ",")),
      paste0("overlappingScore : ", paste(unlist(FPWM@score), collapse = ",")),
      paste0(
        "numberOfBasePairs : ", paste(unlist(FPWM@nSites), collapse = ",")
      ),
      paste0(
        "numberOfOverlappingPeaks : ",
        paste(unlist(FPWM@nPeaks), collapse = ",")
      ),
      paste0("forkPosition : ", FPWM@forkPosition),
      "XX",
      paste("P0", "A", "C", "G", "T", sep = "\t"),
      apply(FPWM@forked, 1, paste, collapse = "\t"),
      "XX",
      "CC transfacFormat: FPWMtransfac format from FPWM",
      paste0("CC matrixFormat: ", matrix_format),
      c("XX", "//")
    )

    writeLines(transfac_vector, fileConn)
    close(fileConn)
    message(paste("FPWM saved as FPWMtransfac file:", fileName))
  }


  if (format == "transfac") {
    if (is.null(fileName)) {
      ynames <- paste(unlist(FPWM@id), collapse = "_&_")
      fileName <- paste0(FPWM@xid, "_+_", ynames, ".tf")
    }

    fileConn <- file(fileName)

    FPWMPO <- FPWM@forked$PO
    from <- min(FPWMPO[duplicated(FPWMPO)])
    to <- max(FPWMPO)
    ix_table <- cbind(which(FPWMPO %in% from), which(FPWMPO %in% to))

    transfac_vector <- character()  # Initialize an empty character vector

    for (ix in seq_along(FPWM@id)) {
      cofactor_id <- paste0(FPWM@xid, "_+_", FPWM@id[ix])
      transfac_vector <- c(
        transfac_vector,
        paste0("AC ", cofactor_id),
        "XX",
        paste0("ID ", cofactor_id),
        "XX",
        paste0("DE ", cofactor_id, " ; from MethMotif"),
        paste("P0", "A", "C", "G", "T", sep = "\t"),
        apply(FPWM@forked[1:FPWM@forkPosition, ], 1, paste, collapse = "\t"),
        apply(
          FPWM@forked[ix_table[ix, 1]:ix_table[ix, 2], ],
          1, paste, collapse = "\t"
        ),
        c("XX", "CC program: forkedTF"),
        paste0("CC numberOfBasePairs: ", FPWM@nPeaks[ix]),
        paste0("CC numberOfOverlappingPeaks: ", FPWM@nSites[ix]),
        paste0("CC matrixFormat: ", matrix_format),
        c("XX", "//")
      )
    }

    writeLines(transfac_vector, fileConn)
    close(fileConn)
    message(paste("FPWM saved as transfac file:", fileName))
  }
}

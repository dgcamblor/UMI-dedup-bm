#' Simple class to work with genetic variants in R.
#' 
#' @slot CHR The chromosome of the variant.
#' @slot POSITION The position of the variant.
#' @slot REF The reference allele.
#' @slot ALT The alternative allele.
#' @slot VAF The variant allele frequency.
#' @slot TYPE The type of variant, dynamically calculated.
#' 
#' @examples
#' # Create a new GeneticVariant object
#' variant <- new("GeneticVariant", CHR = "chr1", POSITION = 12345, REF = "A", ALT = "T", VAF = 0.5)

setClass(
  "GeneticVariant",
  representation(
    CHR = "character",
    POSITION = "numeric",
    REF = "character",
    ALT = "character",
    VAF = "numeric",
    TYPE = "character"
  )
)

get_variant_type <- function(ref, alt) {
  if (nchar(ref) == 1 & nchar(alt) == 1) {
    return("SNV")
  } else if (nchar(ref) > 1 & nchar(alt) > 1) {
    return("MNV")
  } else if (nchar(ref) == 1 & nchar(alt) > 1) {
    return("DEL")
  } else if (nchar(ref) > 1 & nchar(alt) == 1) {
    return("INS")
  } else {
    return("COMPLEX")
  }
}

setMethod(
  "initialize",
  "GeneticVariant",
  function(.Object, CHR, POSITION, REF, ALT, VAF) {
    .Object@CHR <- CHR
    .Object@POSITION <- POSITION
    .Object@REF <- REF
    .Object@ALT <- ALT
    .Object@VAF <- VAF
    .Object@TYPE <- get_variant_type(REF, ALT)
    return(.Object)
  }
)

# Print the variant
setMethod(
  "show",
  "GeneticVariant",
  function(object) {
    cat(paste0(object@CHR, ":", object@POSITION, object@REF, ">", object@ALT, "(", object@VAF, ")"))
  }
)

# Establish equality between two variants
setMethod(
  "==",
  signature(e1 = "GeneticVariant", e2 = "GeneticVariant"),
  function(e1, e2) {
    return(e1@CHR == e2@CHR & e1@POSITION == e2@POSITION & e1@REF == e2@REF & e1@ALT == e2@ALT)
  }
)

# Getters ----------------------------------------------------------------------
setGeneric("get_variant", function(object) {standardGeneric("get_variant")})
setGeneric("get_CHR", function(object) {standardGeneric("get_CHR")})
setGeneric("get_POSITION", function(object) {standardGeneric("get_POSITION")})
setGeneric("get_REF", function(object) {standardGeneric("get_REF")})
setGeneric("get_ALT", function(object) {standardGeneric("get_ALT")})
setGeneric("get_VAF", function(object) {standardGeneric("get_VAF")})
setGeneric("get_TYPE", function(object) {standardGeneric("get_TYPE")})

setMethod(
  "get_variant",
  "GeneticVariant",
  function(object) {
    return(paste0(object@CHR, ":", object@POSITION, object@REF, ">", object@ALT))
  }
)

setMethod(
  "get_CHR",
  "GeneticVariant",
  function(object) {
    return(object@CHR)
  }
)

setMethod(
  "get_POSITION",
  "GeneticVariant",
  function(object) {
    return(object@POSITION)
  }
)

setMethod(
  "get_REF",
  "GeneticVariant",
  function(object) {
    return(object@REF)
  }
)

setMethod(
  "get_ALT",
  "GeneticVariant",
  function(object) {
    return(object@ALT)
  }
)

setMethod(
  "get_VAF",
  "GeneticVariant",
  function(object) {
    return(object@VAF)
  }
)

setMethod(
  "get_TYPE",
  "GeneticVariant",
  function(object) {
    return(object@TYPE)
  }
)
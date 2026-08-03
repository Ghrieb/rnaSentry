#' Lock or unlock an rnaSentry signature
#'
#' Toggles the reproducibility guardrail on an \code{rnaSentry_signature}
#' object. Locking records a timestamp and an audit entry in the signature's
#' \code{flags} table; downstream stages that must not be rerun after
#' publication can then be protected. The function is copy-on-write: the
#' input object is never modified and the updated signature is returned.
#'
#' @param sig An object of class \code{"rnaSentry_signature"}, as returned by
#'   \code{\link{build_signature}}.
#' @param lock Logical. \code{TRUE} (default) locks the signature,
#'   \code{FALSE} unlocks it.
#'
#' @return A new \code{rnaSentry_signature} object with updated
#'   \code{locked}, \code{lock_time} and \code{flags} fields.
#'
#' @examples
#' sig <- structure(list(genes = "gene1", outcome = "overall_survival",
#'                       locked = FALSE),
#'                  class = "rnaSentry_signature")
#' locked <- lock_signature(sig)
#' locked$locked
#' locked$lock_time
#' unlocked <- lock_signature(locked, lock = FALSE)
#' unlocked$locked
#'
#' @export
lock_signature <- function(sig, lock = TRUE) {
  if (!inherits(sig, "rnaSentry_signature")) {
    stop("'sig' must be an rnaSentry_signature object.", call. = FALSE)
  }
  if (!is.logical(lock) || length(lock) != 1 || is.na(lock)) {
    stop("'lock' must be a single TRUE or FALSE.", call. = FALSE)
  }

  out <- sig
  fl <- out$flags
  if (is.null(fl) || !is.data.frame(fl)) {
    fl <- data.frame(check = character(0), severity = character(0),
                     detail = character(0), stringsAsFactors = FALSE)
  }
  add_flag <- function(flags, check, severity, detail) {
    rbind(flags, data.frame(check = check, severity = severity,
                             detail = detail, stringsAsFactors = FALSE))
  }

  if (lock) {
    if (isTRUE(out$locked)) {
      stop("The signature is already locked.", call. = FALSE)
    }
    out$locked <- TRUE
    out$lock_time <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                            " UTC")
    out$flags <- add_flag(fl, "signature_locked", "info",
                          "Signature locked against downstream mutation.")
  } else {
    if (!isTRUE(out$locked)) {
      stop("The signature is not locked.", call. = FALSE)
    }
    out$locked <- FALSE
    out$lock_time <- NULL
    out$flags <- add_flag(fl, "signature_unlocked", "info",
                          "Signature unlocked; downstream stages may be rerun.")
  }
  out
}

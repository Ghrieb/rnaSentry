#' Lock or unlock an rnaSentry signature
#'
#' Toggles the reproducibility guardrail on an \code{rnaSentry_signature}
#' object. Locking records a timestamp and an audit entry in the signature's
#' \code{flags} table; downstream stages that must not be rerun after
#' publication can then be protected. The function is copy-on-write: the
#' input object is never modified and the updated signature is returned.
#'
#' @section Enforceable lock:
#' Locking stores a gene-set fingerprint in a session-level environment.
#' While this fingerprint is set, \code{build_signature()} will refuse to
#' run, preventing silent re-selection of genes after survival analysis.
#' The lock can be reversed by calling \code{lock_signature(sig, lock = FALSE)},
#' but doing so after survival stages have been run is not recommended.
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
    fl <- .new_flags("lock_signature")
  }

  if (lock) {
    if (isTRUE(out$locked)) {
      stop("The signature is already locked.", call. = FALSE)
    }
    out$locked <- TRUE
    out$lock_time <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                            " UTC")
    # Store gene-set fingerprint in session environment for lock enforcement
    locked_env <- get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
    gene_hash <- paste(sort(out$genes), collapse = "|")
    locked_env[[gene_hash]] <- TRUE
    out$flags <- .add_flag(fl, "signature_locked", "info",
                           "Signature locked against downstream mutation.",
                           stage = "lock_signature")
  } else {
    if (!isTRUE(out$locked)) {
      stop("The signature is not locked.", call. = FALSE)
    }
    out$locked <- FALSE
    out$lock_time <- NULL
    # Remove fingerprint from session environment
    locked_env <- get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
    if (length(ls(locked_env)) > 0) {
      rm(list = ls(locked_env), envir = locked_env)
    }
    out$flags <- .add_flag(fl, "signature_unlocked", "info",
                           "Signature unlocked; downstream stages may be rerun.",
                           stage = "lock_signature")
  }
  out
}

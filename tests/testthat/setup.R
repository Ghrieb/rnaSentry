# Clear the lock-enforcement environment before each test file so that locks
# set by one test file do not leak into subsequent files.
if (exists(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))) {
  locked_env <- get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
  if (length(ls(locked_env)) > 0) {
    rm(list = ls(locked_env), envir = locked_env)
  }
}

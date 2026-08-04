library(SummarizedExperiment)

# Clear lock-enforcement environment. Fixture-based lock tests and successful
# run_rnaSentry() calls leave fingerprints in the session environment, so any
# test that asserts on env state (or re-runs the pipeline) must start from a
# clean slate.
.clear_locks <- function() {
  if (exists(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))) {
    locked_env <- get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
    if (length(ls(locked_env)) > 0) {
      rm(list = ls(locked_env), envir = locked_env)
    }
  }
}
.clear_locks()

test_that("lock_signature rejects non-signature input", {
  expect_error(lock_signature(list(genes = "x")), "rnaSentry_signature")
  expect_error(lock_signature(1), "rnaSentry_signature")
})

test_that("lock_signature validates the lock argument", {
  sig <- make_signature_for_testing()
  expect_error(lock_signature(sig, lock = NA), "single TRUE or FALSE")
  expect_error(lock_signature(sig, lock = "yes"), "single TRUE or FALSE")
  expect_error(lock_signature(sig, lock = c(TRUE, FALSE)),
               "single TRUE or FALSE")
})

test_that("locking a signature sets state and records an audit entry", {
  sig <- make_signature_for_testing()
  locked <- lock_signature(sig)
  expect_true(isTRUE(locked$locked))
  expect_true(nzchar(locked$lock_time))
  expect_true(any(locked$flags$check == "signature_locked"))
  expect_false(isTRUE(sig$locked))
  expect_null(sig$lock_time)
  # Cleanup: release the fingerprint this test created.
  lock_signature(locked, lock = FALSE)
})

test_that("locking an already locked signature errors", {
  sig <- make_signature_for_testing(locked = TRUE)
  expect_error(lock_signature(sig), "already locked")
})

test_that("unlocking reverses the state and records an audit entry", {
  sig <- make_signature_for_testing(locked = TRUE)
  unlocked <- lock_signature(sig, lock = FALSE)
  expect_false(isTRUE(unlocked$locked))
  expect_null(unlocked$lock_time)
  expect_true(any(unlocked$flags$check == "signature_unlocked"))
})

test_that("unlocking an unlocked signature errors", {
  sig <- make_signature_for_testing()
  expect_error(lock_signature(sig, lock = FALSE), "not locked")
})

test_that("lock and unlock round-trip", {
  sig <- make_signature_for_testing()
  l1 <- lock_signature(sig)
  expect_error(lock_signature(l1), "already locked")
  u1 <- lock_signature(l1, lock = FALSE)
  expect_false(isTRUE(u1$locked))
  l2 <- lock_signature(u1)
  expect_true(isTRUE(l2$locked))
  expect_true(all(c("signature_locked", "signature_unlocked",
                    "signature_locked") %in% l2$flags$check))
  # Cleanup: release the fingerprint this test created.
  lock_signature(l2, lock = FALSE)
})

test_that("flags carry the stage schema for the audit ledger", {
  sig <- make_signature_for_testing()
  locked <- lock_signature(sig)
  expect_true(all(c("check", "severity", "detail", "stage") %in%
                    colnames(locked$flags)))
  expect_equal(locked$flags$stage[locked$flags$check == "signature_locked"],
               "lock_signature")
  # Cleanup: release the fingerprint this test created.
  lock_signature(locked, lock = FALSE)
})

test_that("print shows the lock state", {
  sig <- make_signature_for_testing()
  expect_output(print(sig), "Not locked")
  locked <- lock_signature(sig)
  expect_output(print(locked), "Locked")
  # Cleanup: release the fingerprint this test created.
  lock_signature(locked, lock = FALSE)
})

test_that("lock enforcement prevents run_rnaSentry when locked", {
  .clear_locks()
  # Simulate the session-environment enforcement: lock a signature, then
  # verify that run_rnaSentry() refuses to re-run until unlocked.
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 3, repeats = 1,
                         folds = 2, seed = 1)
  locked <- lock_signature(sig)
  expect_true(length(ls(
    get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
  )) > 0)
  expect_error(run_rnaSentry(se, time_col = "time", event_col = "event",
                              top_n = 3, repeats = 1, folds = 2, seed = 1,
                              render_report = FALSE),
               "locked signature already exists")
  # Unlock and verify pipeline works again
  unlocked <- lock_signature(locked, lock = FALSE)
  expect_equal(length(ls(
    get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
  )), 0)
  result <- run_rnaSentry(se, time_col = "time", event_col = "event",
                           top_n = 3, repeats = 1, folds = 2, seed = 1,
                           render_report = FALSE)
  expect_true("stages" %in% names(result))
})

test_that("unlocking one signature leaves other locked signatures protected", {
  .clear_locks()
  locked_env <- function() get(".rnaSentry_locked_sigs",
                               envir = asNamespace("rnaSentry"))
  se <- make_survival_se()
  locked_a <- lock_signature(make_signature_for_testing(n_genes = 20))
  locked_b <- lock_signature(make_signature_for_testing(n_genes = 21))
  expect_equal(length(ls(locked_env())), 2)
  # Unlocking B must not release A's fingerprint: the guardrail stays up.
  lock_signature(locked_b, lock = FALSE)
  expect_equal(length(ls(locked_env())), 1)
  expect_error(run_rnaSentry(se, time_col = "time", event_col = "event",
                              top_n = 3, repeats = 1, folds = 2, seed = 1,
                              render_report = FALSE),
               "locked signature already exists")
  # Cleanup: unlock A so later tests start with an empty lock environment.
  lock_signature(locked_a, lock = FALSE)
  expect_equal(length(ls(locked_env())), 0)
})

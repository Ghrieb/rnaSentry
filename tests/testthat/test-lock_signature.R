library(SummarizedExperiment)

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
})

test_that("flags carry the stage schema for the audit ledger", {
  sig <- make_signature_for_testing()
  locked <- lock_signature(sig)
  expect_true(all(c("check", "severity", "detail", "stage") %in%
                    colnames(locked$flags)))
  expect_equal(locked$flags$stage[locked$flags$check == "signature_locked"],
               "lock_signature")
})

test_that("print shows the lock state", {
  sig <- make_signature_for_testing()
  expect_output(print(sig), "Not locked")
  expect_output(print(lock_signature(sig)), "Locked")
})

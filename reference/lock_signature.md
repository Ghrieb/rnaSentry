# Lock or unlock an rnaSentry signature

Toggles the reproducibility guardrail on an `rnaSentry_signature`
object. Locking records a timestamp and an audit entry in the
signature's `flags` table; downstream stages that must not be rerun
after publication can then be protected. The function is copy-on-write:
the input object is never modified and the updated signature is
returned.

## Usage

``` r
lock_signature(sig, lock = TRUE)
```

## Arguments

- sig:

  An object of class `"rnaSentry_signature"`, as returned by
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md).

- lock:

  Logical. `TRUE` (default) locks the signature, `FALSE` unlocks it.

## Value

A new `rnaSentry_signature` object with updated `locked`, `lock_time`
and `flags` fields.

## Enforceable lock

Locking stores a gene-set fingerprint in a session-level environment.
While this fingerprint is set,
[`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
will refuse to run, preventing silent re-selection of genes after
survival analysis. The lock can be reversed by calling
`lock_signature(sig, lock = FALSE)`, but doing so after survival stages
have been run is not recommended.

## Examples

``` r
sig <- structure(list(genes = "gene1", outcome = "overall_survival",
                      locked = FALSE),
                 class = "rnaSentry_signature")
locked <- lock_signature(sig)
locked$locked
#> [1] TRUE
locked$lock_time
#> [1] "2026-08-05 17:02:25 UTC"
unlocked <- lock_signature(locked, lock = FALSE)
unlocked$locked
#> [1] FALSE
```

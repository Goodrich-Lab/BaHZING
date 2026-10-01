## R CMD check results

0 errors | 0 warnings | 2 notes

## Resubmission

This is a resubmission of BaHZING, which was archived on 2026-03-25.

The archived version produced the following NOTE:

  Format_BaHZING: no visible binding for global variable 'id'

This has been corrected by declaring the non-standard evaluation variables
with utils::globalVariables(). The current R CMD check reports no problems
in R code analysis.

## Test environments

* Local: macOS Sequoia 15.7.4, R 4.5.2

## Notes

* The first NOTE is the expected CRAN incoming note that this is a new
  submission of a package that was previously archived, as explained above.
* The second NOTE was emitted while checking for future file timestamps:
  `unable to verify current time`. The local check environment could not
  perform external time verification; no package files were identified as
  having future timestamps.
* JAGS 4.x.y is required, as declared in SystemRequirements.

# Base R only. Invented source-normalization rule, no assay input or helper import.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L || file.exists(args[1])) stop('Expected fresh synthetic receipt path')
Sys.setlocale('LC_COLLATE', 'C')
ids <- c('toy:1:A:C', 'toy:1:A:G', 'toy:1:AA:CC')
stopifnot(identical(order(ids), 1:3))
dp <- c(100, 100, 90)
cao <- c(10, 20, 5)
depth_candidates <- unique(dp)
wt_candidates <- depth_candidates - sum(cao)
wt_keys <- rep('S:1', length(depth_candidates))
assigned <- wt_candidates[match(rep('S:1', 3), wt_keys)]
stopifnot(identical(depth_candidates, c(100,90)), identical(wt_candidates,c(65,55)))
stopifnot(identical(assigned,c(65,65,65)))
# First match is not intrinsically maximum depth.
dp2 <- c(80, 80, 100)
wt2 <- unique(dp2)-sum(cao)
assigned2 <- wt2[match(rep('S:1',3),rep('S:1',length(wt2)))]
stopifnot(identical(wt2,c(45,65)),identical(assigned2,c(45,45,45)))
# Row order at the denominator-table stage is mathematically consequential.
reversed_assignment <- rev(wt_candidates)[match(rep('S:1',3),wt_keys)]
stopifnot(identical(reversed_assignment,c(55,55,55)))
writeLines(c('PASS_BASE_R_SYNTHETIC_FIRST_MATCH',
  'invented_DP=100,100,90; invented_CAO=10,20,5',
  'unique_DP=100,90; candidate_WT=65,55; first_match_WT=65,65,65',
  'counterexample_DP=80,80,100; first_match_WT=45,45,45; max_depth_WT=65',
  'reversing_denominator_rows_selects_WT55; source_sort_order_must_be_retained',
  'No actual CAO, CAF, WT, q or downstream outcome was read'),args[1])

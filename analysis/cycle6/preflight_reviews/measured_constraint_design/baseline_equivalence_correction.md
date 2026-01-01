# Explicit correction: the two additive spans coincide on this FAS domain

This addendum preserves the original report, protocol review and receipt. My earlier message to root calling nucleotide additivity inadequate relative to arbitrary codon effects **on this specific72-domain** was overgeneralized. My mathematical PASS of draft SHA256 `731620b9d8d4057053b80143b2c476b94f8c3feb46c6342c9841ff531ad3393d` also missed its incorrect description of nucleotide additivity as a weaker sensitivity model. That sentence needs correction before freeze. The interval/order theorem, conservative graph bounds and codon-gate counterexample are unchanged.

## Exact span proof

The authenticated factor description supplied by the source lane is two three-state Leu paths and three binary factors. For the first Leu path CTT,CTG,TTG, let a indicate first-base T and b indicate third-base G. Its three observed states are00,01,11. The two nonreference codon indicators satisfy

`1[state1]=b-a`, `1[state2]=a`.

Conversely, a=1[state2] and b=1[state1]+1[state2]. Thus an intercept and the two varying nucleotide indicators represent every possible three-state codon score, and vice versa. The missing fourth corner does not supply an extra constraint on this domain.

For the second path TTG,CTG,CTC, let c indicate first-base T and d indicate third-base G. Its observed states are11,01,00. Then

`1[state1]=d-c`, `1[state2]=1-d`,

with inverse c=1-1[state1]-1[state2] and d=1-1[state2]. Again the spans coincide. Each remaining binary factor has the same one free contrast under either parameterization. Summing over factors gives identical eight-dimensional linear score spaces including the intercept. It follows that allowing any common nondecreasing h also gives identical null classes. Equal ranks alone would not prove equality of two arbitrary spans; these explicit identities do.

The gate-valley example makes the distinction particularly clear. Scores(1,0,1) on the first Leu path equal the nucleotide-additive field1+a-b. Every path inside P3 crosses its low middle state. On the *complete* binary nucleotide square, the unused corner10 would have score2 and supply the compensating route. That corner is not a vertex of the declared Leu path. Hence what fails is transferring a complete-binary-cube connectivity theorem to a coding/design graph with a missing corner, not equality of the two score-model spans here.

## Required consequence

Treat nucleotide-versus-codon parameterizations as an exact implementation equivalence check on this72-domain, not two competing biological models or evidence of extra codon interactions. A disagreement in forced-order feasibility or fitted score span is an implementation/domain error. This correction does not make full72 disconnection diagnostic of interactions: both equivalent additive classes can already show the gate valley. It also does not exclude an unmeasured outside-design synonymous bypass.

No observed outcome was read or fitted. The original files and receipt remain unchanged. A revised protocol can be approved only after this sentence is corrected; production-code verification remains separate.

## Inspected final draft delta

Root's draft SHA256 `eb6e76ed299f27227f68c541b9a0a3af2aee1792ca3efac239e1a46fe211eb51` now explicitly treats the encodings as an exact implementation equivalence check, retains the P3 gate valley and distinguishes the 64 raw-design linear contrasts from statistical residual degrees of freedom for an unknown-monotone-response model. It retains exact primal/dual verification, the explicit mean of nine per-output relative ES values and the open/closed obstruction interval. The mathematical portions of this exact corrected draft pass the bounded review. There is no verdict here on the count reader, production solver or actual observations.

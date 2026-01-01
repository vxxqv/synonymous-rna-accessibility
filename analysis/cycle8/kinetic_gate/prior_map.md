# Bounded primary-source and overlap audit

Checked 25 October 2026. This is an exact two-source comparison requested by the root, not an exhaustive novelty search. Mathbox proof-audit/literature-check conventions were used to separate model algebra, source authentication and assay applicability. The new derivation remains self-reviewed pending a fresh challenge.

## Müller et al. 2024

Müller, J. M., Moos, K., Baar, T., Maier, K. C., Zumer, K., and Tresch, A. (2024). *Nuclear export is a limiting factor in eukaryotic mRNA metabolism*. PLOS Computational Biology, 20(5), e1012059. [Published article](https://doi.org/10.1371/journal.pcbi.1012059).

Authentication: official PLOS published XML, DOI/title/authors checked; SHA256 72a7bf9eb523fc7d63dcbf40ab8f3566f68c630d15d09b964868d7b016bbc97d. Section “A two-compartment model of RNA metabolism,” equations 1 through 4, supplies the RNA subsystem and distinguishes nuclear removal from export. Notation dictionary: mu=s, tau=e, nu=d_N, lambda=d_C. Its fraction-resolved labeling trajectories supply information absent from an equilibrium snapshot. Sections on compartment abundance and metabolic-parameter estimation explicitly treat abundance scaling and combined rates. Their kinetic estimates do not validate COMT proxies or identify our protein module. The RNA equilibrium algebra is known after translation; its full positive fiber is a formal algebraic corollary, not a new kinetic principle.

Application status: verified for the model translation and the distinction between kinetic trajectories and one equilibrium. Inapplicable as evidence that COMT has calibrated compartments, negligible d_N, or known rates. The paper's biological findings are not assumed to transport to this reporter.

## Cenik et al. 2015

Cenik, C., et al. (2015). *Integrative analysis of RNA, translation, and protein levels reveals distinct regulatory variation across humans*. Genome Research, 25(11), 1610-1621. [Published article](https://doi.org/10.1101/gr.193342.115); [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4617958/).

Authentication: final primary XML obtained from Europe PMC, DOI/title/volume/pages checked; SHA256 70611934c91cca147478b70208cf16919f7cf068e6b19c982fb66c1f69c2ebf5. Introduction distinguishes RNA expression, translation and protein stability. Sections “Integrative analysis of RNA expression, ribosome occupancy, and protein levels” and “Gene expression variability between individuals” address integrated measurements and explicitly retain a variable untranslated RNA pool as an alternative to buffering. This already defeats any claim that cross-layer discordance or an untranslated-pool caveat is new here. It does not prove this package's exact six-parameter fiber or authenticate the COMT assay map. Application status: verified adjacent precedent; prediction/association is not identification of rates.

## Increment and limits

The positive fiber, clock quotient, inverse-information examples and proxy counterexamples provide an inspectable **assay-specific identification contract**. They add no measured rate, structure mechanism or broad priority claim. The protein extension and algebraic elimination are elementary consequences of a standard linear production/loss model. No publication headline should be based on calling this a new theorem.

The sources were checked directly in their primary versions, not inferred from citing-paper snippets. No external theorem is needed for the balance-equation proof. Matrix stability uses the elementary triangular eigenvalue criterion. Plugin attempts: Scite search_literature returned an access-plan error; Amass search returned a usage-limit error. Those unsuccessful calls supplied no scientific evidence. Primary web access to Cenik's publisher returned 403 and PMC displayed a browser check; the official Europe PMC fullTextXML route succeeded without bypassing a protection. PLOS official XML also succeeded. No further plugin retries were made. No exhaustive editorial-notice/correction search is claimed; the retrieved published primary versions and their identifiers are pinned, and no correction was imported into the argument.

## Dependency graph

    positive constant model + exact steady-state amounts
        -> three balance equations
        -> positive inverse parameterization and full fiber
        -> nonidentification, clock quotient and exact consistency conditions

    local COMT source contract + published measurement definitions
        -> no authenticated observation map to n,c,p or kinetic flux
        -> no physical kinetic conclusion from a proxy fit

Neither source authentication nor finite synthetic agreement proves that COMT satisfies the model assumptions. The first branch is mathematical; the second is a source-limited applicability decision.

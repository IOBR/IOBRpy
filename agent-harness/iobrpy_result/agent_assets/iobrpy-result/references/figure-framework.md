# Figure Brief And Design Framework

Complete this brief in working notes before plotting. Keep it concise and revise it after
the first data inspection.

Use this framework only after the user explicitly requests creation or modification of a
visualization. Do not use it for interpretation-only requests.

```text
User question:
Result files:
Provenance ledger:
  confirmed native:
  likely or reusable:
  external:
  metadata:
  excluded or unresolved:
Selected backend: Python or R
Result profile:
  method and scale:
  rows / columns:
  sample key:
  metadata coverage:
  fit / uncertainty fields:
Provisional conclusion:
Comparison unit:
Evidence hierarchy:
  hero:
  context:
  robustness:
Main reviewer risk:
Candidate design A:
Candidate design B:
Chosen design and why:
Feature-selection rule:
Ordering rule:
Statistics:
Final dimensions and formats:
```

Do not complete this brief for a figure request until the user has explicitly selected the
backend. The selected backend is exclusive for plotting, previewing, exporting, and visual
QA.

## Design Selection

Compare at least two materially different designs. A stacked bar with a second stacked bar
is not a second design. Consider whether the question is best answered by:

- sample composition versus group effect;
- global structure versus selected biology;
- one method in depth versus cross-method concordance;
- point estimates alone versus estimates paired with fit or uncertainty;
- one integrated figure versus several independent figures.

Choose the design that minimizes likely misinterpretation. A more complex layout is useful
only when each panel contributes non-redundant evidence.

## Result Profiles

Use a short result profile to stop visual design from outrunning semantics:

```text
result class:
measurement unit:
bounded or unbounded:
compositional:
sample orientation:
quality fields:
known transformation:
valid comparison:
invalid comparison:
```

Examples:

- CIBERSORT relative output: bounded, compositional, sample-level fit available, valid
  within-cell-type cohort comparison, invalid cross-method magnitude comparison.
- MCPcounter: unbounded method-specific abundance score, non-compositional, valid relative
  within-feature comparison, invalid percentage interpretation.
- LR score: minimum of paired bulk-expression signals after transformation, valid ranking
  of co-expression-compatible pairs, invalid proof of cell-cell interaction.

## Evidence Integration

When several outputs are available, integrate them only when they answer one shared
question. Useful integration patterns include:

- immune fraction plus immune signature plus ESTIMATE context;
- BayesPrism theta plus theta CV;
- cluster membership plus defining features plus phenotype metadata;
- TRUST4 diversity plus read depth plus top-clone concentration;
- HLA calls plus cohort allele frequency, without inferring antigen presentation.

Do not create a dashboard of every available result. Breadth without a claim is not a
scientific figure.

When a directory mixes IOBRpy and non-IOBRpy files, preserve provenance in every integrated
artifact. Panel titles, legends, captions, and exported source data must distinguish native
IOBRpy, compatible reusable, and external results. Unknown files stay excluded.

## Interpretation Prompts

After rendering, ask:

1. What exact visual pattern is robust to reasonable ordering or transformation choices?
2. Which parts come directly from data and which come from the method model?
3. What alternative technical or biological explanation remains?
4. Does fit, uncertainty, depth, or missingness weaken any highlighted sample?
5. What new evidence would most efficiently distinguish the leading explanations?

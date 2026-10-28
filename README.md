# Remitly 409A — Common Stock Valuation Re-performance (OPM, PWERM, Hybrid, DLOM)

An independent re-performance of three pre-IPO 409A valuations of **Remitly Global, Inc.** (NASDAQ: RELY), built entirely from public SEC filings: the S-1, the charter, the 10-Q/10-K and Remitly's "cheap stock" response letter to the SEC.

## Headline
| Valuation date | Method | Model | Reported by Remitly |
|---|---|---|---|
| 15-Aug-2019 | OPM backsolve to Series E ($5.9566) | $2.52 | $2.51 (+0.3%) |
| 05-Aug-2020 | OPM backsolve to Series F ($9.1456) | $3.57 | $3.54 (+0.7%) |
| 31-Aug-2021 | Hybrid PWERM / OPM | implied IPO value **$34.5–$36.1** | $28.73 vs $38–42 range, $43 IPO |

Both OPM backsolves reproduce the reported values within 1% using a single volatility (54%, from Remitly's own disclosures) and the company's DLOMs. The 2021 hybrid answers the cheap-stock question: $28.73 is consistent with an expected IPO value about 10–14% below the range midpoint set ten days later.

## Method
- **Waterfall:** eight pari passu series of 1x non-participating preferred stock ($399.8M preference), common stock, options, warrants and RSUs. Breakpoints are computed from the charter terms.
- **OPM:** each tranche between breakpoints is a Black-Scholes call spread on total equity. Total equity is backsolved so the latest preferred round prices at its issue price.
- **DLOM:** company DLOMs are tested against the Finnerty (2012) average-strike put and the Chaffe (1993) protective put.
- **Hybrid (2021):** Remitly's disclosed scenario weights (77.5% / 15% / 7.5%) and DLOMs. The expected IPO value is solved to reproduce $28.73. A second calibration also matches the disclosed 14.2% blended DLOM.
- **Bridge:** the $40 midpoint less probability and timing (−$0.94), less DLOM (−$5.75), less the lower expected IPO value (−$4.57), gives **$28.73**.

## Files
| File | What it is |
|---|---|
| `remitly_409a_model.xlsx` | Formula-driven Excel model: inputs, breakpoint OPM for each date, hybrid (Goal Seek on the IPO value), DLOM, sensitivities, results, issues log, sources |
| `remitly_409a_validation.ipynb` | Python validation: cap table checks, waterfall plots, backsolves, (T, σ) frontier, DLOM cross-check, hybrid solve, bridge, two-condition calibration, value path |
| `val409a.py` | Valuation engine (waterfall, OPM, backsolve, DLOM models, hybrid), used by the notebook |
| `remitly_409a_valuation_memo.pdf` | Client-style re-performance memo with a mock auditor Q&A |

## Run
```bash
pip install numpy scipy pandas matplotlib jupyter
jupyter notebook remitly_409a_validation.ipynb
```

## Key findings
- The OPM backsolve is robust: a ±10% band covers most plausible (T, σ) pairs, and the estimated 2019 option strike moves value by only about ±$0.02.
- Company DLOMs (30% and 25%) sit close to the midpoint of the Finnerty–Chaffe range. DLOM is the largest judgement after the backsolve itself.
- Across reasonable discount rates and IPO timing, the implied IPO value stays between $34.1 and $35.6. Matching the disclosed blended DLOM moves it to $36.1.
- Open item: the underwriters' August 2021 valuation indications (not public) are the document that would close the cheap-stock question.

*Prepared by Parnika Sabharwal. Educational project based on public information; not a valuation opinion and not affiliated with Remitly.*

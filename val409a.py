"""409A common stock valuation engine — Remitly Global, Inc. (pre-IPO).
OPM (Black-Scholes breakpoint method) with backsolve to a preferred round, PWERM, hybrid,
DLOM via Finnerty (average-strike put) and Chaffe (protective put). Shares and $ in units, per-share in $."""
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

# ---------- capital structure (Remitly 10-Q Q3-2021 Note 9; charter Ex. 3.1; 10-K FY2021) ----------
SERIES = ["Seed","Seed Prime","A","B","C","D","E","F"]
OIP    = np.array([0.2729, 0.2961, 0.49764, 0.8805, 1.7032, 3.7914, 5.9566, 9.1456])
SH_2020 = np.array([6446322, 8643665, 10359546, 14196476, 25146777, 30331802, 22663933, 9294084], float)

def cap_table(date):
    if date == "2019-08-15":   # before Series F
        pref = SH_2020.copy(); pref[7] = 0
        return dict(pref=pref, common=22425112., options=19045751., opt_k=1.50, warrants=256250., war_k=0.42, rsu=0.)
    if date == "2020-08-05":
        return dict(pref=SH_2020.copy(), common=24289906., options=21034424., opt_k=1.88, warrants=256250., war_k=0.42, rsu=0.)
    if date == "2021-08-31":   # June 30, 2021 counts; extra 328,026 preferred assigned to Seed (immaterial)
        pref = SH_2020.copy(); pref[0] += 127410631 - SH_2020.sum()
        return dict(pref=pref, common=26385643., options=25355906., opt_k=3.13, warrants=256250., war_k=0.42, rsu=617696.)
    raise ValueError(date)

# ---------- waterfall ----------
def _participants(ct, p):
    """value V as a function of the per-share value p received by common-equivalent holders."""
    conv = OIP < p                                   # non-participating: convert when as-converted > preference
    V = np.sum(OIP * ct["pref"] * ~conv) + p * (ct["common"] + ct["rsu"] + np.sum(ct["pref"] * conv))
    V += max(p - ct["war_k"], 0) * ct["warrants"] + max(p - ct["opt_k"], 0) * ct["options"]
    return V

def breakpoints(ct):
    """kinks of the payoff in V: total preference, then each conversion / exercise threshold."""
    LP = np.sum(OIP * ct["pref"])
    ks = sorted(set([0.0] + list(OIP[ct["pref"] > 0]) + [ct["war_k"], ct["opt_k"]]))
    bps = [LP] + [_participants(ct, k) for k in ks if k > 0]
    return np.array(sorted(set(bps))), LP

def payoff(ct, V):
    """payoff to each class at total exit value V: dict of arrays (per class totals)."""
    LP = np.sum(OIP * ct["pref"])
    out = {}
    if V <= LP:
        for i, s in enumerate(SERIES): out[s] = V * OIP[i] * ct["pref"][i] / LP
        out["common"] = 0.0; out["options"] = 0.0; out["warrants"] = 0.0; out["rsu"] = 0.0
        return out
    p = brentq(lambda x: _participants(ct, x) - V, 0, max(V, 1) / 1.0)
    conv = OIP < p
    for i, s in enumerate(SERIES): out[s] = (p if conv[i] else OIP[i]) * ct["pref"][i]
    out["common"] = p * ct["common"]; out["rsu"] = p * ct["rsu"]
    out["warrants"] = max(p - ct["war_k"], 0) * ct["warrants"]; out["options"] = max(p - ct["opt_k"], 0) * ct["options"]
    return out

def bs_call(S, K, T, r, sig):
    if K <= 0: return S
    d1 = (np.log(S / K) + (r + sig**2 / 2) * T) / (sig * np.sqrt(T)); d2 = d1 - sig * np.sqrt(T)
    return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

def opm(ct, V0, T, r, sig):
    """class values via call spreads at breakpoints (payoff piecewise linear in V)."""
    bps, LP = breakpoints(ct)
    grid = np.concatenate([[0.0], bps])
    names = SERIES + ["common", "rsu", "warrants", "options"]
    pay = {n: np.array([payoff(ct, v)[n] for v in grid]) for n in names}
    # slopes on each segment [grid_k, grid_k+1] and beyond last bp
    vals = {}
    calls = np.array([bs_call(V0, k, T, r, sig) for k in grid])
    for n in names:
        y = pay[n]
        slopes = np.diff(y) / np.diff(grid)
        last = payoff(ct, grid[-1] * 1.5)[n]; slopes = np.append(slopes, (last - y[-1]) / (grid[-1] * 0.5))
        ds = np.diff(np.concatenate([[0.0], slopes]))
        vals[n] = float(np.sum(ds * calls))
    return vals, bps

def per_share(ct, vals):
    ps = {s: vals[s] / ct["pref"][i] for i, s in enumerate(SERIES) if ct["pref"][i] > 0}
    ps["common"] = vals["common"] / ct["common"]
    return ps

def backsolve(ct, series, price, T, r, sig):
    i = SERIES.index(series)
    f = lambda V0: opm(ct, V0, T, r, sig)[0][series] / ct["pref"][i] - price
    V0 = brentq(f, 1e8, 5e10)
    vals, bps = opm(ct, V0, T, r, sig)
    return V0, vals, bps

# ---------- DLOM ----------
def finnerty(T, sig):
    """Finnerty (2012) average-strike put, zero dividend, as fraction of value."""
    s2t = sig**2 * T
    v2t = s2t + np.log(2 * (np.exp(s2t) - s2t - 1)) - 2 * np.log(np.exp(s2t) - 1)
    v = np.sqrt(v2t)
    return norm.cdf(v / 2) - norm.cdf(-v / 2)

def chaffe(T, r, sig):
    """Chaffe (1993) at-the-money European protective put / spot."""
    d1 = (r + sig**2 / 2) * T / (sig * np.sqrt(T)); d2 = d1 - sig * np.sqrt(T)
    return np.exp(-r * T) * norm.cdf(-d2) - norm.cdf(-d1)

# ---------- PWERM / hybrid (31-Aug-2021) ----------
HYBRID_2021 = dict(  # weights and scenario DLOMs from Remitly's SEC response letter (15-Sep-2021)
    early=dict(w=0.775, dlom=0.135, t=22/365),      # IPO expected Sept-2021 (priced 22-Sep-2021)
    median=dict(w=0.15, dlom=0.16, t=0.50),         # assumption: ~Q1-2022
    opm=dict(w=0.075, dlom=0.25, T=2.0, sig=0.50),  # delayed exit; σ top of 2021 disclosed option-vol range
    r=0.0042, k=0.20, reported=28.73, midpoint=40.0, ipo_price=43.0)

def equity_from_price(ct, P):
    return P * (ct["common"] + ct["rsu"] + ct["pref"].sum()) + max(P - ct["opt_k"], 0) * ct["options"] + max(P - ct["war_k"], 0) * ct["warrants"]

def hybrid(P, h=HYBRID_2021, ct=None, parts=False, haircut=0.0):
    """haircut = fraction by which the delayed-exit (OPM) equity value sits below the IPO-implied equity value."""
    ct = ct or cap_table("2021-08-31")
    e, m, o = h["early"], h["median"], h["opm"]
    pv_e = P * np.exp(-h["k"] * e["t"]); pv_m = P * np.exp(-h["k"] * m["t"])
    vals, _ = opm(ct, equity_from_price(ct, P) * (1 - haircut), o["T"], h["r"], o["sig"]); pv_o = vals["common"] / ct["common"]
    comp = dict(early=e["w"] * pv_e * (1 - e["dlom"]), median=m["w"] * pv_m * (1 - m["dlom"]), opm=o["w"] * pv_o * (1 - o["dlom"]))
    marketable = e["w"] * pv_e + m["w"] * pv_m + o["w"] * pv_o
    if parts: return sum(comp.values()), comp, marketable
    return sum(comp.values())

def implied_ipo_price(target, **kw):
    return brentq(lambda P: hybrid(P, **kw) - target, 5, 100)

def two_condition_solve(target=28.73, blended=0.142):
    """solve IPO value P and OPM-scenario haircut so BOTH the reported FV and the disclosed blended DLOM are matched."""
    def blend(hc):
        P = implied_ipo_price(target, haircut=hc); v, comp, mkt = hybrid(P, parts=True, haircut=hc); return 1 - v / mkt
    hc = brentq(lambda h: blend(h) - blended, 0.0, 0.95)
    return implied_ipo_price(target, haircut=hc), hc

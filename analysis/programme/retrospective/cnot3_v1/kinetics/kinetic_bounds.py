"""Exact scalar identification sets for a declared finite two-rate model.

No data, biological package or calibrated uncertainty is used. All supplied
numbers are exact Fraction-compatible rates or externally stipulated bounds.
"""
from fractions import Fraction as Q


def bound_pathway(k0_interval, k1_interval, r_interval, delta_interval, physical=True):
    """Sharp projection onto p for independent finite closed parameter boxes.

    Model: k0=b+p, k1=b+delta+r*p, b,p>=0, 0<=r<=1.
    physical additionally requires b+delta>=0; it is not silently assumed.
    """
    l0, u0 = map(Q, k0_interval)
    l1, u1 = map(Q, k1_interval)
    rl, ru = map(Q, r_interval)
    dl, du = map(Q, delta_interval)
    if not (0 <= l0 <= u0 and 0 <= l1 <= u1 and 0 <= rl <= ru <= 1 and dl <= du):
        raise ValueError("Invalid finite closed rate, residual-fraction or shift bounds")
    lower, upper, contradictions = Q(0), u0, []
    # Eliminate the unobserved background and residual pathway contribution.
    constraints = [("lower", 1 - rl, l0 + dl - u1),
                   ("upper", 1 - ru, u0 + du - l1),
                   ("upper", rl, u1 - dl)]
    if physical:
        constraints += [("upper", Q(1), u0 + du), ("upper", rl, u1)]
    for sense, coefficient, rhs in constraints:
        if coefficient == 0:
            if (sense == "upper" and rhs < 0) or (sense == "lower" and rhs > 0):
                contradictions.append((sense, coefficient, rhs))
        elif sense == "lower":
            lower = max(lower, rhs / coefficient)
        else:
            upper = min(upper, rhs / coefficient)
    return {"status": "feasible" if not contradictions and lower <= upper else "infeasible",
            "p_interval": (lower, upper) if not contradictions and lower <= upper else None,
            "constraints": constraints, "physical_background_nonnegative": physical}


def witness(p, k0_interval, k1_interval, r_interval, delta_interval, physical=True):
    """Construct a model witness for any p in the returned interval."""
    p = Q(p)
    l0, u0 = map(Q, k0_interval)
    l1, u1 = map(Q, k1_interval)
    rl, ru = map(Q, r_interval)
    dl, du = map(Q, delta_interval)
    interval = bound_pathway(k0_interval, k1_interval, r_interval, delta_interval, physical)
    if interval["status"] != "feasible" or not interval["p_interval"][0] <= p <= interval["p_interval"][1]:
        raise ValueError("p is outside the identification set")
    bl, bu = max(Q(0), l0 - p), u0 - p
    nl, nu = bl + dl, bu + du
    if physical:
        nl = max(Q(0), nl)
    cl, cu = rl * p, ru * p
    k1 = max(l1, nl + cl)
    c = max(cl, k1 - nu)
    n = k1 - c
    b = max(bl, n - du)
    delta = n - b
    r = c / p if p else rl
    answer = {"k0": b + p, "k1": k1, "p": p, "b": b, "r": r, "delta": delta}
    if not (l0 <= answer["k0"] <= u0 and l1 <= k1 <= u1 and rl <= r <= ru and dl <= delta <= du
            and b >= 0 and p >= 0 and (not physical or b + delta >= 0)
            and answer["k0"] == b + p and k1 == b + delta + r * p):
        raise AssertionError("Constructive witness failed")
    return answer


def fold_and_drop(k0, k1):
    k0, k1 = Q(k0), Q(k1)
    if k0 <= 0 or k1 <= 0:
        raise ValueError("Finite positive half-lives require positive rates")
    return {"half_life_fold_t1_over_t0": k0 / k1, "absolute_net_rate_reduction": k0 - k1,
            "fractional_net_rate_reduction": 1 - k1 / k0}

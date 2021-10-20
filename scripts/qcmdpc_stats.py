#!/usr/bin/python
from math import log2, log1p, expm1, inf
from gmpy2 import bincoef

from statsmodels.stats.proportion import proportion_confint


ALPHA = 0.01


def log_one_minus(p, n):
    if p == 0:
        return inf
    p = float(p)
    return log2(-expm1(n * log1p(-p)))


def multiplicity_distance(r, d, m):
    return r * bincoef(d - 1, d - m - 1) * bincoef(r - d - 1, d - m - 1) / (d - m)


def get_density(n0, r, d, t, weak, error_floor, m, l):
    n = n0 * r
    w = n0 * d
    proba = 0.
    distance = None
    if weak == 1:
        proba += log_one_minus(
            bincoef(r - m, d - m) / bincoef(r, d), 2 * r * (r // 2))
    elif weak == 2:
        p = multiplicity_distance(r, d, m)
        proba += log_one_minus(p / bincoef(r, d), 2 * (r // 2))
    elif weak == 3:
        proba += log_one_minus(
            bincoef(d, m) * bincoef(r - d, d - m) / bincoef(r, d), r)

    if error_floor == 1:
        proba += log_one_minus(
            bincoef(d, l) * bincoef(n - d, t - l) / bincoef(n, t), n)
        distance = d + t - 2 * l
    elif error_floor == 2:
        proba += log_one_minus(
            bincoef(w, l) * bincoef(n - w, t - l) / bincoef(n, t), r * r)
        distance = w + t - 2 * l
    elif error_floor == 3:
        proba += log_one_minus(
            bincoef(w, l) * bincoef(n - w, t - l) / bincoef(n, t), r)
        distance = w + t - 2 * l
    return proba, distance


def cumul_iteration(distr):
    failures = 0 if -1 not in distr else distr[-1]
    total = sum(distr.values())
    cumul = {}

    for it, ct in sorted((distr.items()), reverse=True):
        if it == -1:
            break
        cumul[it] = failures
        failures += ct
    return cumul, total


def get_dfr(distr, alpha=ALPHA):
    def dfr_with_ci(f, n):
        dfr = f / n
        lower, upper = proportion_confint(f, n, alpha, method='beta')
        return log2(lower), log2(dfr), log2(upper)

    dfr = {}

    cumul, total = cumul_iteration(distr)

    for it, failures in sorted((cumul.items()), reverse=True):
        if it == -1:
            break
        dfr[it] = dfr_with_ci(failures, total)

    return dfr


def parse_param(s):
    def int_float_str(x):
        if isinstance(x, (int, float)):
            return x
        try:
            return int(x)
        except (TypeError, ValueError):
            pass
        try:
            return float(x)
        except (TypeError, ValueError):
            return x

    s = s.split()
    l = [p.split('=') for p in s]
    l = list(map(lambda x: (x[0][2:].lower(), int_float_str(x[1])), l))
    return dict(l)


def get_data(filename):
    stats = {}

    s_param = ""
    s_results = ""
    with open(filename, 'r') as file:
        for line in file:
            line = line.rstrip()
            if len(line) == 0:
                continue
            if line[0] == '-':
                s_param = " ".join(sorted(line.split()))
                continue
            s_results = line[:]

    entry = parse_param(s_param)

    results = None
    s_results = s_results.split()
    results = [e.split(":") for e in s_results[1:]]
    if results and results[-1][0][0] == '>':
        results[-1][0] = '-1'
    results = dict(map(lambda e: (int(e[0]), int(e[1])), results))

    if not results:
        return None

    entry['iteration'] = results
    entry['density'], entry['distance'] = get_density(
        *[entry[p] for p in ['index', 'block_length', 'block_weight', 'error_weight', 'weak', 'error_floor', 'weak_p', 'error_floor_p']])

    return entry


if __name__ == '__main__':
    import sys

    if len(sys.argv) <= 1:
        print("{} filename".format(sys.argv[0]))
        sys.exit(1)

    filename = sys.argv[1]
    print(get_data(filename))

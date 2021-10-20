from gmpy2 import mpfr, exp, log, log1p, lngamma
from scipy import integrate

# Default values.
EPSILON = mpfr(-32*log(2))
TOLERANCE = mpfr(1e-3)
STEP = 0.1


def lbeta(x, y):
    x *= 1.
    y *= 1.
    return lngamma(x) + lngamma(y) - lngamma(x+y)


def f(s, t, f1, n1, f2, n2, A):
    """
    Return
    \[
        \frac{1}{K}
        e^{s \frac{f2 + 1}{1 + A}}
        e^{t \frac{f1 + 1}{A} + \frac{f2 + 1}{1 + A}}
        (1 - e^\frac{t}{A})^{n1 - f1}
        (1 - e^\frac{s + t}{1 + A})^{n2 - f2}
    \]
    with
    \[
        K = A (1 + A) B(f1 + 1, n1 - f1 + 1) B(f2 + 1, n2 - f2 + 1) .
    \]
    """
    res = mpfr(0)
    res -= lbeta(f2 + 1, n2 - f2 + 1)
    res -= lbeta(f1 + 1, n1 - f1 + 1)
    res -= log(A * (1 + A))
    res += ((f2 + 1) / (1 + A)) * s
    res += ((f1 + 1) / A + (f2 + 1) / (1 + A)) * t
    res += (n2 - f2) * log1p(-exp((s + t) / (1 + A)))
    res += (n1 - f1) * log1p(-exp(t / A))
    return exp(res)


def get_support(f, s0, t0, epsilon, h):
    """
    Find a rectangular region where the 2-dimensional function 'f' is nonzero
    (greater than 'epsilon').

    Assume that 'f' is unimodal.

    Start from a point (s0, t0) where 'f' is non-zero, then move in one
    direction and naively look for the maximum value in the other direction,
    with steps of size 'h'. When this maximum is less than 'epsilon', the
    boundary of the support has been reached.

    Repeat the process for each direction to find a rectangular region
    including the support of 'f'.
    """
    def naive_search(g, x0, y0, step):
        x, y = x0, y0
        while x + step <= 0:
            best_g, best_y = g(x + step, y), y

            # Arbitrary bound at 100
            for mul in range(100):
                tests = [best_y + delta for delta in [-mul *
                                                      step, mul * step] if best_y + delta <= 0]
                if len(tests) == 0:
                    break

                new_g, new_y = max(map(lambda v: (g(x + step, v), v), tests))
                if new_g < best_g:
                    break
                best_g, best_y = new_g, new_y

            x, y = x + step, best_y
            if log(best_g) < epsilon:
                break
        return x

    smin = naive_search(lambda x, y: f(y, x), t0, s0, -h)
    smax = naive_search(lambda x, y: f(y, x), t0, s0, h)
    tmin = naive_search(lambda x, y: f(x, y), s0, t0, -h)
    tmax = naive_search(lambda x, y: f(x, y), s0, t0, h)

    return smin, smax, tmin, tmax


def conf_int(r3, r1, f1, n1, r2, f2, n2, alpha, epsilon=EPSILON, tolerance=TOLERANCE, step=STEP):
    """
    Return log(p-), log(p), log(p+) where p is the observed failure rate and
    [p-, p+] is the confidence interval.

    Procede using a bisection method.
    """
    A = (r3 - r2) / (r2 - r1)
    p1 = f1 / n1
    p2 = f2 / n2
    p3 = exp(log(p2) * (1 + A) + log(p1) * (-A))
    supp = get_support(lambda s, t: f(s, t, f1, n1, f2, n2, A),
                       log(p3), A*log(p1), epsilon, step)

    def bisect(f_int, a, b, x):
        while True:
            c = (a + b)/2
            I = f_int(c)
            if abs((I - x) / x) < tolerance:
                break
            if I < x:
                a = c
            else:
                b = c
        return c

    def int0(c):
        return integrate.dblquad(lambda s, t: f(s, t, f1, n1, f2, n2, A), supp[0], supp[1], supp[2], c)[0]

    def int1(c):
        return -integrate.dblquad(lambda s, t: f(s, t, f1, n1, f2, n2, A), supp[0], supp[1], c, supp[3])[0]

    c0 = bisect(int0, supp[2], log(p3), alpha / 2)
    c1 = bisect(int1, log(p3), supp[3], -alpha / 2)

    return c0, log(p3), c1


if __name__ == '__main__':
    import sys

    if len(sys.argv) <= 8:
        print("{} r3 r1 f1 n1 r2 f2 n2 alpha [epsilon] [tolerance]".format(
            sys.argv[0]))
        sys.exit(1)
    r3 = int(sys.argv[1])
    r1 = int(sys.argv[2])
    f1 = int(sys.argv[3])
    n1 = int(sys.argv[4])
    r2 = int(sys.argv[5])
    f2 = int(sys.argv[6])
    n2 = int(sys.argv[7])
    alpha = float(sys.argv[8])
    epsilon = float(sys.argv[9]) if len(sys.argv) > 9 else EPSILON
    tolerance = float(sys.argv[10]) if len(sys.argv) > 10 else TOLERANCE

    p0, p, p1 = conf_int(r3, r1, f1, n1, r2, f2, n2, alpha, epsilon, tolerance)

    print(" ".join(["{:.2}"] * 3).format(p0/log(2), p/log(2), p1/log(2)))

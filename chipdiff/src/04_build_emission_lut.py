#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
import numpy as np
from scipy.special import gammaln, betaln

def smart_open(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def read_header_counts_tops(path):
    """
    Read header metadata and compute maxima from a tabular "bin_counts" file.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the input file. The file is opened using smart_open(path, "rt"),
        so both local and remote (e.g. S3) paths supported.

    Returns
    -------
    tuple[int, int, int, int, int, int]
        A 6-tuple: (n1, n2, m, bin_size, max_x1, max_x2)
        - n1, n2, m, bin_size: integers parsed from the header tokens 'n1=',
          'n2=', 'm=' and 'bin=' on the first line (first line must start with '#').
        - max_x1, max_x2: maximum integer values observed in the 4th and 5th
          tab-separated columns (0-based indices 3 and 4) across all data rows.

    Behavior
    --------
    - Expects the first line to be a header starting with '#' and containing tokens
      of the form 'n1=...', 'n2=...', 'm=...', 'bin=...'. Tokens may be separated
      by spaces or tabs.
    - The second line may be an optional column header and is consumed if present;
      if the file has no data rows after the header, the function exits with SystemExit.
    - Subsequent lines are scanned for data rows. Lines that are empty, start with '#',
      or start with 'chrom' (case-insensitive) are skipped.
    - For each data row with at least five tab-separated columns, columns 4 and 5
      (indices 3 and 4) are parsed as integers and used to update max_x1 and max_x2.
    - If any of n1, n2, m, or bin_size cannot be parsed from the header, the
      function exits with SystemExit.

    Exceptions
    ----------
    SystemExit
        Raised when:
        - The first line does not start with '#' or required header tokens are missing.
        - The file contains no data rows after the header.
        - Required header values (n1, n2, m, bin) could not be parsed.

    Notes
    -----
    - The function intentionally ignores malformed data rows with fewer than five columns.
    - Column indices for maxima use zero-based indexing (columns 3 and 4).
    """
    n1 = n2 = m = bin_size = None
    max_x1 = 0
    max_x2 = 0
    with smart_open(path, "rt") as f:
        head = f.readline()
        if not head.startswith("#"):
            raise SystemExit("bin_counts missing header with n1/n2/m/bin.")
        toks = head.strip("# \n").replace("\t", " ").split()
        for t in toks:
            if t.startswith("n1="):
                n1 = int(t.split("=",1)[1])
            elif t.startswith("n2="):
                n2 = int(t.split("=",1)[1])
            elif t.startswith("m="):
                m  = int(t.split("=",1)[1])
            elif t.startswith("bin="):
                bin_size = int(t.split("=",1)[1])

        # consume optional column header
        line2 = f.readline()
        if not line2:
            raise SystemExit("bin_counts has no rows.")

        # scan data rows to get maxima
        for line in f:
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("chrom"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            x1 = int(parts[3])
            x2 = int(parts[4])
            if x1 > max_x1:
                max_x1 = x1
            if x2 > max_x2:
                max_x2 = x2

    if any(v is None for v in (n1,n2,m,bin_size)):
        raise SystemExit("Failed to parse n1/n2/m/bin from header.")
    return n1, n2, m, bin_size, max_x1, max_x2

def log_binom_coeff(n,k):
    """
    Compute the natural logarithm of the binomial coefficient "n choose k".

    Parameters
    ----------
    n : int or array_like
        Number of trials. For combinatorial interpretation `n` should be a
        non-negative integer. Real-valued inputs are accepted because the
        implementation uses the log-gamma function which extends factorials.
    k : int or array_like
        Number of successes. Must be broadcastable with `n`.

    Returns
    -------
    float or ndarray
        The natural logarithm of the binomial coefficient, ln(C(n, k)). If
        `n` and/or `k` are array-like, an array is returned following NumPy
        broadcasting rules.

    Notes
    -----
    This function uses the identity
        ln(C(n, k)) = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    which is numerically stable for large `n`. For integer inputs with
    0 <= k <= n the result equals the logarithm of the usual combinatorial
    binomial coefficient. Inputs outside this range or non-integer inputs
    may produce NaN or values corresponding to the generalized binomial
    coefficient via the gamma function.

    """
    return gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)

def build_log_grid(npts, pmin):
    """
    Build a grid of probability midpoints and corresponding weights that is uniform
    in log-probability space.

    The grid is constructed by partitioning u = log(p) uniformly on the interval
    [log(pmin), log(1 - pmin)] into npts bins, taking the midpoints u_mid of those
    bins, and mapping back to probability space with p = exp(u_mid). The weight
    associated with each midpoint is w = exp(u_mid) * du, which is the
    midpoint-approximation of the integral dp over the corresponding bin.

    Parameters
    ----------
    npts : int
        Number of grid points (bins). Must be >= 1.
    pmin : float
        Lower cutoff for p; the grid covers probabilities in (pmin, 1 - pmin).
        Must satisfy 0 < pmin < 0.5.

    Returns
    -------
    p : numpy.ndarray, shape (npts,)
        Array of probability midpoints (p = exp(u_mid)), each in the open interval
        (pmin, 1 - pmin).
    w : numpy.ndarray, shape (npts,)
        Array of weights for each bin: w_i = exp(u_mid_i) * du_i. These weights
        approximate the integral of dp over each bin; their sum approximates
        (1 - 2 * pmin).

    Notes
    -----
    - The grid is uniform in log-space, so it places more resolution at smaller p.
    - Inputs are not validated by the function; passing invalid values (e.g.
      pmin <= 0, pmin >= 0.5, or npts < 1) will typically raise a ValueError from
      the underlying log or produce an empty/invalid grid.
    - For sufficiently large npts, the midpoint weights provide a good
      approximation to the exact integral of e^u over each bin.
    """
    # uniform in u = log p on [log pmin, log(1-pmin)]
    u_edges = np.linspace(np.log(pmin), np.log(1.0 - pmin), npts + 1)
    u_mid   = 0.5 * (u_edges[:-1] + u_edges[1:])
    du      = (u_edges[1:] - u_edges[:-1])  # (npts,)
    p       = np.exp(u_mid)
    w       = np.exp(u_mid) * du            # dp = e^u du
    return p, w

def emission_logprob_for_pair(x1, x2, n1, n2, alpha, beta, tau,
                              P1, P2, W, state):
    # binomial coefficients and Beta normalizer
    logC1 = log_binom_coeff(n1, x1)
    logC2 = log_binom_coeff(n2, x2)
    logB  = betaln(alpha, beta)  # log Beta(alpha,beta)

    eps = 1e-300
    P1s = np.maximum(P1, eps)
    P2s = np.maximum(P2, eps)

    # ratio region constraints
    if state == "A1":
        region_mask = (P1s / P2s) >= tau
    elif state == "A2":
        region_mask = (P2s / P1s) >= tau
    elif state == "A0":
        r = P2s / P1s
        region_mask = (r >= 1.0/tau) & (r <= tau)
    else:
        raise ValueError("state must be one of 'A0','A1','A2'")

    # log-kernel (excluding constants)
    a1 = x1 + alpha - 1.0
    b1 = (n1 - x1) + beta - 1.0
    a2 = x2 + alpha - 1.0
    b2 = (n2 - x2) + beta - 1.0

    logK = (a1 * np.log(P1s) + b1 * np.log(1.0 - P1s)
           +a2 * np.log(P2s) + b2 * np.log(1.0 - P2s))

    # weighted sum over region: sum( exp(logK) * W )
    neg_inf = -1e300
    masked = np.where(region_mask, logK + np.log(W), neg_inf)

    mval = np.max(masked)
    if mval <= -1e290:
        return -1e300  # effectively zero

    s = np.exp(masked - mval).sum()
    log_integral = mval + np.log(s)

    # combine constants (independent priors => -2*logB)
    log_prob = logC1 + logC2 - 2.0*logB + log_integral
    return log_prob

def main():
    """
    Build and save an emission look-up table (LUT) for a 3-state HMM using a
    paper-faithful numeric quadrature near p ≈ 0.

    This function is intended to be used as the script entry point. It parses
    command-line arguments, reads sequencing depth and binning metadata from a
    counts file header, constructs a 2D numeric integration grid over allele
    frequencies (p1, p2) using a log-spaced 1D grid near zero, computes
    emission log-probabilities for each possible observed read-count pair
    (x1, x2) for three latent states ("A0", "A1", "A2"), and saves the LUT
    together with metadata as a compressed NumPy .npz file.

    Command-line arguments (via argparse)
    - --counts (required): path to counts file (expected header produced by
        read_header_counts_tops). The function calls read_header_counts_tops(args.counts)
        and expects it to return (n1, n2, m, bin_size, max_x1, max_x2).
    - --tau (optional, float, default=3.0): fold-change threshold τ used by
        emission probability calculations.
    - --grid (optional, int, default=81): number of grid points per axis for
        the 1D log-spaced quadrature.
    - --pmin (optional, float, default=None): minimum p for the log grid. If
        not provided the code uses max(1e-12, 1/(10*max(n1,n2))).
    - --out (required): output path for the compressed NPZ file (e.g.
        results/emission_lut.npz). Parent directories will be created if needed.

    Key steps performed
    1. Read header values: n1, n2, m, bin_size, max_x1, max_x2.
    2. Set beta = float(m) and alpha = 1.0 (fixed).
    3. Determine pmin if not specified.
    4. Build a 1D log-spaced grid p and corresponding integration weights w
         via build_log_grid(args.grid, pmin).
    5. Create 2D meshes P1, P2 and combined weights W for quadrature.
    6. Allocate a LUT array of shape (max_x1+1, max_x2+1, 3) initialized to a
         very small log-probability (-1e300) and compute entries by calling
         emission_logprob_for_pair(x1, x2, n1, n2, alpha, beta, tau, P1, P2, W, state)
         for each of the three states "A0","A1","A2".
    7. Periodically log progress to stderr (per 1000 x1 rows).
    8. Save results to a compressed .npz file and write a final message to stderr.

    Output (.npz) contents and types
    - lut: ndarray float64 of shape (max_x1+1, max_x2+1, 3) containing log-prob
        values for states in axis 2 ordered as ["A0","A1","A2"].
    - states: ndarray of unicode strings ["A0","A1","A2"].
    - n1, n2, m, bin_size: int64 metadata from counts header.
    - alpha, beta, tau, pmin: float64 numeric parameters used.
    - grid: int64 grid size used per axis.
    - grid_type: ndarray of unicode strings (e.g. ["log"]).

    Side effects and errors
    - Creates parent directories for the output path if necessary.
    - Writes progress and status messages to stderr.
    - Relies on the helper functions read_header_counts_tops, build_log_grid,
        emission_logprob_for_pair, and the availability of numpy and os modules.
    - May raise exceptions on missing/invalid counts file, permission errors
        when writing output, or failures inside the helper functions.

    Performance notes
    - The LUT computation iterates over all (x1, x2) pairs up to max_x1 and
        max_x2 and calls the emission integrator three times per pair. For large
        maxima this can be CPU- and memory-intensive; consider parallelization
        or reducing the grid size for testing.
    """
    ap = argparse.ArgumentParser(
        description=("Build emission LUT with log-spaced quadrature near p≈0 "
                     "(paper-faithful model; improved numeric grid).")
    )
    ap.add_argument("--counts", required=True, help="results/bin_counts.tsv.gz")
    ap.add_argument("--tau", type=float, default=3.0, help="fold-change threshold τ (default 3.0)")
    ap.add_argument("--grid", type=int, default=81, help="grid points per axis (default 81)")
    ap.add_argument("--pmin", type=float, default=None,
                    help="min p for log grid; default=max(1e-12, 1/(10*max(n1,n2)))")
    ap.add_argument("--out", required=True, help="output NPZ (e.g., results/emission_lut.npz)")
    args = ap.parse_args()

    n1, n2, m, bin_size, max_x1, max_x2 = read_header_counts_tops(args.counts)
    alpha, beta = 1.0, float(m)

    # choose pmin small enough for your depths
    pmin = args.pmin
    if pmin is None:
        pmin = max(1e-12, 1.0 / (10.0 * max(n1, n2)))

    # 1D log grid + weights
    p, w = build_log_grid(args.grid, pmin)

    # 2D meshes for p1,p2 and weights
    P1, P2 = np.meshgrid(p, p, indexing="ij")
    W1, W2 = np.meshgrid(w, w, indexing="ij")
    W = W1 * W2

    # LUT allocation: (max_x1+1, max_x2+1, 3 states [A0,A1,A2])
    shape = (max_x1 + 1, max_x2 + 1, 3)
    lut = np.full(shape, -1e300, dtype=np.float64)
    states = ["A0", "A1", "A2"]

    for x1 in range(max_x1 + 1):
        if x1 % 1000 == 0:
            sys.stderr.write(f"[lut] x1={x1}/{max_x1}\n")
        for x2 in range(max_x2 + 1):
            lut[x1, x2, 0] = emission_logprob_for_pair(
                x1, x2, n1, n2, alpha, beta, args.tau, P1, P2, W, "A0"
            )
            lut[x1, x2, 1] = emission_logprob_for_pair(
                x1, x2, n1, n2, alpha, beta, args.tau, P1, P2, W, "A1"
            )
            lut[x1, x2, 2] = emission_logprob_for_pair(
                x1, x2, n1, n2, alpha, beta, args.tau, P1, P2, W, "A2"
            )

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez_compressed(
        args.out,
        lut=lut,
        states=np.array(states, dtype='U'),
        n1=np.int64(n1),
        n2=np.int64(n2),
        m=np.int64(m),
        bin_size=np.int64(bin_size),
        alpha=np.float64(alpha),
        beta=np.float64(beta),
        tau=np.float64(args.tau),
        grid=np.int64(args.grid),
        pmin=np.float64(pmin),
        grid_type=np.array(["log"], dtype='U')
    )
    sys.stderr.write(f"[lut] saved {args.out} with shape {lut.shape} (x1,x2,states); pmin={pmin}\n")

if __name__ == "__main__":
    main()

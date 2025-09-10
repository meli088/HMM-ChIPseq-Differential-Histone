#!/usr/bin/env python3
import sys
import os
import gzip
import argparse
import numpy as np
from collections import defaultdict

def smart_open(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def parse_header_counts(path):
    """
    Parse the header line of a bin counts file and return the expected integer fields.

    The function opens the file at the given path using smart_open(path, "rt"), reads the
    first line and expects a header beginning with '#' containing key=value tokens.
    Accepted separators in the header line are spaces or tabs. The header must include
    the integer keys 'n1', 'n2', 'm' and 'bin' (for example: "# n1=100 n2=100 m=2 bin=200").
    If present, an additional second line (e.g. a column header) is consumed and ignored.

    Args:
        path (str or os.PathLike): Path to the counts file to parse.

    Returns:
        tuple[int, int, int, int]: A 4-tuple (n1, n2, m, bin) parsed from the header.

    Raises:
        SystemExit: If the first line is missing or does not start with '#' or if any of
                    the required keys ('n1', 'n2', 'm', 'bin') is absent from the header.
        ValueError: If a found key has a non-integer value (propagates from int()).
    """
    with smart_open(path, "rt") as f:
        line = f.readline()
        if not line or not line.startswith("#"):
            raise SystemExit("Header missing in bin_counts (expected: '# n1=.. n2=.. m=.. bin=..').")
        toks = line.strip("# \n").replace("\t", " ").split()
        vals = {}
        for t in toks:
            if "=" in t:
                k,v = t.split("=",1)
                vals[k] = int(v)
        for k in ("n1","n2","m","bin"):
            if k not in vals:
                raise SystemExit(f"Header incomplete: {k} absent.")
        # skip optional column header line
        _line2 = f.readline()
    return vals["n1"], vals["n2"], vals["m"], vals["bin"]

def load_bin_counts(path):
    """
    Load per-bin count data from a tabular counts file.

    This function reads a counts file (opened via smart_open) and returns:
    - a mapping of (chromosome, bin_index) -> (x1, x2) count tuple,
    - a mapping of chromosome -> number_of_bins (computed as max bin_index + 1),
    - and the header metadata values (n1, n2, m, bin_size) as returned by parse_header_counts.

    File format expectations and behavior:
    - The first lines of the file may contain header information; parse_header_counts(path)
        is used to extract n1, n2, m and bin_size before reading the body.
    - The reader skips one initial line (f.readline()) and then ignores any lines that are
        empty, start with '#' (comments), or start with the case-insensitive token "chrom".
    - Each data row is expected to be a tab-separated line with five columns:
        chrom, start, end, x1, x2
        where start, end, x1 and x2 are integers (start and end define the genomic interval).
    - The bin index is computed as start // bin_size (integer floor division).
    - The returned counts dictionary contains an entry for each observed (chrom, bin_index)
        with the integer tuple (x1, x2). If the same (chrom, bin_index) appears multiple times
        in the file, the last occurrence will overwrite earlier ones.
    - The chrom_bins mapping stores, for each chromosome, the number of bins equal to the
        largest observed bin_index + 1. Bins that are not present in the file will not appear
        in the counts mapping.

    Returns:
            tuple:
                    counts (dict): {(chrom (str), bin_index (int)): (x1 (int), x2 (int)), ...}
                    chrom_bins (dict): {chrom (str): num_bins (int), ...}
                    n1 (int): first count/header parameter from parse_header_counts
                    n2 (int): second count/header parameter from parse_header_counts
                    m (int): third header parameter from parse_header_counts
                    bin_size (int): bin size (base pairs) used to compute bin indices

    Raises:
            IOError: if the file cannot be opened by smart_open.
            ValueError: if a data line cannot be parsed into the expected five tab-separated fields
                                    or if numeric conversion to int fails.

    Notes:
            - smart_open allows reading from regular or compressed files (e.g., .gz).
            - The function intentionally does not fill missing bins; only observed bins are returned.
            - parse_header_counts(path) must exist in the same module or be importable and must
                return (n1, n2, m, bin_size).
    """
    n1, n2, m, bin_size = parse_header_counts(path)
    counts = {}
    chrom_bins = defaultdict(int)
    with smart_open(path, "rt") as f:
        # skip header lines
        f.readline()
        for line in f:
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("chrom"):
                continue
            chrom, s_s, e_s, x1_s, x2_s = line.rstrip("\n").split("\t")
            s, _e = int(s_s), int(e_s)
            x1, x2 = int(x1_s), int(x2_s)
            bidx = s // bin_size
            counts[(chrom, bidx)] = (x1, x2)
            if bidx + 1 > chrom_bins[chrom]:
                chrom_bins[chrom] = bidx + 1
    return counts, dict(chrom_bins), n1, n2, m, bin_size

def load_regions(path):
    """Charge les régions BED (chrom, start, end)."""
    regs = []
    with smart_open(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, s_s, e_s = parts[0], parts[1], parts[2]
            regs.append((chrom, int(s_s), int(e_s)))
    return regs

def load_lut(path):
    """
    Load a lookup table (LUT) and associated metadata from a NumPy .npz file.

    Parameters
    ----------
    path : str or os.PathLike
        Path to a NumPy .npz archive produced by the pipeline. The archive is
        expected to contain the following keys:
          - "lut"    : ndarray of shape (max_x1+1, max_x2+1, 3) containing
                       log-probabilities.
          - "states" : sequence/array-like of state labels.
          - "tau"    : scalar numeric value.
          - "grid"   : integer (grid size).

    Returns
    -------
    tuple
        A tuple (lut, states, tau, grid) where:
          - lut (numpy.ndarray): The LUT array of log-probabilities.
          - states (list): A Python list of state labels.
          - tau (float): The tau parameter converted to float.
          - grid (int): The grid parameter converted to int.

    Raises
    ------
    FileNotFoundError
        If the file at `path` does not exist.
    KeyError
        If any required key ("lut", "states", "tau", "grid") is missing from the archive.
    ValueError
        If the loaded values cannot be converted to the expected types
        (for example, non-scalar "tau" or non-integer "grid").

    Notes
    -----
    The .npz file is loaded with numpy.load(..., allow_pickle=True) to permit
    pickled objects for entries such as "states". The function performs
    explicit conversion of `tau` to float and `grid` to int before returning.
    """
    d = np.load(path, allow_pickle=True)
    lut = d["lut"]           # shape (max_x1+1, max_x2+1, 3) in log prob
    states = list(d["states"])
    tau = float(d["tau"])
    grid = int(d["grid"])
    return lut, states, tau, grid

def build_sequences(counts, bin_size, regions, max_regions=10000, rng=None):
    """
    Build a list of binned count sequences for genomic regions.

    For each region (chrom, start, end) this function computes the inclusive
    bin range covered by the region using the provided bin_size, then
    collects the values from the counts mapping for each bin in that range.
    Bins that are not present in counts are skipped. Regions with no
    available bins (after skipping missing bins) are omitted from the result.

    Parameters
    ----------
    counts : Mapping[(str, int), number]
        Mapping from (chromosome, bin_index) to a numerical count (e.g. int or float).
        The function accesses values via counts.get((chrom, bin_index)).
    bin_size : int
        Size of each bin in genomic coordinates. start and end coordinates are
        converted to bin indices by integer division. The last bin index is
        computed as (end - 1) // bin_size so that an interval [start, end) that
        ends at a bin boundary does not include an extra empty bin.
    regions : Sequence[Tuple[str, int, int]]
        Iterable of regions, each given as (chrom, start, end). Regions with
        end <= start are ignored.
    max_regions : int, optional
        If the number of input regions exceeds max_regions, a uniform random
        subset of size max_regions is selected (without replacement).
        Default is 10000.
    rng : numpy.random.Generator or None, optional
        Random number generator used to sample a subset of regions when needed.
        If None, a Generator seeded with 42 is created internally, producing
        deterministic sampling by default.

    Returns
    -------
    List[List[number]]
        A list of sequences, one per region with at least one available bin.
        Each sequence is a list of the count values retrieved from counts for
        consecutive bin indices spanning the region. The ordering of returned
        sequences follows the (possibly subsampled and implicitly shuffled)
        order of the regions after sampling.

    Notes
    -----
    - Missing bins (keys absent in counts) are silently skipped; a region may
      therefore produce a shorter sequence than the full bin span or be
      dropped entirely if no bins are present.
    - When subsampling is triggered, the sampled region order is the order of
      the indices returned by rng.choice (i.e. effectively randomized).
    - The function does not raise on missing keys or invalid inputs; it simply
      ignores invalid regions or missing bins.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    if len(regions) > max_regions:
        idx = rng.choice(len(regions), size=max_regions, replace=False)
        regions = [regions[i] for i in idx]

    sequences = []
    for chrom, start, end in regions:
        if end <= start:
            continue
        b0 = start // bin_size
        b1 = (end - 1) // bin_size
        seq = []
        for b in range(b0, b1 + 1):
            x = counts.get((chrom, b))
            if x is None:
                # bin inexistant (rare si tailles chrom bien matching)
                continue
            seq.append(x)
        if len(seq) > 0:
            sequences.append(seq)
    return sequences

def emission_logprobs_for_seq(seq, lut):
    """
    Compute per-position emission log-probabilities for a sequence of discrete observations.

    Parameters
    ----------
    seq : Sequence[Tuple[int, int]]
        Sequence of length T where each element is a pair of non-negative integers (x1, x2)
        representing the observed counts/indices at a time position. Elements must be indexable
        (e.g., list of tuples, NumPy array with shape (T,2)). Values larger than the LUT bounds
        are clamped to the maximum available index (see Notes).
    lut : numpy.ndarray
        A lookup table of shape (M, N, 3) containing precomputed log-probabilities for each
        possible (x1, x2) pair and for each of the three hidden states. The last dimension
        corresponds to the three states (A0, A1, A2). dtype should be float or float64.

    Returns
    -------
    numpy.ndarray
        Array B of shape (T, 3) and dtype float64 where B[t, :] contains the log-probabilities
        for the three states at time t (in the same order as the last dimension of `lut`).

    Notes
    -----
    - The function assumes exactly three emission states (K == 3). It does not infer K from
      `lut.shape[-1]`, but the LUT must have size 3 in its last dimension.
    - For any observed index x1 > M-1 or x2 > N-1, the function uses the last row/column
      of `lut` (i.e., it clamps indices to M-1 and N-1 respectively). This implements a
      saturation strategy for out-of-range observations instead of raising an error.
    - The returned values are log-probabilities (not normalized probabilities).

    Examples
    --------
    >>> # Given lut with shape (5, 5, 3) and a two-position sequence:
    >>> seq = [(0, 1), (7, 2)]         # second position x1=7 will be clamped to 4 (M-1)
    >>> B = emission_logprobs_for_seq(seq, lut)
    >>> B.shape
    (2, 3)
    """
    T = len(seq)
    K = 3
    B = np.empty((T, K), dtype=np.float64)
    max_x1 = lut.shape[0] - 1
    max_x2 = lut.shape[1] - 1
    for t,(x1,x2) in enumerate(seq):
        i = x1 if x1 <= max_x1 else max_x1
        j = x2 if x2 <= max_x2 else max_x2
        B[t,:] = lut[i, j, :]  # [A0,A1,A2] en log
    return B  # (T,3)

def logsumexp(a, axis=None):
    """
    Compute the logarithm of the sum of exponentials of input elements in a numerically stable way.

    Parameters
    ----------
    a : array_like
        Input array. Values are promoted to an ndarray for computation.
    axis : None, int, or tuple of ints, optional
        Axis or axes over which the sum is taken. If None (default), the sum is taken over all elements.

    Returns
    -------
    res : ndarray or scalar
        The value(s) of log(sum(exp(a))) computed in a numerically stable manner. If `axis` is None a scalar is returned;
        otherwise an array with the specified axis/axes removed.

    Notes
    -----
    The direct computation of log(sum(exp(a))) can overflow/underflow for large/small values of `a`. This function
    subtracts the maximum value along the reduction axis before exponentiation to improve numerical stability:
        log(sum(exp(a))) = max(a) + log(sum(exp(a - max(a)))).
    The implementation uses `keepdims=True` for the intermediate reduction and then squeezes the reduced axis to match
    the conventional reduced shape.

    NaN handling follows NumPy semantics (NaNs propagate).

    Raises
    ------
    ValueError
        May be raised by underlying NumPy operations if the input is empty along the reduction axis.

    """
    amax = np.max(a, axis=axis, keepdims=True)
    s = np.log(np.sum(np.exp(a - amax), axis=axis, keepdims=True)) + amax
    return np.squeeze(s, axis=axis)

def forward_backward_log(B, A, pi):
    """
    Compute the forward-backward posterior quantities for a Hidden Markov Model using log-domain arithmetic.

    Parameters
    ----------
    B : ndarray, shape (T, K)
        Log-emission (or log-likelihood) matrix: B[t, k] = log p(obs_t | z_t = k).
        T is the number of time steps and K is the number of hidden states.
    A : ndarray, shape (K, K)
        Log-transition matrix: A[i, j] = log p(z_{t+1} = j | z_t = i).
        Rows index the "from" state i and columns index the "to" state j.
    pi : ndarray, shape (K,)
        Log initial state distribution: pi[k] = log p(z_0 = k).

    Returns
    -------
    log_gamma : ndarray, shape (T, K)
        log_gamma[t, k] = log p(z_t = k | observations) — the posterior marginal (in log-space)
        for each time t and state k.
    xi_sum : ndarray, shape (K, K)
        Pairwise transition posterior summed over time in log-space:
        xi_sum[i, j] = log( sum_{t=0}^{T-2} p(z_t = i, z_{t+1} = j | observations) ).
        Note: the function accumulates these pairwise posteriors in the linear domain
        (with small regularisation) and stores the result as a log value for numerical
        stability of later parameter updates.
    loglik : float
        Log marginal likelihood of the entire observation sequence: log p(observations).
        Computed as log-sum-exp over the final forward log probabilities.

    Notes
    -----
    - All inputs are expected to be in log-space. This implementation keeps computations
      primarily in the log domain and uses log-sum-exp to avoid underflow/overflow.
    - Time complexity is O(T * K^2) due to the pairwise transition terms.
    - The function relies on a numerically stable log-sum-exp implementation (e.g. scipy.special.logsumexp
      or an equivalent). A tiny constant (1e-300) is used when converting small positive values
      to/from the linear domain while accumulating xi to avoid taking log(0).
    - The returned xi_sum is convenient for M-step updates that require expected transition counts;
      since it is the log of the summed joint posteriors over time, exponentiating xi_sum yields the
      expected counts matrix (up to numerical precision).
    """
    T, K = B.shape
    # Forward
    log_alpha = np.empty((T, K))
    log_alpha[0,:] = pi + B[0,:]
    for t in range(1, T):
        # log_alpha[t,k] = B[t,k] + logsum_j(log_alpha[t-1,j] + logA[j,k])
        log_alpha[t,:] = B[t,:] + logsumexp(log_alpha[t-1,:][:,None] + A, axis=0)

    loglik = logsumexp(log_alpha[-1,:], axis=0)

    # Backward
    log_beta = np.empty((T, K))
    log_beta[-1,:] = 0.0
    for t in range(T-2, -1, -1):
        # log_beta[t,j] = logsum_k( logA[j,k] + B[t+1,k] + log_beta[t+1,k] )
        log_beta[t,:] = logsumexp(A + (B[t+1,:] + log_beta[t+1,:])[None,:], axis=1)

    # Gamma
    log_gamma = log_alpha + log_beta - loglik  # log P(z_t=k | obs)
    # Xi (sum on t)
    # log_xi[t,i,j] ∝ log_alpha[t,i] + logA[i,j] + B[t+1,j] + log_beta[t+1,j]
    xi_sum = np.full_like(A, -np.inf)
    for t in range(T-1):
        M = (log_alpha[t,:][:,None] + A + (B[t+1,:] + log_beta[t+1,:])[None,:]) - loglik
        # sum stable in log: log(exp(x1)+exp(x2)+...)
        # we accumulate in linear but via log-sum-exp pairwise
        # To simplify: we go through the linear space locally (K small)
        xi_lin = np.exp(M)  # (K,K)
        if np.isneginf(xi_sum).all():
            xi_sum = np.log(xi_lin + 1e-300)
        else:
            # accumulation stable: log( exp(old) + xi_lin )
            xi_sum = np.log(np.exp(xi_sum) + xi_lin + 1e-300)
    return log_gamma, xi_sum, loglik

def normalize_rows_log(mat_log):
    """
    Normalize each row of a matrix of log-values using the log-sum-exp trick.

    This function takes a 2-D array-like of log-scores and returns a new array
    of the same shape where each row has been normalized in log-space so that
    the exponentiated values in each row sum to 1 (i.e., the rows represent
    log-probability distributions). The implementation is numerically stable
    because it uses the log-sum-exp operation to avoid overflow/underflow.

    Parameters
    ----------
    mat_log : array-like, shape (n_rows, n_cols)
        Input matrix containing log-values (log-scores, log-likelihoods, or
        unnormalized log-probabilities). Each row is normalized independently.

    Returns
    -------
    normalized : ndarray, shape (n_rows, n_cols)
        Array of log-values where each row has been shifted by the log-sum-exp
        of that row. For each row i:
            exp(normalized[i]).sum() == 1 (up to floating-point precision)

    Notes
    -----
    - The input is not modified in-place; a new array is returned.
    - Works correctly for rows containing large negative values and -inf entries.
    - Expects a 2-D input; behavior with other shapes is unspecified.
    """
    return mat_log - logsumexp(mat_log, axis=1)[:,None]

def train_transitions(counts_path, regions_path, lut_path, out_path,
                      max_iters=20, train_regions=10000, seed=123):
    counts, chrom_bins, n1, n2, m, bin_size = load_bin_counts(counts_path)
    regs = load_regions(regions_path)
    lut, states, tau, grid = load_lut(lut_path)
    K = 3  # A0,A1,A2

    rng = np.random.default_rng(seed)
    seqs = build_sequences(counts, bin_size, regs, max_regions=train_regions, rng=rng)
    if not seqs:
        raise SystemExit("Aucune séquence de formation construite. Vérifie putative_regions et bin_counts.")

    # Emissions (log) par séquence
    seq_B = [emission_logprobs_for_seq(s, lut) for s in seqs]

    # Initialisation des transitions A en log: diagonale favorisée (états persistent)
    A = np.array([
        [0.90, 0.05, 0.05],
        [0.05, 0.90, 0.05],
        [0.05, 0.05, 0.90],
    ], dtype=np.float64)
    A = np.log(A)

    # π fixed at A0 (initial state = no difference)
    pi = np.array([0.0, -np.inf, -np.inf], dtype=np.float64)

    ll_hist = []
    for it in range(1, max_iters+1):
        # E-step: cumulate gamma & xi
        xi_num = np.full((K,K), -np.inf)   # acumulator in log
        gamma_i_sum = np.zeros(K, dtype=np.float64)  # sum_t P(z_t=i | obs) to normalize A
        total_ll = 0.0

        for B in seq_B:
            log_gamma, xi_sum, loglik = forward_backward_log(B, A, pi)
            total_ll += float(loglik)
            # cumulate xi (in linear via log-sum-exp on sequences)
            if np.isneginf(xi_num).all():
                xi_num = xi_sum
            else:
                xi_num = np.log(np.exp(xi_num) + np.exp(xi_sum) + 1e-300)
            # cumulate gamma except the last time (like in Baum-Welch for A)
            gamma_lin = np.exp(log_gamma[:-1,:])  # (T-1,K)
            gamma_i_sum += gamma_lin.sum(axis=0)

        ll_hist.append(total_ll)

        # M-step: A_ij = sum_t xi_ij / sum_t gamma_i(t)
        xi_lin = np.exp(xi_num)  # (K,K)
        # no division by 0
        denom = gamma_i_sum + 1e-300
        A_new = xi_lin / denom[:,None]
        # secure: normalise line
        A_new = A_new / A_new.sum(axis=1, keepdims=True)
        A = np.log(A_new)

        sys.stderr.write(f"[train] iter {it:02d}  loglik={total_ll:.3f}\n")

        # stop criterion (very small growth)
        if it > 1 and abs(ll_hist[-1] - ll_hist[-2]) < 1e-3:
            sys.stderr.write("[train] early stop (Δloglik < 1e-3)\n")
            break

    # save
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savez_compressed(
        out_path,
        A=np.exp(A),                 # transition matrix in prob (not log)
        A_log=A,                     # matrix in log (for reference)
        states=np.array(["A0","A1","A2"], dtype="U"),
        ll_hist=np.array(ll_hist, dtype=np.float64),
        train_regions=np.int64(train_regions),
        seed=np.int64(seed)
    )
    sys.stderr.write(f"[train] saved {out_path}\n")

def main():
    ap = argparse.ArgumentParser(
        description=("Baum–Welch (transitions seules) sur ~10k régions putatives, "
                     "émissions fixées via LUT; état initial A0; papier-fidèle.")
    )
    ap.add_argument("--obs", required=True, help="results/bin_counts.tsv.gz")
    ap.add_argument("--putative", required=True, help="results/putative_regions.bed.gz")
    ap.add_argument("--lut", required=True, help="results/emission_lut.npz")
    ap.add_argument("--out", required=True, help="results/transitions.npz")
    ap.add_argument("--max-iters", type=int, default=20)
    ap.add_argument("--train-regions", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    train_transitions(args.obs, args.putative, args.lut, args.out,
                      max_iters=args.max_iters,
                      train_regions=args.train_regions,
                      seed=args.seed)

if __name__ == "__main__":
    main()

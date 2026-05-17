import numpy as np
import random
import logging
import numpy.typing as npt
from scipy.stats import gamma


def simulate_pores_decay(
    ratio_roi: npt.NDArray,
    ratio_roi_alone: npt.NDArray,
    size_dist: npt.NDArray,
    simulation_time: int,
    initial_pores: int,
    max_pores: int,
    decay: bool,
    debug=False,
    WGS: bool = False,
):

    logging.info(
        f"Enrichment mode simulation to evaluate the optimal reference length:\nLength distribution mean: {np.mean(size_dist)}, initial number of pores: {initial_pores} \n simulation time (hours) {simulation_time / 3600}"
    )

    # Initialize simulation state
    pores = [0] * initial_pores  # Time remaining for each pore (0 if free)

    accepted_molecules = []
    current_time = 0
    num_pores = initial_pores
    p = 0
    ## dead pores
    alive_pores = [1] * initial_pores  # Track which pores are dead

    if WGS:
        p = random.choice(ratio_roi_alone)
    else:
        p = random.choice(ratio_roi)  # probablity of finding a read mapping to the ROI
    total_reads = 0
    while current_time < simulation_time:
        for i in range(num_pores):
            if pores[i] <= 0:  # Pore is free                # Generate a new molecule
                # if the pore is dead and decay is on, it should be not be processing reads
                if decay:
                    if alive_pores[i] == 0:
                        # print(f"Pore number {i} is dead continue to next")
                        continue
                total_reads += 1

                if random.random() < p:  # Molecule is accepted
                    length = random.choice(
                        size_dist
                    )  # get a random size from the distribution
                    pores[i] = length / 400  # Time proportional to length
                    logging.debug(f"Molecule accepted with probability {p}")
                    accepted_molecules.append(length)
                else:  # Molecule is rejected
                    if WGS:  # WGS mode will sequence everything
                        length = random.choice(
                            size_dist
                        )  # get a random size from the distribution
                        pores[i] = length / 400
                        # accepted_molecules.append(length)
                    else:
                        pores[i] = 2  # Fixed rejection time

            else:
                pores[i] -= 1  # Decrement time remaining for this pore

        current_time += 1

    return {
        "total_reads": total_reads,
        "accepted_molecules": accepted_molecules,
        "alive_pores": alive_pores.count(1),
    }


def size_dist_ratio_from_array(
    size,
    k=1000,
    q=0.99,
    roi_count=7,
    roi_length=500,
    m=100,
    _e=500,
    ref_size=0,
    genome_size=4_600_000,
):
    """
    This function takes as input a mean size for the reads and ouputs ratio of kmers
    from the ROI or the optimal Ref with regard to total amount kmers present in 1 single genome

    """

    raw_reads_len = size
    total_kmers = np.int32(genome_size) - raw_reads_len
    if ref_size == 0:
        _roi_size_optimal = (
            2 * gamma_mle_q(raw_reads_len, q)["q_hat"] - roi_length
        )  # full region for a full refere
    else:
        _roi_size_optimal = ref_size

    roi_kmers = roi_count * (
        _roi_size_optimal + raw_reads_len - 2 * m
    )  ## Kmers mapping to the chosen reference (Default: optimal size)

    # TODO add a condition for reads longer than s-2e where s is the size of ref and *e* is the distance of detection (500bp)
    raw_reads_len_alt_mask = raw_reads_len > (_roi_size_optimal - (2 * _e))
    roi_kmers[raw_reads_len_alt_mask] = roi_count * (
        2 * (_e + _roi_size_optimal - 2 * 100)
    )
    ## Kmers mapping to the ROI alone
    roi_kmers_alone = roi_count * (roi_length + raw_reads_len - (2 * m))

    # TODO fix the ratio that shoot over 1 when the size of the Ref decreases
    roi_ratio_alone = roi_kmers_alone / total_kmers

    ratio_roiandoptimal_ref = (
        roi_kmers_alone / roi_kmers
    )  # ratio of reads containing the full ROI out of the reads falling into the optimal ref
    ratio_roiandoptimal_ref[raw_reads_len_alt_mask] = 0

    non_roi_kmers = total_kmers - roi_kmers
    ratio_opt_ref_kmers = roi_kmers / total_kmers  # optimal ref accepting ratio

    return (
        ratio_opt_ref_kmers,
        raw_reads_len,
        ratio_roiandoptimal_ref,
        _roi_size_optimal,
        roi_ratio_alone,
    )


def find_optimal_quantile(
    size,
    time,
    q_start=0.5,
    q_end=0.99,
    steps=10,
    wgs=False,
    genome_size=4_600_000,
    roi_length=500,
    _num_sim=3,
    _roi_count=7,
):
    q_range = np.linspace(q_start, q_end, steps)
    group_id = []
    q_list = []
    _mean_count = []
    _opt_size_list = []
    # The reads are going to be fit to a gamma distribution using MLE and the parameters are going to be used to generate new reads for the simulation.
    gamma_fit = gamma_mle_q(size, q=q_end)
    rng = np.random.default_rng(seed=1234)
    size = np.int64(
        rng.gamma(shape=gamma_fit["shape"], scale=gamma_fit["scale"], size=2000)
    )

    for _q in q_range:
        _ratio_roi, _size_dist, _ratio_roi_opt, _optimal_size, _ratio_roi_alone = (
            size_dist_ratio_from_array(
                size=size,
                roi_count=_roi_count,
                q=_q,
                roi_length=roi_length,
                genome_size=genome_size,
            )
        )
        print(
            f"Average ratio of reads containing roi: {np.mean(_ratio_roi_opt)}\n for an extended ref of {_optimal_size}"
        )
        _results = [
            simulate_pores_decay(
                ratio_roi=tuple(_ratio_roi),
                ratio_roi_alone=_ratio_roi_alone,
                size_dist=tuple(_size_dist),
                simulation_time=time,
                initial_pores=400,
                max_pores=500,
                decay=False,
                WGS=wgs,
            )
            for i in range(_num_sim)
        ]
        group_id.append(_q)
        _l = [len(r.get("accepted_molecules")) * _ratio_roi_opt for r in _results]
        _opt_size_list.append(_optimal_size)
        q_list.append(_q)
        _mean_count.append(np.mean(_l))
    q_max = q_list[np.argmax(_mean_count)]
    _opt_size = _opt_size_list[np.argmax(_mean_count)]
    print(
        f"Max read count {max(_mean_count)}  found at {round(q_max, 2)} with ext size of {int(_opt_size)}"
    )
    return {"q_max": round(float(q_max), 2), "ext_length": (int(_opt_size))}


def gamma_size(mu, sigma, n: int, seed=1234):
    """
    Wrapper function to create gamma distributed reads with
    mu and sigma and n(how many samples) as a an input
    """
    rng = np.random.default_rng(seed=seed)

    a_gen = (mu / sigma) ** 2
    b_gen = mu / sigma**2

    return np.int64(rng.gamma(shape=a_gen, scale=1 / b_gen, size=n))


def gamma_mle_q(x, q=0.9, bootstrap=False, n_boot=1000, ci=0.95, random_state=None):
    """
    Estimate a quantile (default q90) assuming Gamma distribution via MLE.

    Parameters
    ----------
    x : array-like
        Positive observations (e.g., read lengths)
    q : float
        Quantile to estimate (default 0.9)
    bootstrap : bool
        Whether to compute bootstrap CI
    n_boot : int
        Number of bootstrap samples
    ci : float
        Confidence interval level (e.g., 0.95)
    random_state : int or None
        Seed for reproducibility

    Returns
    -------
    dict with:
        - q_hat: estimated quantile
        - shape: k
        - scale: theta
        - (optional) ci_low, ci_high
    """

    x = np.asarray(x)

    # --- basic validation ---
    if np.any(x == 0):
        print(f"{len(x[x == 0])} reads are zero setting them to 1")
        x[x <= 0] = 1

    if np.any(x < 0):
        raise ValueError("All observations must be positive for Gamma MLE.")
    if not (0 < q < 1):
        raise ValueError("q must be in (0,1).")

    # --- fit gamma (MLE) ---
    # force loc=0 → standard gamma
    k, loc, theta = gamma.fit(x, floc=0)

    # --- compute quantile ---
    q_hat = gamma.ppf(q, a=k, scale=theta)

    result = {"q_hat": q_hat, "shape": k, "scale": theta}

    # --- optional bootstrap ---
    if bootstrap:
        rng = np.random.default_rng(random_state)
        boot_q = []

        for _ in range(n_boot):
            sample = rng.choice(x, size=len(x), replace=True)
            k_b, _, theta_b = gamma.fit(sample, floc=0)
            boot_q.append(gamma.ppf(q, a=k_b, scale=theta_b))

        alpha = 1 - ci
        low = np.percentile(boot_q, 100 * (alpha / 2))
        high = np.percentile(boot_q, 100 * (1 - alpha / 2))

        result["ci_low"] = low
        result["ci_high"] = high

    return result

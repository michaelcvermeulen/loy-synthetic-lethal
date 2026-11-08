from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp


@dataclass(frozen=True)
class GaussianMixture1D:
    means: np.ndarray
    sigmas: np.ndarray
    weights: np.ndarray
    log_likelihood: float

    @property
    def separation(self) -> float:
        pooled = np.sqrt(np.mean(np.square(self.sigmas)))
        return float((self.means[1] - self.means[0]) / pooled)


def _log_density(values: np.ndarray, means: np.ndarray, sigmas: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return (
        np.log(weights)[None, :]
        - np.log(sigmas)[None, :]
        - 0.5 * np.log(2.0 * np.pi)
        - 0.5 * np.square((values[:, None] - means[None, :]) / sigmas[None, :])
    )


def fit_two_gaussian(
    values: np.ndarray,
    seed: int = 20260828,
    n_init: int = 8,
    tolerance: float = 1e-8,
    max_iterations: int = 500,
) -> GaussianMixture1D:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.ndim != 1 or data.size < 10:
        raise ValueError("at least ten finite one-dimensional values are required")
    spread = float(np.std(data))
    if spread <= 0:
        raise ValueError("reporter values have no variance")

    rng = np.random.default_rng(seed)
    base = np.quantile(data, [0.25, 0.75])
    best: GaussianMixture1D | None = None
    for _ in range(n_init):
        means = np.sort(base + rng.normal(0.0, 0.05 * spread, size=2))
        sigmas = np.full(2, max(spread, 1e-6), dtype=float)
        weights = np.full(2, 0.5, dtype=float)
        previous = -np.inf
        valid_fit = True
        for _iteration in range(max_iterations):
            log_prob = _log_density(data, means, sigmas, weights)
            normalizer = logsumexp(log_prob, axis=1)
            log_likelihood = float(normalizer.sum())
            responsibilities = np.exp(log_prob - normalizer[:, None])
            effective = responsibilities.sum(axis=0)
            if np.any(effective <= 1e-8):
                valid_fit = False
                break
            weights = effective / data.size
            means = (responsibilities * data[:, None]).sum(axis=0) / effective
            variances = (
                responsibilities * np.square(data[:, None] - means[None, :])
            ).sum(axis=0) / effective
            sigmas = np.sqrt(np.maximum(variances, 1e-12))
            if abs(log_likelihood - previous) <= tolerance * (1.0 + abs(log_likelihood)):
                break
            previous = log_likelihood
        if not valid_fit:
            continue
        order = np.argsort(means)
        candidate = GaussianMixture1D(
            means=np.asarray(means[order], dtype=float),
            sigmas=np.asarray(sigmas[order], dtype=float),
            weights=np.asarray(weights[order], dtype=float),
            log_likelihood=log_likelihood,
        )
        if best is None or candidate.log_likelihood > best.log_likelihood:
            best = candidate
    if best is None:
        raise ValueError("two-component Gaussian fitting failed")
    return best


def equal_posterior_threshold(
    model: GaussianMixture1D,
    minimum_weight: float = 0.01,
    minimum_separation: float = 2.0,
) -> float:
    if np.any(model.weights < minimum_weight):
        raise ValueError("component weight below minimum")
    if model.separation < minimum_separation:
        raise ValueError("insufficient component separation")
    mean1, mean2 = model.means
    sigma1, sigma2 = model.sigmas
    weight1, weight2 = model.weights
    a = -1.0 / (2.0 * sigma1**2) + 1.0 / (2.0 * sigma2**2)
    b = mean1 / sigma1**2 - mean2 / sigma2**2
    c = (
        -mean1**2 / (2.0 * sigma1**2)
        + mean2**2 / (2.0 * sigma2**2)
        + np.log(weight1 / sigma1)
        - np.log(weight2 / sigma2)
    )
    if abs(a) < 1e-12:
        roots = np.array([-c / b])
    else:
        roots = np.roots([a, b, c])
    between = [float(root.real) for root in roots if abs(root.imag) < 1e-8 and mean1 < root.real < mean2]
    if len(between) != 1:
        raise ValueError("no unique equal-posterior threshold between component means")
    return between[0]


def classify_reporters(
    gfp: np.ndarray,
    mcherry: np.ndarray,
    gfp_threshold: float,
    mcherry_threshold: float,
) -> np.ndarray:
    gfp_positive = np.asarray(gfp, dtype=float) >= gfp_threshold
    mcherry_positive = np.asarray(mcherry, dtype=float) >= mcherry_threshold
    if gfp_positive.shape != mcherry_positive.shape:
        raise ValueError("reporter arrays must have identical shapes")
    result = np.full(gfp_positive.shape, "Double-negative", dtype=object)
    result[gfp_positive & ~mcherry_positive] = "GFP-only"
    result[~gfp_positive & mcherry_positive] = "mCherry-only"
    result[gfp_positive & mcherry_positive] = "Double-positive"
    return result


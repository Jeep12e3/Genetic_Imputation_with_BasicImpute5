"""
Simplified IMPUTE5-style genotype imputation demo.

This script is meant for presentation/assignment purposes. It follows the
main implementation idea from the IMPUTE5 paper:

    PBWT-like state selection -> reduced HMM -> forward-backward -> imputation

The real IMPUTE5 implementation is much more optimized and supports large
reference panels, indexed imp5 files, BGEN output, parallelization, and a more
complete Li and Stephens model. This demo focuses on clarity.
"""

from bisect import bisect_left
from dataclasses import dataclass
from typing import Iterable, Optional


Allele = int
MaybeAllele = Optional[int]


@dataclass(frozen=True)
class ImputationResult:
    selected_states: list[int]
    posterior_by_marker: list[list[float]]
    imputed_prob_allele_1: list[Optional[float]]


def reverse_prefix_key(haplotype: list[Allele], positions: list[int], marker: int) -> tuple[int, ...]:
    """Return alleles at observed positions up to marker, read right-to-left.

    PBWT sorts haplotypes by reverse prefixes. For a small teaching example,
    we can construct the reverse-prefix key directly.
    """
    prefix_positions = [pos for pos in positions if pos <= marker]
    return tuple(haplotype[pos] for pos in reversed(prefix_positions))


def target_reverse_prefix_key(
    target: list[MaybeAllele],
    positions: list[int],
    marker: int,
) -> tuple[int, ...]:
    prefix_positions = [pos for pos in positions if pos <= marker]
    return tuple(int(target[pos]) for pos in reversed(prefix_positions))


def choose_neighbours(sorted_items: list[tuple[tuple[int, ...], int]], insert_at: int, count: int) -> list[int]:
    """Select nearby haplotypes around the target position in prefix order."""
    chosen: list[int] = []
    left = insert_at - 1
    right = insert_at

    while len(chosen) < count and (left >= 0 or right < len(sorted_items)):
        if left >= 0:
            chosen.append(sorted_items[left][1])
            left -= 1
            if len(chosen) == count:
                break

        if right < len(sorted_items):
            chosen.append(sorted_items[right][1])
            right += 1

    return chosen


def pbwt_inspired_state_selection(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    observed_positions: list[int],
    neighbours_per_marker: int = 4,
) -> list[int]:
    """Select locally similar reference haplotypes.

    Real IMPUTE5 inserts/locates the target haplotype inside the PBWT of the
    reference panel and selects nearby states at sparse selection markers.
    This readable version explicitly sorts reverse-prefix keys at each observed
    marker and takes neighbouring reference haplotypes.
    """
    selected: list[int] = []

    for marker in observed_positions:
        target_key = target_reverse_prefix_key(target, observed_positions, marker)

        sorted_items = sorted(
            (
                reverse_prefix_key(ref_hap, observed_positions, marker),
                ref_index,
            )
            for ref_index, ref_hap in enumerate(reference_panel)
        )

        keys_only = [item[0] for item in sorted_items]
        insert_at = bisect_left(keys_only, target_key)
        selected.extend(choose_neighbours(sorted_items, insert_at, neighbours_per_marker))

    # Keep order stable while removing duplicates.
    unique_selected: list[int] = []
    seen: set[int] = set()
    for state in selected:
        if state not in seen:
            unique_selected.append(state)
            seen.add(state)

    return unique_selected


def emission_probability(
    observed_allele: MaybeAllele,
    reference_allele: Allele,
    error_rate: float,
) -> float:
    """Probability of observing target allele given copied reference state."""
    if observed_allele is None:
        return 1.0
    return 1.0 - error_rate if observed_allele == reference_allele else error_rate


def transition_probability(
    previous_state: int,
    current_state: int,
    number_of_states: int,
    recombination_rate: float,
) -> float:
    """Simple HMM transition model.

    Staying in the same copying state is likely. Switching states represents a
    recombination event.
    """
    if number_of_states == 1:
        return 1.0
    if previous_state == current_state:
        return 1.0 - recombination_rate
    return recombination_rate / (number_of_states - 1)


def normalize(values: Iterable[float]) -> list[float]:
    values = list(values)
    total = sum(values)
    if total == 0.0:
        return [1.0 / len(values)] * len(values)
    return [value / total for value in values]


def forward_algorithm(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    error_rate: float,
    recombination_rate: float,
) -> list[list[float]]:
    marker_count = len(target)
    state_count = len(selected_states)
    alpha: list[list[float]] = [[0.0] * state_count for _ in range(marker_count)]

    initial_probability = 1.0 / state_count
    alpha[0] = normalize(
        initial_probability
        * emission_probability(target[0], reference_panel[ref_index][0], error_rate)
        for ref_index in selected_states
    )

    for marker in range(1, marker_count):
        row = []
        for current_state_position, ref_index in enumerate(selected_states):
            incoming_probability = sum(
                alpha[marker - 1][previous_state_position]
                * transition_probability(
                    previous_state_position,
                    current_state_position,
                    state_count,
                    recombination_rate,
                )
                for previous_state_position in range(state_count)
            )
            row.append(
                incoming_probability
                * emission_probability(target[marker], reference_panel[ref_index][marker], error_rate)
            )
        alpha[marker] = normalize(row)

    return alpha


def backward_algorithm(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    error_rate: float,
    recombination_rate: float,
) -> list[list[float]]:
    marker_count = len(target)
    state_count = len(selected_states)
    beta: list[list[float]] = [[0.0] * state_count for _ in range(marker_count)]
    beta[-1] = [1.0] * state_count

    for marker in range(marker_count - 2, -1, -1):
        row = []
        for current_state_position in range(state_count):
            probability = sum(
                transition_probability(
                    current_state_position,
                    next_state_position,
                    state_count,
                    recombination_rate,
                )
                * emission_probability(
                    target[marker + 1],
                    reference_panel[next_ref_index][marker + 1],
                    error_rate,
                )
                * beta[marker + 1][next_state_position]
                for next_state_position, next_ref_index in enumerate(selected_states)
            )
            row.append(probability)
        beta[marker] = normalize(row)

    return beta


def posterior_probabilities(alpha: list[list[float]], beta: list[list[float]]) -> list[list[float]]:
    posterior: list[list[float]] = []

    for marker in range(len(alpha)):
        posterior.append(
            normalize(
                alpha[marker][state] * beta[marker][state]
                for state in range(len(alpha[marker]))
            )
        )

    return posterior


def impute_missing_alleles(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    posterior: list[list[float]],
) -> list[Optional[float]]:
    """Return P(allele = 1) for missing markers; observed markers stay None."""
    imputed: list[Optional[float]] = []

    for marker, observed_allele in enumerate(target):
        if observed_allele is not None:
            imputed.append(None)
            continue

        probability_allele_1 = sum(
            posterior[marker][state_position]
            for state_position, ref_index in enumerate(selected_states)
            if reference_panel[ref_index][marker] == 1
        )
        imputed.append(probability_allele_1)

    return imputed


def run_imputation(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    neighbours_per_marker: int = 4,
    error_rate: float = 0.01,
    recombination_rate: float = 0.08,
) -> ImputationResult:
    observed_positions = [index for index, allele in enumerate(target) if allele is not None]

    selected_states = pbwt_inspired_state_selection(
        reference_panel,
        target,
        observed_positions,
        neighbours_per_marker=neighbours_per_marker,
    )

    alpha = forward_algorithm(
        reference_panel,
        target,
        selected_states,
        error_rate=error_rate,
        recombination_rate=recombination_rate,
    )
    beta = backward_algorithm(
        reference_panel,
        target,
        selected_states,
        error_rate=error_rate,
        recombination_rate=recombination_rate,
    )
    posterior = posterior_probabilities(alpha, beta)
    imputed = impute_missing_alleles(reference_panel, target, selected_states, posterior)

    return ImputationResult(
        selected_states=selected_states,
        posterior_by_marker=posterior,
        imputed_prob_allele_1=imputed,
    )


def print_reference_panel(reference_panel: list[list[Allele]]) -> None:
    print("Reference haplotypes:")
    for index, haplotype in enumerate(reference_panel):
        print(f"  h{index}: {haplotype}")


def print_result(target: list[MaybeAllele], result: ImputationResult) -> None:
    print("\nTarget haplotype:")
    print(f"  t : {target}")

    print("\nPBWT-inspired selected copying states:")
    print(f"  {['h' + str(index) for index in result.selected_states]}")

    print("\nImputed missing markers:")
    for marker, probability in enumerate(result.imputed_prob_allele_1):
        if probability is None:
            continue
        predicted = 1 if probability >= 0.5 else 0
        print(
            f"  marker {marker}: P(allele=1) = {probability:.3f}, "
            f"predicted allele = {predicted}"
        )

    completed = [
        target[marker] if target[marker] is not None
        else round(float(result.imputed_prob_allele_1[marker]), 3)
        for marker in range(len(target))
    ]
    print("\nCompleted haplotype, using probabilities at missing markers:")
    print(f"  {completed}")


def main() -> None:
    # Rows are reference haplotypes, columns are genetic markers/SNPs.
    # Alleles are encoded as 0/1 for a biallelic variant.
    reference_panel = [
        [0, 1, 0, 1, 1, 0, 0, 1, 0, 1],
        [0, 1, 0, 1, 0, 0, 0, 1, 0, 1],
        [0, 1, 1, 1, 1, 0, 1, 1, 0, 1],
        [1, 0, 1, 0, 0, 1, 1, 0, 1, 0],
        [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
        [0, 1, 0, 1, 1, 0, 0, 0, 0, 1],
        [1, 1, 0, 1, 1, 0, 0, 1, 1, 1],
        [0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
    ]

    # None means this marker was not observed in the target sample and should
    # be imputed from the reference panel.
    target = [0, 1, None, 1, None, 0, None, 1, None, 1]

    print_reference_panel(reference_panel)
    result = run_imputation(reference_panel, target)
    print_result(target, result)


if __name__ == "__main__":
    main()

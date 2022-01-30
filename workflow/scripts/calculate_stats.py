import pandas as pd
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from multiprocessing import Pool
import numpy as np


def build_truth_trees():
    truth = pd.read_csv(snakemake.input.truth, sep="\t", comment="#")
    trees_truth = defaultdict(IntervalTree)
    for i, row in truth.iterrows():
        tree = trees_truth[row["chrom"]]
        tree.addi(row["start"], row["stop"], snakemake.wildcards.coverage)
    return trees_truth


def build_found_trees(fdr=1.0):
    import numpy as np
    from functools import reduce

    found = pd.read_csv(snakemake.input.table, sep=",", comment="#").sort_values(
        by=["EVENT"]
    )
    threshold = np.log(1.0 - fdr)
    # print(threshold)
    trees_found = defaultdict(IntervalTree)
    num_below_fdr = 0
    num_events = 0
    skipped = []
    for event, group in found.groupby("EVENT"):
        num_events += 1
        logprob_present = np.log(group["PROB_PRESENT"].mean())
        logprob_absent = np.log(group["PROB_ABSENT"].mean())
        logprob_artifact = np.log(group["PROB_ARTIFACT"].mean())
        logprobs = [logprob_present]
        logprob_sum = reduce(np.logaddexp, logprobs)
        if logprob_sum < threshold:
            num_below_fdr += 1
            continue
        if len(group) < 2:
            skipped.append(group)
            continue
        first = group.iloc[0]
        second = group.iloc[1]
        if len(group) != 2:
            if len(set(group["CHROM"]) | set(group["OTHER_CHROM"])) == 1:
                first = group.iloc[0]
                second = group.iloc[-1]
            else:
                # TODO handle non-simple events
                skipped.append(group)
                continue

        assert first["CHROM"] == second["CHROM"]
        start, stop = min(first["start"], second["start"]), max(
            first["stop"], second["stop"]
        )
        start, stop = min(start, stop), max(start, stop)
        chrom = first["CHROM"]
        tree = trees_found[chrom]
        if start == stop:
            stop = start + 1
        tree.addi(
            start,
            stop,
            (
                first["PROB_PRESENT"],
                first["PROB_ABSENT"],
                first["PROB_ARTIFACT"],
            ),
        )
    # print(f"below fdr of {fdr}: {num_below_fdr}")
    # print(f"skipped: {len(skipped)}")
    # print(f"total number of events: {num_events}")
    return trees_found, skipped


def calc_stats(fdr=1.0):
    trees_truth = build_truth_trees()
    trees_found, skipped = build_found_trees(fdr=fdr)
    chroms = sorted(set(trees_truth.keys()) | set(trees_found.keys()))

    entries = []
    empty_probs = [None, None, None]
    num_entries_truth = 0
    num_entries_found = 0
    for chrom in chroms:
        tree_truth, tree_found = trees_truth[chrom], trees_found[chrom]
        for entry in tree_truth:
            num_entries_truth += 1
            loc = [chrom, entry.begin, entry.end]
            length = entry.end - entry.begin
            coverage = entry.data
            matches = tree_found.overlap(entry.begin, entry.end)
            if matches:
                if len(matches) > 1:
                    # print(len(matches))
                    pass
                for match in matches:
                    assert len(match.data) == 3
                    probs = list(match.data)
                    entries.append(loc + [length, True, True, coverage] + probs)
                    break
            else:
                entries.append(loc + [length, True, False, coverage] + empty_probs)
        for entry in tree_found:
            num_entries_found += 1
            loc = [chrom, entry.begin, entry.end]
            length = entry.end - entry.begin
            probs = list(entry.data)
            matches = tree_truth.overlap(entry.begin, entry.end)
            if matches:
                if len(matches) > 1:
                    # print(len(matches))
                    pass
                for match in matches:
                    loc = [chrom, match.begin, match.end]
                    length = match.end - match.begin
                    coverage = match.data
                    entries.append(loc + [length, True, True, coverage] + probs)
                    break
            else:
                entries.append(
                    loc + [length, False, True, snakemake.wildcards.coverage] + probs
                )
    # print(num_entries_truth, num_entries_found)
    for group in skipped:
        loc = ["*", -1, -1]
        entries.append(
            loc
            + [-1, False, True, snakemake.wildcards.coverage]
            + list(group[["PROB_PRESENT", "PROB_ABSENT", "PROB_ARTIFACT"]].mean())
        )

    # print("len(entries):", len(entries))
    stats = pd.DataFrame(
        data=entries,
        columns=[
            "chrom",
            "start",
            "stop",
            "length",
            "in_truth_dataset",
            "in_found_dataset",
            "coverage",
            "PROB_PRESENT",
            "PROB_ABSENT",
            "PROB_ARTIFACT",
        ],
    )
    # stats.set_index(["chrom", "start", "stop"], inplace=True)
    stats.sort_values(by=list(stats.columns), inplace=True)
    stats.drop_duplicates(
        subset=[
            "chrom",
            "start",
            "stop",
            "length",
            "in_truth_dataset",
            "in_found_dataset",
            # "coverage",
        ],
        keep="first",
        inplace=True,
    )

    stats["in_both_datasets"] = stats["in_truth_dataset"] & stats["in_found_dataset"]

    stats["true_positive"] = (
        stats["in_truth_dataset"] & stats["in_found_dataset"]
    ).astype(int)
    stats["true_negative"] = (
        ~stats["in_found_dataset"] & ~stats["in_truth_dataset"]
    ).astype(int)
    stats["false_positive"] = (
        stats["in_found_dataset"] & ~stats["in_truth_dataset"]
    ).astype(int)
    stats["false_negative"] = (
        ~stats["in_found_dataset"] & stats["in_truth_dataset"]
    ).astype(int)

    return stats


global summary  # workaround for multiprocessing


def summary(fdr):
    stats = calc_stats(fdr=fdr)

    grouped = stats.groupby(["coverage"])

    precision_by_coverage = grouped["true_positive"].sum() / (
        grouped["true_positive"].sum() + grouped["false_positive"].sum()
    )
    recall_by_coverage = grouped["true_positive"].sum() / (
        grouped["true_positive"].sum() + grouped["false_negative"].sum()
    )
    specificity_by_coverage = grouped["true_negative"].sum() / (
        grouped["true_negative"].sum() + grouped["false_positive"].sum()
    )

    num_hits_by_coverage = grouped.count()["chrom"]
    stats_by_coverage = pd.concat(
        [
            precision_by_coverage,
            recall_by_coverage,
            specificity_by_coverage,
            num_hits_by_coverage,
        ],
        axis=1,
    )
    stats_by_coverage.columns = [
        "precision",
        "recall",
        "specificity",
        "count",
    ]
    stats_by_coverage["fdr"] = fdr
    return stats_by_coverage


with Pool(22) as pool:
    by_fdr = pd.concat(pool.map(summary, np.linspace(0, 1, 21))).reset_index()

del summary  # workaround for multiprocessing

g = by_fdr.groupby(["fdr"])


def nanmean(a, weights):
    indices = ~np.isnan(a)
    return np.average(a[indices], weights=weights[indices])


precision_by_fdr = g.apply(
    lambda x: nanmean(x["precision"], weights=x["count"])
).rename("precision")
recall_by_fdr = g.apply(lambda x: nanmean(x["recall"], weights=x["count"])).rename(
    "recall"
)

by_fdr.to_csv(snakemake.output.stats, index=False, sep="\t")

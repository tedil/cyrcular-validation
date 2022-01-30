rule precision_recall_plots:
    input:
        stats=expand(
            "results/validation/stats/{model}/{coverage}.tsv",
            allow_missing=True,
            coverage=config["simulation"]["target_coverages"],
        ),
    output:
        plot0="results/validation/plots/{model}.0.{ext}",
        plot1="results/validation/plots/{model}.1.{ext}",
        plot2="results/validation/plots/{model}.2.{ext}",
    log:
        "../logs/plot_precision_recall_{model}_{ext}.log",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_precision_recall.py"


rule stats:
    input:
        truth=config["simulation"]["regions"],
        table="results/calling/tables/{model}_{coverage}.sorted.csv",
    output:
        stats="results/validation/stats/{model}/{coverage}.tsv",
    log:
        "../logs/calculate_stats_{model}_{coverage}.log",
    conda:
        "../envs/calculate_stats.yaml"
    script:
        "../scripts/calculate_stats.py"

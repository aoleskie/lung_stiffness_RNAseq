import pandas as pd
import os

RESULTS   = config["results_dir"]
RESOURCES = config["resources_dir"]

samples = (
    pd.read_csv(config["samples"], sep="\t", comment="#")
      .set_index("sample_id", drop=False)
)
SAMPLES        = samples.index.tolist()
MATRIX_SAMPLES = samples.index[samples["input_type"] == "matrix"].tolist()

# Alignment-lesson targets (STARsolo). Empty list when none configured.
LESSON = config.get("fastq_lesson_samples") or []
LESSON_TARGETS = expand(
    f"{RESULTS}/quant/{{s}}/Solo.out/Gene/filtered/matrix.mtx.gz", s=LESSON
)

def sample_path(wildcards):
    return samples.loc[wildcards.sample, "path"]

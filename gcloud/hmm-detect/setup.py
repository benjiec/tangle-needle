import os
import argparse
from pathlib import Path
from tangle import unique_batch
from heap.gcloud import split_fasta, from_template, GCHMMHelper

PREP_GENERAL_TEMPLATE = "prep-general.sh.template"
PREP_GENERAL_SCRIPT = "prep-general.sh"
PREP_TASK_TEMPLATE = "prep-task.sh.template"
PREP_TASK_SCRIPT = "prep-task.sh"
JOB_TEMPLATE = "job.json.template"
JOB_JSON = "job.json"
INSTRUCTION_TEMPLATE = "instruction.template"
INSTRUCTION_FILE = "README"


if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-accession", required=True)
    ap.add_argument("query_fna")
    ap.add_argument("--entries-per-task", type=int, default=50)
    ap.add_argument("--parallelism", type=int, default=45)
    ap.add_argument("--run-dir-parent", default=".")
    args = ap.parse_args()

    gc = GCHMMHelper(args.run_dir_parent)

    fn_prefix = "input.fna."
    script_dir = Path(__file__).resolve().parent

    gc_input_path_pre_index = f"{gc.gc_run_input_dir}/{fn_prefix}"
    fasta_files = split_fasta(args.query_fna, args.entries_per_task, gc.host_run_input_dir, fn_prefix)

    gc.setenv(
      PARALLELISM=args.parallelism,
      NTASKS=len(fasta_files),
      GC_INPUT_PATH_PRE_INDEX=gc_input_path_pre_index,
      GENOME_ACCESSION=args.genome_accession
    )

    print(gc.host_run_dir)
    gc.instantiate_template(script_dir / PREP_GENERAL_TEMPLATE, PREP_GENERAL_SCRIPT)
    gc.instantiate_template(script_dir / PREP_TASK_TEMPLATE, PREP_TASK_SCRIPT)
    gc.instantiate_template(script_dir / JOB_TEMPLATE, JOB_JSON)
    gc.instantiate_template(script_dir / INSTRUCTION_TEMPLATE, INSTRUCTION_FILE)

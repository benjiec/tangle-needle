import os
import argparse
from pathlib import Path
from tangle import unique_batch
from tangle.sequence import write_fasta_from_dict, read_fasta_as_dict
from heap.gcloud import from_template, GCHMMHelper
from needle.sequence import split_sequence_dictionary

PREP_GENERAL_TEMPLATE = "prep-general.sh.template"
PREP_GENERAL_SCRIPT = "prep-general.sh"
PREP_TASK_TEMPLATE = "prep-task.sh.template"
PREP_TASK_SCRIPT = "prep-task.sh"
JOB_TEMPLATE = "job.json.template"
JOB_JSON = "job.json"
INSTRUCTION_TEMPLATE = "instruction.template"
INSTRUCTION_FILE = "README"

TARGET_SEQLEN = 2000000

if __name__ == "__main__":

    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-accession", required=True)
    ap.add_argument("query_fna")
    ap.add_argument("--parallelism", type=int, default=45)
    ap.add_argument("--run-dir-parent", default=".")
    args = ap.parse_args()

    gc = GCHMMHelper(args.run_dir_parent)

    fn_prefix = "input.fna."
    script_dir = Path(__file__).resolve().parent
    gc_input_path_pre_index = f"{gc.gc_run_input_dir}/{fn_prefix}"

    full_sequence_dict = read_fasta_as_dict(args.query_fna)
    splitted_dicts = split_sequence_dictionary(full_sequence_dict, TARGET_SEQLEN)
    for i,splitted in enumerate(splitted_dicts):
        write_fasta_from_dict(splitted, str(Path(gc.host_run_input_dir) / fn_prefix)+str(i))

    gc.setenv(
      PARALLELISM=args.parallelism,
      NTASKS=len(splitted_dicts),
      GC_INPUT_PATH_PRE_INDEX=gc_input_path_pre_index,
      GENOME_ACCESSION=args.genome_accession
    )

    print(gc.host_run_dir)
    gc.instantiate_template(script_dir / PREP_GENERAL_TEMPLATE, PREP_GENERAL_SCRIPT)
    gc.instantiate_template(script_dir / PREP_TASK_TEMPLATE, PREP_TASK_SCRIPT)
    gc.instantiate_template(script_dir / JOB_TEMPLATE, JOB_JSON)
    gc.instantiate_template(script_dir / INSTRUCTION_TEMPLATE, INSTRUCTION_FILE)

import os
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd
from loguru import logger
from tqdm import tqdm

logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")

info = pd.read_csv("data/info_with_type.csv")
pdbids = info.pdb_id.tolist()

tmp_dir = Path("tmp")
if tmp_dir.exists():
    subprocess.run(["rm", "-rf", str(tmp_dir)])
tmp_dir.mkdir(exist_ok=True)

MMSEQS_PATH = "mmseqs"


def run_command(cmd):
    logger.info(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode == 0:
        return True
    else:
        logger.error(
            f"Command failed: {' '.join(cmd)}. Error: {result.stderr.decode()}"
        )
        return False


os.makedirs("protein_seq_identities", exist_ok=True)

protein_types = info["type"].unique()

for my_type in tqdm(protein_types, desc="Processing protein types"):
    logger.info(f"Processing protein type: {my_type}")
    type_df = info[info["type"] == my_type]

    n = len(type_df)
    if n < 2:
        logger.warning(f"Not enough sequences for type {my_type}. Skipping.")
        continue

    fasta_file = tmp_dir / f"{my_type}.fasta"
    db_name = tmp_dir / f"{my_type}_db"
    result_db = tmp_dir / f"{my_type}_result"
    result_tsv = tmp_dir / f"{my_type}_result.tsv"
    tmp_work_dir = tmp_dir / f"work_{my_type}"
    tmp_work_dir.mkdir(exist_ok=True)

    with open(fasta_file, "w") as f:
        for idx, row in type_df.iterrows():
            seq_id = row["pdb_id"]  # Используем pdbid как идентификатор
            sequence = row["seq"]
            f.write(f">{seq_id}\n{sequence}\n")
    logger.info(f"FASTA file created: {fasta_file}")

    if not run_command([MMSEQS_PATH, "createdb", str(fasta_file), str(db_name)]):
        continue

    if not run_command(
        [
            MMSEQS_PATH,
            "search",
            str(db_name),
            str(db_name),
            str(result_db),
            str(tmp_work_dir),
            "--threads",
            "112",
            "--alignment-mode",
            "3",
        ]
    ):
        continue

    if not run_command(
        [
            MMSEQS_PATH,
            "convertalis",
            str(db_name),
            str(db_name),
            str(result_db),
            str(result_tsv),
            "--format-output",
            "query,target,pident",
        ]
    ):
        continue

    if not result_tsv.exists():
        logger.error(f"Error: {result_tsv} not found.")
        continue

    logger.info(f"Reading results for {my_type}")
    results = pd.read_csv(result_tsv, sep="\t", names=["query", "target", "pident"])

    result_dict = {}
    indices = type_df.pdb_id.tolist()
    index_mapping = {pdbid: idx for idx, pdbid in enumerate(indices)}

    for _, row in results.iterrows():
        i = index_mapping.get(row["query"])
        j = index_mapping.get(row["target"])
        if i is not None and j is not None and i < j:
            identity = row["pident"] / 100.0
            result_dict[(i, j)] = identity

    with open(f"protein_seq_identities/{my_type}.pkl", "wb") as f:
        pickle.dump((indices, result_dict), f)
    logger.info(f"Results saved for {my_type}")

    os.remove(fasta_file)
    subprocess.run(["rm", "-rf", db_name, result_db, tmp_work_dir])

protein_identities = {}
for item in os.listdir("protein_seq_identities"):
    if item.endswith(".pkl"):
        with open(f"protein_seq_identities/{item}", "rb") as f:
            data = pickle.load(f)
        item_name = item.split(".")[0]
        protein_identities[item_name] = data

with open("data/protein_identities_mmseqs.pkl", "wb") as f:
    pickle.dump(protein_identities, f)

logger.info("Finished creating protein_identities_mmseqs.pkl")

if tmp_dir.exists():
    subprocess.run(["rm", "-rf", str(tmp_dir)])
    subprocess.run(["rm", "-rf", "protein_seq_identities"])

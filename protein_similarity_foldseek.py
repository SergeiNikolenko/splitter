import os
import pickle  # nosec
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from loguru import logger
from tqdm import tqdm

logger.add(sys.stdout, format="{message}", level="INFO")

info = pd.read_csv("data/info_with_type.csv", index_col=0)
pdb_dir = Path("/mnt/ligandpro/db/LPCE/final")
tmp_dir = Path("tmp")

if tmp_dir.exists():
    shutil.rmtree(tmp_dir)
tmp_dir.mkdir(exist_ok=True)

FOLDSEEK_PATH = "foldseek"


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
    tqdm.write(f"Processing protein type: {my_type} (Foldseek)")
    type_df = info[info["type"] == my_type]

    pdb_files = []
    for pdbid in type_df.index:
        pdb_file = pdb_dir / f"{pdbid}.pdb"  # Берем файлы только с цепью A
        if (
            pdb_file.exists() and pdb_file.stat().st_size > 0
        ):  # Проверяем, что файл существует и не пуст
            pdb_files.append(pdb_file)
        else:
            logger.warning(f"PDB file {pdb_file} not found or is empty, skipping.")

    if not pdb_files:
        logger.warning(f"No valid PDB files found for type {my_type}. Skipping.")
        continue

    pdb_list_file = tmp_dir / f"{my_type}_pdb_list.txt"
    with open(pdb_list_file, "w") as f:
        for pdb_file in pdb_files:
            f.write(f"{pdb_file}\n")
    tqdm.write(f"PDB list file created: {pdb_list_file}")

    result_tsv = tmp_dir / f"{my_type}_result.m8"
    tmp_work_dir = tmp_dir / f"work_{my_type}"
    tmp_work_dir.mkdir(exist_ok=True)

    if not run_command(
        [
            FOLDSEEK_PATH,
            "easy-search",
            str(pdb_list_file),
            str(pdb_list_file),
            str(result_tsv),
            str(tmp_work_dir),
        ]
    ):
        continue

    if not result_tsv.exists():
        logger.error(f"Error: {result_tsv} not found.")
        continue

    tqdm.write(f"Reading results for {my_type}")

    # Обработка результатов
    results = pd.read_csv(
        result_tsv,
        sep="\t",
        names=[
            "query",
            "target",
            "fident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "evalue",
            "bits",
        ],
    )

    result_dict = {}
    indices = type_df.index.tolist()
    index_mapping = {pdbid: idx for idx, pdbid in enumerate(indices)}

    for _, row in results.iterrows():
        i = index_mapping.get(row["query"])
        j = index_mapping.get(row["target"])
        if i is not None and j is not None and i < j:
            identity = row["fident"] / 100.0
            result_dict[(i, j)] = identity

    with open(f"protein_seq_identities/{my_type}_foldseek.pkl", "wb") as f:
        pickle.dump((indices, result_dict), f)  # nosec
    tqdm.write(f"Results saved for {my_type}")

    os.remove(pdb_list_file)
    subprocess.run(["rm", "-rf", tmp_work_dir])

# Создание итогового файла
protein_identities = {}
for item in os.listdir("protein_seq_identities"):
    if item.endswith("_foldseek.pkl"):
        with open(f"protein_seq_identities/{item}", "rb") as f:
            data = pickle.load(f)  # nosec
        item_name = item.split(".")[0]
        protein_identities[item_name] = data

with open("data/protein_identities_foldseek.pkl", "wb") as f:
    pickle.dump(protein_identities, f)  # nosec

logger.info("Finished creating protein_identities_foldseek.pkl")

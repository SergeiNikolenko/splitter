import functools
import pickle  # nosec
import sys
from typing import Any

import datamol as dm
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from loguru import logger
from rdkit import RDLogger
from scipy.spatial import distance
from sklearn.metrics import pairwise_distances
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")


logger.remove()
logger.add(sys.stdout, format="{message}", level="INFO")


def custom_pdist(
    mols: list[str | dm.Mol],
    n_jobs: int | None = -1,
    squareform: bool = True,
    fold_size: int = 1024,
    metric: str = "jaccard",
    **fp_args: Any,
) -> np.ndarray:
    """Computes the pairwise distance between the fingerprints of all molecules in the
    input set using the specified metric.

    Args:
        mols: List of molecules.
        n_jobs: Number of jobs for parallelization.
                Set to 1 for no parallelization.
                Set to -1 to use all available cores.
        squareform: Whether to return in square form (matrix) or in
                    condensed form (1D vector).
        fold_size: The size to which fingerprints will be folded.
        metric: Metric to use for distance calculation ('jaccard' or 'dice').
        **fp_args: Additional arguments to pass to to_fp().

    Returns:
        dist_mat: Pairwise distance matrix.
    """

    logger.info("Calculating fingerprints for molecules...")
    fps = dm.parallelized(
        functools.partial(dm.to_fp, as_array=True, fold_size=None, **fp_args),
        mols,
        n_jobs=n_jobs,
    )
    fps_array = [fp for fp in fps if fp is not None]
    n = len(fps_array)
    if n == 0:
        logger.warning(
            "No valid fingerprints generated. Returning empty distance matrix."
        )
        return np.array([])
    logger.info(f"Total valid fingerprints: {n}")
    logger.info("Converting fingerprints to numpy array...")
    fps_array = np.array(fps_array, dtype=bool)
    logger.info("Calculating pairwise distances using sklearn's pairwise_distances...")
    if metric not in ["jaccard", "dice"]:
        raise ValueError("Unknown metric. Use 'jaccard' or 'dice'.")
    dist_mat = pairwise_distances(fps_array, metric=metric, n_jobs=n_jobs)
    return dist_mat if squareform else distance.squareform(dist_mat, force="tomatrix")


def calculate_similarity_matrix_with_checks(
    mols: list[dm.Mol], smiles_list: list[str], radius=2, fp_size=1024, metric="jaccard"
) -> np.ndarray:
    """Computes the similarity matrix with checks for identical molecules.

    Args:
        mols: List of molecules.
        radius: Radius for fingerprint calculation.
        fp_size: Size of the fingerprints.
        metric: Metric to use for distance calculation ('jaccard' or 'dice').

    Returns:
        sim_matrix: Similarity matrix.
    """
    logger.info("Calculating distance matrix")
    dist_mat = custom_pdist(
        mols,
        n_jobs=-1,
        squareform=True,
        radius=radius,
        fold_size=fp_size,
        metric=metric,
    )
    if metric in ["jaccard", "dice"]:
        sim_matrix = 1 - dist_mat
    else:
        raise ValueError("Unknown metric. Use 'jaccard' or 'dice'.")
    logger.info("Checking for identical molecules in the similarity matrix")
    for i in tqdm(range(len(mols)), desc="Processing molecules"):
        for j in range(i + 1, len(mols)):
            if sim_matrix[i, j] == 1:
                if smiles_list[i] == smiles_list[j]:
                    sim_matrix[i, j] = 1.0
                else:
                    new_fpsi = dm.to_fp(
                        mols[i], fp_type="ecfp", as_array=True, radius=4, fpSize=fp_size
                    )
                    new_fpsj = dm.to_fp(
                        mols[j], fp_type="ecfp", as_array=True, radius=4, fpSize=fp_size
                    )
                    sim = (
                        1 - distance.jaccard(new_fpsi, new_fpsj)
                        if metric == "jaccard"
                        else 1 - distance.dice(new_fpsi, new_fpsj)
                    )
                    if sim == 1:
                        new_fpsi = dm.to_fp(
                            mols[i],
                            fp_type="ecfp",
                            as_array=True,
                            radius=10,
                            fpSize=fp_size,
                        )
                        new_fpsj = dm.to_fp(
                            mols[j],
                            fp_type="ecfp",
                            as_array=True,
                            radius=10,
                            fpSize=fp_size,
                        )
                        sim = (
                            1 - distance.jaccard(new_fpsi, new_fpsj)
                            if metric == "jaccard"
                            else 1 - distance.dice(new_fpsi, new_fpsj)
                        )
                        if sim == 1 and smiles_list[i] != smiles_list[j]:
                            sim_matrix[i, j] = 0.99
            sim_matrix[j, i] = sim_matrix[i, j]
    logger.info("Similarity matrix calculation complete")
    return sim_matrix


def save_similarity_matrix(
    sim_matrix: np.ndarray, valid_indices: list[int], output_path: str
):
    """Saves the similarity matrix and associated valid indices to a file.

    Args:
        sim_matrix: Similarity matrix to be saved.
        valid_indices: List of valid indices corresponding to the rows/columns in the matrix.
        output_path: Path to the output file.
    """
    logger.info("Saving similarity matrix")
    rounded_matrix = np.round(sim_matrix, 2)
    with open(output_path, "wb") as f:
        pickle.dump({"sim_matrix": rounded_matrix, "valid_indices": valid_indices}, f)  # nosec
    matrix_size_bytes = rounded_matrix.nbytes
    if matrix_size_bytes < 1024:
        matrix_size = f"{matrix_size_bytes} bytes"
    elif matrix_size_bytes < 1024**2:
        matrix_size = f"{matrix_size_bytes / 1024:.2f} KB"
    elif matrix_size_bytes < 1024**3:
        matrix_size = f"{matrix_size_bytes / 1024**2:.2f} MB"
    else:
        matrix_size = f"{matrix_size_bytes / 1024**3:.2f} GB"
    logger.info(
        f"Rounded similarity matrix and indices saved to {output_path}, size: {matrix_size}"
    )


if __name__ == "__main__":
    input_file = "data/info_with_type.csv"
    metric = "dice"  # 'dice' or 'jaccard'
    output_file = f"data/ligand_similarity_{metric}.pkl"

    info = pd.read_csv(input_file)
    info = info.iloc[:1000]
    info.dropna(subset=["smiles"], inplace=True)

    mols = Parallel(n_jobs=-1)(
        delayed(dm.to_mol)(smi)
        for smi in tqdm(info["smiles"], desc="Converting SMILES to molecules")
    )

    mols = [mol for mol in mols if mol is not None]
    valid_indices = info.index[: len(mols)].tolist()

    smiles_list = info["smiles"][: len(mols)].tolist()

    sim_matrix = calculate_similarity_matrix_with_checks(
        mols, smiles_list, metric=metric
    )

    save_similarity_matrix(sim_matrix, valid_indices, output_file)

.PHONY: protein_similarity_mmsec protein_similarity_foldseek ligand_similarity test

protein_similarity_mmsec:
	python protein_similarity_mmsec.py

protein_similarity_foldseek:
	python protein_similarity_foldseek.py

ligand_similarity:
	clear
	python ligand_similarity.py

test:
	clear
	python ligand_similarity_test.py


all: 
	clear
	python 2_data_processing.py
	python 3_ligand_similarity_test.py
	python 4_ligand_similarity.py
	python 5_protein_similarity_foldseek.py
	python 6_protein_similarity_mmsec.py

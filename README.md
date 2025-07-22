# Mirror-Image-Glycans
This repo contains code to predict the interactions of protein and glycans using existing ML models.

To get LectinOracle predictions for the full dataset (100 proteins and 611 glycans), run the following command:
```
python predict.py --config configs/full_data.yaml
```

To run with different glycans, change the 'glycan_path' parameter in the config file.

To run with different proteins, change the 'protein_path' and 'fasta_path' parameters in the config file. For every protein in the 'protein_path' file, the corresponding amino acid sequence will be looked for in the 'fasta_path' file (using its' UniProtID). If the protein cannot be matched to a sequence (its' UniProtID is not provided or cannot be found in the fasta file), then the protein will not be included in predictions.
# Using Multiple Drug Similarity Networks to Promote Adverse Drug Event Detection

### Instruction
- Download the Signal Scores and Drug Similarity Matrices from [here](https://buckeyemailosu-my.sharepoint.com/:f:/g/personal/padhi_3_buckeyemail_osu_edu/EsPkyFORsqNJgrxGmNZWYl0B0Qvmipgh-fCffpoCeRDT6w?e=1EsMaF).
- Copy the downloaed Signal_Scores and Similarity_Matrices directory into the cloned repository.

### Running example
```
python run.py --input Signal_Scores/FAERS_signal_data_BCPNN --method BCPNN --year 22 --quarter 2 --eval_metrics all --similarity go_mf --output output/results_BCPNN/go_mf/results_22Q2.csv --soc_output output/soc_results_BCPNN/go_mf/soc_results_22Q2.csv --split True
```

### Parameters
- --input, input original signal scores files. 
- --method, signal detection algorithm (i.e., MGPS, BCPNN).
- --year, year of the signal data file in 'YY' format.
- --quarter, quarter of signal data file (i.e., 1, 2, 3, 4).
- --eval_metrics, evaluation metrics (i.e., AUC, AUPR, Precision, Recall, etc.) choices=['all', 'specificity-sensitivity'].
- --similarity,  which drug similarity to use choices=['chem', 'atc', 'sw', 'go_bp', 'go_cc', 'go_mf'].
- --output, output file path for the result.
- --soc_output, output file path for the result of MedDRA SOC.
- --split, whether to split entire dataset into validation set and testing set (i.e., True/False). 
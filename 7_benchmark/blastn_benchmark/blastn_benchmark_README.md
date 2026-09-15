Process for making predictions with BLASTn\
Each subdirectory contains scripts to make BLAST predictions under different validation strategies\
All use the same basic scripts\
1: run the `cluster_seqs_to_drop.R` script to generate a list of clusters to drop\
2: use the seqkit utility to run the command in `seqkit_drop_seqs_cmd_blocked.txt` to drop these sequences\
3: generate custom BLAST database and query using the bash script contained in e.g. `generate_blastdb_run_blastn_cmd_blocked.txt` or `generate_blastdb_run_blastn_cmd_loo.txt`}
4: make predictions from BLAST results using the `blastn_host_predict.R` script


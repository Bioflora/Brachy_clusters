## Brachy_clusters

This repository contains code and data for phylogenetic analyses of different Brachypodium species.

### add_gene_seqs2clusters.pl

This script takes as input precomputed clusters of transcripts and produces enlarged clusters with sequences from species or longer gene sequences matching transcripts from previous species.

You must configure variables %blastdirs, %references, @order, %max_references_seqs before analyses, see examples in the script. 

The script is called as: 
```
perl add_gene_seqs2clusters.pl input_cluster_dir output_dir
```

This is how input clusters look like, note the abbreviated species name at the end of the FASTA header:
```
>c19543_g7_i3_Chr09_Bsta
ATGGAGGGCAAGGAGGAGGATGTCCGCCTGGGGGCGAACAAGTTCTCGGAGCGGCAGCCG
>Bradi3g51387.1_Bdis
ATGGAAGGGAAAGAGGAGGACGTGCGGCTGGGGGCGAACAAGTTCTCGGAGCGGCAGCCG
>c37694_g4_i1_chr4_Barb
ATGGAGGGGAAAGAGGAGGACGTGCGGCTGGGGGCGAACAAGTTCTCGGAGCGGCAGCCG
>c29451_g2_i9_chr9_Bpin
ATGGAGGGCAAGGAGGAGGACGTGCGCCTGGGCGCGAACCGCTACTCGGAGCGGCAGCCC
>Brasy-EspV1T35246_chr4_Bsyl
ATGGAGGGCAAGGAGGAGGATGTCCGCCTGGGGGCGAACAAGTTCTCGGAGCGGCAGCCG
>LOC_Os02g44630.1_Osat
ATGGAGGGGAAGGAGGAGGACGTGCGGCTGGGGGCGAACAGGTACTCGGAGAGGCAGCCG
>HORVU6Hr1G064140.2_Hvul
CTCGGAGCGGCAGCCC
>c35008_g2_i7_Bd3_Bhyb
ATGGAAGGGAAAGAGGAGGACGTGCGGCTGGGGGCGAACAAGTTCTCGGAGCGGCAGCCC
>c38520_g4_i1_Bd3_Bboi
CGTCTACTGCACCGCCGGCATCTCAGGAGGACACATCAACCCTGCAGTGACT
>c39171_g2_i3_Chr09_Bmex
ATGGAGGGCAAGGAGGAGGATGTCCGCCTGGGGGCGAACAAGTTCTCTGAGCGGCAGCCC
>c51910_g3_i1_chr4_Bpho
CGTCTACTGCACCGCCGGCATCTCAGGAGGACACATCAACCCTGCAGTGACT
```


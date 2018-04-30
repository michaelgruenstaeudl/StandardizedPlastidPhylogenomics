INP=rpl32-trnL_Alignment_2018.04.24.2100.fasta
OUP=${INP%.fasta*}.trimAl.fasta
trimal -in $INP -out $OUP -gt 0.8 -st 0.001 -cons 60

RepeatMasker -e rmblast -pa 30 -qq -lib TElib.clean.fa SV.fa

perl /home/jianianhua/miniconda3/share/EDTA/bin/count_base.pl INS.fa > INS.stats
perl /home/jianianhua/miniconda3/share/EDTA/bin/count_base.pl DEL.fa > DEL.stats

perl /home/jianianhua/miniconda3/share/EDTA/bin/buildSummary.pl -maxDiv 40 -stats INS.stats INS.out > INS.TE.sum 2>/dev/null
perl /home/jianianhua/miniconda3/share/EDTA/bin/buildSummary.pl -maxDiv 40 -stats DEL.stats DEL.out > DEL.TE.sum 2>/dev/null
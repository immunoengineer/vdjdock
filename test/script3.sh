target="H1N1" #-2007"

#H3N2-1997
for depth in 11 # 2 3 4 5 6 7
do
 
  for file in split-*txt # *top.top.txt *bot.top.txt *bot.tail.txt *top.tail.txt
  do
    bsub "../bin/epimorph.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target --centivax_min_conc=$depth > centivax-${target}-${file}-${depth}.txt"
  done
done

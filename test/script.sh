target="H1N1-2007"

#H3N2-1997

for file in split-*txt # *top.top.txt *bot.top.txt *bot.tail.txt *top.tail.txt
do
  bsub "../bin/epimorph.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target > ${target}-${file}"
done

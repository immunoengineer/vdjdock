# vdjdock
Computational structural immunology tool for analyzing immune interfaces

Usage: /data/apps/vdjdock/bin/epimorph.pl
  Standard options
  --epitopedb=epitopes.txt  database of epitopes, with each line containing
                            space separated list of positions that make up 
                            the epitope.
                            example: A:119 A:122 A:124 A:125 A:127 A:133 A:146 A:149 A:152 A:67
                                     B:125 B:147 B:150 B:151 B:160 B:161 C:4 C:5 D:11 D:9
                                     B:14 B:150 B:20 B:25 B:28 B:30 B:31 B:32 B:33 B:34
  --alignment=MSA.fa        a set of sequences aligned to the structural   
                            reference antigen. It is an MSA with every residue of
                            the reference structure present and no gaps relative to
                            the reference structure, such that MSA length = ref sequence length
                            example: >structure
                                     EDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAG
                                     >H1N1-2009
                                     EDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAG
  --numbering=number.txt    a mapping file assigning a structural number ID
                            to each position in the sequence of the structural
                            reference found in the MSA.fa. Format is three tabs
                            separated seqposition|AA|structuralnumber
 		               example: Residue_Number	Residue_AA	Structure_Number
 					0	E	A:4
 					1	D	A:5
 					2	Q	A:6
 					3	I	A:7
  --antigen=[state]         one of the sequence names as the antigen. then for each epitope
                            code will report all MSA variants that share that epitope with
                            antigen
  --sanity_check=1          [OPTIONAL] if antigen is set to reference structure sequence, it makes sure that
                            the numbering in numbering.txt always matches the alignment or returns
                            an error message. Note --antigen needs to be the structural ref for this!
  --centivax_min_conc=1     [OPTIONAL] minimum number of centivax members to share an epitope here for this epitope
                            to count as above threshold dose. If set to 3, then at least one epitope must be shared
                            between three members of centivax for this site to be counted as present in the mixture
  --forbidden_contact_list  a file that contains forbidden structIDs, one per line. For disallowing certain epitopes
                            (stem binding, etc)


Reference set of human antibody docking templates:
559 antibodies:
1ADQ,1AQK,1BEY,1DFB,1DN0,1DQL,1FGV,1G9M,1HZH,1IGA,1IGM,1IQD,1L7I,1NFD,1NL0,1Q1J,1RHH,1RZ7,1RZF,1RZG,1RZI,1T4K,1TJG,1TZI,1U6A,1W72,2A9N,2AGJ,2AJ3,2B0S,2B2X,2CMR,2D7T,2EIZ,2FGW,2FJH,2FL5,2H9G,2HFF,2HWZ,2J6E,2QQN,2QSC,2R0L,2R56,2UZI,2VXQ,2VXV,2X7L,2XA8,2XQB,2XRA,2XWT,2XZA,2YBR,3AAZ,3BQU,3EO0,3EO9,3EYF,3F12,3FN0,3FZU,3G04,3G6A,3G6J,3GBM,3GHE,3GO1,3H0T,3H42,3HC0,3HI5,3HMW,3IDX,3IYW,3J5M,3J6U,3JCB,3K2U,3KR3,3KYK,3L95,3LMJ,3LRS,3LZF,3M8O,3MLR,3MLW,3MLX,3N85,3NAC,3NGB,3NH7,3P0V,3P30,3PIQ,3Q6G,3QEG,3QEH,3QPX,3SDY,3SE8,3SKJ,3SQO,3T2N,3TCL,3THM,3TNM,3TNN,3TV3,3TWC,3U0T,3U1S,3U30,3U6R,3U7W,3UAJ,3UJI,3UJJ,3ULS,3W9D,3WD5,3WHE,3WLW,3WSQ,3X3F,4C2I,4CNI,4D9L,4D9Q,4DAG,4DGV,4DKF,4EDW,4EOW,4ERS,4F57,4FNL,4FQ1,4FQ2,4FQH,4FQJ,4FQL,4FQQ,4FZ8,4FZE,4G5Z,4G6F,4GLR,4GSD,4GW4,4GXU,4HBC,4HCR,4HF5,4HFU,4HFW,4HG4,4HH9,4HIE,4HIX,4HS6,4HS8,4HT1,4I77,4IDJ,4J4P,4J6R,4JAM,4JB9,4JFZ,4JHA,4JM2,4JO1,4JO3,4JPV,4JPW,4JY6,4K3J,4KRP,4KTD,4KTE,4KVN,4KY1,4LEO,4LMQ,4LRI,4LSP,4LSU,4M5Y,4M6O,4MA3,4MWF,4MXV,4N0Y,4N90,4N9G,4NHH,4NKI,4NM4,4NP4,4NPY,4NRX,4NUJ,4NZU,4O58,4O9H,4OAW,4OCS,4OCW,4OD1,4OD2,4OD3,4OLX,4OQT,4OSU,4PS4,4PTT,4PY7,4Q2Z,4QHL,4R26,4R4B,4R7D,4R8W,4R90,4RAV,4RFE,4RFO,4RX4,4RZC,4S1Q,4S1R,4S1S,4TSA,4U6V,4UAO,4UIF,4UOK,4UT6,4UT7,4UTA,4UU9,4UV4,4V1D,4WUU,4WV1,4XHJ,4XI5,4XMK,4XMP,4XNY,4XNZ,4XVJ,4XVS,4XXD,4Y5V,4Y5Y,4YAQ,4YDI,4YDJ,4YDK,4YDL,4YDV,4YE4,4YFL,4YHZ,4YPG,4YWG,4Z0X,4Z5R,4ZS6,4ZS7,4ZTO,4ZYK,5A3I,5ALB,5ANM,5AWN,5B71,5BK0,5BK3,5BMF,5BQ7,5BZD,5BZW,5C0N,5C6T,5C7X,5CCK,5CD3,5CEX,5CEZ,5CGY,5CHN,5CIL,5CJX,5CZV,5CZX,5D1Q,5D1X,5D1Z,5D6C,5DD6,5DMG,5DR5,5DRW,5DRX,5DRZ,5DSC,5DTF,5DUM,5DUR,5DWU,5E8E,5EA0,5ESV,5EWI,5F6H,5F6I,5F89,5F96,5F9O,5F9W,5FEH,5FGC,5FHA,5FHB,5FUO,5GGQ,5GGT,5GJS,5GMQ,5GS0,5I1E,5I5K,5I8C,5I8K,5I8O,5I9Q,5IBT,5IBU,5IIE,5IJK,5ITB,5JO4,5JO5,5JRP,5JZ7,5K9J,5K9Q,5KAN,5KAQ,5KEM,5KVL,5KW9,5L6Y,5LSP,5N4G,5N4J,5N7W,5NGV,5NYX,5O14,5O4G,5OB5,5OCK,5OTJ,5SX4,5T33,5T3X,5TE4,5TLK,5TPL,5TPN,5TPP,5TQA,5TRP,5TY6,5TZT,5U3J,5U3K,5U3M,5U3P,5U4R,5UBZ,5UD9,5UEK,5UEL,5UEM,5UG0,5UIX,5UKO,5UMI,5USL,5UXQ,5V2A,5V6L,5V6M,5V7R,5V7U,5VAG,5VIC,5VIG,5VK2,5VL7,5VOB,5VOD,5VQM,5W08,5W1G,5W1K,5W42,5W6G,5WB9,5WCA,5WCC,5WCD,5WDF,5WHJ,5WKO,5WNA,5WUV,5X8L,5X8M,5XAJ,5XHV,5XMH,5XWD,5Y2K,5Y9J,5YOY,5YY5,5ZIA,5ZV3,6A4K,6A67,6AL4,6APB,6AVN,6AXK,6AXL,6B08,6B0A,6B0E,6B0G,6B0H,6B3M,6B5L,6B9J,6BA5,6BCK,6BE2,6BF4,6BFQ,6BGT,6BKB,6BKC,6BKD,6BLA,6BLI,6BP2,6BQB,6BTJ,6C6X,6C9U,6CA6,6CA7,6CJK,6CMG,6CNR,6CT7,6CWT,6CYF,6D11,6D2P,6DB5,6DB6,6DB7,6DC3,6DCV,6DEZ,6DF0,6DFI,6DL8,6DLB,6DW2,6DWI,6E3H,6E4X,6E56,6E62,6E63,6EAY,6FG1,6FGB,6FOE,6FY0,6FY3,6GFE,6GG0,6GKU,6GLW,6GLX,6HIG,6HJP,6I04,6I9I,6IAP,6IEC,6IEK,6II4,6II9,6IUT,6IUV,6IVZ,6JEP,6K7O,6MED,6MEE,6MEG,6MHR,6MID,6MJZ,6MNQ,6MNR,6MQC,6MQE,6MQR,6MQS,6MTO,6MTP,6MTQ,6MTR,6MTS,6MTT,6N16,6N7J,6N81,6N8D,6NB3,6NB6,6NC2,6NN3,6NOV,6NZ7,6O39,6OBZ,6OE4,6OGX,6OL5,6OL7,6ORO,6ORP,6RCQ,6RCS,7FAB,8FAB

Examples:

#H3N2-1997

for file in split-*txt # *top.top.txt *bot.top.txt *bot.tail.txt *top.tail.txt
do
  bsub "../bin/epimorph-conservation.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target > ${target}-${file}-conservation"
 
  #bsub "../bin/epimorph.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target > ${target}-${file}"
done
target="H1N1" #-2007"

#H3N2-1997
for depth in 11 # 2 3 4 5 6 7
do
 
  for file in split-*txt # *top.top.txt *bot.top.txt *bot.tail.txt *top.tail.txt
  do
    bsub "../bin/epimorph.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target --centivax_min_conc=$depth > centivax-${target}-${file}-${depth}.txt"
  done
done
target="H1N1-2007"

#H3N2-1997

for file in split-*txt # *top.top.txt *bot.top.txt *bot.tail.txt *top.tail.txt
do
  bsub "../bin/epimorph.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target > ${target}-${file}"
done
target="H3N2" #-2007"

#H3N2-1997
for depth in 2 3 4 5 6 7 # 8 9 10 #11 # 2 3 4 5 6 7
do
  for file in split-*txt # *top.top.txt *bot.top.txt *bot.tail.txt *top.tail.txt
  do
    bsub "../bin/epimorph.pl --epitopedb=$file --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=$target --centivax_min_conc=$depth > centivax-${target}-${file}-${depth}.txt"
  done
done

#../bin/epimorph.pl --epitopedb=10k.sample.txt --alignment=alignment-publication-curated.msa.fa --numbering=numbering.txt --antigen=H1N1 --centivax_min_conc=2


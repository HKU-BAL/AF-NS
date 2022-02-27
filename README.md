# AF-NS
# Introduction
Workflow for identifying long-read novel sequences
![Workflow](http://www.bio8.cs.hku.hk/RNA/tmp/AF-NS.png)<br> 

# Prerequisites 
Make sure you have installed all of the following prerequisites on your machine：<br> 
• python2 <br> 
•	minimap2 <br> 
•	NUCmer <br> 
•	porechop<br> 
•	NanoFilt<br> 
•	bedtools<br> 
•	kraken2<br> 
•	RepeatMasker<br> 

# Run
python AF-NS.py -kraken_db db_folder -i input.fq -r ref.fa -o output_folder <br> 

## kraken2 DB for human
We build kraken2 DB including archaea, bacteria, fungi, plasmid, viral and UniVec datasets, the link is as follows <br> 
https://www.bio8.cs.hku.hk/RNA/kraken2_db.tar

## output
Novel sequences: output_folder/novel.fa





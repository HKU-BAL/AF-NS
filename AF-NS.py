import sys, argparse
import subprocess as sp

if __name__ == '__main__':
    parser = argparse.ArgumentParser('novel seq identification')
    required = parser.add_argument_group('required named arguments')
    parser.add_argument("-i",required=True, help='long-read input fasta')
    parser.add_argument('-r',required=True, help='reference fasta')
    parser.add_argument('-t',type=str, default='24', help='threads')
    parser.add_argument('-o',required=True, help='output folder')
    parser.add_argument('-l',type=str, default='300', help='minimum novel sequence length (default: 300)')
    parser.add_argument('-q',type=str, default='10', help= 'Qscore threshold (default: 10)')
    parser.add_argument('-repeat_pct',type=str, default='0.8' ,help= 'maximum percentage of simple_repeat/low_complexity in a read (default: 0.8); must be in [0,1]')
    parser.add_argument('-kraken_db', required=True, help='Kraken2 DB')
    parser.add_argument('-rp_species',type=str, default='human',help= 'Specify the species when using RepeatMasker (default: human)')
    parser.add_argument('-kraken_conf',type=str, default='0.05',help= 'kraken confidence score (default: 0.05)')
    parser.add_argument('-cov',type=str, default='0.8', help= 'coverage threshold when clustering (default: 0.8); must be in [0,1]')
    args = parser.parse_args()
    args, unknown = parser.parse_known_args()
   # procedure 1
    ref_index_cmd = ['minimap2','-d',args.o + '/ref.min', args.r]
    sp.call(ref_index_cmd)
    mkdir_cmd = ['mkdir' ,  args.o + '/tmp']
    sp.call(mkdir_cmd ) 
    alignment_cmd = 'minimap2 -t ' + args.t + ' -ax map-ont ' + args.o + '/ref.min ' + args.i + ' | samtools view -Sb -o ' + args.o + '/tmp/read2ref.bam - && htsbox samview -p ' + args.o + '/tmp/read2ref.bam > ' + args.o + '/tmp/read2ref.paf'    
    sp.call(alignment_cmd, shell=True)
    read_cmd = '''awk  'OFS="\t"{print $1, 0, $2}' ''' +  args.o + '/tmp/read2ref.paf' + ' | sort | uniq | bedtools sort -i - > ' + args.o + '/tmp/reads.bed'
    sp.call(read_cmd, shell=True)
    alignbed_cmd = " grep -v '*' " + args.o + "/tmp/read2ref.paf | " + '''awk 'OFS="\t"{print $1, $3, $4}' ''' + ' > ' + args.o + '/tmp/maps.bed' 
    sp.call(alignbed_cmd, shell=True)
    filter_cmd = 'bedtools sort -i ' + args.o+'/tmp/maps.bed | bedtools merge -i - | bedtools subtract -a ' + args.o + "/tmp/reads.bed -b - | awk '{if($3-$2>=" + args.l + ") print $1}' - | sort |uniq | seqtk subseq " + args.i + ' - | NanoFilt -q ' + args.q + " > " + args.o+'/tmp/highQ.fa' + ' && porechop -t ' + args.t + ' -i ' +  args.o+'/tmp/highQ.fa' + ' -o ' + args.o+"/tmp/trim.fq  && sed -n '1~4s/^@/>/p;2~4p' " + args.o+'/tmp/trim.fq  > ' + args.o + "/tmp/trim.fa"
    sp.call(filter_cmd, shell=True)

    # procedure 2
    round1_align = 'minimap2 -t ' + args.t + ' -x map-ont -o ' + args.o + '/tmp/1st.paf ' + args.o + '/ref.min ' + args.o + '/tmp/trim.fa'
    sp.call(round1_align, shell=True)
    sp.call([sys.executable, 'bin/reads.py',args.o + "/tmp/trim.fa", args.o + "/tmp/trim.bed"])
    bed_cmd = '''awk 'OFS="\t"{print $1, $3, $4}' ''' + args.o+'/tmp/1st.paf | bedtools sort -i - | bedtools subtract -a ' +  args.o + "/tmp/trim.bed -b - > " + args.o+'/tmp/1st_remain.bed'   
    sp.call(bed_cmd,shell=True)
    sp.call([sys.executable, 'bin/filter.py',args.o +'/tmp/', "trim.fa", "1st_remain.bed", 'novel_1st', args.l])       
    round2 = 'nucmer -t ' + args.t + ' -l 15 -c 31 -p ' + args.o + '/tmp/2nd ' + args.r + ' ' + args.o +'/tmp/novel_1st.fa && show-coords -q -o -T -H  -l ' + args.o + '/tmp/2nd.delta > ' +  args.o + '/tmp/2nd.tab && rm ' + args.o + '/tmp/2nd.delta  && ' + '''awk 'OFS="\t"{if($3<$4) print $11, $3-1, $4; else print $11, $4-1, $3}' ''' + args.o + '/tmp/2nd.tab |bedtools sort -i - |bedtools merge -i - | bedtools subtract -a ' + args.o + '/tmp/novel_1st.bed -b -  > ' + args.o + '/tmp/2nd_remain.bed'
    sp.call(round2, shell=True)
    sp.call([sys.executable, 'bin/filter.py',args.o +'/tmp/', "novel_1st.fa", "2nd_remain.bed", 'novel_2nd', args.l]) 
    round3_align = 'minimap2 -t ' + args.t + ' -x map-ont -o ' + args.o + '/tmp/3rd.paf ' + args.o + '/ref.min ' + args.o + '/tmp/novel_2nd.fa'
    sp.call(round3_align, shell=True)

    # procedure 3
    kraken = 'kraken2 --db ' + args.kraken_db + ' --threads ' + args.t + ' --confidence ' + args.kraken_conf + ' ' + args.o + '/tmp/novel_2nd.fa  --output ' + args.o + '/tmp/kraken --report ' + args.o + '/tmp/kraken_report --unclassified-out ' + args.o + '/tmp/clean.fa'   
    sp.call(kraken, shell=True)
    sp.call([sys.executable, 'bin/novel_cand.py', args.o +'/tmp/', args.l])
    cluster_cmd = 'minimap2 -t ' + args.t + ' -x ava-ont -o ' + args.o + '/tmp/overlap.paf ' + args.o + '/tmp/novel_tmp.fa ' + args.o + '/tmp/novel_tmp.fa'
    sp.call(cluster_cmd, shell=True)
    sp.call([sys.executable, 'bin/cluster.py', args.o +'/tmp/', args.cov])
    repeat_cmd = 'RepeatMasker -pa ' + args.t + ' -s  -species ' + args.rp_species ' + args.o + '/tmp/cluster.fa  -dir ' + args.o + '/tmp/repeatmasker'
    sp.call(repeat_cmd, shell=True)
    repeat_bed = ''' awk 'OFS="\t"{if($11=="Simple_repeat" || $11=="Low_complexity") print $5, $6-1, $7}' ''' + args.o + '/tmp/repeatmasker/cluster.fa.out |bedtools sort -i - |bedtools merge -i -  > ' + args.o + '/tmp/repeat.bed'
    sp.call(repeat_bed, shell=True)
    sp.call([sys.executable, 'bin/final.py', args.o,args.repeat_pct])
    rm_cmd = 'rm -r ' + args.o + '/tmp'
    sp.call(rm_cmd, shell=True)





#    bed_cmd = '''awk 'OFS="\t"{print $1, $3, $4}' ''' + args.o+'/tmp/3rd.paf | bedtools sort -i - | bedtools subtract -a ' +  args.o + "/tmp/novel_2nd.bed -b - > " + args.o+'/tmp/3rd_remain.bed'   
 #   sp.call(bed_cmd,shell=True)
  #  sp.call([sys.executable, 'bin/filter.py',args.o +'/tmp/', "novel_2nd.fa", "3rd_remain.bed", 'novel_3rd', args.l])   
  #  read_cmd = ["awk '{OFS='\t'}{print $1, 0, $2}' " +  args.o + '/tmp/read2ref.paf' + '| sort | uniq | bedtools sort -i - ']
 #   subprocess.call(read_cmd, stdout=open(args.tmp+'/reads.bed','w')
#    alignbed_cmd = [" grep -v '*' " + args.o + "/tmp/read2ref.paf |awk '{OFS='\t'}{print $1, $3, $4}' " ]
#    alignbed_cmd = [ "awk '{OFS='\t'}{print $1, $3, $4}' "  + args.o + "/tmp/read2ref.paf" ]

#    print alignbed_cmd
 #   sp.call(alignbed_cmd,stdout=open(args.o+'/tmp/maps.bed','w'), shell=True)
 #   select_cmd = ['bedtools sort -i' + args.tmp+'/tmp/maps.bed |bedtools merge -i - |bedtools subtract -a ' + args.tmp+'/reads.bed' + " -b - " ]
  #  subprocess.call(select_cmd, stdout=open(args.tmp+'/tmp




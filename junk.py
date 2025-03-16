

data='/exports/ana-scarlab/xqnguyen/course/sequencing/data/ABE_demux'
#### Find input fastq files ###

fq1s = sorted(glob.glob('/exports/ana-scarlab/xqnguyen/course/sequencing/ABE_demux_R1.fastq.gz')) 
fq1s = sorted(glob.glob(data + '*_R1*.fastq.gz')) 
fq2s = sorted(glob.glob(data + '*_R2*.fastq.gz')) 
print(fq1s, fq2s)
asdasd
if len(fq1s) != len(fq2s):
    sys.exit("Please, different number of input and output fastq files")

if len(fq1s) == len(fq2s) == 0:
    sys.exit('fastq files not found')

#### Read barcodes and expand set according to input hamming distance ####
if not os.path.isfile(args.cbcfile):
    sys.exit("Barcode file not found")


print('check barcode data frame')
bc_df = read_csv(args.cbcfile, sep = ',', names = ['bc','cellID'], index_col = 0)
print(bc_df.head())
print('check compatible_bcs')
bc_df['compatible_bcs'] = bc_df.apply(lambda x: find_compatible_barcodes(x.name, args.cbchd), axis = 1)
print('check cnt_allbcs')
cnt_allbcs = Counter([x for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs']])
print('check cellID')
allbc_df = pd.DataFrame({x: {'cellID': bc_df.loc[idx,'cellID'], 'original': idx} for idx in bc_df.index for x in bc_df.loc[idx, 'compatible_bcs'] if cnt_allbcs[x]==1}).T

### Create output directory if it does not exist ####
os.system('mkdir -p '+args.outdir)

#### Read fastq files and assign cell barcode and UMI ####
print('read fastqc')
if not args.demux:
    fout1 = open(args.outdir + '/' + args.fqo + '_R1.fastq', 'w')
    fout2 = open(args.outdir + '/' + args.fqo + '_R2.fastq', 'w')
else:
    fout1 = {idx:  open(args.outdir + '/' + args.fqo + '_' + str(bc_df.loc[idx, 'cellID']).zfill(3) + '_R1.fastq', 'w') for idx in bc_df.index}
    fout2 = {idx:  open(args.outdir + '/' + args.fqo + '_' + str(bc_df.loc[idx, 'cellID']).zfill(3) + '_R2.fastq', 'w') for idx in bc_df.index}

print('read fw_primers')
fw_primers = find_compatible_barcodes(args.fwdprimer, HDmax = 2)
print('read rv_primers')
rv_primers = find_compatible_barcodes(args.rvsprimer, HDmax = 2)

ns = 0; nt = 0
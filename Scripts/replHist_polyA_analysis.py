
## Examples of alignments
# A00514:356:HMHFVDRXX:1:1101:1136:1423   83      16      680508  60      98M3S   =       680403  -203    GCCGTGCCGCGCGGCGCCATGAAGGGCAAGGAGGAGAAGGAGGGCGGCGCACGTCTGGGCGCTGGCGGCGGAAGCCCCGAGAAGAGCCCGAGCGCGCAAGG    FFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:F,F,FFFFFFFFFF:FFFFFFF,:FFFFFFFFFF:FFFFFFFF:FFFFF,FFFFFFFFFFFFFFF:FFF,    AS:i:-6 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:53G44      YS:i:-7 YT:Z:CP  NH:i:1  MQ:i:60 MC:Z:2S99M      ms:i:3321
# A00514:356:HMHFVDRXX:1:1101:1136:1423   163     16      680403  60      2S99M   =       680508  203     GGCGGCGGCGGAGCTGGGCCGGGCCCGAGCGGATCGCGGGCTCGGGTTGCGGGGCTCCGGCTGCGGGCGCTGGGCCGCGAGGCGCGGAGCTTGGGAGCGGA    F::FFFFFFFFFFFFFFFF,FFFFF:F:FFFFFFFFFFF:FFFFFFF:FFFFFFFFF,F:FFFFF:FFFFF:F,,FFFFFFFFF,,FFFFFF,F:FFFF,F    AS:i:-7 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:44C54      YS:i:-6 YT:Z:CP  NH:i:1  MQ:i:60 MC:Z:98M3S      ms:i:3468
# A00514:356:HMHFVDRXX:1:1101:1136:25692  99      11      123059615       60      4S97M   =       123059760       243     TGTTTAGTTTGGCATCTCGAAGGGCTTTCTCTACTGGGTCCAGGGTGCCACGGAACAGGTCAGCATTCAGTTCTTCAAATCGGGCACGGGTAATGGAGGTA    ,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:::FF:FFFFFFFF:FFFFFFFF:FFFF:FFFFFFFFFFFFFF,FFFFF,FFFF,F:FF    AS:i:-4 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:97 YS:i:-3  YT:Z:CP NH:i:1  MQ:i:60 MC:Z:98M3S      ms:i:3554
# A00514:356:HMHFVDRXX:1:1101:1136:25692  147     11      123059760       60      98M3S   =       123059615       -243    GTGCTGGAAGAGAGGGTACGCTTAGCACGTTCACAAGCAGTACGGAGGCGTCTTACAGCTCTCTTGTTCTCACTGATGTCCTTCTTATGCTTGCGCTTCGC    FFFFFFFFFFFFFFFFFFFF,FFFFFFFF,FFFFFFFFFFFF:FFFFFFF:FFFFFFFF:FFFF:FFFF:FFFFFFFFFF:FFFFFF,FFFFFFFFFFFFF    AS:i:-3 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:98 YS:i:-4  YT:Z:CP NH:i:1  MQ:i:60 MC:Z:4S97M

## Example of use for this script
# samtools view /Volumes/Seb_GA/replHist_polyA_temp/D356T5.nsorted.bam | head -n 1000 | python3 ./Scripts/ASF1_chaperone/replHist_polyA_analysis_230622.py ./Data/ASF1_chaperone/annotations/replicative_histone_genes.gff | less -S


import sys
import re
import pandas as pd

def dec_to_bin(val, size=12):
    out = [0] * size
    #
    iii = 0
    while val:
        # print(val)
        out[iii] = val % 2
        val = val // 2
        iii += 1
    #
    return(out)


def gff_attr_proc(attr_str):
    # attr_str='''gene_id "ENSG00000273213"; gene_version "4"; gene_name "H3-2"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; replication "dependent";'''
    #
    attrs = attr_str.split(';')[0:-1]
    #
    ## get all values
    outd = {}
    if '=' in attrs[0]:
        for iii, attr in enumerate(attrs):
            ## the value
            val_srch = re.search('=[^=]*$', attr)
            startp = val_srch.start() + 1
            endp = val_srch.end()
            value = attr[startp:endp]
            #
            ## the attribute id
            attr_srch = re.search('^[^=]*=', attr)
            attr_id = attr[0:(attr_srch.end() - 1)].strip()
            #
            ## record the attribute
            outd[attr_id] = value
    else:
        for iii, attr in enumerate(attrs):
            ## the value
            val_srch = re.search('\"[^\"][^\"]*\"', attr)
            startp = val_srch.start() + 1
            endp = val_srch.end() - 1
            value = attr[startp:endp]
            #
            ## the attribute id
            attr_srch = re.search('^[^\"]*\"', attr)
            attr_id = attr[0:(attr_srch.end() - 1)].strip()
            #
            ## record the attribute
            outd[attr_id] = value
    #
    return(outd)


class Ann:
    def __init__(self, line):
        parts=line.strip().split('\t')
        # print(parts)
        self.ref=parts[0]
        self.start=int(parts[3]) - 1
        self.end=int(parts[4]) - 1
        self.strand=parts[6]
        self.len=self.end - self.start + 1
        self.attrs=parts[-1]
        #
        ## parse the attributes
        self.attrs = gff_attr_proc(self.attrs)
    #
    def print(self):
        print(self.ref, self.start, self.end, self.strand)
    #
    def print_bed(self):
        return('\t'.join([self.ref, str(self.start + 1), str(self.end + 1), '.', '0', self.strand]))
    #
    def idx(self):
        return('_'.join([self.ref, str(self.start), str(self.end), self.strand]))
    #
    def extend_threep(elen):
        if self.strand == '-':
            self.start -= elen
        else:
            self.end += elen
        #
        self.len=self.end - self.start + 1
        #
        return(None)


def load_gff(gffp):
    ## retrieve the annotations
    anns = []
    if re.search('\\.gz', gffp):
        with gzip.open(gffp, 'rb') as fich:
            for line in fich:
                if line[0] == '#':
                    continue
                #
                anns.append(Ann(line))
    else:
        with open(gffp, 'r') as fich:
            for line in fich:
                if line[0] == '#':
                    continue
                #
                anns.append(Ann(line))
    #
    annsd = {}
    for ann in anns:
        try:
            annsd[ann.ref].append(ann)
        except:
            annsd[ann.ref] = [ann]
    #
    ## order by starting position, for each ref
    for ref in annsd:
        start_vec = [ann.start for ann in annsd[ref]]
        ord_idx = [i for _,i in sorted(zip(start_vec, list(range(len(start_vec)))))]
        annsd[ref] = [annsd[ref][iii] for iii in ord_idx]
    #
    return(annsd)


def print_all_anns(anns):
    '''Print in table shape all the annotations from a list of 'Ann' objects.'''
    #
    for ref in anns:
        annu = anns[ref]
        #
        for ann in annu:
            print('\t'.join([ann.ref, str(ann.start), str(ann.end), ann.strand] + [attr + ':' + ann.attrs[attr] for attr in ann.attrs]))


def loci_in_any_ann(ref, pos, pair_strand, anns, on_rev_strand={'+': '-', '-': '+'}):
    if ref in anns:
        annu = anns[ref]
        #
        for ann in annu:
            # if pos >= ann.start and pos <= ann.end and (pair_strand in ['*', ann.strand]):
            if pos >= ann.start and pos <= ann.end and (pair_strand in ['*', on_rev_strand[ann.strand]]):
                # print(True)
                return(ann)
    else:
        return(False)
    #
    return(False)


def read_line_proc(readline):
    fields = readline.strip().split('\t')
    mq = int(fields[4]) ## mapping quality
    sf = dec_to_bin(int(fields[1])) ## binary sam flag
    ref = fields[2] ## mapping reference (chrmosome, ...)
    pos = int(fields[3]) ## mapping position
    strand = ['+', '-'][sf[4]] ## strand of the read itself
    pair_strand = pair_strandness(sf) ## strand of the pair to which the read belongs to
    #
    return((fields, mq, sf, ref, pos, strand, pair_strand))


def pair_strandness(sam_flag_bin):
    if sam_flag_bin[2]: ## read is not mapped, thus no strandness available
        return(None)
    #
    pair_strand = ['+', '-'][(sam_flag_bin[4] + sam_flag_bin[7]) % 2] ## strand of the pair to which the read belongs to
    #
    return(pair_strand)



## parameter
# mq_thres = 2

## read from stdin (~in pipeline)
fich = sys.stdin
ann_path = sys.argv[1]

## sam flags
sam_flags = [
    'read paired',
    'read mapped in proper pair',
    'read unmapped',
    'mate unmapped',
    'read reverse strand',
    'mate reverse strand',
    'first in pair',
    'second in pair',
    'not primary alignment',
    'read fails platform/vendor quality checks',
    'read is PCR or optical duplicate',
    'supplementary alignment'
]

## load the gene annotations of human genome (GRCh38)
gff = load_gff(ann_path)

for reada in fich:
    # print(reada, end='') # temp
    in_annot = False
    #
    ## parse BAM read 1 ('a')
    fieldsa, mqa, sfa, refa, posa, stranda, pair_stranda = read_line_proc(reada)
    #
    ## technical check (is read 1 primary alignment)
    while sfa[8]:
        reada = fich.readline()
        fieldsa, mqa, sfa, refa, posa, stranda, pair_stranda = read_line_proc(reada)
    #
    ## technical check (in read pair)
    # if not (sfa[0] and sfa[6]):
    if not sfa[0]:
        continue
    #
    ## same parsing for read 2 ('b')
    readb = fich.readline()
    fieldsb, mqb, sfb, refb, posb, strandb, pair_strandb = read_line_proc(readb)
    #
    ## technical check (is read 2 primary alignment)
    while sfb[8]:
        readb = fich.readline()
        fieldsb, mqb, sfb, refb, posb, strandb, pair_strandb = read_line_proc(readb)
    #
    ## technical check (real read1 and real read2)
    if (not sfa[6] and not sfb[7]): ## reads are inverted
        fieldsa, mqa, sfa, refa, posa, stranda, pair_stranda = read_line_proc(readb)
        fieldsb, mqb, sfb, refb, posb, strandb, pair_strandb = read_line_proc(reada)
    elif (sfa[6] and sfb[7]): ## the two reads are reported in right order
        pass
    else: ## something weird, do not process and go to next read pair
        continue
    #
    ## case where both reads are mapped (in proper pair) or case where only read 2 is mapped (not proper pair, read 1 unmapped, mate of read 2 unmapped)
    if (sfa[1] and (not sfa[2])) or ((not sfa[1]) and sfa[2] and sfb[3]): # (1100xx10 and 1100xx01) or (1010xx10 and 1001xx01)
        in_annot = loci_in_any_ann(ref=refb, pos=posb, pair_strand=pair_strandb, anns=gff)
    #
    if in_annot:
        ## filter for pairs with read 1 containing stretch of Ts on 5' end
        if bool(re.search('(^T{22})', fieldsa[9])) and bool(re.search('(^[0-9]*S|\*)', fieldsa[5])):
            ann = in_annot
            gene_id = ann.attrs['gene_id']
            reada_mod = '%s GI:%s\n' % (reada.strip(), gene_id)
            #
            # print(reada, end='')
            # print(readb, end='')
            print(reada_mod, end='')
            # print('\t'.join([ann.ref, str(ann.start), str(ann.end), ann.strand] + [attr + ':' + ann.attrs[attr] for attr in ann.attrs]))

####

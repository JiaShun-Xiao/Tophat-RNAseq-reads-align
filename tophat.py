# -*- coding: utf-8 -*-
__author__ = 'jiashun'

sam = open('/2_disk/xiaojs/RNAseq/cut_align/SRR1793917_1_cut.sam', 'r')
sam_split = open('/2_disk/xiaojs/RNAseq/mapping/SRR1793917_1_split.sam', 'w')
change = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'M', 'X', 'Y']
chan = {'1':1, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '2':2, '20':20, '21':21, '22':22, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, 'M':23, 'X':24, 'Y':25}

chroms = []
for c in change:
    chrom = open('/2_disk/xiaojs/python/blast/chrom_cal/chromosome'+str(c)+".txt", 'r')
    chroms.append(chrom.read().strip())


def alignment(qu, re):
    if len(qu) > len(re):
        return 0
    for i in xrange(len(qu)):
        if qu[i] == 'N' or re[i] == 'N':
            continue
        else:
            if qu[i] != re[i]:
                return i
    if i == len(qu)-1:
        return i+1


def complement(query):
    que_ = list(query)
    que_.reverse()
    com = []
    for q in que_:
        if q == 'A':
            com.append('T')
        elif q == 'C':
            com.append('G')
        elif q == 'G':
            com.append('C')
        elif q == 'T':
            com.append('A')
        else:
            com.append(q)
    return ''.join(com)


def revers_alignment(qu, re):
    if len(qu) > len(re):
        return 0
    qu = list(qu)
    qu.reverse()
    qu = ''.join(qu)
    re = list(re)
    re.reverse()
    re = ''.join(re)
    for i in xrange(len(qu)):
        if qu[i] == 'N' or re[i] == 'N':
            continue
        else:
            if qu[i] != re[i]:
                return i
    if i == len(qu)-1:
        return i+1


def find_all(arr, item):
    return [i for i,a in enumerate(arr) if item == a]


for i in xrange(27):
    line = sam.readline()
    sam_split.write(line)

nomis = 0
onemis = 0
twomis = 0
threemis = 0

while line:
    line = sam.readline().strip()
    if line == '':
        break
    line_s = line.split()
    reads = line_s[0][0:-2]
    line_nu = line_s[1]
    line_ch = line_s[2]
    line_lo = line_s[3]
    line_mq = line_s[4]
    line_seq = line_s[9]
    line_qual = line_s[10]
    l2 = sam.readline().strip()
    l2_s = l2.split()
    l2_nu = l2_s[1]
    l2_ch = l2_s[2]
    l2_lo = l2_s[3]
    l2_mq = l2_s[4]
    l2_seq = l2_s[9]
    l2_qual = l2_s[10]
    l3 = sam.readline().strip()
    l3_s = l3.split()
    l3_nu = l3_s[1]
    l3_ch = l3_s[2]
    l3_lo = l3_s[3]
    l3_mq = l3_s[4]
    l3_seq = l3_s[9]
    l3_qual = l3_s[10]
    l4 = sam.readline().strip()
    l4_s = l4.split()
    l4_nu = l4_s[1]
    l4_ch = l4_s[2]
    l4_lo = l4_s[3]
    l4_mq = l4_s[4]
    l4_seq = l4_s[9]
    l4_qual = l4_s[10]
    all_nu = [line_nu, l2_nu, l3_nu, l4_nu]
    if '4' in all_nu:
        all_nu.remove('4')
    all_nu = list(set(all_nu))
    if len(all_nu) != 1:
        continue
    all_ch = [line_ch, l2_ch, l3_ch, l4_ch]
    if '*' in all_ch:
        all_ch.remove('*')
    all_ch = list(set(all_ch))
    if len(all_ch) != 1:
        continue
    lo1 = int(line_lo)
    lo2 = int(l2_lo)
    lo3 = int(l3_lo)
    lo4 = int(l4_lo)
    test = [line_nu, l2_nu, l3_nu, l4_nu]
    index = find_all(test, '4')
    if len(index) == 0:
        nomis += 1
        ch = line_ch
        if line_nu == '0':
            loo1 = lo1
            if lo1 != lo2 - 25:
                if line_mq == '42':
                    seq1 = line_seq
                    qual1 = line_qual
                else:
                    seq1 = line_seq[0:len(line_seq)-1]
                    qual1 = line_qual[0:len(line_qual)-1]
                if l2_mq == '42':
                    loo2 = lo2
                    seq2 = l2_seq + l3_seq + l4_seq
                    qual2 = l2_qual +l3_qual + l4_qual
                else:
                    loo2 = lo2+1
                    seq2 = l2_seq[1:len(l2_seq)] + l3_seq + l4_seq
                    qual2 = l2_qual[1:len(l2_qual)] + l3_qual + l4_qual
            elif lo1 == lo2 - 25 != lo3 - 50:
                if l2_mq == '42':
                    seq1 = line_seq + l2_seq
                    qual1 = line_qual + l2_qual
                else:
                    seq1 = line_seq + l2_seq[0:len(l2_seq)-1]
                    qual1 = line_qual + l2_qual[0:len(l2_qual)-1]
                if l3_mq == '42':
                    loo2 = lo3
                    seq2 = l3_seq + l4_seq
                    qual2 = l3_qual + l4_qual
                else:
                    loo2 = lo3+1
                    seq2 = l3_seq[1:len(l3_seq)] + l4_seq
                    qual2 = l3_qual[1:len(l3_qual)] + l4_qual
            else:
                if l3_mq == '42':
                    seq1 = line_seq + l2_seq + l3_seq
                    qual1 = line_qual + l2_qual + l3_qual
                else:
                    seq1 = line_seq + l2_seq + l3_seq[0:len(l3_seq)-1]
                    qual1 = line_qual + l2_qual + l3_qual[0:len(l3_qual)-1]
                if l4_mq == '42':
                    loo2 = lo4
                    seq2 = l4_seq
                    qual2 = l4_qual
                else:
                    loo2 = lo4+1
                    seq2 = l4_seq[1:len(l4_seq)]
                    qual2 = l4_qual[1:len(l4_qual)]
            sam_split.write(reads+'_1'+'\t'+str(0)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
            sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
        elif line_nu == '16':
            loo2 = lo4
            if lo1 != lo2 + 25:
                if line_mq == '42':
                    loo1 = lo1
                    seq1 = line_seq
                    qual1 = line_qual
                else:
                    loo1 = lo1 + 1
                    seq1 = line_seq[1:len(line_seq)]
                    qual1 = line_qual[1:len(line_qual)]
                if l2_mq == '42':
                    seq2 = l4_seq + l3_seq + l2_seq
                    qual2 = l4_qual + l3_qual + l2_qual
                else:
                    seq2 = l4_seq + l3_seq + l2_seq[0:len(l2_seq)-1]
                    qual2 = l4_qual + l3_qual + l2_qual[0:len(l2_qual)-1]
            elif lo1 == lo2 + 25 != lo3 + 50:
                if l2_mq == '42':
                    loo1 = lo2
                    seq1 = l2_seq + line_seq
                    qual1 = l2_qual + line_qual
                else:
                    loo1 = lo2+1
                    seq1 = l2_seq[1:len(l2_seq)] + line_seq
                    qual1 = l2_qual[1:len(l2_qual)] + line_qual
                if l3_mq == '42':
                    seq2 = l4_seq + l3_seq
                    qual2 = l4_qual + l3_qual
                else:
                    seq2 = l4_seq + l3_seq[0:len(l3_seq)-1]
                    qual2 = l4_qual + l3_qual[0:len(l3_qual)-1]
            else:
                if l3_mq == '42':
                    loo1 = lo3
                    seq1 = l3_seq + l2_seq + line_seq
                    qual1 = l3_qual + l2_qual + line_qual
                else:
                    loo1 = lo3+1
                    seq1 = l3_seq[1:len(l3_seq)] + l2_seq + line_seq
                    qual1 = l3_qual[1:len(l3_qual)] + l2_qual + line_qual
                if l4_mq == '42':
                    seq2 = l4_seq
                    qual2 = l4_qual
                else:
                    seq2 = l4_seq[0:len(l4_seq)-1]
                    qual2 = l4_qual[0:len(l4_qual)-1]
            sam_split.write(reads+'_1'+'\t'+str(16)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
            sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
    elif len(index) == 1:
        onemis += 1
        if len(find_all(test, '0')) > 0:
            if index[0] == 0:
                ch = l2_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                ref = chromo[lo2-50000:lo2-1]
                gap = revers_alignment(line_seq, ref)
                loo2 = lo2 - gap
                seq1 = line_seq[0:len(line_seq)-gap]
                qual1 = line_qual[0:len(line_seq)-gap]
                #print seq1
                loo1 = ref.find(seq1)
                loo3 = -1
                if lo2+50 == lo3 + 25 == lo4:
                    seq2 = line_seq[len(line_seq)-gap:len(line_seq)] + l2_seq + l3_seq + l4_seq
                    qual2 = line_qual[len(line_seq)-gap:len(line_seq)] + l2_qual + l3_qual + l4_qual
                elif lo2+50 == lo3 + 25 != lo4:
                    seq2 = line_seq[len(line_seq)-gap:len(line_seq)] + l2_seq + l3_seq
                    qual2 = line_qual[len(line_seq)-gap:len(line_seq)] + l2_qual + l3_qual
                    loo3 = lo4
                    seq3 = l4_seq
                    qual3 = l4_qual
                elif lo2+50 != lo3 + 25 == lo4:
                    seq2 = line_seq[len(line_seq)-gap:len(line_seq)] + l2_seq
                    qual2 = line_qual[len(line_seq)-gap:len(line_seq)] + l2_qual
                    loo3 = lo3
                    seq3 = l3_seq + l4_seq
                    qual3 = l3_qual + l4_qual
                else:
                    continue
                if loo1 > 0:
                    loo1 += lo2 - 50000 + 1
                    sam_split.write(reads+'_1'+'\t'+str(0)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                    sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                    if loo3 != -1:
                        sam_split.write(reads+'_3'+'\t'+str(0)+'\t'+ch+'\t'+str(loo3)+'\t'+str(42)+'\t'+str(len(seq3))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq3+'\t'+qual3+'\n')
                else:
                    sam_split.write(reads+'_1'+'\t'+str(0)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                    if loo3 != -1:
                        sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo3)+'\t'+str(42)+'\t'+str(len(seq3))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq3+'\t'+qual3+'\n')
            if index[0] == 1:
                ch = line_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                loo1 = lo1
                ref = chromo[lo1+25-1:lo1+50]
                gap = alignment(l2_seq, ref)
                seq1 = line_seq+l2_seq[0:gap]
                qual1 = line_qual+l2_qual[0:gap]
                ref = chromo[lo3-26:lo3-1]
                '''
                print '-----------'
                print l2_seq
                print ref
                '''
                gap = revers_alignment(l2_seq, ref)
                loo2 = lo3 - gap
                seq2 = l2_seq[25-gap:25] + l3_seq + l4_seq
                qual2 = l2_qual[25-gap:25] + l3_qual + l4_qual
                #print seq2
                #print chromo[loo2-1:loo2+75]
                #print '--------------'
                sam_split.write(reads+'_1'+'\t'+str(0)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                '''
                print lo1, lo2, lo3, lo4
                print line_seq+' '+l2_seq+' '+l3_seq+' '+ l4_seq
                print seq1, len(seq1)
                print chromo[loo1-1:loo1+100]
                print '--------------'
                print seq2, len(seq2)
                print chromo[loo2-1:loo2+100]
                print '--------------'
                '''
            if index[0] == 2:
                ch = line_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                loo1 = lo1
                ref = chromo[lo2+25-1:lo2+50]
                gap = alignment(l3_seq, ref)
                seq1 = line_seq+l2_seq+l3_seq[0:gap]
                qual1 = line_qual+l3_qual+l3_qual[0:gap]
                ref = chromo[lo4-26:lo4-1]
                gap = revers_alignment(l3_seq, ref)
                loo2 = lo4 - gap
                seq2 = l3_seq[25-gap:25] + l4_seq
                qual2 = l3_qual[25-gap:25] + l4_qual
                sam_split.write(reads+'_1'+'\t'+str(0)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
            if index[0] == 3:
                ch = line_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                loo1 = lo1
                ref = chromo[lo3+25-1:lo3+50000]
                '''
                print '-----------'
                print l4_seq
                print ref[0:25]
                '''
                gap = alignment(l4_seq, ref)
                seq_end = l4_seq[gap:25]
                qual_end = l4_seq[gap:25]
                loo_end = ref.find(seq_end)
                loo2 = -1
                if lo1+50 == lo2 + 25 == lo3:
                    seq1 = line_seq+l2_seq+l3_seq+l4_seq[0:gap]
                    qual1 = line_qual+l3_qual+l3_qual+l4_qual[0:gap]
                elif lo1+50 == lo2 + 25 != lo3:
                    seq1 = line_seq+l2_seq
                    qual1 = line_qual+l3_qual
                    loo2 = lo3
                    seq2 = l3_seq+l4_seq[0:gap]
                    qual2 = l3_qual+l4_qual[0:gap]
                elif lo1+50 != lo2 + 25 == lo3:
                    seq1 = line_seq
                    qual1 = line_qual
                    seq2 = l2_seq+l3_seq+l4_seq[0:gap]
                    qual2 = l2_qual+l3_qual+l4_qual[0:gap]
                    loo2 = lo2
                sam_split.write(reads+'_1'+'\t'+str(0)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                if loo2 != -1:
                    sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                    if loo_end > 0:
                        loo_end += lo3 + 25
                        sam_split.write(reads+'_3'+'\t'+str(0)+'\t'+ch+'\t'+str(loo_end)+'\t'+str(42)+'\t'+str(len(seq_end))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq_end+'\t'+qual_end+'\n')
                else:
                    if loo_end > 0:
                        loo_end += lo3 + 25
                        sam_split.write(reads+'_2'+'\t'+str(0)+'\t'+ch+'\t'+str(loo_end)+'\t'+str(42)+'\t'+str(len(seq_end))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq_end+'\t'+qual_end+'\n')
                '''
                print lo1, lo2, lo3, lo4
                print seq1, len(seq1)
                print chromo[loo1-1:loo1+100]
                print '--------------'
                if loo2 != -1:
                    print seq2, len(seq2)
                    print chromo[loo2-1:loo2+100]
                    print '--------------'
                if loo_end != -1:
                    print seq_end, len(seq_end)
                    print chromo[loo_end-1:loo_end+100]
                    print '--------------'
                '''
        elif len(find_all(test, '16')) > 0:
            if index[0] == 0:
                ch = l2_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                ref = chromo[lo2+25-1:lo2+50000]
                gap = alignment(complement(line_seq), ref)
                seq1 = complement(line_seq[0:len(line_seq)-gap])
                qual1 = line_qual[0:len(line_seq)-gap]
                #print seq1
                loo1 = ref.find(seq1)
                loo3 = -1
                if lo2-50 == lo3 - 25 == lo4:
                    loo2 = lo4
                    seq2 = l4_seq + l3_seq + l2_seq + complement(line_seq[len(line_seq)-gap:len(line_seq)])
                    qual2 = l4_qual + l3_qual + l2_qual + line_qual[len(line_seq)-gap:len(line_seq)]
                elif lo2 - 50 == lo3 - 25 != lo4:
                    loo2 = lo3
                    seq2 = l3_seq + l2_seq + complement(line_seq[len(line_seq)-gap:len(line_seq)])
                    qual2 = l3_qual + l2_qual + line_qual[len(line_seq)-gap:len(line_seq)]
                    loo3 = lo4
                    seq3 = l4_seq
                    qual3 = l4_qual
                elif lo2 - 50 != lo3 - 25 == lo4:
                    loo2 = lo2
                    seq2 = l2_seq + complement(line_seq[len(line_seq)-gap:len(line_seq)])
                    qual2 = l2_qual + line_qual[len(line_seq)-gap:len(line_seq)]
                    loo3 = lo4
                    seq3 = l4_seq + l3_seq
                    qual3 = l4_qual + l3_qual
                else:
                    continue
                if loo1 > 0:
                    loo1 += lo2+25
                    sam_split.write(reads+'_1'+'\t'+str(16)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                    sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                    if loo3 != -1:
                        sam_split.write(reads+'_3'+'\t'+str(16)+'\t'+ch+'\t'+str(loo3)+'\t'+str(42)+'\t'+str(len(seq3))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq3+'\t'+qual3+'\n')
                else:
                    sam_split.write(reads+'_1'+'\t'+str(16)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                    if loo3 != -1:
                        sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo3)+'\t'+str(42)+'\t'+str(len(seq3))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq3+'\t'+qual3+'\n')
                '''
                print lo1, lo2, lo3, lo4
                if loo1 != -1:
                    print seq1
                    print chromo[loo1-1:loo1+100]
                    print '--------------'
                if loo2 != -1:
                    print seq2
                    print chromo[loo2-1:loo2+100]
                    print '--------------'
                if loo3 != -1:
                    print seq3
                    print chromo[loo3-1:loo3+100]
                    print '--------------'
                '''
            if index[0] == 1:
                ch = line_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                loo2 = lo4
                ref = chromo[lo3+25-1:lo3+50]
                gap = alignment(complement(l2_seq), ref)
                seq2 = l4_seq+l3_seq+complement(l2_seq)[0:gap]
                qual2 = l4_qual+l3_qual+l2_qual[0:gap]
                ref = chromo[lo1-26:lo1-1]
                gap = revers_alignment(complement(l2_seq), ref)
                loo1 = lo1 - gap
                seq1 = complement(l2_seq)[25-gap:25] + line_seq
                qual1 = l2_qual[25-gap:25] + line_qual
                sam_split.write(reads+'_1'+'\t'+str(16)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
            if index[0] == 2:
                ch = line_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                loo2 = lo4
                ref = chromo[lo4+25-1:lo4+50]
                gap = alignment(complement(l3_seq), ref)
                seq2 = l4_seq+complement(l3_seq)[0:gap]
                qual2 = l4_qual+l3_qual[0:gap]
                ref = chromo[lo2-26:lo2-1]
                gap = revers_alignment(complement(l3_seq), ref)
                loo1 = lo2 - gap
                seq1 = complement(l3_seq)[25-gap:25] + l2_seq + line_seq
                qual1 = l3_qual[25-gap:25] + l2_qual + line_qual
                sam_split.write(reads+'_1'+'\t'+str(16)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                '''
                print lo1, lo2, lo3, lo4
                print line_seq+' '+l2_seq+' '+l3_seq+' '+ l4_seq
                print seq1, len(seq1)
                print chromo[loo1-1:loo1+100]
                print '--------------'
                print seq2, len(seq2)
                print chromo[loo2-1:loo2+100]
                print '--------------'
                '''
            if index[0] == 3:
                ch = l2_ch
                chr = ch[3:len(ch)]
                chromo = chroms[chan[chr]-1]
                ref = chromo[lo3-50000:lo3-1]
                gap = revers_alignment(complement(l4_seq), ref)
                seq_end = complement(l4_seq[0:len(l4_seq)-gap])
                qual_end = line_qual[0:len(l4_seq)-gap]
                #print seq1
                loo_end = ref.find(seq_end)
                loo2 = -1
                if lo1-50 == lo2 - 25 == lo3:
                    loo1 = lo3
                    seq1 = complement(l4_seq)[len(l4_seq)-gap:len(l4_seq)] + l3_seq + l2_seq + line_seq
                    qual1 = l4_qual[len(l4_seq)-gap:len(l4_seq)] + l3_qual + l2_qual + line_qual
                elif lo1-50 == lo2 - 25 != lo3:
                    loo1 = lo2
                    seq1 = l2_seq + line_seq
                    qual1 = l3_qual + line_qual
                    loo2 = lo3 - gap
                    seq2 = complement(l4_seq)[len(l4_seq)-gap:len(l4_seq)] + l3_seq
                    qual2 = l4_qual[len(l4_seq)-gap:len(l4_seq)] + l3_qual
                elif lo1-50 != lo2 - 25 == lo3:
                    loo1 = lo1
                    seq1 = line_seq
                    qual1 = line_qual
                    loo2 = lo3 - gap
                    seq2 = complement(l4_seq)[len(l4_seq)-gap:len(l4_seq)] + l3_seq +l2_seq
                    qual2 = l4_qual[len(l4_seq)-gap:len(l4_seq)] + l3_qual + l2_qual
                else:
                    continue
                sam_split.write(reads+'_1'+'\t'+str(16)+'\t'+ch+'\t'+str(loo1)+'\t'+str(42)+'\t'+str(len(seq1))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq1+'\t'+qual1+'\n')
                if loo2 != -1:
                    sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo2)+'\t'+str(42)+'\t'+str(len(seq2))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq2+'\t'+qual2+'\n')
                    if loo_end > 0:
                        loo_end += lo3-50000+1
                        sam_split.write(reads+'_3'+'\t'+str(16)+'\t'+ch+'\t'+str(loo_end)+'\t'+str(42)+'\t'+str(len(seq_end))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq_end+'\t'+qual_end+'\n')
                else:
                    if loo_end > 0:
                        loo_end += lo3-50000+1
                        sam_split.write(reads+'_2'+'\t'+str(16)+'\t'+ch+'\t'+str(loo_end)+'\t'+str(42)+'\t'+str(len(seq_end))+'M'+'\t'+'*'+'\t'+'0'+'\t'+'0'+'\t'+seq_end+'\t'+qual_end+'\n')
    elif len(index) == 2:
        twomis += 1
    elif len(index) == 3:
        threemis += 1

print nomis, onemis, twomis, threemis
sam.close()
sam_split.close()

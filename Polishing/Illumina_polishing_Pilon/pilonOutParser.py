#!/usr/bin/env python3
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Parsing Pilon polishing standart output.')
#parser.add_argument('-i', '--input', type=argparse.FileType('r'))
#parser.add_argument('-o', '--output', type=argparse.FileType('w'))

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-i', '--input', type=argparse.FileType('r'), help='Input file name', required=True)
requiredNamed.add_argument('-o', '--output', type=argparse.FileType('w'), help='Output file name', required=True)

args = parser.parse_args()
infile = args.input
outfile = args.output

scaffolds = []
s_len = []
total_reads = []
coverage = []
confirmed = []
total_bases = []
snps = []
insertions = []
ins_bases = []
deletions = []
del_bases = []
gaps = []
breaks = []

for line in infile:
    if line[-5:] == "log:\n":
        line = line[:-6]
        scaffolds.append(line.split(':')[0])
        s_len.append(int(line.split(':')[1].split('-')[1]) -
                     int(line.split(':')[1].split('-')[0]) + 1)
    elif len(line.split())>2 and line.split()[2] == 'coverage':
        coverage.append(int(line.split()[-1]))
    elif line[:3] == 'Tot':
        total_reads.append(int(line.split(',')[0].split()[-1]))
    elif line[:4] == 'Conf':
        confirmed.append(int(line.split()[1]))
        total_bases.append(int(line.split()[3]))
    elif line[:4] == 'Corr':
        snps.append(int(line.split(';')[0].split()[1]))
        insertions.append(int(line.split(';')[1].split()[1]))
        ins_bases.append(int(line.split(';')[1].split()[5]))
        deletions.append(int(line.split(';')[1].split()[7]))
        del_bases.append(int(line.split(';')[1].split()[11]))
    elif line[:8] == 'fix gap:':
        splitted = line.split()
        gaps.append((splitted[2].split(';')[0], 
                     int(splitted[4][1:]),
                     int(splitted[5][1:]),
                     splitted[6]))
    elif line[:10] == 'fix break:':
        splitted = line.split()
        breaks.append((splitted[2].split(';')[0], 
                     int(splitted[4][1:]),
                     int(splitted[5][1:]),
                     splitted[6]))

infile.close()

df = pd.DataFrame({'Scaffold': scaffolds, 'Length': s_len, 'total_reads': total_reads,
                   'Coverage': coverage, 'Confirmed': confirmed, 'Total_bases': total_bases,
                   'snps': snps, 'ins': insertions, 'ins_bases': ins_bases,
                   'dels': deletions, 'del_bases': del_bases}).groupby('Scaffold').agg('sum')

print(f'total bases\t{df.Length.sum()}', file=outfile)
print(f'confirmed\t{df.Confirmed.sum()} ({df.Confirmed.sum()/df.Length.sum()*100:.2f}%)', file=outfile)
print(f'snps\t{df.snps.sum()} ({df.snps.sum()/df.Length.sum()*100:.2}%)', file=outfile)
print(f'insertions\t{df.ins.sum()}\ninsert_bases\t{df.ins_bases.sum()} ({df.ins_bases.sum()/df.Length.sum()*100:.2}%)', file=outfile)
print(f'deletions\t{df.dels.sum()}\ndels_bases\t{df.del_bases.sum()} ({df.del_bases.sum()/df.Length.sum()*100:.2}%)', file=outfile)
print(file=outfile)

gaps_df = pd.DataFrame(gaps, columns=['seq_id', 'deleted', 'added', 'type'])

print(f'Gaps filled\t{len(gaps_df)}', file=outfile)
print(f"Gaps closed\t{len(gaps_df[gaps_df['type'] == 'ClosedGap'])}", file=outfile)
print(f"Partially filled\t{len(gaps_df[gaps_df['type'] == 'PartialFill'])}", file=outfile)
print(f'sequences with gaps filled\t{len(gaps_df.seq_id.unique())}', file=outfile)
print(f'deleted\t{gaps_df.deleted.sum()}bp ({gaps_df.deleted.sum()/df.Length.sum()*100:.2}%)', file=outfile)
print(f'added\t{gaps_df.added.sum()}bp ({gaps_df.added.sum()/df.Length.sum()*100:.2}%)', file=outfile)
print(file=outfile)

break_df = pd.DataFrame(breaks, columns=['seq_id', 'deleted', 'added', 'type'])

print(f'Breaks fixed\t{len(break_df)}', file=outfile)
print(f'sequences with break fix\t{len(break_df.seq_id.unique())}', file=outfile)
print(f'deleted\t{break_df.deleted.sum()}bp ({break_df.deleted.sum()/df.Length.sum()*100:.2}%)', file=outfile)
print(f'added\t{break_df.added.sum()}bp ({break_df.added.sum()/df.Length.sum()*100:.2}%)', file=outfile)

outfile.close()
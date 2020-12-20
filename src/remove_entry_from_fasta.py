#!/use/bin/env python3
import argparse
parser = argparse.ArgumentParser(description='Remove sequences from a fasta')
parser.add_argument('--infasta', type=str,help='fasta to remove from ')

parser.add_argument('--txToRemove', type=str, help='entries to remove')

parser.add_argument('--outfasta', type=str,help='file for filtered fasta ')
args=parser.parse_args()

with open(args.infasta) as infasta, open(args.txToRemove) as bad_tx, open(args.outfasta,'w+') as outfasta:
    names=set()
    for line in bad_tx:
        names.add('>'+line.strip())
    oldline=infasta.readline().strip()
    while oldline:
        if oldline not in names and '>' in oldline:
            write=True
        elif oldline in names and '>' in oldline:
            write=False
        if write:
            outfasta.write(oldline+'\n')
        oldline=infasta.readline().strip()

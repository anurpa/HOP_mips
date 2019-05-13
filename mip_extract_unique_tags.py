
#!/usr/bin/env python

__author__      = "Rebecca Barnard"


import argparse, gzip, string, sys
from numpy  import *

def pick_single_reads(sam):
  topreads = {}
  count = 0
  for read in sam:
    count += 1
    #if count > 10: break
    if read.startswith('#'): continue

    read = read.rstrip('\n').split('\t')
    name = read[0].split('-')
    # exclusion conditions: tags with Ns, unmapped reads
    if 'N' in name[-1]: continue 
    if read[1] not in ['0','16']: continue

    qual = sum(map(ord,list(read[10]))) # summarize read quality
    key = ':'.join(read[1:4]) + ':' + name[-1]  # key on flag, pos, tag

    if key in topreads:
      if topreads[key][1] < qual:
        topreads[key] = [read[0],qual]
    else:
      topreads[key] = [read[0],qual]
  # done walking samfile
  return topreads
# done with single-end read picking

def pick_paired_reads(sam):
  topreads = {}
  firstreads = {}
  for read in sam:
    if read.startswith('#'): continue
    read = read.rstrip('\n').split('\t')
    name = read[0].split('-')
    # exclusion condition: tags with Ns
    if 'N' in name[-1]: continue 

    qual = sum(map(ord,list(read[10]))) # summarize read quality

    if read[1] in ['83','99']: # mapped, properly paired, first in pair
      key = ':'.join(read[1:4]) + ':' + name[-1]  # key on flag, pos, tag
      firstreads[read[0]] = [key, name[-1], qual] # store key, tag, qual sum
    else:
      if read[0] in firstreads:
        key = firstreads[read[0]][0]
        qual_pair = firstreads[read[0]][2] + qual
        if key in topreads:
          if topreads[key][1] < qual_pair:
            topreads[key] = [read[0], qual_pair]
        else:
          topreads[key] = [read[0], qual_pair]
        del firstreads[read[0]]
  # done walking samfile
  return topreads
# done with paired-end read picking

def main():
  p = argparse.ArgumentParser(description='Extract unique tags from MIP alignments')

  p.add_argument('samfile', action='store', type=str,  
    help='SAM format MIP alignment')
  p.add_argument('-p', '--paired', action='store_true',
    help='reads are paired-end (default False)')

  args = p.parse_args()

  sam = open(args.samfile, 'r', 100000)
  if args.paired:
    topreads = pick_paired_reads(sam)
  else:
    topreads = pick_single_reads(sam)

  mynames={}
  for r in topreads:
    mynames[topreads[r][0]]=0

  #print len(mynames)

  #sam = open(args.samfile, 'r', 100000)
  sam.seek(0)

  count = 0

  for line in sam:
    count += 1
    #if count > 10: break
    if line.startswith('#'):
      print (line)
    else:
      line=line.rstrip('\n')
      temp=line.split('\t')
      if temp[0] in mynames:
        print (line)
    
  sam.close()
# end main function

if __name__=='__main__':
  main()



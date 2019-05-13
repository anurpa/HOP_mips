#!/home/groups/oroaklab/src/pyenv/versions/2.7.9/bin/python
__author__=" Rebecca Barnard "

import argparse, sys, gzip

# removes tags from fastq seuqence and appends to barcode

p = argparse.ArgumentParser(description='Split tags from FASTQ sequence')

p.add_argument('fq1', action='store', type=str, metavar='FILE',
  help='input fastq file for read1, default=None')
p.add_argument('fq2', action='store', type=str, metavar='FILE',
  help='input fastq file for read2, default=None')
p.add_argument('out1', action='store', type=str, metavar='FILE',
  help='output fastq file for read1, default=None')
p.add_argument('out2', action='store', type=str, metavar='FILE',
  help='output fastq file for read2, default=None')
p.add_argument('-l1', '--length1', action='store', type=int, default=3,
  help='total length of Read 1 MIP tag (default 3)')
p.add_argument('-l2', '--length2', action='store', type=int, default=3,
  help='total length of Read 2 MIP tag (default 3)')

args = p.parse_args()


if args.fq1.endswith(".gz"):
  r1_in = gzip.open(args.fq1, 'rb')
else:
  r1_in = open(args.fq1, 'r')


if args.fq2.endswith(".gz"):
  r2_in = gzip.open(args.fq2, 'rb')
else:
  r2_in = open(args.fq2, 'r')


r1_out = open(args.out1, 'w')
r2_out = open(args.out2, 'w')

while True:
  read1 = []
  read2 = []
  try:
    for i in range(0,4):
      read1.append(r1_in.next())
      read2.append(r2_in.next())
  except StopIteration:
    break

  tag = read1[1][0:args.length1] + read2[1][0:args.length2]

  name1 = read1[0].split("/")
  name2 = read2[0].split("/")
  read1[0] = name1[0] + "-" + tag + "/" + name1[1] # add tag to barcode in read name
  read2[0] = name2[0] + "-" + tag + "/" + name2[1]

  read1[1] = read1[1][args.length1:] # trim SEQ lines
  read2[1] = read2[1][args.length2:]
  read1[3] = read1[3][args.length1:] # trim QUAL lines
  read2[3] = read2[3][args.length2:]


  for r in read1: r1_out.write(r)
  for r in read2: r2_out.write(r)

# done processing fastq lines



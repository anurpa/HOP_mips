#!/usr/bin/env python

__author__=" Rebecca Barnard "

# input: arm info file, single-end SAM alignment
# output: modified single-end SAM alignment

# arm info file
# '-' mips get: chromosome, Start is ligation arm start, the 0 designation.  
#    len(ligation_arm), len(extention_arm) (Note sequenced on forward strand)
# '+' mips get: chromosome, Start is extension start, the 16 designation.  
#    len(ligation_arm), len(extention_arm) (Note sequenced on reverse strand)
# start is always most upstream coordinate

# mip arms take two
# single-end SAM input


import argparse, re, sys

def import_mips(filein):
  mip_dict = {}
  with open(filein) as mip_file:
    for mip in mip_file:
      if mip.startswith("CHR"): continue  # skip header
      mip = mip.rstrip("\n").split("\t")  # chrom pos strand arm_lig arm_ext
      key = mip[0] + ":" + mip[1] + ":" + mip[2]  # key on chrom:pos:strand
      #mip[1] = int(mip[1])  # pos to int
      mip[3] = int(mip[3])  # lig arm size to int
      mip[4] = int(mip[4])  # ext arm size to int

      # retain longest trim regions
      if key in mip_dict.keys():
        if mip_dict[key][0]<mip[3]: mip_dict[key][0] = mip[3]
        if mip_dict[key][1]<mip[4]: mip_dict[key][1] = mip[4]
      else:
        mip_dict[key] = [mip[3],mip[4]]  # lig arm, ext arm lengths
  # done parsing MIP file
  return mip_dict
# end MIP import function


# NOTE: beware of off-by-one...
def mip_info_adjust(mip_dict, outfile):
  with open(outfile,'w') as out:
    for key in mip_dict:
      line = key.split(":")
      line.extend(['0','0'])
      if line[2] == '0':
        line[1] = str( int(line[1]) + mip_dict[key][0] )
      else:
        line[1] = str( int(line[1]) + mip_dict[key][1] )
      out.write("\t".join(line) + "\n")
    # done processing MIPs
  # close outfile
# end MIP info output function


def trim_se_reads(filein, mip_dict, isize, mapq):
  if filein is not None:
    samfile = open(filein)
  else:
    samfile = sys.stdin

  for read in samfile:
    if read.startswith("@"): continue   # skip header lines
    read = read.rstrip('\n').split('\t')
    name = read[0]
    flag = read[1]

    # exclusion criteria
    if flag not in ['0','16']: continue  # single-end reads only
    if read[4]<mapq: continue            # skip low-MAPQ reads
    if "S" in read[5]: continue          # skip soft-clipped reads
    # there don't seem to be many of those
    key = read[2] + ":" + read[3] + ":" + flag
    if key not in mip_dict: continue     # skip reads not matching MIPs

    # pull trim sizes
    if flag == '0':  # forward read, ligation arm first
      trim_for = mip_dict[key][0]
      trim_rev = mip_dict[key][1]
    else:            # reverse read, extension arm first
      trim_for = mip_dict[key][1]
      trim_rev = mip_dict[key][0]

    # break out cigar vals and notation
    cigar_vals = [int(v) for v in re.split('[A-Z]', read[5]) if v != '']
    cigar_types = re.findall('[A-Z]', read[5])
    cigar_list = []
    for i in range(len(cigar_vals)):
      cigar_list.append( (cigar_vals[i], cigar_types[i]) )

    # skip reads with indels in the arms
    # they make life difficult, and there aren't very many anyway
    if cigar_list[0][1] != "M" or cigar_list[0][0]<trim_for:
      continue
    if cigar_list[-1][1] != "M" or cigar_list[-1][0]<trim_rev:
      continue

    # also exclude odd-sized DNA fragments
    if isize != 0:
      alnsize = sum( [v for v,t in cigar_list if t != "I"] )
      if alnsize != isize: continue

    # now do actual trimming
    read[3] = str(int(read[3]) + trim_for)           # move alignment start
    read[9] = read[9][trim_for:len(read[9])-trim_rev]    # clip sequence
    read[10] = read[10][trim_for:len(read[10])-trim_rev] # clip quality scores

    # and fix the cigar string
    # note we've guaranteed the first match is >= than the trim distance

    if cigar_list[0][0] > trim_for:
      cigar_list[0] = ( cigar_list[0][0] - trim_for, cigar_list[0][1] )
    elif cigar_list[0][0] == trim_for:
      cigar_list = cigar_list[1:]

    if cigar_list[-1][0] > trim_rev:
      cigar_list[-1] = ( cigar_list[-1][0]-trim_rev, cigar_list[-1][1] )
    elif cigar_list[-1][0] == trim_rev:
      cigar_list.pop()

    # update cigar string in read
    read[5] = "".join( [str(v)+t for v,t in cigar_list] )

    sys.stdout.write("\t".join(read) + "\n")
    # tmp pass-through

# done processing SAM reads


def main():
  p = argparse.ArgumentParser(description='Trims MIP arm sequences.')

  p.add_argument('mips', action='store', type=str,  
    help='tab-delimited MIP arm info file')
  p.add_argument('-s','--sam', action='store', type=str, default=None,
    help='SAM format aligned reads; default stdin')
  p.add_argument('-i','--insert', action='store', type=int, default=94,
    help='expected insert size; default 94; 0 to disable filter')
  p.add_argument('-m','--mip_out', action='store', type=str, default=None,
    help='write trimmed MIP positions to specified file')
  p.add_argument('-q','--mapq', action='store', type=int, default=0,
    help='require min MAPQ to output read (default 0)')

  args = p.parse_args()

  mip_dict = import_mips(args.mips)

  if args.mip_out is not None:
    mip_info_adjust(mip_dict,args.mip_out)

  trim_se_reads(args.sam, mip_dict, args.insert, args.mapq)


# done with main function


if __name__=='__main__':
  main()



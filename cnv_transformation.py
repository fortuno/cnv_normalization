import argparse
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cnv_input', type=str, required=True,
      help='Segment CNV input file')
    parser.add_argument('--gencode', type=str, required=True,
      help='Gencode v22 reference')
    parser.add_argument('--cnv_output', type=str, required=True,
      help='Output file with weighted mean CNV per gene')            
    args = parser.parse_args()

    return args

def parse_cnv(cnv_input):

    key_data = []
    delimiter = '\t'
    header = None

    with open(cnv_input,'rt') as f:
      for line in f:
          if not header:
             header = line.strip('\n').split(delimiter)
          else:
             line_data = dict(zip(header, line.strip('\n').split(delimiter)))
             key_data.append(line_data)

    return key_data

def calculate_gene_count(cnv_data, gencode):

    delimiter = '\t'
    header = None

    count = dict()
    gene_len   = 0
    gene_nocov = 0
    gene_id    = None
    outfile = open(args.cnv_output, 'wb')
    with open(gencode,'rt') as f:
      for line in f:                 
          line_data  = line.strip('\n').split(delimiter)
          
          # Calculate weighted CNV before changing gene_id
          if gene_id and gene_id != line_data[0]:
                count[gene_id] = np.log2((count[gene_id] + gene_nocov) / gene_len)
                outfile.write("%s\t%.4f\n" % (gene_id, count[gene_id]))       
         
          # Get info from genecode aggregated exomes
          gene_id    = line_data[0] 
          chrom      = line_data[1].replace('chr', '')
          exon_start = int(line_data[2])
          exon_end   = int(line_data[3])
          if not gene_id in count:
              count[gene_id] = 0
              gene_len = 0  
              gene_nocov = 0
          gene_len   += (exon_end - exon_start) 
          gene_nocov += (exon_end - exon_start) 

          # Calculate weight according to overlapping between exon and segments
          for seg in cnv_data:
              if chrom == str(seg['Chromosome']):
                 weight = 0
                 seg_start = int(seg['Start'])
                 seg_end   = int(seg['End'])

                 # Case 1: exon is completely contained in segment
                 if seg_start < exon_start and seg_end > exon_end:   
                    weight = exon_end - exon_start   
                 if seg_start > exon_start and seg_start < exon_end:
                    # Case 2: segment is completely contained in exon
                    if seg_end > exon_start and seg_end < exon_end:
                       weight = seg_end - seg_start
                    # Case 3: exon is partially contained in segment (exon starts first)   
                    else:
                       weight = exon_end - seg_start
                 # Case 4: exon is partially contained in segment (segment start firts)      
                 elif seg_end > exon_start and seg_end < exon_end: 
                       weight = seg_end - exon_start                 
                 count[gene_id] += (2**float(seg['Segment_Mean']))*weight
                 gene_nocov -= weight 

      # Calculated weighted CNV for final gene
      count[gene_id] = np.log2((count[gene_id] + gene_nocov) / gene_len)
      outfile.write("%s\t%.4f\n" % (gene_id, count[gene_id])) 
      outfile.close()


if __name__ == '__main__':
    args = parse_args()

    cnv_data = parse_cnv(args.cnv_input)
    calculate_gene_count(cnv_data, args.gencode)


import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genmodel', type=str, required=True,
      help='Gen model file')    
    parser.add_argument('--output', type=str, required=True,
      help='File with aggregated exons')     
    args = parser.parse_args()

    return args

def aggregate_exons(exons, start_values, gene_id, fout):
    aggr_chr = None
    aggr_end = 0
    aggr_start = 0
    for i in sorted(enumerate(start_values), key=lambda x:x[1]):                    
       index = i[0] 
       exon = exons[index]
       exon_chr   = exon[0]
       exon_start = exon[3]
       exon_end   = exon[4]        
       if not aggr_chr:
          aggr_chr   = exon_chr
          aggr_start = exon_start
          aggr_end   = exon_end
       elif exon_start > aggr_end:
          fout.write("{0}\t{1}\t{2}\t{3}\n".format(gene_id, aggr_chr, aggr_start, aggr_end))  
          aggr_chr   = exon_chr
          aggr_start = exon_start
          aggr_end   = exon_end      
       else:
          if exon_end > aggr_end:
             aggr_end = exon_end

    fout.write("{0}\t{1}\t{2}\t{3}\n".format(gene_id, aggr_chr, aggr_start, aggr_end))  

def parse_gencode(gencode_file, output):

    delimiter = '\t'
    prev_gene = None
    exons = []
    start_values = []
    fout = open(output, 'w')    
    with open(gencode_file,'rt') as f:
      for line in f:
          if line[0:2] == '##':
             print line.strip('\n')
          else:
             line_data = line.strip('\n').split(delimiter)
             if line_data[2] == 'exon':
                 gene_id = line_data[-1].split(';')[0].split('"')[1]    
                 if prev_gene and gene_id != prev_gene:
                    # Aggregate exons in previous gene
                    print "Aggregating " + prev_gene + "..."
                    aggregate_exons(exons, start_values, prev_gene, fout)                  
                    
                    # Start new gene
                    prev_gene = gene_id
                    exons = []
                    start_values = []
                    exons.append(line_data)    
                    start_values.append(line_data[3])                     
                 else:
                    exons.append(line_data)    
                    start_values.append(line_data[3]) 
                    prev_gene = gene_id
  
      print "Aggregating " + gene_id + "..."
      aggregate_exons(exons, start_values, gene_id, fout) 
      fout.close()


def parse_gaf(gaf_file, output):

    delimiter = '\t'
    prev_gene = None
    exons = []
    start_values = []
    fout = open(output, 'w')    
    with open(gaf_file,'rt') as f:
      for line in f:
          if line[0] == '#':
             headers = line.strip('\n').split(delimiter)
             geneidx = headers.index('FeatureID')
             typeidx = headers.index('FeatureType')
             exonidx = headers.index('CompositeCoordinates')
          else:
             line_data = line.strip('\n').split(delimiter)
             if line_data[typeidx] == 'gene':
                 gene_id = line_data[geneidx] #.split('|')[0]   
                 exons = line_data[exonidx].split(',')
                 first = exons[0].split(':')
                 exons[0] = first[1]
                 chrom = first[0]
                 for exon in exons:
                    pos = exon.split('-')
                    pos[1] = pos[1].replace(':+', '')
                    pos[1] = pos[1].replace(':', '')                   
                    fout.write("%s\t%s\t%s\t%s\n" % (gene_id, chrom, pos[0], pos[1]))
    fout.close()


if __name__ == '__main__':
    args = parse_args()

    if os.path.splitext(args.genmodel)[1] == '.gtf':
        parse_gencode(args.genmodel, args.output)
    else:
        parse_gaf(args.genmodel, args.output)    


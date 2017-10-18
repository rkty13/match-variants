import glob
import gzip
import json
from copy import deepcopy
from pprint import pprint

from vcf2clinvar import *
from vcf2clinvar.common import *

def _next_line(filebuffer):
    try:
        next_line = filebuffer.readline()
    except AttributeError:
        next_line = filebuffer.next()
    try:
        next_line = next_line.decode('utf-8')
        return next_line
    except AttributeError:
        return next_line

def parse_clinvar_data(filepath):
    print('Parsing clinvar data...')
    clinvar_file = gzip.open(filepath)
    clinvar_cur_line = _next_line(clinvar_file)
    while clinvar_cur_line.startswith('#'):
        clinvar_cur_line = _next_line(clinvar_file)
    clinvar_data = []
    while clinvar_cur_line:
        clinvar_vcf_line = ClinVarVCFLine(vcf_line=clinvar_cur_line)
        clinvar_data.append(clinvar_vcf_line)
        clinvar_cur_line = _next_line(clinvar_file)
    print('Done parsing clinvar data')
    return clinvar_data

def aggregate_variant_data(variant_dir):
    print('Parsing variant data...')
    variant_paths = glob.glob(variant_dir + '/*.json')
    variant_data = []
    for variant_path in variant_paths:
        with open(variant_path) as variant_file:
            data = json.load(variant_file)
        variant_data.append(data)
    print('Done parsing variant data')
    return variant_data

def _determine_clinvar_good(row):
    for allele in row.alleles:
        if hasattr(allele, 'records') and \
                len([r.sig for r in allele.records \
                if r.sig not in ('0', '1', '2', '3', '255')]):
            return True
    return False

def _generate_gennotes_url(row):
    return 'https://gennotes.herokuapp.com/genevieve-edit/?build=b37&amp;chrom={0}&amp;pos={1}&amp;ref_allele={2}&amp;var_allele={3}'.format(row.chrom, row.start, row.ref_allele, row.alt_alleles)

def _match_variant_clinvar(variant, row):
   pass 

def find_clinvar_in_variant(variant_data):
    data = []
    for i in range(len(variant_data)):
        if 'Summary' in variant_data[i] and \
                'clinvar_data' not in variant_data[i] and \
                '/clinvar/RCV' in variant_data[i]['Summary']:
            data.append(variant_data[i])
    return data

def associate_variant_clinvar(clinvar_data, variant_data):
    print('Associating variant data...')
    modified_variant_data = deepcopy(variant_data)
    total_num = 0
    num_matched = 0
    for i in range(len(modified_variant_data)):
        total_num += 1
        variant = modified_variant_data[i]
        if ('Build 37 Chromosome' in variant and \
                'Build 37 Position' in variant and \
                'Build 37 Variant Allele' in variant) or \
                'dbSNP IDs' in variant:
            print('-------------------------- VARIANT VALID: -------------------------')
            for row in clinvar_data:
                if ('Build 37 Chromosome' in variant and \
                        'Build 37 Position' in variant and \
                        'Build 37 Variant Allele' in variant and \
                        CHROM_INDEX[str(variant['Build 37 Chromosome'])] == CHROM_INDEX[str(row.chrom)] and \
                        int(variant['Build 37 Position']) == row.start and \
                        variant['Build 37 Variant Allele'] in row.alt_alleles) or \
                    ('dbSNP IDs' in variant and \
                        row.dbsnp_id in variant['dbSNP IDs']):
                    num_matched += 1
                    print('It\'s a match!')
                    variant['clinvar_data'] = deepcopy(row).as_dict()
                    variant['clinvar_good'] = _determine_clinvar_good(row)
                    variant['gennotes_url'] = _generate_gennotes_url(row)
                    break
    print('Matched: ', num_matched, '/', total_num)
    return modified_variant_data

if __name__ == '__main__':
    # Load clinvar data and variant data into memory
    clinvar_data = parse_clinvar_data('./clinvar.vcf.gz')
    variant_data = aggregate_variant_data('./variant_data')
    associated_variant_data = associate_variant_clinvar(clinvar_data, variant_data)
    with open('final_variant_clinvar_data.json', 'w') as f:
        json.dump(associated_variant_data, f, indent=4)
    variants = find_clinvar_in_variant(associated_variant_data)
    print(len(variants))
    with open('clinvar_in_variant.json' ,'w') as f:
        json.dump(variants, f, indent=4)

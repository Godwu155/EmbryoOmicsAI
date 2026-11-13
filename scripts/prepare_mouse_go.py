"""Prepare direct mouse BP annotation sets, without ontology propagation.

This deliberately does not combine a GAF with a mismatched ontology snapshot.
"""
import argparse
from collections import defaultdict
import gzip
from pathlib import Path
import yaml

from embryoomics.cellxgene import sha256


def prepare(gaf, out):
    out = Path(out)
    if out.exists():
        raise ValueError('Choose a new output directory')
    sets, ids = defaultdict(set), defaultdict(set)
    headers = []
    with gzip.open(gaf, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('!'):
                headers.append(line.rstrip())
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 15:
                raise ValueError('Malformed GAF row')
            if fields[0] != 'MGI' or fields[8] != 'P' or fields[12].split('|')[0] not in ('taxon:10090', 'NCBITaxon:10090'):
                continue
            if 'NOT' in fields[3].split('|') or fields[6] == 'ND':
                continue
            ids[fields[2]].add(fields[1])
            sets[fields[4]].add(fields[2])
    ambiguous = {symbol for symbol, values in ids.items() if len(values) > 1}
    out.mkdir(parents=True)
    gmt = out / 'mouse_GO_BP_direct.gmt'
    gmt.write_text(''.join(term + '\thttps://amigo.geneontology.org/amigo/term/' + term + '\t' +
                   '\t'.join(sorted(genes - ambiguous)) + '\n' for term, genes in sorted(sets.items()) if genes - ambiguous), encoding='utf-8')
    metadata = {'species': 'mus_musculus', 'id_type': 'gene_symbol',
        'source_url': 'https://current.geneontology.org/annotations/gaf/MOUSE-mod.gaf.gz',
        'version': '; '.join(h for h in headers if h.startswith(('!date-generated:', '!go-version:'))),
        'license': 'CC-BY-4.0; https://www.geneontology.org/page/use-and-license',
        'gaf_sha256': sha256(gaf), 'gmt_sha256': sha256(gmt), 'headers': headers,
        'method': 'Direct biological-process annotations only; mouse MGI; exclude NOT and ND; include other evidence codes; NO ancestor propagation',
        'excluded_ambiguous_gaf_symbols': sorted(ambiguous)}
    (out / 'metadata.yaml').write_text(yaml.safe_dump(metadata, allow_unicode=True, sort_keys=False), encoding='utf-8')
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gaf', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    print(yaml.safe_dump(prepare(args.gaf, args.out), sort_keys=False))

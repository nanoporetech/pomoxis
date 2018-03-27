import argparse
import os
import pysam
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser('Trim alignments in multiple bams to common overlap window.')
    parser.add_argument('bams', nargs='+', help='input bam files')
    parser.add_argument('-r', '--ref_name',
                        help='Reference to process, only needed if bams contain multiple references.')
    parser.add_argument('-o', '--output', default='trimmed.fasta',
                        help='Output fasta file.')

    args = parser.parse_args()

    refs = set()
    for bam in args.bams:
        with pysam.AlignmentFile(bam) as b:
            refs.update(b.references)

    if args.ref_name is None:
        if len(refs) > 1:
            raise RuntimeError('Bams contain multiple references, '
                               +'use the --ref_name argument to specify one of {}.'.format(refs))
        else:
            args.ref_name = list(refs)[0]
    elif args.ref_name not in refs:
        raise KeyError('Ref {} not in bam refs {}'.format(args.ref_name, refs))

    bam_files = { bam: pysam.AlignmentFile(bam) for bam in args.bams}
    reads = {bn: list(bam.fetch(args.ref_name)) for bn,bam in bam_files.items()}
    if not all([ len(l) == 1 for l in reads.values()]):
        raise RuntimeError('Expected just one read per bam, do not chunk and filter to primary alignments')

    reads = {bn: l[0] for bn,l in reads.items()}
    print('Initial alignments:')
    for bn, read in reads.items():
        print('{} {} {} {}'.format(bn, read.query_name, read.reference_start, read.reference_end))
    start = max([r.reference_start for r in reads.values()])
    end = min([r.reference_end for r in reads.values()])
    print('Initial ref start {} ref end {}'.format(start, end))

    # trim back to point where each read is mapped to a common ref
    ref_pos = { k: read.get_reference_positions(full_length=True) for k, read in reads.items()}
    while True:
        start_in_ref_pos = [start in p for p in ref_pos.values()]
        if all(start_in_ref_pos):
            break
        else:
            start += 1
            print('Shifing forward start')

    while True:
        end_in_ref_pos = [end - 1 in p for p in ref_pos.values()]
        if all(end_in_ref_pos):
            break
        else:
            end -= 1
            print('Shifing back start')

    print('Final ref start {} ref end {}'.format(start, end))

    seq_objs = []
    for bn, read in reads.items():
        ref_pos = read.get_reference_positions(full_length=True)
        q_start = ref_pos.index(start)
        q_end = ref_pos.index(end - 1) + 1
        seq = read.query_sequence[q_start:q_end]
        prefix = os.path.splitext(os.path.basename(bn))[0]
        read_id='{}_{}_{}_{}'.format(prefix, read.query_name, q_start, q_end)
        seq_ob = SeqIO.SeqRecord(Seq(seq), id=read_id)
        seq_objs.append(seq_ob)

    with open(args.output, 'w') as fh:
        SeqIO.write(seq_objs, fh, 'fasta')

    for b in bam_files.values():
        b.close()


if __name__ == '__main__':
    main()

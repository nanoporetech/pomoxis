import argparse
import os
import pysam
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser(
        prog='trim_alignments',
        description='Trim alignments in multiple bams to common overlap window.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bams', nargs='+', help='input bam files')
    parser.add_argument('-r', '--ref_name',
                        help='Reference to process, only needed if bams contain multiple references.')
    parser.add_argument('-o', '--output_prefix', default='trimmed',
                        help='Prefix for outputs.')
    parser.add_argument('-f', '--reference_fasta', default=None,
                        help='Reference fasta to trim to alignment window.')

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

    if len(args.bams) > 1:  # expect just 1 read per bam
        reads = {bn: list(bam.fetch(args.ref_name)) for bn,bam in bam_files.items()}
        if not all([ len(l) == 1 for l in reads.values()]):
            raise RuntimeError('Expected just one read per bam, do not chunk and filter to primary alignments')
        reads = {bn: l[0] for bn,l in reads.items()}
    else: # expect more than 1 read
        reads = {r.query_name: r for r in bam_files[args.bams[0]].fetch()}
        if not len(reads.values()) > 1:
            raise RuntimeError('The bam {} contained just 1 read.'.format(args.bam[0]))

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

    print('Final trimmed region: {}:{}-{}'.format(args.ref_name, start, end))

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

    output = '{}_queries.fasta'.format(args.output_prefix)
    with open(output, 'w') as fh:
        SeqIO.write(seq_objs, fh, 'fasta')

    for b in bam_files.values():
        b.close()

    if args.reference_fasta is not None:  # output trimmed reference
        output_ref = '{}_reference.fasta'.format(args.output_prefix)
        ndx = SeqIO.index(args.reference_fasta, 'fasta')
        if not args.ref_name in ndx.keys():
            raise KeyError('Reference {} not in {}'.format(args.ref_name, args.fasta))

        with open(output_ref, 'w') as fh:
            ref_id='{}_{}_{}'.format(args.ref_name, start, end)
            trimmed_ref = SeqIO.SeqRecord(ndx[args.ref_name].seq[start: end], ref_id)
            SeqIO.write([trimmed_ref], fh, 'fasta')
        ndx.close()


if __name__ == '__main__':
    main()

import argparse

_field_name_map_ = {
    'del': 'Q(del)',
    'err_ont': 'Q(acc)',
    'idel': 'Q(iden)',
    'ins': 'Q(ins)',
    'iden': 'Q(iden)',
}

_cols_ = ['name', 'ref', 'ref_cover', 'Q(acc)', 'Q(iden)', 'Q(del)', 'Q(ins)']

_format_str_ = '\t'.join(['{' + c + '}' for c in _cols_])

_median_col_ind_ = {True: 3, False: 1}


def get_ref_cover(lines, ref):
    for l in lines:
        s = l.split(' ')
        if s[0] == ref:
            return s[1]


def extract_vals(vals, name_map, col_ind):
    results = {}
    for line in vals:
        s = line.split()
        name = s[0]
        val = s[col_ind]
        if name in name_map:
            results[name_map[name]] = val
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Extract Q scores from summary_from_stats output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('summaries', nargs='+', help='*summ.txt created by summary_from_stats')
    parser.add_argument('--median', help='Use median. If false, use mean.', action='store_true', default=False)
    parser.add_argument('--ref', default=None, help='process single ref, rather than overall result')
    args = parser.parse_args()
    qscore_from_summary(args.summaries, args.median, args.ref)


def qscore_from_summary(summaries, median, ref):
    if ref is None:
        match = '#  Q Scores'
    else:
        match = '# {} Q Scores'.format(ref)

    col_ind = _median_col_ind_[median]

    print('\t'.join(_cols_))

    results = []
    for f in sorted(summaries):
        with open(f) as fh:
            lines = [l.strip() for l in fh.readlines()]
        if match not in lines:
            continue
        ind = lines.index(match)
        vals = lines[ind + 2: ind + 7]
        data = extract_vals(vals, _field_name_map_, col_ind)
        data['name'] = f
        if ref is None:
            data['ref_cover'] = 'NA'
            data['ref'] = 'all'
        else:
            data['ref_cover'] = get_ref_cover(lines, ref)
            data['ref'] = ref
        print(_format_str_.format(**data))
        results.append(data)

    return results


if __name__ == "__main__":
    main()

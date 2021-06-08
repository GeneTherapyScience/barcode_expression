#!/user/bin/env python3
from merge_readerror import *
from alignment import *
from random import randrange
from collections import Counter, defaultdict
from tqdm import tqdm, trange

Nsampling = 10
Npair=0
distance_bound = 30

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nsampling', type=int, default=Nsampling,
                        help='give Nsampling.')
    parser.add_argument('-p', '--pairs', type=int, default=Npair,
                        help='give Nsampling.')
    parser.add_argument('-t', '--template', default=None,
                        help='give template to compare with.')
    parser.add_argument('-d', '--detail', default=None,
                        help='give template to compare with.')
    parser.add_argument('--noheader', action='store_true',
                        help='the input does not include header.')
    args = parser.parse_args()
    Nsampling = args.nsampling
    Npair = args.pairs

    header, data = inputdata(has_header=(not args.noheader))
    N = len(data)
    if args.template:
        count = defaultdict(int)
        detail = dict()
        for i in trange(N):
            seq, readnum = data[i][:2]
            if seq in detail:
                detail[seq][2] += readnum
            else:
                align, mut = seq_alignment(args.template, seq)
                detail[seq] = [align, mut, readnum]
            count[mut] += readnum
        print('score', 'reads', sep='\t')
        for k in sorted(count.keys()):
            print(k, count[k], sep='\t')
        if args.detail:
            with open(args.detail, 'w') as f:
                print('template:', args.template, file=f)
                print('sequence', 'alignment', 'score', 'reads', sep='\t', file=f)
                for k, v in sorted(detail.items(), key=lambda x: (x[1][1],-x[1][2],x[0])):
                    print(k, *v, sep='\t', file=f)
    elif Npair > 0:
        count = defaultdict(int)
        for _ in trange(Npair):
            while True:
                m, n = randrange(N), randrange(N)
                if m != n:
                    break
            count[levenshtein_distance(data[m][0], data[n][0])] += 1
        print(sorted(count.items()), sep='\t')
    else:
        results = []
        interest_samples = [randrange(N) for _ in range(Nsampling)]
        for i in trange(Nsampling):
            m = interest_samples[i]
            barcode, readnum, _ = data[m]
            count = defaultdict(int)
            for n in range(N):
                if n != m:
                    count[levenshtein_distance(data[n][0], barcode)] += 1
            results.append((barcode, count))

        for barcode, count in results:
            print(barcode, sorted(count.items()), sep='\t')

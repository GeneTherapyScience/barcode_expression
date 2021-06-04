import numpy as np

DNA_match_pairs = set()
for b in 'ATGC':
    DNA_match_pairs.add(b+b)
for n in 'ATGCWSN':
    DNA_match_pairs |= {'N'+n, n+'N'}
for w in 'ATW':
    DNA_match_pairs |= {'W'+w, w+'W'}
for s in 'GCS':
    DNA_match_pairs |= {'S'+s, s+'S'}

match_letters = 'MIXD'
def seq_alignment(base, target, gap=2.5, extend=0.5, substitution=1, bound=None):
    M, N = len(base), len(target)
    gap = int(gap*100)
    extend = int(extend*100)
    substitution = int(substitution*100)
    if bound is None:
        bound = (substitution+gap+extend)*(M+N+1)
    prev = np.zeros((N+1,3), dtype=int)
    cur = np.array([[i, bound, bound] for i in range(1,N+1)] + [(0,0,0)], dtype=int)
    # (normal, in-insertion, in-deletion)
    parent = np.zeros((M,N,3,3), dtype=int)
    for m in range(M):
        k = cur.min()
        if k > bound:
            return k
        cur, prev = prev, cur
        cur[:] = 0
        cur[-1] = (m+1, bound, bound)
        for n in range(N):
            candidates = [
                prev[n-1] + int(base[m]+target[n] not in DNA_match_pairs)*substitution,
                np.array([
                    cur[n-1][0] + gap, # new insertion
                    cur[n-1][1] + extend, # continue insertion
                    cur[n-1][2] + gap, # new insetion
                ], dtype=int),
                np.array([
                    prev[n][0] + gap, # new deletion
                    prev[n][1] + gap, # continue deletion
                    prev[n][2] + extend, # new deletion
                ], dtype=int),
            ]
            k = [candidates[i].argmin() for i in range(3)]
            parent[m,n] = np.array([[m-1,n-1,k[0]], [m,n-1,k[1]], [m-1,n,k[2]]], dtype=int)
            cur[n] = [candidates[i][k[i]] for i in range(3)]

        m, n = M-1, N-1
        i = cur[N-1].argmin()
        distance = cur[N-1][i]/substitution
        backarr = []
        while m >= 0:
            if i == 0:
                if base[m]+target[n] in DNA_match_pairs:
                    match = 'M'
                else:
                    match = 'X'
            elif i == 1:
                match = 'I'
            else:
                match = 'D'
            backarr.append(match)
            m, n, i = parent[m,n,i]
        matchseq = ''.join(backarr[::-1])

    return matchseq, distance


if __name__ == '__main__':
    base = 'GGTGGCTTTACCAACAGTAC'
    target = 'GATTCATCTCATCTATCAGAAAATAAATAAA'
    alignment, score = seq_alignment(base, target)
    print('base:', base)
    print('target:', target)
    print('alignment:', alignment)
    print('score:', score)

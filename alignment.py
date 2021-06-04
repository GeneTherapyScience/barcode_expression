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

match_letters = 'MXID'
def seq_alignment(base, target, gap=2.5, extend=0.5, substitution=1, bound=None):
    M, N = len(base), len(target)
    unit = 10000
    epsilon = 1
    gap = int(gap*unit)
    extend = int(extend*unit) - epsilon
    substitution = int(substitution*unit)
    if bound is None:
        bound = (substitution+gap+extend)*(M+N+1)
    prev = np.zeros((N+1,4), dtype=int) # first cur
    cur = np.ones((N+1,4), dtype=int) * bound # first prev
    cur[-1][1] = 0
    # (normal, in-insertion, in-deletion)
    parent = np.zeros((M,N,4,3), dtype=int)
    for m in range(M):
        k = cur.min()
        if k > bound:
            return k
        cur, prev = prev, cur
        cur[-1] = bound
        for n in range(N):
            base_match = int(base[m]+target[n] in DNA_match_pairs)
            candidates = [
                prev[n-1] + (1-base_match)*bound,
                prev[n-1] + substitution + base_match*bound,
                cur[n-1] + gap, # new insertion
                prev[n] + gap, # new deletion
            ]
            candidates[0][0] -= epsilon
            candidates[2][2] += extend - gap # continue insertion
            candidates[3][3] += extend - gap # delete insertion
            k = [candidates[i].argmin() for i in range(4)]
            parent[m,n] = np.array([[m-1,n-1,k[0]], [m-1,n-1,k[1]],[m,n-1,k[2]], [m-1,n,k[3]]], dtype=int)
            cur[n] = [candidates[i][k[i]] for i in range(4)]

    m, n = M-1, N-1
    i = cur[N-1].argmin()
    print(cur[N-1][i]/substitution)
    distance = round(cur[N-1][i]/substitution*10)/10
    backarr = []
    while m >= 0:
        backarr.append(match_letters[i])
        m, n, i = parent[m,n,i]
    matchseq = ''.join(backarr[::-1])

    return matchseq, distance


if __name__ == '__main__':
    base = 'GGTGGCTTTACCAACAGTAC'
    # target = 'GATTCATCTCATCTATCAGAAAATAAATAAA'
    target = 'GGCTTTACCAACAGTAC'
    # target = 'GGTGGCTTTACCAACAGTAC'
    alignment, score = seq_alignment(base, target)
    print('base:', base)
    print('target:', target)
    print('alignment:', alignment)
    print('score:', score)

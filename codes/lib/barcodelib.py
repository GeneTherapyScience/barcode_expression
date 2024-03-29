import os
import sys
from collections import defaultdict
import shutil
import itertools

ATGC = sorted('ATGC')
CGTA = ATGC[::-1]

def levenshtein_neighbors(barcode, distance, inserts=ATGC, min_distance=0):
    if distance == 0 and min_distance <= 0:
        yield barcode
    else:
        prev = list()
        done = set()
        if min_distance > 0:
            for b in levenshtein_neighbors(barcode, min_distance-1):
                done.add(b)
        for b in levenshtein_neighbors(barcode, distance-1):
            if b not in done:
                prev.append(b)
                done.add(b)
                yield b
        for b in prev:
            N = len(b)
            for i in range(N):
                c = b[:i] + b[i+1:]
                if not c in done:
                    done.add(c)
                    yield c
                if b[i] == 'N':
                    continue
                for n in inserts:
                    if n != b[i]:
                        c = b[:i] + n + b[i+1:]
                        if not c in done:
                            done.add(c)
                            yield c
            for i in range(N+1):
                for n in inserts:
                    c = b[:i] + n + b[i:]
                    if not c in done:
                        done.add(c)
                        yield c

def combinations(N,k):
    if 0 <= k <= N:
        c = [0]*k
        pool = [(-1,-1)]
        while pool:
            i, n = pool.pop()
            if i >= 0:
                c[i] = n
            if i >= k-1:
                yield tuple(c)
            else:
                for m in range(i+N-k+1,n,-1):
                    pool.append((i+1, m))

def hamming_neighbors(barcode, distance, inserts=ATGC, min_distance=0):
    N = len(barcode)
    for k in range(min_distance, distance+1):
        for c in combinations(N,k):
            if {barcode[c[i]] for i in range(k)} - set(inserts):
                continue
            arr = list(barcode)
            for r in itertools.product(ATGC, repeat=k):
                for i in range(k):
                    if barcode[c[i]] == r[i]:
                        break
                    else:
                        arr[c[i]] = r[i]
                else:
                    yield ''.join(arr)

DNA_match_pairs = set()
for b in 'ATGC':
    DNA_match_pairs.add(b+b)
for n in 'ATGCWSN':
    DNA_match_pairs |= {'N'+n, n+'N'}
for w in 'ATW':
    DNA_match_pairs |= {'W'+w, w+'W'}
for s in 'GCS':
    DNA_match_pairs |= {'S'+s, s+'S'}

def levenshtein_distance(barcode0, barcode1, bound=None):
    M, N = len(barcode0), len(barcode1)
    if bound is None:
        bound = M + N
    if M < N:
        M, N, barcode0, barcode1 = N, M, barcode1, barcode0
    cur = [bound]*N + [0]
    for m in range(M):
        k = min(cur[:-1])
        if k > bound:
            return k
        prev = cur
        cur = [0]*(N+1)
        cur[-1] = bound
        for n in range(N):
            cur[n] = min(
                prev[n-1] + int(barcode0[m]+barcode1[n] not in DNA_match_pairs),
                cur[n-1] + 1,
                prev[n] + 1,
                )
    return cur[N-1]

def hamming_distance(barcode0, barcode1, bound=None):
    M, N = len(barcode0), len(barcode1)
    if bound is None:
        bound = max(M,N)
    if M != N:
        return bound
    d = 0
    for i in range(M):
        d += int(barcode0[i]+barcode1[i] not in DNA_match_pairs)
        if d > bound:
            break
    return d

def N_candidates(barcode):
    k = barcode.count('N')
    if k == 0:
        yield barcode
    else:
        barcode_list = list(barcode)
        Npos = [i for i, x in enumerate(barcode) if x == 'N']
        pool = [(n, 0) for n in CGTA]
        while pool:
            c, i = pool.pop()
            barcode_list[Npos[i]] = c
            if i == k-1:
                yield str(barcode_list)
            else:
                for n in CGTA:
                    pool.append((n, i+1))

def inputs(f=sys.stdin):
    while True:
        line = f.readline().rstrip('\n')
        if not line: # EOF
            break
        else:
            yield line

def readref(referencefile, has_header=True, column=0):
    ret = set()
    if referencefile:
        with open(referencefile) as f:
            if has_header:
                f.readline() # drop header
            for line in inputs(f):
                line = line.strip()
                if line[0] == '#' or len(line)==0:
                    continue
                else:
                    barcode = line.split()[column]
                    ret.add(barcode)
    return ret

def inputdata(f=sys.stdin, has_header=True, reads_sort=True):
    if has_header:
        header = f.readline().rstrip()
    else:
        header = None
    data = []
    for line in inputs(f):
        line = line.strip()
        if not line:
            continue
        record = line.split()
        if len(record) >= 3:
            barcode, readnum, mutations = record[:3]
        elif len(record) == 2:
            barcode, readnum = record
            mutations = 0
        else:
            barcode = record[0]
            readnum = 1
            mutations = 0

        readnum = 0 if readnum=='NA' else int(readnum)
        mutations = 0 if mutations=='NA' else int(mutations)
        data.append((barcode, readnum, mutations))
    if reads_sort:
        data.sort(key=lambda x: (x[0].count('N'), -x[1], x[0]))
    return header, data

def load_halfway(loadfile):
    halfway_barcodes = set()
    halfway_readnum = defaultdict(int)
    halfway_mutations = defaultdict(int)
    start_i = 0
    if loadfile:
        with open(loadfile) as f:
            start_i = int(f.readline())
            for line in inputs(f):
                barcode, readnum, mutations = line.split()
                halfway_barcodes.add(barcode)
                halfway_readnum[barcode] = int(readnum)
                halfway_mutations[barcode] = int(mutations)
    return start_i, halfway_barcodes, halfway_readnum, halfway_mutations

def save_halfway(savefile, merged_barcodes, merged_readnum, merged_mutations, step):
    if not savefile:
        return False
    tmpfile = savefile+'.tmp'
    with open(tmpfile, 'w') as f:
        print(step, file=f)
        for barcode in sorted(merged_barcodes):
            readnum, mutations = merged_readnum[barcode], merged_mutations[barcode]
            print(barcode, readnum, mutations, sep='\t', file=f)
    shutil.move(tmpfile, savefile)
    return True

class UnionFind:
    def __init__(self, N=None):
        if N is None or N < 1:
            self.parent = defaultdict(lambda: -1)
        else:
            self.parent = [-1]*int(N)

    def root(self, n):
        stack = []
        while not (isinstance(n, int) and n < 0):
            stack.append(n)
            n = self.parent[n]
        m = stack.pop()
        while stack:
            self.parent[stack.pop()] = m
        return m

    def merge(self, m, n):
        rm = self.root(m)
        rn = self.root(n)
        if rm != rn:
            if -self.parent[rm] < -self.parent[rn] and rm != m: # if m is root, it remains root after merge.
                rm, rn = rn, rm
            self.parent[rm] += self.parent[rn]
            self.parent[rn] = rm

    def size(self, n):
        return -self.parent[self.root(n)]
    
    def setsize(self, n, s):
        rn = self.root(n)
        if rn == n:
            self.parent[n] = -s
        else:
            self.parent[self.root(n)] -= (s-1)

    def connected(self, m, n):
        return self.root(m) == self.root(n)

    def groups(self):
        if isinstance(self.parent,list):
            return list(filter(lambda i: self.parent[i]<0, range(len(self.parent))))
        else: # self.parent: defaultdict
            return list(filter(lambda i: self.parent[i]<0, self.parent.keys()))
 
    def groups_num(self):
        return len(self.groups())

    def elements(self):
        if isinstance(self.parent,list):
            return range(len(self.parent))
        else:
            return self.parent.keys()

def get_distance_dictionary(dict_filename):
    distance_dict = dict()
    with open(dict_filename) as f:
        for line in inputs(f):
            stg, mixd, d = line.split()
            d = float(d)
            distance_dict[stg] = d
    return distance_dict

def get_mixd_dictionary(dict_filename):
    mixd_dict = dict()
    with open(dict_filename) as f:
        for line in inputs(f):
            stg, mixd, d = line.split()
            # d = float(d)
            mixd_dict[stg] = mixd
    return mixd_dict

mixd_distances = dict()
def mixd_distance(mixd):
    if mixd not in mixd_distances:
        di = dx = dd = 0
        prev = ''
        for s in mixd:
            if s == 'X':
                dx += 1
            elif s == 'I':
                if prev == 'I':
                    di += 0.5
                else:
                    di += 2.5
            elif s == 'D':
                if prev =='D':
                    dd += 0.5
                else:
                    dd += 2.5
            prev = s
        d = di + dx + dd
        mixd_distances[mixd] = (d,di,dx,dd)
    return mixd_distances[mixd]

def get_distance_list(dict_filename):
    distance_list = []
    with open(dict_filename) as f:
        for line in inputs(f):
            stg, mixd, d = line.split()
            distance_list.append(float(d))
    return distance_list

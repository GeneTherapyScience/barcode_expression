from scipy import stats
import sys
n = int(sys.argv[1])
x = [float(input()) for _ in range(n)]
y = [float(line)for line in sys.stdin.readlines()]
print(stats.mannwhitneyu(x, y))

# ============================================================
# Riduzione della congettura a un enunciato combinatorio.
# Per grafi G su p = 6, 7 nodi con molti gemelli e per ogni riorientamento G'
# dentro le classi di gemelli NON raggiungibile con cricche-modulo, si cerca un
# insieme K di al piu' 5 nodi con G'[K] fuori dalla chiusura di G[K]. Se esiste
# sempre, la congettura segue da Amendola et al. (Cor. 4.9: G ~ G' implica
# G[K] ~ G'[K]) e dalla classificazione esatta per p <= 5 (11_exact_census.py).
# Solo combinatoria, nessun calcolo numerico.
# Uso:  python3 inst/simulations/identifiability/12_local_witness.py <p> <n_grafi>
# ============================================================
import sys, itertools, random

def skel(d):
    n = len(d); return [[1 if (d[i][j] or d[j][i]) else 0 for j in range(n)] for i in range(n)]

def reach(d):
    n = len(d); R = [[(u == v) or d[v][u] == 1 for v in range(n)] for u in range(n)]
    for k in range(n):
        for u in range(n):
            if R[u][k]:
                for v in range(n):
                    if R[k][v]: R[u][v] = True
    return R

def cliques(s):
    n = len(s); out = []
    def ext(Q, cand):
        for v in cand:
            Q2 = Q + [v]
            if len(Q2) >= 2: out.append(Q2)
            ext(Q2, [u for u in cand if u > v and s[v][u]])
    ext([], list(range(n)))
    return out

def orient_all(d, T):
    pr = [(a, b) for i, a in enumerate(T) for b in T[i + 1:]]
    for o in range(2 ** len(pr)):
        x = [row[:] for row in d]
        for t, (a, b) in enumerate(pr):
            x[a][b] = x[b][a] = 0
            if (o >> t) & 1: x[b][a] = 1
            else: x[a][b] = 1
        yield x

def module_moves(d):
    n = len(d); s = skel(d); R = reach(d); out = []
    for T in cliques(s):
        Ts = set(T)
        pa = [frozenset(u for u in range(n) if d[v][u] and u not in Ts) for v in T]
        ch = [frozenset(u for u in range(n) if d[u][v] and u not in Ts) for v in T]
        if len(set(pa)) > 1 or len(set(ch)) > 1: continue
        P_, C_ = pa[0], ch[0]
        if any(not d[c][u] for u in P_ for c in C_): continue
        oth = [k for k in range(n) if k not in Ts and k not in P_ and k not in C_]
        def trek(k): return any(R[z][k] and R[z][v] for v in T for z in range(n))
        if any(not all(s[k][v] for v in T) and trek(k) for k in oth): continue
        out.extend(orient_all(d, T))
    return out

def key(d): return tuple(tuple(r) for r in d)

def closure(d):
    seen = {key(d)}; queue = [d]; i = 0
    while i < len(queue):
        y = queue[i]; i += 1
        for z in module_moves(y):
            kz = key(z)
            if kz not in seen: seen.add(kz); queue.append(z)
    return seen

def twin_classes(d):
    n = len(d); s = skel(d); N = [frozenset([i] + [j for j in range(n) if s[i][j]]) for i in range(n)]
    cl = {}
    for i in range(n): cl.setdefault(N[i], []).append(i)
    return [c for c in cl.values() if len(c) >= 2]

def candidates(d):
    out = [d]
    for W in twin_classes(d): out = [x for y in out for x in orient_all(y, W)]
    return out

def induced(d, K): return [[d[j][i] for i in K] for j in K]   # stessa convenzione: d[j][i] = 1 se i -> j

def random_graph_with_twins(n, rng):
    while True:
        d = [[0] * n for _ in range(n)]
        # una o due cricche (gemelli potenziali) e archi casuali; i vicini esterni sono
        # condivisi dai membri della cricca con orientamenti casuali (anche misti)
        nodes = list(range(n)); rng.shuffle(nodes)
        t = rng.choice([2, 3, 3, 4]); T = nodes[:t]; rest = nodes[t:]
        for a, b in itertools.combinations(T, 2):
            if rng.random() < .5: d[b][a] = 1
            else: d[a][b] = 1
        O = [u for u in rest if rng.random() < .5]
        for u in O:
            for v in T:
                if rng.random() < .5: d[v][u] = 1
                else: d[u][v] = 1
        for a, b in itertools.combinations(rest, 2):
            if rng.random() < .35:
                if rng.random() < .5: d[b][a] = 1
                else: d[a][b] = 1
        if twin_classes(d) and max(len(c) for c in cliques(skel(d))) <= 5: return d

p, NG = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(1200 + p)
cache = {}
def closure_small(x):
    k = key(x)
    if k not in cache: cache[k] = closure(x)
    return cache[k]
n_cand = n_out = n_w4 = n_w5 = n_none = 0; examples = []
for g in range(NG):
    d = random_graph_with_twins(p, rng); cl = closure(d)
    for x in candidates(d):
        kx = key(x)
        if kx == key(d): continue
        n_cand += 1
        if kx in cl: continue
        n_out += 1
        De = [(a, b) for a in range(p) for b in range(a + 1, p) if d[a][b] != x[a][b] or d[b][a] != x[b][a]]
        found = None
        for size in (3, 4, 5):
            for K in itertools.combinations(range(p), size):
                if not any(a in K and b in K for a, b in De): continue
                if key(induced(x, K)) not in closure_small(induced(d, K)): found = size; break
            if found: break
        if found is None:
            n_none += 1
            if len(examples) < 5: examples.append((d, x))
        elif found <= 4: n_w4 += 1
        else: n_w5 += 1
print(f"p = {p}: {NG} grafi, {n_cand} riorientamenti dentro i gemelli, {n_out} fuori dalla chiusura")
print(f"  testimone con al piu' 4 nodi: {n_w4}; con 5 nodi: {n_w5}; nessun testimone: {n_none}")
for d, x in examples:
    lab = lambda m: " ".join(f"{i+1}->{j+1}" for i in range(p) for j in range(p) if m[j][i])
    print("  G :", lab(d)); print("  G':", lab(x))

# ============================================================
# Classificazione esatta (aritmetica razionale) delle classi di equivalenza dei
# grafi semplici con p <= 5 nodi: dimostrazione assistita dal calcolatore che la
# classe di ogni grafo e' la chiusura sotto il riorientamento delle cricche-modulo.
#
# - Necessita' gia' dimostrata (Amendola et al. 2025, Cor. 4.9, e Teorema 3 delle
#   note): due grafi equivalenti hanno lo stesso scheletro e differiscono solo su
#   archi fra gemelli (N[i] = N[j]). Candidati = riorientamenti dentro le classi
#   di gemelli.
# - Candidati fuori dalla chiusura: si costruisce B intero e stabile su G, si
#   calcola Sigma esatta e si mostra che A_{G'}(Sigma) vec(B') = -vech(I) non ha
#   soluzione (rango esatto). Poiche' A_{G'}(Sigma) ha rango pieno (Amendola et
#   al., Teorema 2.4), Sigma non sta in M_{G'}: i modelli sono diversi.
# - Candidati nella chiusura: per controllo, si verifica esattamente che il
#   sistema ha soluzione (Teorema 1 delle note).
# Uso:  python3 inst/simulations/identifiability/11_exact_census.py <p>
# ============================================================
import sys, itertools
from fractions import Fraction
import numpy as np

p = int(sys.argv[1])
PAIRS = [(r, c) for r in range(p) for c in range(r + 1, p)]
NP = len(PAIRS)

def to_mat(s):                     # s[q] in {0,1,2}: 1 = r -> c, 2 = c -> r; d[j][i] = 1 se i -> j
    d = [[0] * p for _ in range(p)]
    for q, v in enumerate(s):
        r, c = PAIRS[q]
        if v == 1: d[c][r] = 1
        elif v == 2: d[r][c] = 1
    return d

def key(d): return tuple(tuple(row) for row in d)

def skel(d): return [[1 if (d[i][j] or d[j][i]) else 0 for j in range(p)] for i in range(p)]

def reach(d):                       # R[u][v]: u raggiunge v (u = v incluso)
    R = [[(u == v) or d[v][u] == 1 for v in range(p)] for u in range(p)]
    for k in range(p):
        for u in range(p):
            if R[u][k]:
                for v in range(p):
                    if R[k][v]: R[u][v] = True
    return R

def cliques(s):
    out = []
    def ext(Q, cand):
        for v in cand:
            Q2 = Q + [v]
            if len(Q2) >= 2: out.append(Q2)
            ext(Q2, [u for u in cand if u > v and s[v][u]])
    ext([], list(range(p)))
    return out

def orient_all(d, T):              # tutti gli orientamenti semplici degli archi interni a T
    pr = [(a, b) for i, a in enumerate(T) for b in T[i + 1:]]
    for o in range(2 ** len(pr)):
        x = [row[:] for row in d]
        for t, (a, b) in enumerate(pr):
            x[a][b] = x[b][a] = 0
            if (o >> t) & 1: x[b][a] = 1
            else: x[a][b] = 1
        yield x

def module_moves(d):
    s = skel(d); R = reach(d); out = []
    for T in cliques(s):
        Ts = set(T)
        pa = [frozenset(u for u in range(p) if d[v][u] and u not in Ts) for v in T]
        ch = [frozenset(u for u in range(p) if d[u][v] and u not in Ts) for v in T]
        if len(set(pa)) > 1 or len(set(ch)) > 1: continue
        P_, C_ = pa[0], ch[0]
        if any(not d[c][u] for u in P_ for c in C_): continue
        oth = [k for k in range(p) if k not in Ts and k not in P_ and k not in C_]
        def trek(k): return any(R[z][k] and R[z][v] for v in T for z in range(p))
        if any(not all(s[k][v] for v in T) and trek(k) for k in oth): continue
        out.extend(orient_all(d, T))
    return out

def closure(d):
    seen = {key(d)}; queue = [d]; i = 0
    while i < len(queue):
        y = queue[i]; i += 1
        for z in module_moves(y):
            kz = key(z)
            if kz not in seen: seen.add(kz); queue.append(z)
    return seen

def twin_classes(d):
    s = skel(d); N = [frozenset([i] + [j for j in range(p) if s[i][j]]) for i in range(p)]
    cl = {}
    for i in range(p): cl.setdefault(N[i], []).append(i)
    return [c for c in cl.values() if len(c) >= 2]

def candidates(d):                  # riorientamenti dentro le classi di gemelli
    out = [d]
    for W in twin_classes(d):
        out = [x for y in out for x in orient_all(y, W)]
    return out

# ---- aritmetica esatta ----
def rank(M):
    M = [row[:] for row in M]; r = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

def solve(A, b):                     # sistema quadrato non singolare
    n = len(A); M = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c] / M[c][c]
                M[i] = [a - f * bb for a, bb in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]

VECH = [(k, l) for k in range(p) for l in range(k, p)]
def lyap_exact(B):                   # B Sigma + Sigma B^T = -I, incognite vech(Sigma)
    idx = {kl: t for t, kl in enumerate(VECH)}
    def sid(a, b): return idx[(min(a, b), max(a, b))]
    A = []; rhs = []
    for (k, l) in VECH:
        row = [Fraction(0)] * len(VECH)
        for m in range(p):
            row[sid(m, l)] += B[k][m]; row[sid(k, m)] += B[l][m]
        A.append(row); rhs.append(Fraction(-1 if k == l else 0))
    v = solve(A, rhs)
    return [[v[sid(a, b)] for b in range(p)] for a in range(p)]

def in_model(S, d2):                 # esiste B' su G' (+ diagonale) con B'S + SB'^T = -I ?
    cols = [(i, i) for i in range(p)] + [(j, i) for i in range(p) for j in range(p) if d2[j][i]]
    A = []; rhs = []
    for (k, l) in VECH:
        row = []
        for (a, m) in cols:        # contributo di B'[a][m]
            v = Fraction(0)
            if a == k: v += S[m][l]
            if a == l: v += S[k][m]
            row.append(v)
        A.append(row); rhs.append(Fraction(-1 if k == l else 0))
    rA = rank(A)
    assert rA == len(cols)           # rango pieno (Amendola et al., Teorema 2.4)
    return rank([A[i] + [rhs[i]] for i in range(len(A))]) == rA

VALS = [1, -2, 3, 1, -1, 2, -3, 2, 1, -2, 3, -1, 2, 1, -3, 2, -1, 3, 1, -2]
def int_B(d, shift):                 # B intero, diagonale dominante (stabile)
    B = [[Fraction(0)] * p for _ in range(p)]; t = shift
    for j in range(p):
        for i in range(p):
            if d[j][i]: B[j][i] = Fraction(VALS[t % len(VALS)]); t += 1
    for j in range(p): B[j][j] = -(sum(abs(x) for x in B[j]) + 1)
    return B

# ---- rappresentanti delle classi di isomorfismo ----
states = np.array(list(itertools.product(range(3), repeat=NP)), dtype=np.int8)
perms = list(itertools.permutations(range(p)))
pos = {pr: q for q, pr in enumerate(PAIRS)}
def perm_states(pm):                 # stato della coppia q dopo la permutazione
    src = np.zeros(NP, dtype=int); flip = np.zeros(NP, dtype=bool)
    for q, (r, c) in enumerate(PAIRS):
        a, b = pm[r], pm[c]
        if a < b: src[pos[(a, b)]] = q
        else: src[pos[(b, a)]] = q; flip[pos[(b, a)]] = True
    return src, flip
w3 = 3 ** np.arange(NP)
best = None
for pm in perms:
    src, flip = perm_states(pm)
    st = states[:, src].copy()
    st[:, flip] = np.where(st[:, flip] == 0, 0, 3 - st[:, flip])
    code = st.astype(np.int64) @ w3
    best = code if best is None else np.minimum(best, code)
_, reps = np.unique(best, return_index=True)

n_twins = n_cand = n_cert = n_retry = n_eq_ok = 0
for g in reps:
    d = to_mat(states[g])
    if not twin_classes(d): continue
    n_twins += 1
    cl = closure(d); B = int_B(d, 0); S = lyap_exact(B)
    for x in candidates(d):
        kx = key(x)
        if kx == key(d): continue
        n_cand += 1
        if kx in cl:
            assert in_model(S, x), "membro della chiusura non equivalente"
            n_eq_ok += 1
        else:
            shift = 0; Sx = S
            while in_model(Sx, x):    # B non generico per questo candidato: si cambia B
                shift += 1; n_retry += 1; assert shift < 20
                Sx = lyap_exact(int_B(d, shift))
            n_cert += 1
print(f"p = {p}: {len(reps)} classi di isomorfismo, {n_twins} con gemelli")
print(f"  candidati (riorientamenti dentro i gemelli): {n_cand}")
print(f"  nella chiusura delle cricche-modulo, equivalenza verificata esattamente: {n_eq_ok}")
print(f"  fuori dalla chiusura, non equivalenza certificata esattamente: {n_cert} (B cambiato {n_retry} volte)")

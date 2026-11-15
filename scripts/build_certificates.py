import itertools
import json
import random
from pathlib import Path


def all_cuts_masks(n):
    k = n // 2
    return [sum(1 << i for i in (0,) + r) for r in itertools.combinations(range(1, n), k - 1)]


def precompute_g(n, masks):
    M = len(masks)
    G = [[0] * M for _ in range(M)]
    for i in range(M):
        mi = masks[i]
        for j in range(i + 1, M):
            G[i][j] = G[j][i] = 4 * (mi & masks[j]).bit_count() - n
    return G


def tabu_search(M, G, m, W, steps=250, restarts=4, seed=0, cand_limit=100, tabu_tenure=9):
    best_sel, best_val = None, None
    for r in range(restarts):
        rng = random.Random(seed + r * 73)
        sel = rng.sample(range(M), m)
        in_sel = [False] * M
        for x in sel:
            in_sel[x] = True
        out_pool = [v for v in range(M) if not in_sel[v]]
        contrib = [sum(W[v][u] for u in sel) for v in range(M)]
        cur = sum(contrib[u] for u in sel) // 2
        bv, bs = cur, sel[:]
        tabu = {}
        for step in range(steps):
            bd, bp = None, None
            cand_v = out_pool if len(out_pool) <= cand_limit else rng.sample(out_pool, cand_limit)
            for u in sel:
                cu, Wu = contrib[u], W[u]
                for v in cand_v:
                    d = contrib[v] - Wu[v] - cu
                    if tabu.get(v, -1) > step and cur + d >= bv:
                        continue
                    if bd is None or d < bd:
                        bd, bp = d, (u, v)
            if bp is None:
                break
            u, v = bp
            sel.remove(u)
            sel.append(v)
            in_sel[u] = False
            in_sel[v] = True
            out_pool.remove(v)
            out_pool.append(u)
            cur += bd
            Wu, Wv = W[u], W[v]
            for w in range(M):
                contrib[w] += Wv[w] - Wu[w]
            tabu[u] = step + tabu_tenure
            tabu[v] = step + tabu_tenure
            if cur < bv:
                bv, bs = cur, sel[:]
        if best_val is None or bv < best_val:
            best_val, best_sel = bv, bs
    return sorted(best_sel)


def build_witt_and_disjoint():
    # Extended ternary Golay code G_12 -> S(5,6,12)/Z_2
    P = [[0, 1, 1, 1, 1, 1]]
    c = [0, 1, 2, 2, 1]
    for i in range(5):
        P.append([1] + [c[(j - i) % 5] for j in range(5)])
    G_mat = [[1 if i == j else 0 for j in range(6)] + P[i] for i in range(6)]
    witt1 = sorted({
        frozenset(i for i, x in enumerate(cw) if x != 0)
        for coeffs in itertools.product([0, 1, 2], repeat=6)
        for cw in [[sum(coeffs[r] * G_mat[r][col] for r in range(6)) % 3 for col in range(12)]]
        if sum(x != 0 for x in cw) == 6 and 0 in {i for i, x in enumerate(cw) if x != 0}
    }, key=sorted)
    masks1 = [sum(1 << i for i in s) for s in witt1]
    set1 = set(masks1)
    rng = random.Random(42)
    while True:
        perm = list(range(12))
        rng.shuffle(perm)
        cand = []
        for s in witt1:
            ps = frozenset(perm[i] for i in s)
            if 0 not in ps:
                ps = frozenset(range(12)) - ps
            cand.append(sum(1 << i for i in ps))
        if set(cand).isdisjoint(set1):
            masks2 = sorted(cand)
            break
    # Find one Hadamard 11-clique inside witt1 (all blocks orthogonal to each other, containing 0 and say 1? Wait: blocks intersecting in 3)
    adj0 = [[i != j and (masks1[i] & masks1[j]).bit_count() == 3 for j in range(66)] for i in range(66)]
    clique11 = None
    def bt(cand, cur):
        nonlocal clique11
        if clique11 is not None:
            return
        if len(cur) == 11:
            clique11 = [masks1[i] for i in cur]
            return
        need = 11 - len(cur)
        for idx, v in enumerate(cand):
            if len(cand) - idx < need:
                break
            nc = [u for u in cand[idx + 1:] if adj0[v][u]]
            if len(nc) >= need - 1:
                bt(nc, cur + [v])
    bt(list(range(66)), [])
    return sorted(masks1), sorted(masks2), sorted(clique11)


def main():
    n12 = 12
    masks12 = all_cuts_masks(n12)
    M12 = len(masks12)
    G12 = precompute_g(n12, masks12)
    WQ12 = [[G12[i][j] ** 2 if i != j else 0 for j in range(M12)] for i in range(M12)]
    h12_2 = {g: 2 ** ((12 - g) // 2) + 2 ** ((12 + g) // 2) for g in range(-8, 12, 4)}
    W2_12 = [[h12_2[G12[i][j]] if i != j else 0 for j in range(M12)] for i in range(M12)]
    Wno8_Q12 = [[(100000 + G12[i][j] ** 2) if abs(G12[i][j]) == 8 else (G12[i][j] ** 2 if i != j else 0) for j in range(M12)] for i in range(M12)]

    witt1, witt2, had11 = build_witt_and_disjoint()
    # m=22 optimal {0:132, 4:99} (Q=1584, c8=0)
    sel22 = None
    for r22 in range(10):
        cand22 = tabu_search(M12, G12, 22, Wno8_Q12, steps=800, restarts=1, seed=222 + r22 * 19, cand_limit=500, tabu_tenure=8)
        q22 = sum(Wno8_Q12[cand22[i]][cand22[j]] for i in range(22) for j in range(i + 1, 22))
        if q22 == 1584:
            sel22 = cand22
            break
    assert sel22 is not None
    masks22 = [masks12[i] for i in sel22]

    # Separation pairs at n=12: m in [20, 40, 50, 60, 140, 220]
    sep12 = []
    for m in [20, 40, 50, 60, 140, 220]:
        sQ = tabu_search(M12, G12, m, WQ12, steps=300, restarts=5, seed=100 + m)
        if m <= 50:
            sH = tabu_search(M12, G12, m, Wno8_Q12, steps=350, restarts=6, seed=200 + m)
        elif m == 60:
            # Take S(5,6,12)/Z_2 minus 6 mutually orthogonal cuts from had11!
            rem = set(had11[:6])
            sH_masks = sorted(set(witt1) - rem)
            sH = [masks12.index(x) for x in sH_masks]
        else:
            # Minimize d=3 score
            h12_3 = {g: 3 ** ((12 - g) // 2) + 3 ** ((12 + g) // 2) for g in range(-8, 12, 4)}
            W3_12 = [[h12_3[G12[i][j]] if i != j else 0 for j in range(M12)] for i in range(M12)]
            sH = tabu_search(M12, G12, m, W3_12, steps=250, restarts=4, seed=300 + m)
        sep12.append({
            "m": m,
            "Q_masks": [masks12[i] for i in sQ],
            "H_masks": [masks12[i] for i in sH],
        })

    # Separation pairs at n=14: m in [40, 80, 140]
    n14 = 14
    masks14 = all_cuts_masks(n14)
    M14 = len(masks14)
    G14 = precompute_g(n14, masks14)
    WQ14 = [[G14[i][j] ** 2 if i != j else 0 for j in range(M14)] for i in range(M14)]
    h14_3 = {g: 3 ** ((14 - g) // 2) + 3 ** ((14 + g) // 2) for g in range(-10, 14, 4)}
    W3_14 = [[h14_3[G14[i][j]] if i != j else 0 for j in range(M14)] for i in range(M14)]
    sep14 = []
    for m in [40, 80, 140]:
        sQ = tabu_search(M14, G14, m, WQ14, steps=200, restarts=3, seed=400 + m)
        sH = tabu_search(M14, G14, m, W3_14, steps=200, restarts=3, seed=500 + m)
        sep14.append({
            "m": m,
            "Q_masks": [masks14[i] for i in sQ],
            "H_masks": [masks14[i] for i in sH],
        })

    out = {
        "n12": {
            "witt1_masks": witt1,
            "witt2_masks": witt2,
            "hadamard11_in_witt1_masks": had11,
            "m22_masks": masks22,
            "separation_pairs": sep12,
        },
        "n14": {
            "separation_pairs": sep14,
        },
    }
    Path(r"c:\Users\Thijs\Desktop\Projecten\HaarPurityDesigns\scripts\certificates_data.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print("Wrote certificates_data.json successfully!")


if __name__ == "__main__":
    main()

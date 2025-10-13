import json
import numpy as np

def grid_edges(rows=3, cols=3):
    def idx(r, c): return r*cols + c
    E = []
    for r in range(rows):
        for c in range(cols):
            u = idx(r, c)
            if c+1 < cols: E.append([u, idx(r, c+1)])
            if r+1 < rows: E.append([u, idx(r+1, c)])
    return E

def gaussian(t, t0, sigma):
    return np.exp(-0.5*((t - t0)/sigma)**2)

# time base
dt = 0.5e-9
N  = 400                 # 200 ns total
t  = np.arange(N)*dt

rows, cols = 3, 3
nq = rows*cols
edges = grid_edges(rows, cols)

# init envelopes
I = [np.zeros(N) for _ in range(nq)]
Q = [np.zeros(N) for _ in range(nq)]
S = [np.zeros(N) for _ in range(len(edges))]

# example: two single-qubit pulses (I only), two exchange pulses
I[0] = gaussian(t, 40e-9, 6e-9)          # q0
I[8] = gaussian(t, 80e-9, 6e-9)          # q8

def edge_index(edges, a, b):
    pair = [min(a,b), max(a,b)]
    for i,e in enumerate(edges):
        if e == pair: return i
    return None

e01 = edge_index(edges, 0, 1)
e45 = edge_index(edges, 4, 5)
if e01 is not None:
    S[e01] = gaussian(t, 120e-9, 8e-9)
if e45 is not None:
    S[e45] = gaussian(t, 150e-9, 8e-9)

# clip to [0,1]
for k in range(nq):
    I[k] = np.clip(I[k], 0.0, 1.0)
    Q[k] = np.clip(Q[k], 0.0, 1.0)
for i in range(len(edges)):
    S[i] = np.clip(S[i], 0.0, 1.0)

# pack JSON
out = {
    "dt_s": dt,
    "num_samples": N,
    "qubits": [
        {"index": k,
         "I_envelope": list(map(float, I[k])),
         "Q_envelope": list(map(float, Q[k]))}
        for k in range(nq)
    ],
    "edges": [
        {"pair": edges[i],
         "s_envelope": list(map(float, S[i]))}
        for i in range(len(edges))
    ]
}

with open("pulses_envelopes_9q.json", "w") as f:
    json.dump(out, f, indent=2)

print("Wrote pulses_envelopes_9q.json")
print(f"samples: {N}, duration [ns]: {N*dt*1e9:.1f}, dt [ns]: {dt*1e9:.1f}")
print(f"qubits: {nq}, edges: {len(edges)}")


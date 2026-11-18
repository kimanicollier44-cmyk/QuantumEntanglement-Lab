#!/usr/bin/env python3
"""QuantumEntanglement-Lab - educational entanglement-measures lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * This lab computes entanglement MEASURES of small density
    matrices, exactly, from scratch: Wootters concurrence and the
    entanglement of formation built on it, negativity from the
    partial transpose, von Neumann entropy of the reduced state,
    mutual information, the Horodecki CHSH maximum from the
    correlation tensor, and the CKW monogamy balance for the GHZ
    and W three-qubit states.
  * A measure is a number about a state ON PAPER. No photon was
    entangled, no device was benchmarked, and no number here is a
    specification for any hardware, network, or protocol.
  * Correlation is not entanglement: the classically correlated
    mixture in the zoo carries a full bit of mutual information
    with concurrence exactly zero, and the lab prints both numbers
    side by side so the difference cannot be blurred.

Every state in the zoo, every sweep point, and every threshold is
computed from the matrices by the lab's own linear algebra (a
cyclic Jacobi eigensolver under a real embedding for Hermitian
spectra) - nothing is quoted. Hardware, network, cryptographic,
metrology, and quantum-advantage uses are permanent HOLD gates -
with Gemini Argon 4 Extended Thinking as DECIDER (names-only,
tighten-only), Muse as advisor (never overrides HOLD), and
Google Flow plan-only seats.
"""
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"
SQRT2 = math.sqrt(2.0)
CHSH_CLASSICAL = 2.0
CHSH_TSIRLESON = 2.0 * SQRT2

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "quantum_network_deployment": "HOLD - no deploying a quantum network or repeater chain from lab measures",
    "real_entanglement_distribution": "HOLD - no distributing real entanglement; these states exist only as matrices",
    "hardware_benchmark_claim": "HOLD - a concurrence value is not a benchmark of any physical device",
    "quantum_advantage_claim": "HOLD - no quantum-advantage or supremacy claim follows from a paper measure",
    "qkd_security_claim": "HOLD - no QKD security guarantee; lab CHSH values certify no real key",
    "teleportation_service_claim": "HOLD - no teleportation service or protocol deployment claim",
    "vendor_device_comparison": "HOLD - no comparing or ranking real vendor hardware with these numbers",
    "cryptography_strength_claim": "HOLD - no cryptographic strength claim for any real system",
    "device_certification": "HOLD - no certifying real quantum hardware, sources, or detectors",
    "sensing_metrology_claim": "HOLD - no sensing or metrology performance claim from a state on paper",
    "loophole_free_bell_claim": "HOLD - a computed CHSH maximum is not a loophole-free Bell test",
    "investment_claim": "HOLD - no investment claim; a teaching lab is not a prospectus",
    "production_system_claim": "HOLD - no production system may be built or run from this lab",
    "measure_as_physical_proof_claim": "HOLD - a measure is a number about a state on paper, never proof about a physical system",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Linear algebra
# Density matrices are lists of lists of complex numbers, written
# from scratch. Hermitian spectra come from a cyclic Jacobi
# eigensolver run on the real embedding [[Re, -Im], [Im, Re]] of
# the complex matrix: every eigenvalue of the Hermitian matrix
# appears exactly twice in the embedding's spectrum, and the
# complex eigenvectors are recovered by keeping one representative
# of each conjugate pair (greedy independence selection).

def _zeros(n, m=None):
    m = n if m is None else m
    return [[0j] * m for _ in range(n)]


def _matmul(a, b):
    n, k, m = len(a), len(b), len(b[0])
    out = _zeros(n, m)
    for i in range(n):
        for kk in range(k):
            aik = a[i][kk]
            if aik == 0:
                continue
            for j in range(m):
                out[i][j] += aik * b[kk][j]
    return out


def _dagger(a):
    return [[a[j][i].conjugate() for j in range(len(a))]
            for i in range(len(a[0]))]


def _conjugate(a):
    return [[x.conjugate() for x in row] for row in a]


def _trace(a):
    return sum(a[i][i] for i in range(len(a)))


def _kron(a, b):
    n, m, p, q = len(a), len(a[0]), len(b), len(b[0])
    out = _zeros(n * p, m * q)
    for i in range(n):
        for j in range(m):
            for u in range(p):
                for v in range(q):
                    out[i * p + u][j * q + v] = a[i][j] * b[u][v]
    return out


def _expect(rho, op):
    """Tr(rho @ op) - real for Hermitian op, as every use here is."""
    n = len(rho)
    total = 0j
    for i in range(n):
        for j in range(n):
            total += rho[i][j] * op[j][i]
    return total.real


def _jacobi_eigh(a):
    """Eigendecomposition of a real symmetric matrix (cyclic Jacobi).

    Returns (eigenvalues, eigenvectors) sorted by eigenvalue
    descending; eigenvectors are the COLUMNS of the second result:
    a @ V[:, i] == eigenvalues[i] * V[:, i]. Input is not modified.
    """
    n = len(a)
    m = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _sweep in range(100):
        off = math.sqrt(sum(m[p][q] ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < 1e-13:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = m[p][q]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = m[p][p], m[q][q]
                theta = (aqq - app) / (2.0 * apq)
                if theta >= 0:
                    t = 1.0 / (theta + math.sqrt(theta * theta + 1.0))
                else:
                    t = -1.0 / (-theta + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    mkp, mkq = m[k][p], m[k][q]
                    m[k][p] = c * mkp - s * mkq
                    m[k][q] = s * mkp + c * mkq
                for k in range(n):
                    mpk, mqk = m[p][k], m[q][k]
                    m[p][k] = c * mpk - s * mqk
                    m[q][k] = s * mpk + c * mqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    vals = [m[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: -vals[i])
    vals = [vals[i] for i in order]
    vecs = [[v[row][i] for i in order] for row in range(n)]
    return vals, vecs


def hermitian_eigh(h):
    """Full eigendecomposition of a complex Hermitian matrix.

    Returns (eigenvalues, eigenvectors) sorted descending, via the
    real embedding. Eigenvalues of the embedding pair up exactly;
    the greedy independence pass keeps one complex eigenvector per
    pair, including inside degenerate subspaces.
    """
    n = len(h)
    size = 2 * n
    r = [[0.0] * size for _ in range(size)]
    for i in range(n):
        for j in range(n):
            r[i][j] = h[i][j].real
            r[i][j + n] = -h[i][j].imag
            r[i + n][j] = h[i][j].imag
            r[i + n][j + n] = h[i][j].real
    vals, vecs = _jacobi_eigh(r)
    kept = []
    for col in range(size):
        w = [complex(vecs[row][col], vecs[row + n][col]) for row in range(n)]
        for _, kv in kept:
            ip = sum(kv[i].conjugate() * w[i] for i in range(n))
            w = [w[i] - ip * kv[i] for i in range(n)]
        norm = math.sqrt(sum(abs(x) ** 2 for x in w))
        if norm > 1e-6:
            kept.append((vals[col], [x / norm for x in w]))
        if len(kept) == n:
            break
    kept.sort(key=lambda kv: -kv[0])
    return [k[0] for k in kept], [k[1] for k in kept]


def hermitian_eigenvalues(h):
    """Eigenvalues of a complex Hermitian matrix, descending."""
    n = len(h)
    size = 2 * n
    r = [[0.0] * size for _ in range(size)]
    for i in range(n):
        for j in range(n):
            r[i][j] = h[i][j].real
            r[i][j + n] = -h[i][j].imag
            r[i + n][j] = h[i][j].imag
            r[i + n][j + n] = h[i][j].real
    vals, _ = _jacobi_eigh(r)
    return [vals[2 * i] for i in range(n)]


def _sqrt_psd(rho):
    """Matrix square root of a positive semidefinite Hermitian matrix."""
    evals, evecs = hermitian_eigh(rho)
    n = len(rho)
    out = _zeros(n)
    for lam, vec in zip(evals, evecs):
        s = math.sqrt(max(lam, 0.0))
        for i in range(n):
            for j in range(n):
                out[i][j] += s * vec[i] * vec[j].conjugate()
    return out


# ================================================================ 2. States
# Basis order for two qubits is |00>, |01>, |10>, |11>; for three
# qubits index = q0*4 + q1*2 + q2. Every state below is built from
# its amplitudes or its defining mixture - nothing is measured.

SIGMA_X = [[0j, 1 + 0j], [1 + 0j, 0j]]
SIGMA_Y = [[0j, -1j], [1j, 0j]]
SIGMA_Z = [[1 + 0j, 0j], [0j, -1 + 0j]]
PAULIS = (SIGMA_X, SIGMA_Y, SIGMA_Z)
SPIN_FLIP = _kron(SIGMA_Y, SIGMA_Y)


def ket_density(amps):
    """Density matrix |psi><psi| from a list of amplitudes."""
    return [[a * b.conjugate() for b in amps] for a in amps]


def bell_state(name):
    v = [0j] * 4
    s = 1.0 / math.sqrt(2.0)
    if name == "phi+":
        v[0], v[3] = s, s
    elif name == "phi-":
        v[0], v[3] = s, -s
    elif name == "psi+":
        v[1], v[2] = s, s
    elif name == "psi-":
        v[1], v[2] = s, -s
    else:
        raise ValueError(f"unknown Bell state: {name}")
    return ket_density(v)


def product_state():
    """|00> - the plainest separable state."""
    return ket_density([1 + 0j, 0j, 0j, 0j])


def product_plus_state():
    """|+>|+> - separable, but with coherence in every local basis."""
    return ket_density([0.5, 0.5, 0.5, 0.5])


def pure_family(theta):
    """cos(theta)|00> + sin(theta)|11> - the pure sweep family."""
    return ket_density([math.cos(theta), 0j, 0j, math.sin(theta)])


def werner_state(w):
    """w * |Phi+><Phi+| + (1 - w) * I/4 - Bell state in white noise."""
    b = bell_state("phi+")
    return [[w * b[i][j] + (1.0 - w) * (0.25 if i == j else 0.0)
             for j in range(4)] for i in range(4)]


def classical_mixture():
    """(|00><00| + |11><11|)/2 - correlated, but with no entanglement."""
    rho = _zeros(4)
    rho[0][0] = 0.5
    rho[3][3] = 0.5
    return rho


def maximally_mixed():
    return [[0.25 if i == j else 0j for j in range(4)] for i in range(4)]


def ghz_state():
    v = [0j] * 8
    v[0] = v[7] = 1.0 / math.sqrt(2.0)
    return ket_density(v)


def w_state():
    v = [0j] * 8
    for idx in (1, 2, 4):
        v[idx] = 1.0 / math.sqrt(3.0)
    return ket_density(v)


def partial_trace(rho, n_qubits, trace_out):
    """Trace out the given qubit positions (0 = most significant)."""
    keep = [q for q in range(n_qubits) if q not in trace_out]
    dim_f, dim_r = 1 << n_qubits, 1 << len(keep)
    red = _zeros(dim_r)

    def bit(idx, q):
        return (idx >> (n_qubits - 1 - q)) & 1

    for i in range(dim_f):
        for j in range(dim_f):
            if all(bit(i, q) == bit(j, q) for q in trace_out):
                ri = rj = 0
                for q in keep:
                    ri = (ri << 1) | bit(i, q)
                    rj = (rj << 1) | bit(j, q)
                red[ri][rj] += rho[i][j]
    return red


def partial_transpose_b(rho):
    """Partial transpose on the second qubit of a two-qubit state."""
    out = _zeros(4)
    for a in range(2):
        for b in range(2):
            for ap in range(2):
                for bp in range(2):
                    out[2 * a + b][2 * ap + bp] = rho[2 * a + bp][2 * ap + b]
    return out


# ================================================================ 3. Measures
# Every measure is computed from the matrix in front of it. The
# analytic values they are checked against live in run_checks(),
# never inside the measures themselves.

def concurrence(rho):
    """Wootters concurrence of a two-qubit density matrix.

    R = rho (YxY) rho* (YxY); its eigenvalues come from the
    Hermitian matrix sqrt(rho) R' sqrt(rho) with the same spectrum,
    and C = max(0, l1 - l2 - l3 - l4) over their square roots.
    """
    rho_tilde = _matmul(_matmul(SPIN_FLIP, _conjugate(rho)), SPIN_FLIP)
    sq = _sqrt_psd(rho)
    m = _matmul(_matmul(sq, rho_tilde), sq)
    evals = sorted((max(x, 0.0) for x in hermitian_eigenvalues(m)),
                   reverse=True)
    lam = [math.sqrt(x) for x in evals]
    return max(0.0, lam[0] - lam[1] - lam[2] - lam[3])


def binary_entropy(x):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def entanglement_of_formation(rho):
    """EoF = h((1 + sqrt(1 - C^2))/2), Wootters' closed form."""
    c = concurrence(rho)
    return binary_entropy((1.0 + math.sqrt(max(0.0, 1.0 - c * c))) / 2.0)


def negativity(rho):
    """N = (||rho^Gamma||_1 - 1)/2 from the partial-transpose spectrum."""
    evals = hermitian_eigenvalues(partial_transpose_b(rho))
    return (sum(abs(x) for x in evals) - 1.0) / 2.0


def von_neumann_entropy(rho):
    return -sum(x * math.log2(x)
                for x in hermitian_eigenvalues(rho) if x > 1e-12)


def mutual_information(rho):
    """I(A:B) = S(A) + S(B) - S(AB), in bits."""
    ra = partial_trace(rho, 2, (1,))
    rb = partial_trace(rho, 2, (0,))
    return (von_neumann_entropy(ra) + von_neumann_entropy(rb)
            - von_neumann_entropy(rho))


def correlation_tensor(rho):
    """T[i][j] = Tr(rho sigma_i x sigma_j) - real 3x3."""
    return [[_expect(rho, _kron(PAULIS[i], PAULIS[j]))
             for j in range(3)] for i in range(3)]


def chsh_max(rho):
    """Horodecki maximum: 2*sqrt(m1 + m2), m = two largest
    eigenvalues of T^T T. 2 is the classical bound, 2*sqrt(2) the
    Tsirelson bound."""
    t = correlation_tensor(rho)
    m = [[sum(t[k][i] * t[k][j] for k in range(3)) for j in range(3)]
         for i in range(3)]
    vals, _ = _jacobi_eigh(m)
    return 2.0 * math.sqrt(max(vals[0] + vals[1], 0.0))


def _r6(x):
    return round(float(x), 6) + 0.0


def state_measures(rho):
    """The full measure card for a two-qubit state, rounded for print."""
    ra = partial_trace(rho, 2, (1,))
    return {
        "concurrence": _r6(concurrence(rho)),
        "eof": _r6(entanglement_of_formation(rho)),
        "negativity": _r6(negativity(rho)),
        "entropy_reduced": _r6(von_neumann_entropy(ra)),
        "entropy_joint": _r6(von_neumann_entropy(rho)),
        "mutual_information": _r6(mutual_information(rho)),
        "chsh_max": _r6(chsh_max(rho)),
    }


# ================================================================ 4. Zoo + sweeps
def build_zoo():
    """The state zoo: every entry measured by the same code path."""
    return {
        "bell_phi_plus": bell_state("phi+"),
        "bell_phi_minus": bell_state("phi-"),
        "bell_psi_plus": bell_state("psi+"),
        "bell_psi_minus": bell_state("psi-"),
        "product_00": product_state(),
        "product_plus_plus": product_plus_state(),
        "pure_theta_pi_6": pure_family(math.pi / 6.0),
        "werner_w_0.50": werner_state(0.5),
        "werner_w_0.80": werner_state(0.8),
        "classical_mixture": classical_mixture(),
        "maximally_mixed": maximally_mixed(),
    }


def zoo_report():
    return {name: state_measures(rho) for name, rho in build_zoo().items()}


def pure_sweep():
    """The pure family at theta = k*pi/24, with the exact CHSH law
    CHSH_max = 2*sqrt(1 + C^2) printed next to the computed value."""
    rows = []
    for k in range(13):
        theta = k * math.pi / 24.0
        rho = pure_family(theta)
        c = concurrence(rho)
        ch = chsh_max(rho)
        rows.append({
            "theta": _r6(theta),
            "concurrence": _r6(c),
            "negativity": _r6(negativity(rho)),
            "eof": _r6(entanglement_of_formation(rho)),
            "chsh_max": _r6(ch),
            "chsh_law_2sqrt_1+C2": _r6(2.0 * math.sqrt(1.0 + c * c)),
            "abs_error": _r6(abs(ch - 2.0 * math.sqrt(1.0 + c * c))),
        })
    return rows


def chsh_relation_report():
    rows = pure_sweep()
    return {
        "law": "CHSH_max = 2*sqrt(1 + C^2) for pure two-qubit states",
        "points": len(rows),
        "worst_abs_error": _r6(max(r["abs_error"] for r in rows)),
        "classical_bound": CHSH_CLASSICAL,
        "tsirelson_bound": _r6(CHSH_TSIRLESON),
    }


def _bisect_threshold(fn, lo, hi, iterations=80):
    for _ in range(iterations):
        mid = (lo + hi) / 2.0
        if fn(mid):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2.0


def werner_report():
    """Werner family: measures vs the closed forms, and the two
    thresholds (separability at w = 1/3, CHSH violation at
    w = 1/sqrt(2)) located by bisection on the computed measures."""
    rows = []
    for w in (0.0, 0.2, 1.0 / 3.0, 0.34, 0.5, 0.7,
              1.0 / SQRT2, 0.71, 0.8, 1.0):
        rho = werner_state(w)
        rows.append({
            "w": _r6(w),
            "concurrence": _r6(concurrence(rho)),
            "concurrence_formula": _r6(max(0.0, (3.0 * w - 1.0) / 2.0)),
            "negativity": _r6(negativity(rho)),
            "negativity_formula": _r6(max(0.0, (3.0 * w - 1.0) / 4.0)),
            "chsh_max": _r6(chsh_max(rho)),
            "chsh_formula_2sqrt2_w": _r6(2.0 * SQRT2 * w),
        })
    w_sep = _bisect_threshold(
        lambda w: concurrence(werner_state(w)) > 1e-9, 0.0, 1.0)
    w_chsh = _bisect_threshold(
        lambda w: chsh_max(werner_state(w)) > CHSH_CLASSICAL, 0.5, 1.0)
    return {
        "rows": rows,
        "separability_threshold_w": _r6(w_sep),
        "separability_exact_1_over_3": _r6(1.0 / 3.0),
        "chsh_threshold_w": _r6(w_chsh),
        "chsh_exact_1_over_sqrt2": _r6(1.0 / SQRT2),
        "lesson": "between w = 1/3 and w = 1/sqrt(2) the Werner state is "
                  "entangled but violates no CHSH inequality - "
                  "entanglement and nonlocality are different resources",
    }


def _pair_concurrences(rho3):
    out = {}
    for name, out_q in (("AB", 2), ("AC", 1), ("BC", 0)):
        out[name] = concurrence(partial_trace(rho3, 3, (out_q,)))
    return out


def _bipartition_concurrence(rho3, qubit=0):
    """C(A|BC) = sqrt(2 (1 - Tr rho_A^2)) for a pure three-qubit state."""
    others = tuple(q for q in range(3) if q != qubit)
    ra = partial_trace(rho3, 3, others)
    purity = _trace(_matmul(ra, ra)).real
    return math.sqrt(max(0.0, 2.0 * (1.0 - purity)))


def monogamy_report():
    """CKW monogamy: C^2(A|BC) >= C^2(AB) + C^2(AC), computed for
    GHZ (all entanglement in the triple) and W (shared pairwise,
    saturating the inequality)."""
    report = {}
    for name, rho3 in (("ghz", ghz_state()), ("w", w_state())):
        pairs = _pair_concurrences(rho3)
        cab = _bipartition_concurrence(rho3, 0)
        report[name] = {
            "pair_concurrences": {k: _r6(v) for k, v in pairs.items()},
            "c_a_bc": _r6(cab),
            "c2_a_bc": _r6(cab * cab),
            "c2_ab_plus_c2_ac": _r6(pairs["AB"] ** 2 + pairs["AC"] ** 2),
            "ckw_slack": _r6(cab * cab - pairs["AB"] ** 2 - pairs["AC"] ** 2),
            "ckw_holds": cab * cab >= pairs["AB"] ** 2 + pairs["AC"] ** 2
                         - 1e-9,
        }
    return report


# ================================================================ 5. kNN classifier
# Features: concurrence, negativity, CHSH maximum. The training
# grid is labelled by the lab's own verdict policy (policy_label
# below) - the classifier interpolates the policy across the
# measure space, it does not replace the computed measures.

def policy_label(concurrence_value, negativity_value, chsh_value):
    """The lab's verdict policy over the measure card."""
    if concurrence_value <= 1e-6:
        return "separable"
    if concurrence_value >= 0.98 and chsh_value >= 2.80:
        return "maximally-entangled"
    if chsh_value > CHSH_CLASSICAL:
        return "nonlocal"
    return "entangled-local"


def _training_grid():
    rows = []
    for ci in range(0, 11):
        c = ci / 10.0
        for ni in range(0, 6):
            n = ni / 10.0
            for ch in (0.0, 0.5, 1.0, 1.5, 2.0, 2.25, 2.5, CHSH_TSIRLESON):
                rows.append(((c, n, ch), policy_label(c, n, ch)))
    return rows


TRAINING_GRID = _training_grid()


def _scaled(features):
    c, n, ch = features
    return (c, n / 0.5, ch / CHSH_TSIRLESON)


def classify_measures(concurrence_value, negativity_value, chsh_value, k=5):
    """Distance-weighted kNN over the policy-labelled grid,
    deterministic tie-breaking. Weighting by 1/(d + eps) keeps the
    thin maximally-entangled corner of the policy (C >= 0.98 and
    CHSH >= 2.80) from being outvoted by its nonlocal neighbours."""
    target = _scaled((concurrence_value, negativity_value, chsh_value))
    ranked = []
    for features, label in TRAINING_GRID:
        point = _scaled(features)
        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(target, point)))
        ranked.append((dist, features, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-6)
    priority = {"maximally-entangled": 0, "nonlocal": 1,
                "entangled-local": 2, "separable": 3}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


def classifier_report():
    zoo = build_zoo()
    report = {}
    for name, rho in zoo.items():
        c = concurrence(rho)
        n = negativity(rho)
        ch = chsh_max(rho)
        report[name] = {
            "concurrence": _r6(c),
            "negativity": _r6(n),
            "chsh_max": _r6(ch),
            "policy": policy_label(c, n, ch),
            "knn": classify_measures(c, n, ch),
        }
    return report


# ================================================================ 6. Advisory seats
def _trace_hash(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="quantumentanglement-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all hardware and network gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace_hash(payload),
    }


def claude_advisor(topic="quantumentanglement-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the classical mixture next to the Bell states in every report: one bit of mutual information with concurrence exactly zero is the cleanest proof that correlation is not entanglement.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["two qubits on a bench, a thread of light tying them together",
                  "the thread thinning as white noise is poured over the state",
                  "a balance scale weighing pairwise entanglement against the triple",
                  "the GHZ state as a three-way knot no pair can untie",
                  "the W state sharing its thread evenly between every pair"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    def close(a, b, tol=1e-6):
        return abs(a - b) <= tol

    # 1. eigensolver honesty - spectra and reconstruction
    y_evals = sorted(hermitian_eigenvalues(SIGMA_Y))
    chk("eigensolver_sigma_y_spectrum",
        close(y_evals[0], -1.0, 1e-9) and close(y_evals[1], 1.0, 1e-9))
    psi = bell_state("psi+")
    evals, evecs = hermitian_eigh(psi)
    chk("eigensolver_bell_spectrum",
        close(evals[0], 1.0, 1e-9)
        and all(close(x, 0.0, 1e-9) for x in evals[1:]))
    rec = _zeros(4)
    for lam, vec in zip(evals, evecs):
        for i in range(4):
            for j in range(4):
                rec[i][j] += lam * vec[i] * vec[j].conjugate()
    chk("eigensolver_reconstruction",
        all(abs(rec[i][j] - psi[i][j]) < 1e-9
            for i in range(4) for j in range(4)))
    chk("density_trace_one",
        all(close(_trace(rho).real, 1.0, 1e-9)
            for rho in build_zoo().values()))

    # 2. Bell states - maximal on every measure
    for name in ("phi+", "phi-", "psi+", "psi-"):
        rho = bell_state(name)
        chk(f"bell_{name}_concurrence_1", close(concurrence(rho), 1.0, 1e-9))
        chk(f"bell_{name}_negativity_half",
            close(negativity(rho), 0.5, 1e-9))
        chk(f"bell_{name}_eof_1",
            close(entanglement_of_formation(rho), 1.0, 1e-9))
        chk(f"bell_{name}_chsh_tsirelson",
            close(chsh_max(rho), CHSH_TSIRLESON, 1e-9))
        chk(f"bell_{name}_mutual_information_2",
            close(mutual_information(rho), 2.0, 1e-9))
        chk(f"bell_{name}_reduced_entropy_1",
            close(von_neumann_entropy(partial_trace(rho, 2, (1,))), 1.0, 1e-9))

    # 3. separable states - zero entanglement, honest correlations
    chk("product_concurrence_0", concurrence(product_state()) < 1e-9)
    chk("product_negativity_0", abs(negativity(product_state())) < 1e-9)
    chk("product_mutual_information_0",
        abs(mutual_information(product_state())) < 1e-9)
    chk("product_chsh_at_classical_bound",
        close(chsh_max(product_state()), 2.0, 1e-9))
    chk("product_plus_concurrence_0", concurrence(product_plus_state()) < 1e-9)
    chk("maximally_mixed_all_zero",
        concurrence(maximally_mixed()) < 1e-9
        and abs(negativity(maximally_mixed())) < 1e-9
        and chsh_max(maximally_mixed()) < 1e-9
        and abs(mutual_information(maximally_mixed())) < 1e-9)
    cm = classical_mixture()
    chk("classical_mixture_concurrence_0", concurrence(cm) < 1e-9)
    chk("classical_mixture_negativity_0", abs(negativity(cm)) < 1e-9)
    chk("classical_mixture_mutual_information_1",
        close(mutual_information(cm), 1.0, 1e-9))
    chk("classical_mixture_chsh_no_violation",
        close(chsh_max(cm), 2.0, 1e-9))

    # 4. pure family - closed forms across the sweep
    worst = 0.0
    for k in range(13):
        theta = k * math.pi / 24.0
        rho = pure_family(theta)
        worst = max(worst,
                    abs(chsh_max(rho)
                        - 2.0 * math.sqrt(1.0 + concurrence(rho) ** 2)))
    chk("pure_family_chsh_law_worst_error", worst < 1e-6)
    for theta in (math.pi / 12, math.pi / 6, math.pi / 4, math.pi / 3):
        rho = pure_family(theta)
        chk(f"pure_family_concurrence_sin2theta_{theta:.3f}",
            close(concurrence(rho), abs(math.sin(2 * theta)), 1e-9))
        chk(f"pure_family_negativity_half_concurrence_{theta:.3f}",
            close(negativity(rho), concurrence(rho) / 2.0, 1e-9))
        chk(f"pure_family_eof_equals_reduced_entropy_{theta:.3f}",
            close(entanglement_of_formation(rho),
                  von_neumann_entropy(partial_trace(rho, 2, (1,))), 1e-9))
    chk("pure_theta_pi_6_values",
        close(concurrence(pure_family(math.pi / 6)), 0.866025, 1e-6)
        and close(chsh_max(pure_family(math.pi / 6)), 2.645751, 1e-6))

    # 5. Werner family - formulas and both thresholds
    for w in (0.2, 0.34, 0.5, 0.7, 0.8, 1.0):
        rho = werner_state(w)
        chk(f"werner_concurrence_formula_{w}",
            close(concurrence(rho), max(0.0, (3 * w - 1) / 2), 1e-9))
        chk(f"werner_negativity_formula_{w}",
            close(negativity(rho), max(0.0, (3 * w - 1) / 4), 1e-9))
        chk(f"werner_chsh_formula_{w}",
            close(chsh_max(rho), 2 * SQRT2 * w, 1e-9))
    wrep = werner_report()
    chk("werner_separability_threshold",
        abs(wrep["separability_threshold_w"] - 1.0 / 3.0) < 1e-6)
    chk("werner_chsh_threshold",
        abs(wrep["chsh_threshold_w"] - 1.0 / SQRT2) < 1e-6)
    chk("werner_below_threshold_separable",
        concurrence(werner_state(0.30)) < 1e-9
        and concurrence(werner_state(0.40)) > 0.01)
    chk("werner_entangled_but_local_band",
        policy_label(concurrence(werner_state(0.5)),
                     negativity(werner_state(0.5)),
                     chsh_max(werner_state(0.5))) == "entangled-local")

    # 6. monogamy - GHZ vs W (unrounded values, straight from the states)
    mono = monogamy_report()
    ghz_pairs = _pair_concurrences(ghz_state())
    w_pairs = _pair_concurrences(w_state())
    w_cab = _bipartition_concurrence(w_state(), 0)
    chk("ghz_pairwise_concurrence_zero",
        all(v < 1e-9 for v in ghz_pairs.values()))
    chk("ghz_bipartition_concurrence_1",
        close(_bipartition_concurrence(ghz_state(), 0), 1.0, 1e-9))
    chk("w_pairwise_concurrence_two_thirds",
        all(close(v, 2.0 / 3.0, 1e-9) for v in w_pairs.values()))
    chk("w_bipartition_concurrence_squared_8_9",
        close(w_cab * w_cab, 8.0 / 9.0, 1e-9))
    chk("w_ckw_saturated",
        abs(w_cab * w_cab - w_pairs["AB"] ** 2 - w_pairs["AC"] ** 2) < 1e-9)
    chk("ckw_holds_both",
        mono["ghz"]["ckw_holds"] and mono["w"]["ckw_holds"])

    # 7. classifier - kNN interpolates the lab's own policy
    rep = classifier_report()
    chk("classifier_bell_maximally_entangled",
        rep["bell_phi_plus"]["knn"] == "maximally-entangled")
    chk("classifier_product_separable",
        rep["product_00"]["knn"] == "separable")
    chk("classifier_classical_mixture_separable",
        rep["classical_mixture"]["knn"] == "separable")
    chk("classifier_werner_050_entangled_local",
        rep["werner_w_0.50"]["knn"] == "entangled-local")
    chk("classifier_werner_080_nonlocal",
        rep["werner_w_0.80"]["knn"] == "nonlocal")
    chk("classifier_pure_pi_6_nonlocal",
        rep["pure_theta_pi_6"]["knn"] == "nonlocal")
    chk("classifier_matches_policy_on_zoo",
        all(row["knn"] == row["policy"] for row in rep.values()))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_network_hold", gate("quantum_network_deployment") == "HOLD")
    chk("gate_benchmark_hold", gate("hardware_benchmark_claim") == "HOLD")
    chk("gate_advantage_hold", gate("quantum_advantage_claim") == "HOLD")
    chk("gate_bell_test_hold", gate("loophole_free_bell_claim") == "HOLD")
    chk("gate_measure_proof_hold",
        gate("measure_as_physical_proof_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("QuantumEntanglement-Lab - entanglement measures lab")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - what a measure is and is not")
    print("    This lab computes entanglement measures of small density")
    print("    matrices, exactly, from scratch: concurrence, entanglement")
    print("    of formation, negativity, reduced entropy, mutual")
    print("    information, the Horodecki CHSH maximum, and CKW monogamy.")
    print("    A measure is a number about a state ON PAPER. No photon was")
    print("    entangled and no device was benchmarked to produce it.")
    zoo = zoo_report()
    print("\n[2] The state zoo - every state, every measure, one code path")
    print(f"    {'state':20s} {'C':>8s} {'EoF':>8s} {'N':>8s} "
          f"{'S_red':>8s} {'MI':>8s} {'CHSH':>8s}")
    for name, ms in zoo.items():
        print(f"    {name:20s} {ms['concurrence']:8.6f} {ms['eof']:8.6f} "
              f"{ms['negativity']:8.6f} {ms['entropy_reduced']:8.6f} "
              f"{ms['mutual_information']:8.6f} {ms['chsh_max']:8.6f}")
    print("    -> the four Bell states max every measure. The classical")
    print("       mixture carries MI = 1.000000 with C = 0.000000:")
    print("       correlation is not entanglement, in print.")
    print("\n[3] Pure family cos(t)|00> + sin(t)|11> - the CHSH law")
    rel = chsh_relation_report()
    for row in pure_sweep()[::3]:
        print(f"    theta {row['theta']:.4f} | C {row['concurrence']:.6f} | "
              f"CHSH {row['chsh_max']:.6f} | law {row['chsh_law_2sqrt_1+C2']:.6f}")
    print(f"    -> CHSH_max = 2*sqrt(1 + C^2) holds across all "
          f"{rel['points']} sweep points; worst error "
          f"{rel['worst_abs_error']:.6f}.")
    print("\n[4] Werner states - two thresholds, one honest band")
    wrep = werner_report()
    for row in wrep["rows"]:
        print(f"    w {row['w']:.6f} | C {row['concurrence']:.6f} | "
              f"N {row['negativity']:.6f} | CHSH {row['chsh_max']:.6f}")
    print(f"    -> separable up to w = {wrep['separability_threshold_w']:.6f} "
          f"(exact 1/3 = {wrep['separability_exact_1_over_3']:.6f}); CHSH "
          f"violated only past w = {wrep['chsh_threshold_w']:.6f} "
          f"(exact 1/sqrt(2) = {wrep['chsh_exact_1_over_sqrt2']:.6f}).")
    print("       Between them: entangled, but local. Different resources.")
    print("\n[5] Monogamy - entanglement cannot be shared freely (CKW)")
    mono = monogamy_report()
    for name in ("ghz", "w"):
        row = mono[name]
        pairs = row["pair_concurrences"]
        print(f"    {name.upper():3s} pairs AB {pairs['AB']:.6f} "
              f"AC {pairs['AC']:.6f} BC {pairs['BC']:.6f} | "
              f"C(A|BC) {row['c_a_bc']:.6f} | slack {row['ckw_slack']:.6f}")
    print("    -> GHZ: no pair holds anything; the triple holds it all.")
    print("       W: every pair holds 2/3 and the CKW budget is spent")
    print("       exactly - the inequality saturates.")
    print("\n[6] State classifier (TinyML kNN, labels from the lab's policy)")
    for name, row in classifier_report().items():
        print(f"    {name:20s} -> {row['knn']}")
    print("\n[7] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[8] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "zoo":
        print(json.dumps(zoo_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "pure":
        print(json.dumps(pure_sweep(), indent=2, sort_keys=True))
        return 0
    if cmd == "werner":
        print(json.dumps(werner_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "chsh":
        print(json.dumps(chsh_relation_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "monogamy":
        print(json.dumps(monogamy_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "quantum_advantage_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())

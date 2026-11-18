"""Tests for QuantumEntanglement-Lab: eigensolver, measures, zoo,
sweeps, thresholds, monogamy, classifier, gates, determinism."""
import math

import main as m

TOL = 1e-9
SQRT2 = math.sqrt(2.0)


# ---- eigensolver
def test_sigma_y_spectrum():
    vals = sorted(m.hermitian_eigenvalues(m.SIGMA_Y))
    assert abs(vals[0] - (-1.0)) < TOL
    assert abs(vals[1] - 1.0) < TOL


def test_bell_spectrum_and_reconstruction():
    rho = m.bell_state("psi+")
    evals, evecs = m.hermitian_eigh(rho)
    assert abs(evals[0] - 1.0) < TOL
    assert all(abs(x) < TOL for x in evals[1:])
    for i in range(4):
        for j in range(4):
            rec = sum(lam * vec[i] * vec[j].conjugate()
                      for lam, vec in zip(evals, evecs))
            assert abs(rec - rho[i][j]) < TOL


def test_zoo_density_matrices_trace_one():
    for rho in m.build_zoo().values():
        assert abs(m._trace(rho).real - 1.0) < TOL


# ---- Bell states
def test_bell_states_maximal():
    for name in ("phi+", "phi-", "psi+", "psi-"):
        rho = m.bell_state(name)
        assert abs(m.concurrence(rho) - 1.0) < TOL
        assert abs(m.negativity(rho) - 0.5) < TOL
        assert abs(m.entanglement_of_formation(rho) - 1.0) < TOL
        assert abs(m.chsh_max(rho) - 2 * SQRT2) < TOL
        assert abs(m.mutual_information(rho) - 2.0) < TOL


def test_bell_reduced_entropy_one():
    rho = m.bell_state("phi-")
    ra = m.partial_trace(rho, 2, (1,))
    assert abs(m.von_neumann_entropy(ra) - 1.0) < TOL


# ---- separable states
def test_product_state_zero_entanglement():
    rho = m.product_state()
    assert m.concurrence(rho) < TOL
    assert abs(m.negativity(rho)) < TOL
    assert abs(m.mutual_information(rho)) < TOL
    assert abs(m.chsh_max(rho) - 2.0) < TOL


def test_product_plus_state_separable():
    assert m.concurrence(m.product_plus_state()) < TOL


def test_maximally_mixed_all_zero():
    rho = m.maximally_mixed()
    assert m.concurrence(rho) < TOL
    assert abs(m.negativity(rho)) < TOL
    assert m.chsh_max(rho) < TOL
    assert abs(m.mutual_information(rho)) < TOL


def test_classical_mixture_correlation_is_not_entanglement():
    rho = m.classical_mixture()
    assert m.concurrence(rho) < TOL
    assert abs(m.negativity(rho)) < TOL
    assert abs(m.mutual_information(rho) - 1.0) < TOL
    assert abs(m.chsh_max(rho) - 2.0) < TOL


# ---- pure family
def test_pure_family_concurrence_sin_2theta():
    for theta in (math.pi / 12, math.pi / 6, math.pi / 4, math.pi / 3):
        rho = m.pure_family(theta)
        assert abs(m.concurrence(rho) - abs(math.sin(2 * theta))) < TOL


def test_pure_family_negativity_is_half_concurrence():
    for theta in (math.pi / 12, math.pi / 6, math.pi / 3):
        rho = m.pure_family(theta)
        assert abs(m.negativity(rho) - m.concurrence(rho) / 2) < TOL


def test_pure_family_eof_equals_reduced_entropy():
    for theta in (math.pi / 12, math.pi / 6, math.pi / 3):
        rho = m.pure_family(theta)
        s_red = m.von_neumann_entropy(m.partial_trace(rho, 2, (1,)))
        assert abs(m.entanglement_of_formation(rho) - s_red) < TOL


def test_pure_family_chsh_law():
    for k in range(13):
        rho = m.pure_family(k * math.pi / 24)
        law = 2.0 * math.sqrt(1.0 + m.concurrence(rho) ** 2)
        assert abs(m.chsh_max(rho) - law) < 1e-6


def test_pure_theta_pi_6_headline_values():
    rho = m.pure_family(math.pi / 6)
    assert abs(m.concurrence(rho) - 0.866025) < 1e-6
    assert abs(m.chsh_max(rho) - 2.645751) < 1e-6
    assert abs(m.entanglement_of_formation(rho) - 0.811278) < 1e-6


# ---- Werner family
def test_werner_formulas():
    for w in (0.2, 0.34, 0.5, 0.7, 0.8, 1.0):
        rho = m.werner_state(w)
        assert abs(m.concurrence(rho) - max(0.0, (3 * w - 1) / 2)) < TOL
        assert abs(m.negativity(rho) - max(0.0, (3 * w - 1) / 4)) < TOL
        assert abs(m.chsh_max(rho) - 2 * SQRT2 * w) < TOL


def test_werner_thresholds_located():
    rep = m.werner_report()
    assert abs(rep["separability_threshold_w"] - 1 / 3) < 1e-6
    assert abs(rep["chsh_threshold_w"] - 1 / SQRT2) < 1e-6


def test_werner_entangled_but_local_band():
    rho = m.werner_state(0.5)
    label = m.policy_label(m.concurrence(rho), m.negativity(rho),
                           m.chsh_max(rho))
    assert label == "entangled-local"


# ---- monogamy
def test_ghz_monogamy():
    pairs = m._pair_concurrences(m.ghz_state())
    assert all(v < TOL for v in pairs.values())
    assert abs(m._bipartition_concurrence(m.ghz_state(), 0) - 1.0) < TOL
    assert m.monogamy_report()["ghz"]["ckw_holds"]


def test_w_monogamy_saturates_ckw():
    pairs = m._pair_concurrences(m.w_state())
    assert all(abs(v - 2 / 3) < TOL for v in pairs.values())
    cab = m._bipartition_concurrence(m.w_state(), 0)
    assert abs(cab * cab - 8 / 9) < TOL
    assert abs(cab * cab - pairs["AB"] ** 2 - pairs["AC"] ** 2) < TOL
    assert m.monogamy_report()["w"]["ckw_holds"]


# ---- classifier
def test_classifier_zoo_verdicts():
    rep = m.classifier_report()
    assert rep["bell_phi_plus"]["knn"] == "maximally-entangled"
    assert rep["product_00"]["knn"] == "separable"
    assert rep["classical_mixture"]["knn"] == "separable"
    assert rep["werner_w_0.50"]["knn"] == "entangled-local"
    assert rep["werner_w_0.80"]["knn"] == "nonlocal"
    assert rep["pure_theta_pi_6"]["knn"] == "nonlocal"


def test_classifier_matches_policy_on_zoo():
    rep = m.classifier_report()
    assert all(row["knn"] == row["policy"] for row in rep.values())


# ---- gates + seats
def test_hold_gates_14():
    assert len(m.HOLD_GATES) == 14
    assert all(v.startswith("HOLD") for v in m.HOLD_GATES.values())


def test_specific_gates_hold():
    assert m.gate("quantum_network_deployment") == "HOLD"
    assert m.gate("hardware_benchmark_claim") == "HOLD"
    assert m.gate("quantum_advantage_claim") == "HOLD"
    assert m.gate("loophole_free_bell_claim") == "HOLD"
    assert m.gate("measure_as_physical_proof_claim") == "HOLD"
    assert m.gate("run_simulation") == "PROCEED"


def test_seats_never_override():
    brains = m.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["seat"] == "Gemini Argon 4 Extended Thinking"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_hold_demo_exit_code():
    assert m.main(["hold-demo"]) == 2
    assert m.main(["hold-demo", "qkd_security_claim"]) == 2
    assert m.main(["hold-demo", "run_simulation"]) == 0


# ---- checks + determinism
def test_run_checks_all_pass():
    checks = m.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_zoo_report_deterministic():
    assert m.zoo_report() == m.zoo_report()


def test_demo_runs():
    assert m.main(["demo"]) == 0

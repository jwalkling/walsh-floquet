"""
Walsh–Fourier Floquet utilities
================================

Core library for constructing and analysing periodically driven (Floquet)
spin systems using both Walsh and Fourier bases.

This module provides:

- Walsh basis utilities (sequency ordering, product tables, generators).
- Extended-space Quasienergy operators in Walsh and Fourier bases.
- Drive constructors (Walsh / Fourier delta-kick and square-wave drives).
- Static many-body Ising Hamiltonians.
- Numerically exact Floquet unitaries (no truncation) for simple drives.
- Time-resolved evolution for the kicked drive.
- Photon participation entropy (localisation) and Walsh vs Fourier truncation error calculations.
"""

import numpy as np
import matplotlib.pyplot as plt  # kept for possible downstream use
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap

from scipy.linalg import hadamard
from scipy.linalg import logm, expm
import scipy.linalg as sp

from scipy.integrate import solve_ivp, quad
from scipy.stats import linregress
from scipy.sparse import kron, identity, csr_matrix  # for many-body Hamiltonians


# --------------------------------------------------------------------------- #
# Foundational Mathematical Operations
# --------------------------------------------------------------------------- #

# Pauli spin matrices defined globally
PAULI_X = np.array([[0, 1],
                    [1, 0]])

PAULI_Y = np.array([[0, -1j],
                    [1j, 0]])

PAULI_Z = np.array([[1, 0],
                    [0, -1]])


def safe_log10(arr):
    """
    Safely compute base-10 logarithm, ignoring non-positive entries.

    Parameters
    ----------
    arr : array_like
        Input array.

    Returns
    -------
    logvals : ndarray
        `log10(arr[mask])` for entries where `arr > 0`.
    mask : ndarray of bool
        Boolean mask with True where `arr > 0`.

    Notes
    -----
    This is useful for plotting log-scaled data that may contain zeros or
    negative values due to numerical noise.
    """
    arr = np.asarray(arr, dtype=float)
    mask = arr > 0
    return np.log10(arr[mask]), mask


# --------------------------------------------------------------------------- #
# Walsh basis functions and related operations
# --------------------------------------------------------------------------- #

def walsh_generator(N):
    """
    Construct the generator of cyclic translations in the Walsh basis.

    Parameters
    ----------
    N : int
        Dimension of the Walsh/Hadamard system (must be a power of two).

    Returns
    -------
    G : ndarray (N, N), complex
        Generator matrix `G` such that the translation operator
        satisfies `T = exp(2π G / N)` in the Walsh basis.

    Notes
    -----
    The construction proceeds by comparing a shifted Hadamard matrix
    to the original to infer the translation operator.
    """
    walsh = sp.hadamard(N)
    translated = np.roll(walsh, 1, axis=1)

    # Translation operator in Walsh basis
    translation = translated @ walsh.T / N

    # Generator G such that translation = exp(2π G / N)
    return (N / (2 * np.pi)) * sp.logm(translation)


def nat_to_seq_list(size):
    """
    Compute the sequency (number of sign changes) for each Walsh function.

    Parameters
    ----------
    size : int
        Size of the Hadamard matrix (power of two).

    Returns
    -------
    ordering : ndarray of int, shape (size,)
        `ordering[i]` gives the number of sign changes in row `i` of
        the Hadamard matrix (i.e. the sequency of Walsh function W_i).
    """
    H = sp.hadamard(size)
    diffs = np.diff(H, axis=1)          # 0 or ±2 where sign changes
    ordering = (np.abs(diffs) // 2).sum(axis=1).astype(int)
    return ordering


def seq_to_nat(M):
    """
    Convert a matrix from sequency-ordered Walsh basis to natural-ordered Walsh basis.

    Parameters
    ----------
    M : ndarray, shape (N, N)
        Matrix expressed in sequency-ordered Walsh basis; N must be 2^n.

    Returns
    -------
    M_nat : ndarray, shape (N, N)
        Matrix expressed in the natural ordering of Walsh/Hadamard rows.

    Notes
    -----
    If `v_seq` are coefficients in sequency order and `T` is the permutation
    matrix returned implicitly here, then

        v_nat = T @ v_seq

    and this function effectively computes `T @ M @ T^T`.
    """
    size = M.shape[0]
    ordering = nat_to_seq_list(size)
    U = np.zeros((size, size), dtype=int)
    U[np.arange(size), ordering] = 1
    return U @ M @ U.T


def walsh_products(N):
    """
    Construct the full Walsh product table W_i * W_j = W_k.

    Parameters
    ----------
    N : int
        Number of Walsh functions (must be a power of two).

    Returns
    -------
    M : ndarray of int, shape (N, N)
        Product table such that `M[i, j] = k` where
        W_i(t) * W_j(t) = W_k(t) for all t.

    Notes
    -----
    Implemented via a vectorised comparison of all row products against
    the Hadamard matrix.
    """
    H = sp.hadamard(N)                          # (N, N)

    # products[i, j, :] = H[i, :] * H[j, :]
    products = H[:, None, :] * H[None, :, :]    # (N, N, N)

    # eq[i, j, k, l] tests equality of product(i,j) and row k at index l
    eq = products[:, :, None, :] == H[None, None, :, :]   # (N, N, N, N)

    # matches[i, j, k] True if product(i,j) == H[k,:] at all positions
    matches = np.all(eq, axis=-1)                          # (N, N, N)

    # M[i, j] = k where the match is True
    M = np.argmax(matches, axis=2)                        # (N, N)

    return M


def compute_seq_order_walsh(N=16):
    """
    Generate sequency-ordered Walsh functions.

    Parameters
    ----------
    N : int, optional
        Size of the Hadamard matrix, by default 16.

    Returns
    -------
    walsh_rows : list of (int, list)
        Each entry is `(sequency, row_values)` where `row_values` is the
        Walsh function sampled on the Hadamard grid.
    """
    H = sp.hadamard(N)
    walsh_rows = []
    for row in H:
        sign_changes = np.sum(np.diff(np.sign(row)) != 0)
        walsh_rows.append((sign_changes, row.tolist()))
    walsh_rows.sort(key=lambda x: x[0])
    return walsh_rows


def compute_walsh_modes(func, N):
    """
    Compute Walsh coefficients of a function on [0, 1).

    Parameters
    ----------
    func : callable
        Function f(x) defined on [0, 1). Must accept a NumPy array.
    N : int
        Number of Walsh samples; must be power of two for consistency.

    Returns
    -------
    walsh_coeffs : ndarray, shape (N,)
        Walsh coefficients in natural Walsh ordering.
    """
    xvals = np.linspace(0, 1 - 1 / N, N)
    fvals = func(xvals)
    walsh_coeffs = hadamard(N) @ fvals / N
    return walsh_coeffs


def compute_exps(N):
    """
    Construct an unnormalised discrete Fourier transform matrix.

    Parameters
    ----------
    N : int
        Matrix size.

    Returns
    -------
    E : ndarray, shape (N, N), complex
        Matrix with entries exp(i * 2π * a * (N/2 - 1 - b) / N),
        where a is the row index and b the column index.
    """
    a, b = np.ogrid[:N, :N]
    return np.exp(a * 1j * 2 * np.pi / N * (N // 2 - 1 - b))


def walsh_to_fourier(M):
    """
    Convert a matrix from Walsh basis to (discrete) Fourier basis.

    Parameters
    ----------
    M : ndarray, shape (N, N)
        Operator expressed in the Walsh basis.

    Returns
    -------
    M_ft : ndarray, shape (N, N)
        Operator expressed in discrete Fourier basis.

    """
    N = M.shape[0]
    U = hadamard(N) @ compute_exps(N) / N
    return U.conj().T @ M @ U


# --------------------------------------------------------------------------- #
# Quasienergy operators in extended space
# --------------------------------------------------------------------------- #

def Q_walsh(H0, H1, HDrivePattern, omega, n):
    """
    Construct the extended-space quasienergy operator in the Walsh basis.

    Parameters
    ----------
    H0 : ndarray, shape (d, d)
        Time-independent part of the Hamiltonian.
    H1 : ndarray, shape (d, d)
        Driven part of the Hamiltonian, scaled by the drive pattern.
    HDrivePattern : ndarray, shape (N, N)
        Walsh-basis representation of the drive coupling between
        different time/Floquet sectors.
    omega : float
        Drive frequency.
    n : int
        Walsh depth such that N = 2^n is the number of time sectors.

    Returns
    -------
    Q : ndarray, shape (N*d, N*d)
        Extended-space Hamiltonian acting on time ⊗ physical Hilbert space.
    """
    Nwidth = 2 ** n
    hilbert_d = H0.shape[0]

    H0extended = np.kron(np.eye(Nwidth), H0)
    H1extended = np.kron(HDrivePattern, H1)
    freq_extended = -1j * omega * np.kron(walsh_generator(Nwidth),
                                           np.eye(hilbert_d))

    Hextended = H0extended + H1extended + freq_extended
    return Hextended #Return extended space Hamiltonian 


def Q_fourier(H0, H1, HDrivePattern, omega, MaxMode):
    """
    Construct the extended-space quasienergy operator in the Fourier basis.

    Parameters
    ----------
    H0 : ndarray, shape (d, d)
        Time-independent part of the Hamiltonian.
    H1 : ndarray, shape (d, d)
        Driven part of the Hamiltonian, scaled by the Fourier drive pattern.
    HDrivePattern : ndarray, shape (N, N)
        Fourier-basis representation of the drive, where N = 2*MaxMode + 1.
    omega : float
        Drive frequency.
    MaxMode : int
        Maximum Fourier harmonic index retained in the truncation.

    Returns
    -------
    Q : ndarray, shape (N*d, N*d)
        Extended-space Hamiltonian in Fourier basis.
    """
    Nwidth = 2 * MaxMode + 1
    hilbert_d = H0.shape[0]

    H0extended = np.kron(np.eye(Nwidth), H0)
    Hoffextended = np.kron(HDrivePattern, H1)

    exp_deriv_matrix = np.diag(np.arange(-MaxMode, MaxMode + 1, 1))
    freq_extended = -omega * np.kron(exp_deriv_matrix, np.eye(hilbert_d))

    Hextended = H0extended + Hoffextended + freq_extended
    return Hextended #Return extended space Hamiltonian


# --------------------------------------------------------------------------- #
# Drive functions (Walsh and Fourier bases)
# --------------------------------------------------------------------------- #

def compute_kick_walsh(n):
    """
    Construct Walsh-basis representation of an "up-down" delta kick drive.
    Positive pulse at t = 0, negative pulse at t = T/2.

    Parameters
    ----------
    n : int
        Walsh depth, N = 2^n.

    Returns
    -------
    drive_W : ndarray, shape (N, N)
        Drive coupling matrix in Walsh basis for an up–down kick protocol.
    """
    N = 2 ** n
    matrix = walsh_products(N)

    # Keep only entries that are 2 * odd => exactly one sign flip structure
    mask = ~((matrix % 2 == 0) & (matrix % 4 != 0))
    matrix[mask] = 0
    matrix[matrix != 0] = 1

    return seq_to_nat(matrix)


def compute_sqwave_walsh(n):
    """
    Construct Walsh-basis representation of an odd square-wave drive.
    Same as the Walsh function with sequency=1.

    Parameters
    ----------
    n : int
        Walsh depth, N = 2^n.

    Returns
    -------
    drive_W : ndarray, shape (N, N)
        Matrix coupling photon sectors for a square-wave drive.
    """
    half = 2 ** (n - 1)
    return np.diag(np.ones(half), half) + np.diag(np.ones(half), -half)


def compute_genericdrive_walsh(func, n):
    """
    Construct a Walsh-basis drive from an arbitrary scalar function.

    Parameters
    ----------
    func : callable
        Function f(x) defined on [0, 1). Samples determine Walsh coefficients.
    n : int
        Walsh depth, N = 2^n.

    Returns
    -------
    drive_W : ndarray, shape (N, N)
        Drive matrix constructed from Walsh coefficients of `func`.
    """
    N = 2 ** n
    coeffs = compute_walsh_modes(func, N)
    return coeffs[walsh_products(N)]


def compute_kick_fourier(N):
    """
    Construct Fourier-basis representation of an "up-down" delta kick drive.
    Positive pulse at t = 0, negative pulse at t = T/2.
    Parameters
    ----------
    N : int
        Number of Fourier modes.

    Returns
    -------
    drive_F : ndarray, shape (N, N)
        Binary matrix with ones on odd off-diagonals and zeros elsewhere.
    """
    matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if abs(i - j) % 2 == 1:
                matrix[i, j] = 1
    return matrix


def compute_sqwave_fourier(N):
    """
    Construct Fourier-basis representation of an odd square-wave drive.
    Same as the Walsh function with sequency=1.

    Parameters
    ----------
    N : int
        Number of Fourier modes.

    Returns
    -------
    drive_F : ndarray, shape (N, N), complex
        Hermitian matrix encoding the square wave via truncated Fourier series.
    """
    Hfourier = np.zeros(N, dtype=complex)
    for i in range(N):
        if i % 2 != 0:
            Hfourier[i] = 2 / (np.pi * i * 1j)

    Hoffdiag = np.zeros((N, N), dtype=complex)
    for i in np.arange(1, N):
        if i % 2 != 0:
            Hoffdiag += (2 / (np.pi * i * 1j)) * np.diag(np.ones(N - i), i)
    Hoffdiag = Hoffdiag + Hoffdiag.T.conj()
    return Hoffdiag


# --------------------------------------------------------------------------- #
# Static many-body Hamiltonians (Ising model)
# --------------------------------------------------------------------------- #

def H_ising_zz_zfield(J, hz, L):
    """
    Construct a nearest-neighbour Ising ZZ chain with longitudinal Z field.

    Parameters
    ----------
    J : float
        ZZ coupling strength.
    hz : float
        Longitudinal field strength along Z.
    L : int
        Number of spins.

    Returns
    -------
    H : ndarray, shape (2^L, 2^L)
        Many-body Hamiltonian matrix (dense).
    """
    H = csr_matrix((2 ** L, 2 ** L), dtype=complex)
    identity2 = np.eye(2)

    # Nearest-neighbour ZZ terms
    for i in range(L - 1):
        term = 1
        for j in range(L):
            if j == i:
                term = kron(term, PAULI_Z, format="csr")
            elif j == i + 1:
                term = kron(term, PAULI_Z, format="csr")
            else:
                term = kron(term, identity2, format="csr")
        H += -J * term

    # Longitudinal field terms
    for i in range(L):
        term = 1
        for j in range(L):
            if j == i:
                term = kron(term, PAULI_Z, format="csr")
            else:
                term = kron(term, identity2, format="csr")
        H += -hz * term

    return H.toarray()


def H_ising_xfield(hx, L):
    """
    Construct a transverse-field X Hamiltonian on L spins.

    Parameters
    ----------
    hx : float
        Transverse field strength along X.
    L : int
        Number of spins.

    Returns
    -------
    H : ndarray, shape (2^L, 2^L)
        Many-body Hamiltonian matrix (dense).
    """
    H = csr_matrix((2 ** L, 2 ** L), dtype=complex)
    identity2 = np.eye(2)

    for i in range(L):
        term = 1
        for j in range(L):
            if j == i:
                term = kron(term, PAULI_X, format="csr")
            else:
                term = kron(term, identity2, format="csr")
        H += -hx * term

    return H.toarray()


# --------------------------------------------------------------------------- #
# Explicit Floquet unitaries for simple drives
# --------------------------------------------------------------------------- #
'''
In contrast to the extended-space methods above, these functions compute the Floquet unitary
explicitly without truncation, by directly exponentiating the time-evolution operator over one period.
These are the "exact" numerical solutions for Floquet unitaries for simple drives.
'''

def UF_sqwave(omega, H0, H1):
    """
    Compute Floquet quasienergy phases for a square-wave driven system.

    Parameters
    ----------
    omega : float
        Drive frequency.
    H0 : ndarray, shape (d, d)
        Static Hamiltonian.
    H1 : ndarray, shape (d, d)
        Drive term (switched ±H1 over half periods).

    Returns
    -------
    phases : ndarray, shape (d,)
        Quasienergy phases in [-π, π) (up to branch choices).

    Notes
    -----
    The square-wave protocol is:
        H(t) = H0 + H1  for 0 <= t < T/2
             = H0 - H1  for T/2 <= t < T
    """
    T = 2 * np.pi / omega
    Ha = H0 + H1
    Hb = H0 - H1

    UT = sp.expm(-1j * Hb * T / 2) @ sp.expm(-1j * Ha * T / 2)
    return np.exp(-1j * np.linalg.eigh(1j * sp.logm(UT))[0])


def UF_kick(omega, H0, H1):
    """
    Compute Floquet quasienergy phases for an "up-down" delta kick drive.

    Parameters
    ----------
    omega : float
        Drive frequency.
    H0 : ndarray, shape (d, d)
        Static Hamiltonian.
    H1 : ndarray, shape (d, d)
        Kick Hamiltonian.

    Returns
    -------
    phases : ndarray, shape (d,)
        Quasienergy phases e^{-i ε T}, where ε are quasienergies.

    Notes
    -----
    Protocol:
        U(T) = e^{-i H0 T/2} e^{+i H1} e^{-i H0 T/2} e^{-i H1}.
    """
    T = 2 * np.pi / omega
    UT = (sp.expm(-1j * H0 * T / 2)
          @ sp.expm(1j * H1)
          @ sp.expm(-1j * H0 * T / 2)
          @ sp.expm(-1j * H1))
    return np.exp(-1j * np.linalg.eigh(1j * sp.logm(UT))[0])


def UF_W2_and_W13(hx, hy, hz, omega):
    """
    Compute unitary evolution over one period for a bimodal Walsh drive.

    Parameters
    ----------
    hx, hy, hz : float
        Field strengths in X, Y, Z directions.
    omega : float
        Drive frequency.

    Returns
    -------
    U : ndarray, shape (2, 2)
        Single-period time evolution operator (approximate, via Trotter).
    """
    T = 2 * np.pi / omega
    N = 16
    dt = T / N

    U = np.eye(2, dtype=complex)

    W2 = compute_seq_order_walsh()[2][1]
    W13 = compute_seq_order_walsh()[13][1]

    for n in range(N):
        W2val = W2[n]
        W13val = W13[n]
        H = hz * PAULI_Z + W2val * hx * PAULI_X + W13val * hy * PAULI_Y
        U = expm(-1j * H * dt) @ U

    return U


# --------------------------------------------------------------------------- #
# Time-resolved Floquet evolution (kicked drive)
# --------------------------------------------------------------------------- #

def U_kick_t(t, omega, H0, H1):
    """
    Compute time-evolution operator U(t) for the up–down kick drive.

    Parameters
    ----------
    t : float
        Time at which to evaluate the evolution (t >= 0).
    omega : float
        Drive frequency.
    H0 : ndarray, shape (d, d)
        Static Hamiltonian.
    H1 : ndarray, shape (d, d)
        Kick Hamiltonian.

    Returns
    -------
    U : ndarray, shape (d, d)
        Time-evolution operator at time t.

    Notes
    -----
    The evolution is decomposed as t = n T + t_rem and the protocol
    for the kicks is applied piecewise over the residual time t_rem.
    """
    T = 2 * np.pi / omega
    UT = (sp.expm(-1j * H0 * T / 2)
          @ sp.expm(1j * H1)
          @ sp.expm(-1j * H0 * T / 2)
          @ sp.expm(-1j * H1))

    n_full = int(t // T)
    t_rem = t % T

    U = np.linalg.matrix_power(UT, n_full)

    if t_rem == 0:
        return U
    if t_rem < T / 2:
        U = sp.expm(-1j * H0 * t_rem) @ sp.expm(-1j * H1) @ U
    else:
        U = (sp.expm(-1j * H0 * (t_rem - T / 2))
             @ sp.expm(1j * H1)
             @ sp.expm(-1j * H0 * T / 2)
             @ sp.expm(-1j * H1)
             @ U)
    return U


def U_kick_tt0(t, t0, omega, H0, H1):
    """
    Compute evolution operator U(t, t0) from time t0 to t.

    Parameters
    ----------
    t : float
        Final time.
    t0 : float
        Initial time.
    omega : float
        Drive frequency.
    H0 : ndarray
        Static Hamiltonian.
    H1 : ndarray
        Kick Hamiltonian.

    Returns
    -------
    U : ndarray
        Evolution operator from t0 to t.
    """
    return U_kick_t(t, omega, H0, H1) @ U_kick_t(t0, omega, H0, H1).conj().T


def P_kick_t(t, omega, H0, H1):
    """
    Compute micromotion operator P(t) for "up-down" delta kick drive.

    Parameters
    ----------
    t : float
        Time at which to evaluate P(t).
    omega : float
        Drive frequency.
    H0 : ndarray
        Static Hamiltonian.
    H1 : ndarray
        Kick Hamiltonian.

    Returns
    -------
    P : ndarray
        Micromotion operator P(t) such that U(t) ≈ P(t) e^{-i H_F t}.
    """
    T = 2 * np.pi / omega
    UT = U_kick_tt0(T, 0, omega, H0, H1)
    HF = 1j * sp.logm(UT) / T
    return U_kick_tt0(t, 0, omega, H0, H1) @ sp.expm(1j * HF * t)


def u_kick_t(t, omega, H0, H1, u0):
    """
    Evolve a Floquet mode, u0, under the up–down kick drive.

    Parameters
    ----------
    t : float
        Time at which to evaluate the state.
    omega : float
        Drive frequency.
    H0 : ndarray
        Static Hamiltonian.
    H1 : ndarray
        Kick Hamiltonian.
    u0 : ndarray, shape (d,)
        Initial state vector.

    Returns
    -------
    ut : ndarray, shape (d,)
        Evolved state at time t.
    """
    return P_kick_t(t, omega, H0, H1) @ u0


def HFt0(t0, omega, H0, H1):
    """
    Compute Floquet Hamiltonian H_F[t0] at an arbitrary initial time t0.

    Parameters
    ----------
    t0 : float
        Initial time.
    omega : float
        Drive frequency.
    H0 : ndarray
        Static Hamiltonian.
    H1 : ndarray
        Kick Hamiltonian.

    Returns
    -------
    HF : ndarray
        Floquet Hamiltonian defined via U(T + t0, t0) = exp(-i H_F T).
    """
    T = 2 * np.pi / omega
    return 1j * sp.logm(U_kick_tt0(T + t0, t0, omega, H0, H1)) / T


def compute_Floquet_mode_t(n_points, omega, hz, hx):
    """
    Compute the Floquet mode u(t) over one period for an "up-down" kicked single spin.

    Parameters
    ----------
    n_points : int
        Number of time samples (2^n_points total).
    omega : float
        Drive frequency.
    hz : float
        Static Z field.
    hx : float
        Kick strength along X.

    Returns
    -------
    FloquetModes : ndarray, shape (2, 2^n_points), complex
        Time-dependent Floquet mode components for spin-up and spin-down.
    """
    T = 2 * np.pi / omega
    H0 = hz * PAULI_Z
    H1 = hx * PAULI_X

    UT = U_kick_t(T, omega, H0, H1)
    _HF = 1j * sp.logm(UT) / T  # kept for reference, not explicitly used

    # Pick one eigenvector of U(T)
    u0 = np.linalg.eig(UT)[1][:, 0]

    s = 2 ** n_points
    FloquetModes = np.zeros((2, s), dtype=complex)

    for i in range(s):
        t = T * i / s
        ut = u_kick_t(t, omega, H0, H1, u0)
        FloquetModes[:, i] = ut

    return FloquetModes


def compute_Floquet_mode_t_up(n_points, omega, hz, hx):
    """
    Compute spin-up Floquet mode u_up(t) and its Walsh coefficients 
    for an "up-down" delta kick drive.

    Parameters
    ----------
    n_points : int
        Power-of-two exponent for number of time samples, 2^n_points.
    omega : float
        Drive frequency.
    hz : float
        Static Z field.
    hx : float
        Kick strength along X.

    Returns
    -------
    FM : ndarray, shape (2^n_points,), complex
        Time series of the spin-up component.
    FM_walsh : ndarray, shape (2^n_points,), complex
        Walsh coefficients of the spin-up component.
    """
    FM = compute_Floquet_mode_t(n_points, omega, hz, hx)[0, :]
    FM_walsh = hadamard(2 ** n_points) @ FM / (2 ** n_points)
    return FM, FM_walsh


def compute_Floquet_mode_t_down(n_points, omega, hz, hx):
    """
    Compute spin-down Floquet mode u_down(t) and its Walsh coefficients
    for an "up-down" delta kick drive.

    Parameters
    ----------
    n_points : int
        Power-of-two exponent for number of time samples, 2^n_points.
    omega : float
        Drive frequency.
    hz : float
        Static Z field.
    hx : float
        Kick strength along X.

    Returns
    -------
    FM : ndarray, shape (2^n_points,), complex
        Time series of the spin-down component.
    FM_walsh : ndarray, shape (2^n_points,), complex
        Walsh coefficients of the spin-down component.
    """
    FM = compute_Floquet_mode_t(n_points, omega, hz, hx)[1, :]
    FM_walsh = hadamard(2 ** n_points) @ FM / (2 ** n_points)
    return FM, FM_walsh


# --------------------------------------------------------------------------- #
# Measuring error and localisation
# --------------------------------------------------------------------------- #

def compute_PPE(H, Hilbertd):
    """
    Compute photon participation entropy (PPE) of a Floquet Hamiltonian.
    Measure of localisation where we trace out the physical degrees of freedom.

    Parameters
    ----------
    H : ndarray, shape (M, M)
        Extended-space Hamiltonian, with photon ⊗ spin structure.
    Hilbertd : int
        Dimension of the physical (spin) Hilbert space.

    Returns
    -------
    S_mean : float
        Mean photon participation entropy averaged over all eigenstates.

    Notes
    -----
    The eigenvectors are reshaped into (N_photon, Hilbertd, ...), and
    probabilities are summed over the spin index to obtain photon-sector
    distributions.
    """
    M = H.shape[0]
    N = M // Hilbertd
    vals = []

    eigvals, eigvecs = np.linalg.eigh(H)
    for i in range(M):
        P = np.abs(eigvecs[:, i]) ** 2
        P_photon = P.reshape(N, Hilbertd, -1).sum(axis=1)
        PlogP = -P_photon * np.log(np.clip(P_photon, 1e-15, None))
        vals.append(PlogP.sum())

    return float(np.mean(vals))


def compute_MFIM_kick_errors(omega, J, hzs, hxs, L, n):
    """
    Compare Walsh vs Fourier Floquet errors for a multi-spin Ising model 
    for an "up-down" delta kick drive.

    Parameters
    ----------
    omega : float
        Drive frequency.
    J : float
        Ising ZZ coupling strength.
    hzs : array_like
        Grid of longitudinal fields h_z to scan.
    hxs : array_like
        Grid of transverse fields h_x to scan.
    L : int
        Number of spins.
    n : int
        Walsh depth such that N = 2^n is number of Walsh modes.

    Returns
    -------
    Errors_F : ndarray, shape (len(hzs), len(hxs))
        Phase errors (Fourier basis) vs exact Floquet spectrum.
    Errors_W : ndarray, shape (len(hzs), len(hxs))
        Phase errors (Walsh basis) vs exact Floquet spectrum.

    Notes
    -----
    - Builds extended-space quasienergy operators Q_walsh and Q_fourier.
    - Extracts the physical first Floquet Brillouin zone.
    - Compares phases against the exact two-kick unitary UF_kick.
    - Uses a rolling scheme in phase space to handle branch cuts and
      zone-boundary crossings robustly. Picks best possible alignment.
    """

    N = 2 ** n
    maxmode = int((N - 2) / 2)
    M = 2 ** L * N

    T = 2 * np.pi / omega

    drive_fourier = omega / np.pi * compute_kick_fourier(2 * maxmode + 1)
    drive_walsh = omega / np.pi * compute_kick_walsh(n)

    m1 = len(hzs)
    m2 = len(hxs)

    Errors_F = np.empty((m1, m2))
    Errors_W = np.empty((m1, m2))

    for i, hz in enumerate(hzs):
        for j, hx in enumerate(hxs):
            H0 = H_ising_zz_zfield(J, hz, L)
            H1 = H_ising_xfield(hx, L)

            # Construct extended-space quasienergy operators
            Q_F = Q_fourier(H0, H1, drive_fourier, omega, maxmode)
            Q_W = Q_walsh(H0, H1, drive_walsh, omega, n)

            sol_F = np.linalg.eigh(Q_F)
            sol_W = np.linalg.eigh(Q_W)

            # Exact numerical solution for comparison
            numsol = 1j * np.log(UF_kick(omega, H0, H1))
            conj_numsol = np.conjugate(np.exp(-1j * np.sort(numsol)))

            E_W = sol_W[0][int(M / 2) - 2 ** L:int(M / 2)] * T
            E_W_phases = np.exp(-1j * np.sort(E_W))

            E_F = sol_F[0][maxmode * 2 ** L:(maxmode + 1) * 2 ** L] * T
            E_F_phases = np.exp(-1j * np.sort(E_F))

            Error_W = np.average(np.abs(np.angle(E_W_phases * conj_numsol)))
            Error_F = np.average(np.abs(np.angle(E_F_phases * conj_numsol)))

            # Rolling correction around branch cuts
            for roll in range(-8, 9):
                EF_roll = np.roll(E_F_phases, roll)
                EW_roll = np.roll(E_W_phases, roll)

                trialF = np.average(np.abs(np.angle(EF_roll * conj_numsol)))
                trialW = np.average(np.abs(np.angle(EW_roll * conj_numsol)))

                # Pick best alignment
                if trialF < Error_F:
                    Error_F = trialF
                if trialW < Error_W:
                    Error_W = trialW

            Errors_W[i, j] = Error_W
            Errors_F[i, j] = Error_F

    return Errors_F, Errors_W

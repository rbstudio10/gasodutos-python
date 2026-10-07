# -*- coding: utf-8 -*-
"""
gasutil.py - Constants and helper functions shared by the gas pipeline programs.

Original author (HP PPL): Antonio Ricardo Andrade Bozolla
Theoretical basis: Stuckenbruck, S. - Escoamento em Dutos (Pipe Flow), Volume B
(PUC-Rio, 2014).
HP PPL -> Python 3 translation (standard library only: math).
"""
import math

# ----------------------------------------------------------------------
# Constants (cf. Table 8.3)
# ----------------------------------------------------------------------
G = 9.806            # m/s2
R_AR = 287.0         # J/(kg.K)
P_STD = 101325.0     # Pa
T_STD = 293.15       # K
Z_STD = 1.0
RHO_AR_STD = 1.2043  # kg/m3 (air at standard conditions)
C1_TEORICO = 13.305
C1_AGA = 13.303
TOL = 1e-6           # relative tolerance on f (fixed-point loop)
MAXIT = 30           # maximum iterations of the fixed-point loop
N_COLEBROOK = 10     # ITERATE(...) iterations of the PPL

# Empirical models: (name, C, exponent of lambda, exponent of the radicand,
#                    exponent of D)
#   Q = C * n * (T_std/P_std) * [(P1^2-P2^2-Grav)/(lam^a * L * Zm * Tm)]^b * D^c
MODELOS_EMPIRICOS = (
    ("WEYMOUTH",    137.32, 1.0,    0.5,    2.6667),
    ("PANHANDLE-A",  99.51, 0.8539, 0.5394, 2.6182),
    ("PANHANDLE-B", 137.24, 0.9608, 0.51,   2.53),
    ("IGT",          88.06, 0.8,    0.5555, 2.6667),
    ("MUELLER",      87.51, 0.74,   0.5747, 2.724),
    ("FRITZSCHE",    94.26, 0.8580, 0.5382, 2.6911),
)


class GasPipelineError(ValueError):
    """Equivalent to PRINT(message) + BREAK in the PPL code."""


def exigir(condicao_ok, mensagem):
    if not condicao_ok:
        raise GasPipelineError(mensagem)


def emitir(avisos, estrito, *linhas):
    """Diagnostic warning. estrito=True aborts (like the PPL BREAK);
    estrito=False only records the warning and carries on."""
    if estrito:
        raise GasPipelineError("\n".join(linhas))
    avisos.extend(linhas)


# ----------------------------------------------------------------------
# Validation of the common inputs (SI: L in m, P in Pa, Tm in K)
# ----------------------------------------------------------------------
def validar_comuns(L, eps, lam, mu, k, P1, P2, Zm, Tm, n, Ca, h):
    exigir(L >= 10e3, "The length of the duct must be >= 10_km")
    exigir(L <= 350e3, "The length of the duct must be <= 350_km")
    exigir(eps >= 0, "The absolute roughness must be >= 0_mm")
    exigir(eps <= 3, "The absolute roughness must be <= 3_mm")
    exigir(lam > 0, "The density must be > 0")
    exigir(lam <= 0.85, "The density must be <= 0.85")
    exigir(mu > 0, "The dynamic viscosity must be > 0_Pa*s")
    exigir(mu <= 1e-4, "The dynamic viscosity must be <= 1E-4_Pa*s")
    exigir(k >= 1.05, "The isentropic coeff. must be >= 1.05")
    exigir(k <= 1.8, "The isentropic coeff. must be <= 1.8")
    exigir(P1 <= 20.9e6, "P1 must be <= 20.9_MPa")
    exigir(P1 >= 0.5e6, "P1 must be >= 0.5_MPa")
    exigir(P2 < P1, "P2 must be < P1")   # the PPL tested P2>P1 (P2=P1 would give KER=0)
    exigir(P2 > 0, "P2 must be > 0_MPa")
    exigir(Zm >= 0.7, "The coeff. compressibility must be >= 0.7")
    exigir(Zm < 1, "The coeff. compressibility must be < 1")
    exigir(Tm >= 243.15, "The average temperature must be >= 243.15_K")
    exigir(Tm <= 393.15, "The average temperature must be <= 393.15_K")
    exigir(n > 0.80, "The duct efficiency must be > 0.8")
    exigir(n <= 1, "The duct efficiency must be <= 1")
    exigir(Ca > 0.9, "The coefficient for AGA-A must be > 0.9")
    exigir(Ca < 1, "The coefficient for AGA-A must be < 1")
    exigir(h <= 2e3, "The height difference (h) must be <= 2E3_m")
    exigir(h >= -1e3, "The height difference (h) must be >= -1E3_m")


# ----------------------------------------------------------------------
# Base quantities (Rg, rho_std, Grav, KER)
# ----------------------------------------------------------------------
class Base:
    pass


def estado_base(L, P1, P2, Tm, lam, Zm, h):
    b = Base()
    b.L, b.P1, b.P2, b.Tm = L, P1, P2, Tm
    b.Rg = R_AR / lam
    b.Pm = (2.0 / 3.0) * ((-P1 * P2 / (P1 + P2)) + P1 + P2)
    b.rho_s = RHO_AR_STD * lam
    b.grav = (2 * lam * b.Pm ** 2 / (Zm * R_AR * Tm)) * G * h
    b.num = P1 ** 2 - P2 ** 2 - b.grav
    exigir(b.num > 0, "P1^2 - P2^2 - Grav <= 0: no flow possible with these data")
    b.ker = (b.num / (lam * L * Zm * Tm)) ** 0.5
    return b


# ----------------------------------------------------------------------
# Hydraulics
# ----------------------------------------------------------------------
def colebrook(f0, eps_mm, D, Re, n=N_COLEBROOK):
    """ITERATE(1/(-2 log10(eps/(3.7D) + 2.51/(Re sqrt(f))))^2, f, f0, n)."""
    rr = eps_mm * 1e-3 / (3.7 * D)
    f = f0
    for _ in range(n):
        f = (1.0 / (-2.0 * math.log10(rr + 2.51 / (Re * math.sqrt(f))))) ** 2
    return f


def reynolds(Q, D, rho_s, mu):
    return 4 * rho_s * Q / (math.pi * mu * D)


def velocidade(Q, D, P2, Zm, Tm):
    """In-situ velocity at the outlet (P2)."""
    return (Zm * Tm / P2) * (P_STD / (Z_STD * T_STD)) * (Q / (math.pi * D ** 2 / 4))


def mach(V, k, Rg, Tm):
    return V / math.sqrt(k * Rg * Tm)


def velocidade_erosao(P2, Zm, Rg, Tm):
    return 120.0 / math.sqrt(P2 / (Zm * Rg * Tm))


def vazao_empirica(C, a, b, c, n, base, lam, Zm, D):
    return (C * n * (T_STD / P_STD)
            * (base.num / (lam ** a * base.L * Zm * base.Tm)) ** b * D ** c)


def diametro_empirico(C, a, b, c, n, base, lam, Zm, Q):
    fator = (C * n * (T_STD / P_STD)
             * (base.num / (lam ** a * base.L * Zm * base.Tm)) ** b)
    return (Q / fator) ** (1.0 / c)


# ----------------------------------------------------------------------
# Diagnostics (messages identical to the PPL ones)
# ----------------------------------------------------------------------
def verificar_reynolds(Re, avisos, estrito):
    if Re < 30e3:
        emitir(avisos, estrito, "Check dP, very low speed (Re<30E3)")
    elif Re > 1e8:
        emitir(avisos, estrito, "Check dP, very high speed (Re>1E8)")


def verificar_rugosidade(eps_mm, D, avisos, estrito):
    rr = eps_mm / (D * 1e3 / 2)
    if rr <= 2e-8:
        emitir(avisos, estrito, "Check absolute roughness value", "Very smooth tube")
    elif rr >= 0.1:
        emitir(avisos, estrito, "Check absolute roughness value", "Very rough tube")


# ----------------------------------------------------------------------
# Output formatting
# ----------------------------------------------------------------------
def formatar_bloco(i, nome, r, linha_final):
    f = f"{r['f']:.4e}" if r["f_sci"] else f"{r['f']:.0f}"
    return [
        f" {i}. {nome}:",
        f"  friction = {f}",
        f"  Reynolds = {r['Re']:.4e}",
        f"  Mach = {r['Ma']:.3f}",
        f"  Vel. = {r['V']:.1f} m/s",
        f"  {linha_final}",
        "",
    ]


def formatar_rodape(Ve, Vm, dPdx):
    return [
        " EROSION VELOCITY & MAX. RECOMMENDED:",
        f"  Erosion ~ {Ve:.1f} m/s",
        f"  Max. Recommended ~ {Vm:.1f} m/s",
        "",
        " AVERAGE PRESSURE GRADIENT:",
        f"  dP/dx = {dPdx:.1f} kPa/km",
    ]

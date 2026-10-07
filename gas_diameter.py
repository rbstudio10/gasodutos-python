# -*- coding: utf-8 -*-
"""
gas_diameter.py - Gas pipeline diameter (compressible flow)

Author: Antonio Ricardo Andrade Bozolla (revision R1 - 21/07/26)
Theoretical basis: Stuckenbruck, S. - Escoamento em Dutos (Pipe Flow), Volume B
(PUC-Rio, 2014)
Translation of the HP PPL program GasDiameter to Python 3.

Models: Theoretical (Colebrook), Weymouth, Panhandle A/B, IGT, Mueller,
Fritzsche, AGA-A e AGA-B.

Inputs (same units as the PPL program):
    L9   length [km]                 Q9   standard flow rate [Nm3/d]
    eps  absolute roughness [mm]     lam  gas relative density (air = 1) [-]
    mu   dynamic viscosity [Pa.s]    k    isentropic coefficient [-]
    P19  inlet pressure [bar abs]    P29  outlet pressure [bar abs]
    Zm   average compressibility factor [-]
    Tm9  average temperature [degC]  n    pipe efficiency [-]
    Ca   AGA-A coefficient [-]       h    elevation difference [m]
"""
import math
from gasutil import *


def gas_diameter(L9, Q9, eps, lam, mu, k, P19, P29, Zm, Tm9, n, Ca, h,
                 estrito=True):
    avisos = []
    L = L9 * 1000.0
    Q = Q9 / 86400.0
    P1 = P19 * 100e3
    P2 = P29 * 100e3
    Tm = Tm9 + 273.15

    exigir(L >= 10e3, "The length of the duct must be >= 10_km")
    exigir(Q >= 0.5, "The standard flow must be > 0.5_m^3/s")
    validar_comuns(L, eps, lam, mu, k, P1, P2, Zm, Tm, n, Ca, h)

    b = estado_base(L, P1, P2, Tm, lam, Zm, h)
    res = {}

    # ---------------- THEORETICAL (convergence loop on f) ----------------
    f2, err, it = 0.02, 1.0, 0
    while err > TOL and it < MAXIT:
        C2 = 1.0 / math.sqrt(f2)
        Db = (Q / (n * C1_TEORICO * C2 * (T_STD / P_STD) * b.ker)) ** (1 / 2.5)
        R2 = reynolds(Q, Db, b.rho_s, mu)
        fx = colebrook(f2, eps, Db, R2)
        err = abs(fx - f2) / f2
        f2 = fx
        it += 1
    if err > TOL:
        avisos.append("THEORIC: fixed-point loop did not converge (maxit)")

    verificar_rugosidade(eps, Db, avisos, estrito)
    verificar_reynolds(R2, avisos, estrito)
    if Db > 0.889:
        emitir(avisos, estrito, "Review input data",
               "Maximum diameter reached (35_in)")
    V2b = velocidade(Q, Db, P2, Zm, Tm)
    res["THEORIC"] = dict(f=f2, Re=R2, Ma=mach(V2b, k, b.Rg, Tm), V=V2b,
                          D=Db, f_sci=True)

    # ---------------- EMPIRICAL MODELS ----------------
    for nome, C, a, bb, c in MODELOS_EMPIRICOS:
        D = diametro_empirico(C, a, bb, c, n, b, lam, Zm, Q)
        V = velocidade(Q, D, P2, Zm, Tm)
        R = reynolds(Q, D, b.rho_s, mu)
        verificar_reynolds(R, avisos, estrito)
        res[nome] = dict(f=1.0, Re=R, Ma=mach(V, k, b.Rg, Tm), V=V, D=D,
                         f_sci=False)

    # ---------------- AGA-A (convergence loop on f) ----------------
    f10a, Ra1, err, it = f2, R2, 1.0, 0
    while err > TOL and it < MAXIT:
        C2d = 2 * Ca * math.log10((Ra1 * math.sqrt(f10a)) / 2.51)
        Da1 = (Q / (n * C1_AGA * C2d * (T_STD / P_STD) * b.ker)) ** (1 / 2.5)
        Ra1 = reynolds(Q, Da1, b.rho_s, mu)
        fx = colebrook(f10a, eps, Da1, Ra1)
        err = abs(fx - f10a) / f10a
        f10a = fx
        it += 1
    if err > TOL:
        avisos.append("AGA-A: fixed-point loop did not converge (maxit)")
    verificar_reynolds(Ra1, avisos, estrito)
    Va1 = velocidade(Q, Da1, P2, Zm, Tm)
    res["AGA-A"] = dict(f=f10a, Re=Ra1, Ma=mach(Va1, k, b.Rg, Tm), V=Va1,
                        D=Da1, f_sci=True)

    # ---------------- AGA-B (convergence loop on f) ----------------
    exigir(eps > 0, "AGA-B requires absolute roughness > 0_mm")
    Db1, fb, err, it = Db, f2, 1.0, 0
    while err > TOL and it < MAXIT:
        C2e = 2 * math.log10(3.7 * Db1 / (eps * 1e-3))
        Db1 = (Q / (n * C1_AGA * C2e * (T_STD / P_STD) * b.ker)) ** (1 / 2.5)
        Rb = reynolds(Q, Db1, b.rho_s, mu)
        fx = colebrook(fb, eps, Db1, Rb)
        err = abs(fx - fb) / fb
        fb = fx
        it += 1
    if err > TOL:
        avisos.append("AGA-B: fixed-point loop did not converge (maxit)")
    verificar_reynolds(Rb, avisos, estrito)
    Vb = velocidade(Q, Db1, P2, Zm, Tm)
    res["AGA-B"] = dict(f=fb, Re=Rb, Ma=mach(Vb, k, b.Rg, Tm), V=Vb,
                        D=Db1, f_sci=True)

    Ve = velocidade_erosao(P2, Zm, b.Rg, Tm)
    return dict(modelos=res, Ve=Ve, Vm=0.4 * Ve, dPdx=(P1 - P2) / L,
                avisos=avisos)


def formatar(r):
    linhas = [" == GAS PIPELINE DIAMETERS =="]
    for a in r["avisos"]:
        linhas.append(" ! " + a)
    for i, (nome, m) in enumerate(r["modelos"].items(), 1):
        linhas += formatar_bloco(i, nome, m, f"ID = {m['D'] * 39.37:.2f} in.")
    linhas += formatar_rodape(r["Ve"], r["Vm"], r["dPdx"])
    return "\n".join(linhas)


if __name__ == "__main__":
    # Example: edit the values below and run.
    try:
        r = gas_diameter(L9=100, Q9=3.0e6, eps=0.046, lam=0.6, mu=1.1e-5,
                         k=1.3, P19=70, P29=40, Zm=0.9, Tm9=25, n=0.92,
                         Ca=0.95, h=0)
        print(formatar(r))
    except GasPipelineError as e:
        print("INPUT ERROR:", e)

# -*- coding: utf-8 -*-
"""
gas_flow_rate.py - Calculo de Vazao para Gasodutos (Escoamento Compressivel)

Autor: Antonio Ricardo Andrade Bozolla (revisao R1 - 21/07/26)
Base teorica: Stuckenbruck, S. - Escoamento em Dutos, Volume B (PUC-Rio, 2014)
Traducao de HP PPL (GasFlowRate) para Python 3.

Entradas (mesmas unidades do programa PPL):
    L9   comprimento [km]            D    diametro interno [m]
    eps  rugosidade absoluta [mm]    lam  densidade relativa ao ar [-]
    mu   viscosidade dinamica [Pa.s] k    coef. isentropico [-]
    P19  pressao de entrada [bar abs]  P29 pressao de saida [bar abs]
    Zm   fator de compressibilidade medio [-]
    Tm9  temperatura media [degC]    n    eficiencia do duto [-]
    Ca   coeficiente do AGA-A [-]    h    desnivel [m]

AGA-B: com D dado, C2 = 2*log10(3.7*D/eps) e explicito (passe unico).
"""
import math
from gasutil import *


def gas_flow_rate(L9, D, eps, lam, mu, k, P19, P29, Zm, Tm9, n, Ca, h,
                  estrito=True):
    avisos = []
    L = L9 * 1000.0
    P1 = P19 * 100e3
    P2 = P29 * 100e3
    Tm = Tm9 + 273.15

    exigir(L >= 10e3, "The length of the duct must be >= 10_km")
    exigir(D >= 0.1016, "The minimum diameter must be >= 0.1016_m")
    exigir(D <= 0.8001, "The maximum diameter must be <= 0.8001_m")
    validar_comuns(L, eps, lam, mu, k, P1, P2, Zm, Tm, n, Ca, h)

    b = estado_base(L, P1, P2, Tm, lam, Zm, h)
    verificar_rugosidade(eps, D, avisos, estrito)
    res = {}

    # ---------------- TEORICO (laco de convergencia em f) ----------------
    f3, err, it = 0.02, 1.0, 0
    while err > TOL and it < MAXIT:
        C2 = 1.0 / math.sqrt(f3)
        Q2 = C1_TEORICO * C2 * n * (T_STD / P_STD) * b.ker * D ** 2.5
        R3 = reynolds(Q2, D, b.rho_s, mu)
        fx = colebrook(f3, eps, D, R3)
        err = abs(fx - f3) / f3
        f3 = fx
        it += 1
    if err > TOL:
        avisos.append("THEORIC: fixed-point loop did not converge (maxit)")
    verificar_reynolds(R3, avisos, estrito)
    V2b = velocidade(Q2, D, P2, Zm, Tm)
    res["THEORIC"] = dict(f=f3, Re=R3, Ma=mach(V2b, k, b.Rg, Tm), V=V2b,
                          Q=Q2, f_sci=True)

    # ---------------- MODELOS EMPIRICOS ----------------
    for nome, C, a, bb, c in MODELOS_EMPIRICOS:
        Q = vazao_empirica(C, a, bb, c, n, b, lam, Zm, D)
        V = velocidade(Q, D, P2, Zm, Tm)
        R = reynolds(Q, D, b.rho_s, mu)
        verificar_reynolds(R, avisos, estrito)
        res[nome] = dict(f=1.0, Re=R, Ma=mach(V, k, b.Rg, Tm), V=V, Q=Q,
                         f_sci=False)

    # ---------------- AGA-A (laco de convergencia em f) ----------------
    fa, Ra1, err, it = f3, R3, 1.0, 0
    while err > TOL and it < MAXIT:
        C2d = 2 * Ca * math.log10((Ra1 * math.sqrt(fa)) / 2.51)
        Qa1 = C1_AGA * C2d * n * (T_STD / P_STD) * b.ker * D ** 2.5
        Ra1 = reynolds(Qa1, D, b.rho_s, mu)
        fx = colebrook(fa, eps, D, Ra1)
        err = abs(fx - fa) / fa
        fa = fx
        it += 1
    if err > TOL:
        avisos.append("AGA-A: fixed-point loop did not converge (maxit)")
    verificar_reynolds(Ra1, avisos, estrito)
    Va1 = velocidade(Qa1, D, P2, Zm, Tm)
    res["AGA-A"] = dict(f=fa, Re=Ra1, Ma=mach(Va1, k, b.Rg, Tm), V=Va1,
                        Q=Qa1, f_sci=True)

    # ---------------- AGA-B (passe unico) ----------------
    exigir(eps > 0, "AGA-B requires absolute roughness > 0_mm")
    C2e = 2 * math.log10(3.7 * D / (eps * 1e-3))
    Qb = C1_AGA * C2e * n * (T_STD / P_STD) * b.ker * D ** 2.5
    Vb = velocidade(Qb, D, P2, Zm, Tm)
    Rb = reynolds(Qb, D, b.rho_s, mu)
    verificar_reynolds(Rb, avisos, estrito)
    fb = colebrook(0.02, eps, D, Rb)   # apenas para exibicao
    res["AGA-B"] = dict(f=fb, Re=Rb, Ma=mach(Vb, k, b.Rg, Tm), V=Vb, Q=Qb,
                        f_sci=True)

    Ve = velocidade_erosao(P2, Zm, b.Rg, Tm)
    return dict(modelos=res, Ve=Ve, Vm=0.4 * Ve, dPdx=(P1 - P2) / L,
                avisos=avisos)


def formatar(r):
    linhas = [" == GAS PIPELINE FLOW RATE =="]
    for a in r["avisos"]:
        linhas.append(" ! " + a)
    for i, (nome, m) in enumerate(r["modelos"].items(), 1):
        linhas += formatar_bloco(i, nome, m,
                                 f"Qstd = {m['Q'] * 86400:.0f} Nm3/d")
    linhas += formatar_rodape(r["Ve"], r["Vm"], r["dPdx"])
    return "\n".join(linhas)


if __name__ == "__main__":
    # Exemplo: edite os valores abaixo e execute.
    try:
        r = gas_flow_rate(L9=14, D=0.15, eps=0.020, lam=0.72, mu=1.03e-5,
                          k=1.46, P19=9.2, P29=6, Zm=0.985, Tm9=20.5, n=1,
                          Ca=0.97, h=910)
        print(formatar(r))
    except GasPipelineError as e:
        print("ERRO DE ENTRADA:", e)

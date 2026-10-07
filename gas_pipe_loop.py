# -*- coding: utf-8 -*-
"""
gas_pipe_loop.py - Escoamento Compressivel em Gasoduto com Loop

Autor: Antonio Ricardo Andrade Bozolla (versao aprimorada 14/03/26)
Base teorica: Stuckenbruck, S. - Escoamento em Dutos, Volume B (PUC-Rio, 2014)
Traducao de HP PPL (GasPipeLoop) para Python 3.

Topologia: A --(1)--> B --(2: B-C-E)--+--> E --(4)--> F
                                  \\--(3: B-D-E)--/

Entradas (mesmas unidades do programa PPL):
    P1        pressao em A [bar abs]       lam   densidade relativa [-]
    L1..L4    comprimentos [km]            D1..D4 diametros internos [in]
    Zm        compressibilidade media [-]  Tm    temperatura media [K]  (!)
    f         fator de atrito de Darcy [-] n     eficiencia do duto [-]
    Q         vazao std total [Nm3/s]
    Za,Zb,Ze,Zf  altitudes dos pontos A, B, E, F [m]
    n_iter    iteracoes do ponto fixo de pressao (ITERATE do PPL = 10)

Atencao: aqui Tm e dado em kelvin (como no PPL), diferente de
gas_diameter/gas_flow_rate, que recebem graus Celsius.
"""
import math
from gasutil import *

C1 = 13.305
EXP_C = 2.5      # expoente de D
EXP_B = 0.5      # expoente de Q


def _pm(pa, pb):
    return (2.0 / 3.0) * (pa + pb - (pa * pb / (pa + pb)))


def gas_pipe_loop(P1, lam, L1, L2, L3, L4, D1, D2, D3, D4, Zm, Tm, f, n, Q,
                  Za=749.0, Zb=1017.0, Ze=1095.0, Zf=1172.0, n_iter=10):
    # ---------------- clausulas iniciais ----------------
    exigir(P1 >= 10, "P1 must be >= 10 bar")
    exigir(P1 <= 200, "P1 must be <= 200 bar")
    exigir(lam > 0.3, "lambda must be > 0.3")
    exigir(lam <= 0.85, "lambda must be <= 0.85")
    for nome, L in (("L1", L1), ("L2", L2), ("L3", L3), ("L4", L4)):
        exigir(L >= 1, f"{nome} must be >= 1 km")
        exigir(L <= 100, f"{nome} must be <= 100 km")
    exigir(D1 >= 4.016, "D1 must be >= 4.016 in")
    exigir(D1 >= D2, "D1 must be >= D2")
    exigir(D1 >= D3, "D1 must be >= D3")
    exigir(D1 <= 36, "D1 must be <= 36 in")
    exigir(D2 >= 4.016, "D2 must be >= 4.016 in")
    exigir(D2 <= 36, "D2 must be <= 36 in")
    exigir(D3 >= 4.016, "D3 must be >= 4.016 in")
    exigir(D3 <= D1, "D3 must be <= D1")
    exigir(D3 <= 36, "D3 must be <= 36 in")
    exigir(D4 >= 4.016, "D4 must be >= 4.016 in")
    exigir(D4 <= 36, "D4 must be <= 36 in")
    exigir(Q >= 1, "The flow rate (std) must be > 1 Nm3/s")
    exigir(f > 0, "The friction factor must be > 0")
    exigir(f < 0.045, "The friction factor must be < 0.045")
    exigir(n > 0.807, "The efficiency of the duct must be > 0.807")
    exigir(n <= 1, "The efficiency of the duct must be <= 1")
    exigir(Zm >= 0.7, "The compressibility should be >= 0.7")
    exigir(Zm <= 1, "The compressibility should be <= 1")
    exigir(Tm >= 243.15, "The avg temperature should be >= 243.15 K")
    exigir(Tm <= 393.15, "The avg temperature should be <= 393.15 K")

    # ---------------- conversoes ----------------
    p1 = P1 * 1e5
    l1, l2, l3, l4 = (x * 1e3 for x in (L1, L2, L3, L4))
    d1, d2, d3, d4 = (x * 0.0254 for x in (D1, D2, D3, D4))

    C2 = 1.0 / math.sqrt(f)
    Rg = R_AR / lam
    fator = n * C1 * C2 * T_STD / P_STD
    inv_b = 1.0 / EXP_B

    def Kt(L, D):
        return lam * L * Zm * Tm / (fator * D ** EXP_C) ** inv_b

    Kt_ab, Kt_bce = Kt(l1, d1), Kt(l2, d2)
    Kt_bde, Kt_ef = Kt(l3, d3), Kt(l4, d4)
    KT = (Kt_bce ** -0.5 + Kt_bde ** -0.5) ** -2

    def coef_grav(za, zb):
        return -2 * lam * G * (zb - za) / (R_AR * Zm * Tm)

    def resolver(pa, kt_eq, za, zb):
        """ITERATE(sqrt(coefg*Pm(pa,p)^2 - Kt*Q^(1/b) + pa^2), p, 0.6*pa, n_iter)"""
        p = 0.6 * pa
        cg = coef_grav(za, zb)
        for _ in range(n_iter):
            arg = cg * _pm(pa, p) ** 2 - kt_eq * Q ** inv_b + pa ** 2
            exigir(arg > 0, "No solution: flow rate too high for the given "
                            "pressure/diameters")
            p = math.sqrt(arg)
        return p

    p2 = resolver(p1, Kt_ab, Za, Zb)    # ponto B
    p3 = resolver(p2, KT, Zb, Ze)       # ponto E
    p4 = resolver(p3, Kt_ef, Ze, Zf)    # ponto F

    # ---------------- vazoes nos ramos do loop ----------------
    grav_c = (2 * lam * G * (Ze - Zb) / (R_AR * Zm * Tm)) * _pm(p2, p3) ** 2
    Q_bce = math.sqrt((p2 ** 2 - p3 ** 2 - grav_c) / Kt_bce)
    Q_bde = Q - Q_bce

    # ---------------- velocidades in situ e de erosao ----------------
    V_b = velocidade(Q, d1, p2, Zm, Tm)
    V_bce = velocidade(Q_bce, d2, p3, Zm, Tm)
    V_bde = velocidade(Q_bde, d3, p3, Zm, Tm)
    V_f = velocidade(Q, d4, p4, Zm, Tm)
    Vers_b = velocidade_erosao(p2, Zm, Rg, Tm)
    Vers_e = velocidade_erosao(p3, Zm, Rg, Tm)   # vale p/ bce e bde (mesmo P3)
    Vers_f = velocidade_erosao(p4, Zm, Rg, Tm)

    Q_bde_alt = math.sqrt(max(p2 ** 2 - p3 ** 2 - grav_c, 0.0) / Kt_bde)

    return dict(
        Kt=dict(ab=Kt_ab, bce=Kt_bce, bde=Kt_bde, ef=Kt_ef),
        P_barg=dict(B=p2 * 1e-5 - 1.01325, E=p3 * 1e-5 - 1.01325,
                    F=p4 * 1e-5 - 1.01325),
        V=dict(ab=V_b, bce=V_bce, bde=V_bde, ef=V_f),
        Vers=dict(ab=Vers_b, bce=Vers_e, bde=Vers_e, ef=Vers_f),
        Q_bce=Q_bce, Q_bde=Q_bde,
        fechamento_Q_bde=Q_bde_alt / Q_bde - 1.0 if Q_bde else float("nan"),
        dPdx=dict(ab=(p1 - p2) / l1, bce=(p2 - p3) / l2,
                  bde=(p2 - p3) / l3, ef=(p3 - p4) / l4),
    )


def formatar(r):
    L = [" == RESULTS FOR GAS PIPELINE LOOP ==",
         " 1. SEGMENT RESIST. COEFFICIENT (Kt):"]
    for s in ("ab", "bce", "bde", "ef"):
        L.append(f"  {s} = {r['Kt'][s]:.4e}")
    L += ["", " 2. PRESSURES POINTS B, E, AND F:"]
    for p in ("B", "E", "F"):
        L.append(f"  P{p} = {r['P_barg'][p]:.1f} bar(g)")
    L += ["", " 3. DOWNSTREAM VELOCITIES SEGMENTS:"]
    for s in ("ab", "bce", "bde", "ef"):
        L.append(f"  {s} = {r['V'][s]:.1f} m/s -> V(ers) = {r['Vers'][s]:.1f} m/s")
    L += ["", " 4. VOLUMETRIC FLOW RATES DUCTS C AND D:",
          f"  Qstd_(c) = {r['Q_bce'] * 86400:.1f} Nm3/d",
          f"  Qstd_(d) = {r['Q_bde'] * 86400:.1f} Nm3/d",
          "", " AVG. PRESSURE GRADIENT (dP/dx) SEGMENTS:"]
    for s in ("ab", "bce", "bde", "ef"):
        L.append(f"  {s} = {r['dPdx'][s]:.1f} kPa/km")
    return "\n".join(L)


if __name__ == "__main__":
    # Exemplo: edite os valores abaixo e execute.
    try:
        r = gas_pipe_loop(P1=83, lam=0.62, L1=56, L2=71, L3=71, L4=89,
                          D1=14, D2=10, D3=12, D4=12, Zm=0.91, Tm=299.15,
                          f=0.011, n=0.95, Q=25.4629,
                          Za=749, Zb=1017, Ze=1095, Zf=1172)
        print(formatar(r))
    except GasPipelineError as e:
        print("ERRO DE ENTRADA:", e)

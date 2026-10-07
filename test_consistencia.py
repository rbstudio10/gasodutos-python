# -*- coding: utf-8 -*-
"""
test_consistencia.py - Testes de consistencia interna (sem depender da HP Prime).

1) Ida e volta: gas_diameter(Q) -> D_modelo ; gas_flow_rate(D_modelo) -> Q.
   Cada modelo deve devolver a vazao original (tolerancia 1e-4).
2) Loop: Q_bce + Q_bde = Q e Q_bde recalculado pela equacao do ramo D
   deve coincidir com Q - Q_bce (fechamento), e sensibilidade a n_iter.
3) Casos de rejeicao de entrada (equivalentes aos BREAK do PPL).

Execute:  python test_consistencia.py
"""
from gas_diameter import gas_diameter
from gas_flow_rate import gas_flow_rate
from gas_pipe_loop import gas_pipe_loop
from gasutil import GasPipelineError

TOL_REL = 1e-4
falhas = 0


def ok(cond, msg):
    global falhas
    print(("  OK   " if cond else "  FALHA ") + msg)
    if not cond:
        falhas += 1


casos = [
    dict(L9=100, Q9=3.0e6, eps=0.046, lam=0.60, mu=1.1e-5, k=1.30,
         P19=70, P29=40, Zm=0.90, Tm9=25, n=0.92, Ca=0.95, h=0),
    dict(L9=180, Q9=5.0e6, eps=0.05, lam=0.65, mu=1.2e-5, k=1.28,
         P19=100, P29=60, Zm=0.88, Tm9=30, n=0.95, Ca=0.96, h=200),
    dict(L9=60, Q9=2.0e6, eps=0.03, lam=0.58, mu=1.0e-5, k=1.31,
         P19=50, P29=30, Zm=0.92, Tm9=20, n=0.90, Ca=0.94, h=-100),
]

print("1) Ida e volta diametro -> vazao")
for i, c in enumerate(casos, 1):
    rd = gas_diameter(**c)
    Q_alvo = c["Q9"] / 86400.0
    for nome, m in rd["modelos"].items():
        D = m["D"]
        if not (0.1016 <= D <= 0.8001):
            print(f"  (caso {i}, {nome}: D={D:.4f} m fora da faixa de gas_flow_rate; ignorado)")
            continue
        args = {k: v for k, v in c.items() if k != "Q9"}
        rf = gas_flow_rate(D=D, **args)
        Q_calc = rf["modelos"][nome]["Q"]
        erro = abs(Q_calc / Q_alvo - 1.0)
        ok(erro < TOL_REL, f"caso {i} {nome:12s} D={D*39.37:6.2f} in  erro rel. Q = {erro:.2e}")

print("\n2) Loop")
base = dict(P1=90, lam=0.6, L1=30, L2=25, L3=25, L4=20, D1=20, D2=12,
            D3=12, D4=20, Zm=0.9, Tm=298.15, f=0.012, n=0.92, Q=15.0)
r10 = gas_pipe_loop(**base)
r100 = gas_pipe_loop(**base, n_iter=100)
ok(abs(r10["Q_bce"] + r10["Q_bde"] - base["Q"]) < 1e-9, "Q_bce + Q_bde = Q")
ok(abs(r10["fechamento_Q_bde"]) < 1e-3,
   f"fechamento do ramo bde = {r10['fechamento_Q_bde']:.2e}")
for ponto in ("B", "E", "F"):
    d = abs(r10["P_barg"][ponto] - r100["P_barg"][ponto])
    ok(d < 1e-3, f"P_{ponto}: 10 vs 100 iteracoes difere {d:.2e} bar")
# ramos com diametros diferentes: bce (12") e bde (10") -> bde deve levar menos
assim = dict(base, D3=10)
ra = gas_pipe_loop(**assim)
ok(ra["Q_bde"] < ra["Q_bce"], "ramo de menor diametro leva menos vazao")
ok(abs(ra["fechamento_Q_bde"]) < 1e-3,
   f"fechamento (D3=10 in) = {ra['fechamento_Q_bde']:.2e}")

print("\n3) Rejeicao de entradas")
rejeicoes = [
    (lambda: gas_diameter(**dict(casos[0], L9=5)), "L < 10 km"),
    (lambda: gas_diameter(**dict(casos[0], P29=80)), "P2 > P1"),
    (lambda: gas_flow_rate(D=0.05, **{k: v for k, v in casos[0].items() if k != 'Q9'}), "D < 0.1016 m"),
    (lambda: gas_pipe_loop(**dict(base, D3=24)), "D3 > D1"),
    (lambda: gas_pipe_loop(**dict(base, Q=500.0)), "vazao excessiva"),
]
for fn, rotulo in rejeicoes:
    try:
        fn()
        ok(False, f"{rotulo}: deveria rejeitar")
    except GasPipelineError as e:
        ok(True, f"{rotulo}: rejeitado ({str(e).splitlines()[0]})")

print("\nRESULTADO:", "TUDO OK" if falhas == 0 else f"{falhas} FALHA(S)")

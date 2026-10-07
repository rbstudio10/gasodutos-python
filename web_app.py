# -*- coding: utf-8 -*-
"""
web_app.py - Local web interface (dark theme) for the gas pipeline programs:

    /diameter   gas_diameter.py   Diameter (known flow rate)
    /flow       gas_flow_rate.py  Flow rate (known diameter)
    /loop       gas_pipe_loop.py  Looped gas pipeline (pressures, flows, velocities)

How to use (Pydroid 3 on Android, or any Python 3 on a PC):
  1. Keep this file in the SAME folder as gas_diameter.py, gas_flow_rate.py,
     gas_pipe_loop.py and gasutil.py.
  2. Open web_app.py and run it.
  3. Open the browser at  http://127.0.0.1:8000
     (if the port is busy, the program uses the next one and prints the address).
  4. To quit, stop the program (in Pydroid: back to the editor and tap the stop button).

Standard library only. The server accepts connections from the device itself
only (127.0.0.1).

Calculation author: Antonio Ricardo Andrade Bozolla
Theoretical basis: Stuckenbruck, S. - Escoamento em Dutos (Pipe Flow), Volume B
(PUC-Rio, 2014)
"""
import html
import math
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from gas_diameter import gas_diameter
from gas_flow_rate import gas_flow_rate
from gas_pipe_loop import gas_pipe_loop
from gasutil import GasPipelineError

HOST = "127.0.0.1"
PORTA = 8000
POL = 39.37          # same factor as the PPL program (m -> in)
POL_M = 0.0254       # m per inch (exact input conversion)
P_ATM = 1.01325      # bar


# ----------------------------------------------------------------------
# Form fields
# ----------------------------------------------------------------------
def _c(nome, sim, rotulo, unidade, padrao, lo=None, hi=None,
       lo_ok=True, hi_ok=True, largo=False, dica="", rel=None, faixa=None,
       mil=False):
    """rel = ("<" | "<=", other_field_name): this field must be less than
    (or less than or equal to) the other one."""
    return dict(nome=nome, sim=sim, rotulo=rotulo, unidade=unidade,
                padrao=padrao, lo=lo, hi=hi, lo_ok=lo_ok, hi_ok=hi_ok,
                largo=largo, dica=dica, rel=rel, faixa=faixa, mil=mil)


def _compartilhados():
    return dict(
        L9=_c("L9", "L", "Length", "km", "100", 10, 350),
        eps=_c("eps", "ε", "Absolute roughness", "mm", "0.046", 0, 3, lo_ok=False),
        h=_c("h", "h", "Elevation difference", "m", "0", -1000, 2000,
             dica="Positive = uphill in the flow direction."),
        n=_c("n", "n", "Efficiency", "–", "0.92", 0.8, 1, lo_ok=False),
        Ca=_c("Ca", "Ca", "AGA-A coeff.", "–", "0.95", 0.9, 1, lo_ok=False, hi_ok=False),
        lam=_c("lam", "λ", "Relative density", "ar = 1", "0.6", 0, 0.85, lo_ok=False),
        mu=_c("mu", "μ", "Dynamic viscosity", "×10⁻⁵ Pa·s", "1.1", 0, 10,
              lo_ok=False, dica="1.1 means 1.1×10⁻⁵ Pa·s."),
        k=_c("k", "k", "Exponent k", "–", "1.3", 1.05, 1.8),
        Zm=_c("Zm", "Z", "Avg. Z factor", "–", "0.9", 0.7, 1, hi_ok=False),
        Tm9=_c("Tm9", "T", "Average temperature", "°C", "25", -30, 120),
        P19=_c("P19", "P₁", "Inlet pressure", "bar", "70", 5, 209),
        P29=_c("P29", "P₂", "Outlet pressure", "bar", "40", 0, None, lo_ok=False,
               dica="Must be lower than P₁.", rel=("<", "P19"), faixa="0 < P₂ < P₁"),
    )


def _padroes(grupos, novos):
    """Apply the initial values (same as the HP Prime webscreen)."""
    for _, cs in grupos:
        for c in cs:
            if c["nome"] in novos:
                c["padrao"] = novos[c["nome"]]
    return grupos


PADRAO_DIAMETRO = dict(L9="185", eps="0.018", Q9="2500000", h="405", n="0.95",
                       Ca="0.94", lam="0.56", mu="1.31", k="1.299", Zm="0.91",
                       Tm9="22", P19="98", P29="45")
PADRAO_VAZAO = dict(L9="14", D="5.9055", eps="0.02", h="910", n="1", Ca="0.97",
                    lam="0.72", mu="1.03", k="1.46", Zm="0.985", Tm9="20.5",
                    P19="9.213", P29="6.013")
PADRAO_LOOP = dict(Q9="2200000", P1="83", n="0.95", f="0.011", lam="0.62",
                   Zm="0.91", Tm9="26", L1="56", L2="71", L3="71", L4="89",
                   D1="14", D2="10", D3="12", D4="12")


def _grupos_diametro():
    b = _compartilhados()
    return _padroes([
        ("Pipe and flow", [
            b["L9"], b["eps"],
            _c("Q9", "Q", "Standard flow rate", "Nm³/d", "3000000", 43200, None,
               largo=True, mil=True, dica="Accepts 3000000 or 3,000,000."),
            b["h"], b["n"], b["Ca"]]),
        ("Gas", [b["lam"], b["mu"], b["k"], b["Zm"], b["Tm9"]]),
        ("Absolute pressures", [b["P19"], b["P29"]]),
    ], PADRAO_DIAMETRO)


def _grupos_vazao():
    b = _compartilhados()
    return _padroes([
        ("Pipe", [
            b["L9"], b["eps"],
            _c("D", "D", "Inside diameter", "in", "20", 4, 31.5, largo=True),
            b["h"], b["n"], b["Ca"]]),
        ("Gas", [b["lam"], b["mu"], b["k"], b["Zm"], b["Tm9"]]),
        ("Absolute pressures", [b["P19"], b["P29"]]),
    ], PADRAO_VAZAO)


def _grupos_loop():
    b = _compartilhados()
    alt = dict(lo=-500, hi=5000)
    return _padroes([
        ("Network and flow", [
            _c("Q9", "Q", "Total standard flow rate", "Nm³/d", "2200000", 86400, None,
               largo=True, mil=True, dica="2,200,000 Nm³/d ≈ 25.46 Nm³/s."),
            _c("P1", "P₁", "Pressure at A (abs.)", "bar", "90", 10, 200),
            _c("n", "n", "Efficiency", "–", "0.92", 0.807, 1, lo_ok=False),
            _c("f", "f", "Friction factor", "–", "0.012", 0, 0.045,
               lo_ok=False, hi_ok=False)]),
        ("Gas", [
            _c("lam", "λ", "Relative density", "ar = 1", "0.6", 0.3, 0.85, lo_ok=False),
            _c("Zm", "Z", "Avg. Z factor", "–", "0.9", 0.7, 1), b["Tm9"]]),
        ("Segment lengths", [
            _c("L1", "L₁", "Segment 1 (A–B)", "km", "30", 1, 100),
            _c("L2", "L₂", "Segment 2 (B–C–E)", "km", "25", 1, 100),
            _c("L3", "L₃", "Segment 3 (B–D–E)", "km", "25", 1, 100),
            _c("L4", "L₄", "Segment 4 (E–F)", "km", "20", 1, 100)]),
        ("Inside diameters", [
            _c("D1", "D₁", "Segment 1", "in", "20", 4.016, 36),
            _c("D2", "D₂", "Segment 2", "in", "12", 4.016, 36, rel=("<=", "D1"),
               faixa="4.016 ≤ D₂ ≤ D₁"),
            _c("D3", "D₃", "Segment 3", "in", "12", 4.016, 36, rel=("<=", "D1"),
               faixa="4.016 ≤ D₃ ≤ D₁"),
            _c("D4", "D₄", "Segment 4", "in", "20", 4.016, 36)]),
        ("Elevations", [
            _c("Za", "Zₐ", "Point A", "m", "749", **alt),
            _c("Zb", "Zᵦ", "Point B", "m", "1017", **alt),
            _c("Ze", "Zₑ", "Point E", "m", "1095", **alt),
            _c("Zf", "Zf", "Point F", "m", "1172", **alt)]),
    ], PADRAO_LOOP)


NOMES = {
    "THEORIC": ("Theoretical", "Theoretical (Colebrook)", "teo"),
    "WEYMOUTH": ("Weymouth", "Weymouth", "emp"),
    "PANHANDLE-A": ("Panhandle A", "Panhandle A", "emp"),
    "PANHANDLE-B": ("Panhandle B", "Panhandle B", "emp"),
    "IGT": ("IGT", "IGT", "emp"),
    "MUELLER": ("Mueller", "Mueller", "emp"),
    "FRITZSCHE": ("Fritzsche", "Fritzsche", "emp"),
    "AGA-A": ("AGA-A", "AGA-A", "aga"),
    "AGA-B": ("AGA-B", "AGA-B", "aga"),
}
CATEGORIAS = {"teo": "Theoretical", "emp": "Empirical", "aga": "AGA"}

AVISOS_PT = {
    "Check absolute roughness value": "Check the absolute roughness value.",
    "Very smooth tube": "Tube is too smooth for the method's range.",
    "Very rough tube": "Tube is too rough for the method's range.",
    "Check dP, very low speed (Re<30E3)":
        "Check ΔP: very low velocity (Re < 30,000) in at least one model.",
    "Check dP, very high speed (Re>1E8)":
        "Check ΔP: very high velocity (Re > 10⁸) in at least one model.",
    "Review input data": "Review the input data.",
    "Maximum diameter reached (35_in)": "Maximum diameter of the method reached (35 in).",
}


# ----------------------------------------------------------------------
# Numbers: parsing (accepts thousands separators) and en-US formatting
# ----------------------------------------------------------------------
def ler_numero(s):
    """Parse a number typed by the user.

    Accepts 3000000, 3,000,000, 3.000.000 and 1,5 (decimal comma) as well as
    the usual 0.046. A comma followed by exactly three digits is read as a
    thousands separator (1,500 = 1500), except after a leading zero
    (0,046 = 0.046)."""
    s = (s or "").strip().replace(" ", "")
    if not s:
        return None
    if "," in s and "." in s:
        if s.rfind(",") > s.rfind("."):          # 1.234,5 -> decimal comma
            s = s.replace(".", "").replace(",", ".")
        else:                                     # 1,234.5
            s = s.replace(",", "")
    elif "," in s:
        partes = s.lstrip("+-").split(",")
        if len(partes) > 2 or (len(partes) == 2 and len(partes[1]) == 3
                               and 1 <= len(partes[0]) <= 3
                               and not partes[0].startswith("0")):
            s = s.replace(",", "")                # 3,000,000 / 1,500
        else:
            s = s.replace(",", ".")               # 1,5
    elif s.count(".") > 1:
        s = s.replace(".", "")                    # 3.000.000
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def fnum(x, nd=2):
    return f"{x:.{nd}f}"


def fmil(x):
    return f"{x:,.0f}"


_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def fsci(x, nd=2):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    return f"{fnum(m, nd)}×10{str(e).translate(_SUP)}"


def faixa_txt(c):
    if c["faixa"]:
        return c["faixa"]
    lo, hi, s = c["lo"], c["hi"], c["sim"]
    f = fmil if c["mil"] else (lambda v: fnum(v, 3).rstrip("0").rstrip("."))
    a = f"{f(lo)} {'≤' if c['lo_ok'] else '<'} " if lo is not None else ""
    b = f" {'≤' if c['hi_ok'] else '<'} {f(hi)}" if hi is not None else ""
    if lo is not None and hi is None:
        return f"{s} {'≥' if c['lo_ok'] else '>'} {f(lo)}"
    return f"{a}{s}{b}"


# ----------------------------------------------------------------------
# Validation and calculation
# ----------------------------------------------------------------------
def validar(brutos, campos):
    vals, erros = {}, {}
    por = {c["nome"]: c for c in campos}
    for c in campos:
        v = ler_numero(brutos.get(c["nome"], ""))
        if v is None:
            erros[c["nome"]] = "Enter a number."
            continue
        lo, hi = c["lo"], c["hi"]
        if lo is not None and (v < lo if c["lo_ok"] else v <= lo):
            erros[c["nome"]] = "Out of range: " + faixa_txt(c)
        elif hi is not None and (v > hi if c["hi_ok"] else v >= hi):
            erros[c["nome"]] = "Out of range: " + faixa_txt(c)
        vals[c["nome"]] = v
    for c in campos:
        n = c["nome"]
        if c["rel"] and n in vals and n not in erros and c["rel"][1] in vals:
            op, o = c["rel"]
            a, b = vals[n], vals[o]
            if (a >= b) if op == "<" else (a > b):
                erros[n] = "%s must be %s %s." % (
                    c["sim"], "lower than" if op == "<" else "lower than or equal to",
                    por[o]["sim"])
    return vals, erros


def _calc_diametro(v):
    return gas_diameter(L9=v["L9"], Q9=v["Q9"], eps=v["eps"], lam=v["lam"],
                        mu=v["mu"] * 1e-5, k=v["k"], P19=v["P19"], P29=v["P29"],
                        Zm=v["Zm"], Tm9=v["Tm9"], n=v["n"], Ca=v["Ca"],
                        h=v["h"], estrito=False)


def _calc_vazao(v):
    return gas_flow_rate(L9=v["L9"], D=round(v["D"] * POL_M, 6), eps=v["eps"],
                         lam=v["lam"], mu=v["mu"] * 1e-5, k=v["k"],
                         P19=v["P19"], P29=v["P29"], Zm=v["Zm"], Tm9=v["Tm9"],
                         n=v["n"], Ca=v["Ca"], h=v["h"], estrito=False)


def _calc_loop(v):
    r = gas_pipe_loop(P1=v["P1"], lam=v["lam"], L1=v["L1"], L2=v["L2"],
                      L3=v["L3"], L4=v["L4"], D1=v["D1"], D2=v["D2"],
                      D3=v["D3"], D4=v["D4"], Zm=v["Zm"],
                      Tm=round(v["Tm9"] + 273.15, 6), f=v["f"], n=v["n"],
                      Q=v["Q9"] / 86400.0, Za=v["Za"], Zb=v["Zb"], Ze=v["Ze"],
                      Zf=v["Zf"])
    r["entrada"] = dict(v)
    r["avisos"] = []
    return r


def traduz_erro(msg):
    msg = msg.replace("\n", " ")
    if msg.startswith("No solution"):
        return ("no solution: the flow rate is too high for the given pressure "
                "and diameters (the pressure drops to zero in some segment).")
    return msg.replace("_", " ")


def calcular(prog, brutos):
    """Returns (result, global_error, per_field_errors)."""
    vals, erros = validar(brutos, prog["campos"])
    if erros:
        return None, "Fix the highlighted fields.", erros
    try:
        r = prog["calc"](vals)
    except GasPipelineError as e:
        return None, "Input rejected by the calculation: " + traduz_erro(str(e)), {}
    except (ArithmeticError, ValueError) as e:
        return None, f"The calculation failed with these data ({e}).", {}
    vistos, avisos = set(), []
    for a in r.get("avisos", []):
        t = AVISOS_PT.get(a)
        if t is None:
            t = ("Model %s: the iterative loop did not converge in 30 iterations."
                 % a.split(":")[0].title()) if "fixed-point" in a else a
        if t not in vistos:
            vistos.add(t)
            avisos.append(t)
    r["avisos_pt"] = avisos
    return r, "", {}


# ----------------------------------------------------------------------
# SVG charts
# ----------------------------------------------------------------------
def _ticks(lo, hi, n=4):
    span = hi - lo
    raw = span / n
    mag = 10 ** math.floor(math.log10(raw))
    step = mag
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if span / step <= n:
            break
    t = math.ceil(lo / step - 1e-9) * step
    out = []
    while t <= hi + 1e-9:
        out.append(round(t, 10))
        t += step
    nd = 0
    while abs(step * 10 ** nd - round(step * 10 ** nd)) > 1e-9 and nd < 4:
        nd += 1
    return out, nd, step


def _esc(s):
    return html.escape(str(s), quote=True)


def _grafico_pontos(modelos, val, aria, tip):
    """Dot plot, one point per model. val(m) -> number; tip(key, m) -> text."""
    linhas = list(modelos.items())
    n = len(linhas)
    W, X0, X1, RH, TOP, BOT = 360, 112, 316, 30, 8, 30
    H = TOP + n * RH + BOT
    vals = [val(m) for _, m in linhas]
    vmin, vmax = min(vals), max(vals)
    span = vmax - vmin
    pad = span * 0.14 if span > 0 else max(abs(vmax) * 0.02, 0.5)
    lo, hi = vmin - pad, vmax + pad
    ticks, nd, _ = _ticks(lo, hi, 4)

    def X(v):
        return X0 + (v - lo) / (hi - lo) * (X1 - X0)

    yb = TOP + n * RH
    p = [f'<svg class="chart" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="{_esc(aria)} from {fnum(vmin)} to {fnum(vmax)}. '
         f'A table with the values is provided below.">']
    for t in ticks:
        x = X(t)
        p.append(f'<line class="gl" x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{yb + 4}"/>')
        p.append(f'<text class="ax" x="{x:.1f}" y="{yb + 19}" text-anchor="middle">{fnum(t, nd)}</text>')
    p.append(f'<line class="bl" x1="{X0}" y1="{yb}" x2="{X1}" y2="{yb}"/>')
    xt = X(vals[0])
    p.append(f'<line class="ref" x1="{xt:.1f}" y1="{TOP}" x2="{xt:.1f}" y2="{yb}"/>')
    for i, (chave, m) in enumerate(linhas):
        curto, _, cat = NOMES[chave]
        cy = TOP + i * RH + RH / 2
        x = X(vals[i])
        p.append(f'<text class="rl" x="{X0 - 10}" y="{cy + 4:.1f}" text-anchor="end">{_esc(curto)}</text>')
        p.append(f'<circle class="dot c-{cat}" cx="{x:.1f}" cy="{cy:.1f}" r="6"/>')
        if vals[i] in (vmin, vmax):
            p.append(f'<text class="vl" x="{x + 11:.1f}" y="{cy + 4:.1f}">{fnum(vals[i])}</text>')
        p.append(f'<rect class="hit" x="0" y="{TOP + i * RH}" width="{W}" height="{RH}" '
                 f'tabindex="0" data-tip="{_esc(tip(chave, m))}"/>')
    p.append("</svg>")
    return "".join(p)


def _grafico_barras(modelos, Vm, Ve):
    linhas = list(modelos.items())
    n = len(linhas)
    W, X0, X1, RH, TOP, BOT, BH = 360, 112, 316, 28, 36, 30, 14
    H = TOP + n * RH + BOT
    vmax = max(max(m["V"] for _, m in linhas), Ve)
    ticks, nd, step = _ticks(0, vmax * 1.06, 4)
    xmax = math.ceil(vmax * 1.06 / step - 1e-9) * step
    ticks = [t for t in ticks if t <= xmax + 1e-9]

    def X(v):
        return X0 + v / xmax * (X1 - X0)

    yb = TOP + n * RH
    p = [f'<svg class="chart" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="Outlet velocity by model, in meters per second, '
         f'with the maximum recommended velocity of {fnum(Vm, 1)} and the erosion '
         f'velocity of {fnum(Ve, 1)}. A table with the values is provided below.">']
    for t in ticks:
        x = X(t)
        p.append(f'<line class="gl" x1="{x:.1f}" y1="{TOP - 4}" x2="{x:.1f}" y2="{yb}"/>')
        p.append(f'<text class="ax" x="{x:.1f}" y="{yb + 19}" text-anchor="middle">{fnum(t, nd)}</text>')
    p.append(f'<line class="bl" x1="{X0}" y1="{TOP - 4}" x2="{X0}" y2="{yb}"/>')
    for i, (chave, m) in enumerate(linhas):
        curto, _, cat = NOMES[chave]
        y = TOP + i * RH + (RH - BH) / 2
        x1 = X(m["V"])
        p.append(_barra(X0, x1, y, BH, f"bar c-{cat}"))
        ym = y + BH / 2 + 4
        p.append(f'<text class="rl" x="{X0 - 10}" y="{ym:.1f}" text-anchor="end">{_esc(curto)}</text>')
        p.append(f'<text class="vl" x="{x1 + 7:.1f}" y="{ym:.1f}">{fnum(m["V"], 1)}</text>')
        tip = f"{NOMES[chave][1]}: V = {fnum(m['V'], 1)} m/s · Mach {fnum(m['Ma'], 3)}"
        p.append(f'<rect class="hit" x="0" y="{TOP + i * RH}" width="{W}" height="{RH}" '
                 f'tabindex="0" data-tip="{_esc(tip)}"/>')
    xm, xe = X(Vm), X(Ve)
    p.append(f'<line class="th warn" x1="{xm:.1f}" y1="14" x2="{xm:.1f}" y2="{yb}"/>')
    p.append(f'<line class="th crit" x1="{xe:.1f}" y1="26" x2="{xe:.1f}" y2="{yb}"/>')
    p.append(f'<text class="tl" x="{xm + 5:.1f}" y="18">max. recommended {fnum(Vm, 1)}</text>')
    p.append(f'<text class="tl" x="{xe - 5:.1f}" y="30" text-anchor="end">erosion {fnum(Ve, 1)}</text>')
    p.append("</svg>")
    return "".join(p)


def _barra(X0, x1, y, BH, cls):
    """Horizontal bar with a rounded data end (4 px)."""
    if x1 - X0 > 8:
        d = (f"M{X0},{y:.1f}H{x1 - 4:.1f}Q{x1:.1f},{y:.1f} {x1:.1f},{y + 4:.1f}"
             f"V{y + BH - 4:.1f}Q{x1:.1f},{y + BH:.1f} {x1 - 4:.1f},{y + BH:.1f}H{X0}Z")
        return f'<path class="{cls}" d="{d}"/>'
    return (f'<rect class="{cls}" x="{X0}" y="{y:.1f}" width="{max(x1 - X0, 1):.1f}" '
            f'height="{BH}"/>')


def _grafico_perfil(pontos):
    """Pressure profile: pontos = [(name, km, bar_g)]. Single series."""
    W, X0, X1, TOP, BOT, H = 360, 46, 334, 26, 40, 200
    xs = [p[1] for p in pontos]
    ps = [p[2] for p in pontos]
    pmin, pmax = min(ps), max(ps)
    pad = (pmax - pmin) * 0.18 if pmax > pmin else 1.0
    lo, hi = pmin - pad, pmax + pad
    ticks, nd, _ = _ticks(lo, hi, 4)
    xmax = max(xs) or 1.0
    yb = H - BOT

    def X(v):
        return X0 + 10 + v / xmax * (X1 - X0 - 20)

    def Y(v):
        return TOP + (hi - v) / (hi - lo) * (yb - TOP)

    p = [f'<svg class="chart" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="Gauge pressure profile along the pipeline, from '
         f'{fnum(pmax, 1)} bar at point A to {fnum(pmin, 1)} bar at point F. '
         f'A table with the values is provided below.">']
    for t in ticks:
        y = Y(t)
        p.append(f'<line class="gl" x1="{X0}" y1="{y:.1f}" x2="{X1}" y2="{y:.1f}"/>')
        p.append(f'<text class="ax" x="{X0 - 8}" y="{y + 4:.1f}" text-anchor="end">{fnum(t, nd)}</text>')
    p.append(f'<line class="bl" x1="{X0}" y1="{yb}" x2="{X1}" y2="{yb}"/>')
    pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for _, a, b in pontos)
    p.append(f'<polyline class="ln" points="{pts}"/>')
    ult = len(pontos) - 1
    for i, (nome, km, pg) in enumerate(pontos):
        x, y = X(km), Y(pg)
        p.append(f'<circle class="dot c-teo" cx="{x:.1f}" cy="{y:.1f}" r="5"/>')
        anc = "start" if i == 0 else ("end" if i == ult else "middle")
        dx = -4 if i == 0 else (4 if i == ult else 0)
        p.append(f'<text class="vl" x="{x - dx:.1f}" y="{y - 11:.1f}" text-anchor="{anc}">{fnum(pg, 1)}</text>')
        p.append(f'<text class="rl" x="{x:.1f}" y="{yb + 17}" text-anchor="middle">{nome}</text>')
        ok_km = all(abs(X(km) - X(o)) > 34 for j, (_, o, _) in enumerate(pontos) if j != i)
        if ok_km:
            p.append(f'<text class="ax" x="{x:.1f}" y="{yb + 31}" text-anchor="middle">{fnum(km, 0)} km</text>')
        tip = f"Point {nome}: {fnum(pg, 1)} bar(g), {fnum(km, 1)} km from A"
        p.append(f'<rect class="hit" x="{x - 28:.1f}" y="{TOP - 10}" width="56" height="{yb - TOP + 40}" '
                 f'tabindex="0" data-tip="{_esc(tip)}"/>')
    p.append("</svg>")
    return "".join(p)


def _grafico_trechos(segs):
    """Velocity bars per segment with Vm and Ve marks on each row.
    segs = [(rotulo, V, Vers)]"""
    n = len(segs)
    W, X0, X1, RH, TOP, BOT, BH = 360, 62, 300, 38, 8, 30, 14
    H = TOP + n * RH + BOT
    vmax = max(max(s[1], s[2]) for s in segs)
    ticks, nd, step = _ticks(0, vmax * 1.06, 4)
    xmax = math.ceil(vmax * 1.06 / step - 1e-9) * step
    ticks = [t for t in ticks if t <= xmax + 1e-9]

    def X(v):
        return X0 + v / xmax * (X1 - X0)

    yb = TOP + n * RH
    p = [f'<svg class="chart" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="Velocity at the end of each segment, in meters per second, '
         f'with marks for the maximum recommended and the erosion velocity. '
         f'A table is provided below.">']
    for t in ticks:
        x = X(t)
        p.append(f'<line class="gl" x1="{x:.1f}" y1="{TOP}" x2="{x:.1f}" y2="{yb}"/>')
        p.append(f'<text class="ax" x="{x:.1f}" y="{yb + 19}" text-anchor="middle">{fnum(t, nd)}</text>')
    p.append(f'<line class="bl" x1="{X0}" y1="{TOP}" x2="{X0}" y2="{yb}"/>')
    for i, (rot, V, Ve) in enumerate(segs):
        y = TOP + i * RH + (RH - BH) / 2
        ym = y + BH / 2 + 4
        x1 = X(V)
        p.append(_barra(X0, x1, y, BH, "bar c-teo"))
        p.append(f'<text class="rl" x="{X0 - 8}" y="{ym:.1f}" text-anchor="end">{_esc(rot)}</text>')
        xm, xe = X(0.4 * Ve), X(Ve)
        y0, y1 = TOP + i * RH + 3, TOP + (i + 1) * RH - 3
        p.append(f'<line class="th warn" x1="{xm:.1f}" y1="{y0}" x2="{xm:.1f}" y2="{y1}"/>')
        p.append(f'<line class="th crit" x1="{xe:.1f}" y1="{y0}" x2="{xe:.1f}" y2="{y1}"/>')
        p.append(f'<text class="vl" x="{max(x1, xe) + 8:.1f}" y="{ym:.1f}">{fnum(V, 1)}</text>')
        tip = (f"Segment {rot}: V = {fnum(V, 1)} m/s · max. recommended "
               f"{fnum(0.4 * Ve, 1)} · erosion {fnum(Ve, 1)}")
        p.append(f'<rect class="hit" x="0" y="{TOP + i * RH}" width="{W}" height="{RH}" '
                 f'tabindex="0" data-tip="{_esc(tip)}"/>')
    p.append("</svg>")
    return "".join(p)


# ----------------------------------------------------------------------
# HTML pieces
# ----------------------------------------------------------------------
ICONES = {
    "ok": '<svg class="ic ok" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7"/><path d="M4.6 8.3l2.3 2.3 4.5-4.9"/></svg>',
    "warn": '<svg class="ic warn" viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.8l6.6 11.6H1.4z"/><path d="M8 6.4v3.2M8 11.4v.4"/></svg>',
    "crit": '<svg class="ic crit" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7"/><path d="M5.4 5.4l5.2 5.2M10.6 5.4l-5.2 5.2"/></svg>',
}


def situacao(V, Vm, Ve):
    if V <= Vm:
        return "ok", "Within recommended"
    if V <= Ve:
        return "warn", "Above recommended"
    return "crit", "Above erosion velocity"


def _campo(c, valor, erro):
    n = c["nome"]
    lo = "" if c["lo"] is None else repr(c["lo"])
    hi = "" if c["hi"] is None else repr(c["hi"])
    cls = "f" + (" largo" if c["largo"] else "") + (" bad" if erro else "")
    dica = f'<span class="dica">{_esc(c["dica"])}</span>' if c["dica"] else ""
    rel = ""
    if c["rel"]:
        rel = (f' data-rel="{c["rel"][0]}:{c["rel"][1]}" data-rel-sym="{_esc(c["sim"])}"')
    return (
        f'<div class="{cls}">'
        f'<label for="i-{n}"><span>{_esc(c["rotulo"])}</span>'
        f'<em>{_esc(faixa_txt(c))}</em></label>'
        f'<span class="in"><input id="i-{n}" name="{n}" type="text" inputmode="decimal" '
        f'autocomplete="off" value="{_esc(valor)}" data-sym="{_esc(c["sim"])}" '
        f'data-lo="{lo}" data-hi="{hi}" data-lo-ok="{int(c["lo_ok"])}" '
        f'data-hi-ok="{int(c["hi_ok"])}" data-faixa="{_esc(faixa_txt(c))}"{rel} '
        f'aria-describedby="m-{n}">'
        f'<span class="u">{_esc(c["unidade"])}</span></span>'
        f'{dica}<span class="msg" id="m-{n}" role="alert">{_esc(erro or "")}</span></div>'
    )


def _formulario(prog, brutos, erros):
    partes = [f'<form id="form" class="noprint" method="post" action="/{prog["slug"]}" novalidate>']
    for titulo, campos in prog["grupos"]:
        partes.append(f'<fieldset><legend>{_esc(titulo)}</legend><div class="grid">')
        for c in campos:
            partes.append(_campo(c, brutos.get(c["nome"], c["padrao"]), erros.get(c["nome"])))
        partes.append("</div></fieldset>")
    partes.append(f'<div class="acoes"><button type="submit" class="prim">{_esc(prog["botao"])}</button>'
                  f'<a class="sec" href="/{prog["slug"]}">Reset to example</a></div></form>')
    return "".join(partes)


def _avisos_html(avisos):
    if not avisos:
        return ""
    itens = "".join(f"<li>{_esc(a)}</li>" for a in avisos)
    return f'<div class="aviso">{ICONES["warn"]}<div><strong>Warning</strong><ul>{itens}</ul></div></div>'


def _secao(prog, brutos, r, hero, tiles, legenda, paineis, tit_cards, cartoes,
           tabelas):
    """Common skeleton of the results area."""
    t = "".join(f'<div><span>{_esc(a)}</span><b>{b}</b> {_esc(u)}</div>' for a, b, u in tiles)
    leg = ""
    if legenda:
        leg = f'<div class="legenda" aria-label="Legend">{legenda}</div>'
    pn = "".join(
        f'<figure class="painel"><figcaption>{_esc(a)}<small>{_esc(b)}</small></figcaption>{svg}</figure>'
        for a, b, svg in paineis)
    entradas = "".join(
        f"<div><dt>{_esc(c['rotulo'])} ({_esc(c['sim'])})</dt>"
        f"<dd>{_esc(brutos.get(c['nome'], ''))} {_esc(c['unidade'])}</dd></div>"
        for _, cs in prog["grupos"] for c in cs)
    return f"""
<section id="resultado" aria-labelledby="t-res">
  <h2 id="t-res" class="sr">Result</h2>
  {_avisos_html(r.get("avisos_pt"))}
  {hero}
  <div class="tiles">{t}</div>
  {leg}
  {pn}
  <h2 class="tit">{_esc(tit_cards)}</h2>
  <div class="cards">{cartoes}</div>
  <details class="tab"><summary>View as table</summary>
    <div class="rolar" tabindex="0" role="region" aria-label="Results table">{tabelas}</div>
  </details>
  <details class="tab" id="entradas"><summary>Input data</summary><dl class="kv2">{entradas}</dl></details>
  <div class="acoes noprint"><button type="button" class="sec" id="imprimir">Print or save as PDF</button></div>
  <p class="nota">Calculation aid. Results must be verified by a qualified
  professional before any use in a design. Theoretical basis: Stuckenbruck,
  <i>Escoamento em Dutos</i> (Pipe Flow), Vol. B (PUC-Rio, 2014).</p>
</section>"""


def _hero(rotulo, numero, unidade, sub):
    return (f'<div class="hero"><p class="rot">{_esc(rotulo)}</p>'
            f'<p class="num">{numero}<span> {_esc(unidade)}</span></p>'
            f'<p class="sub">{_esc(sub)}</p></div>')


def _legenda_modelos():
    return "".join(f'<span class="lg"><i class="sw c-{k}"></i>{CATEGORIAS[k]}</span>'
                   for k in ("teo", "emp", "aga"))


def _resultado_modelos(prog, brutos, r, modo):
    """Results of the diameter (modo='D') and flow-rate (modo='Q') programs."""
    m = r["modelos"]
    Vm, Ve = r["Vm"], r["Ve"]
    teo = m["THEORIC"]
    if modo == "D":
        val = lambda x: x["D"] * POL
        un, un_h = "in", "in"
        rot = "Inside diameter, theoretical model (Colebrook)"
        sub_u = "across the nine models: %s to %s in"
        sub_extra = f"{fnum(teo['D'] * 1000, 0)} mm"
        t1 = ("Inside diameter by model (in)", "The vertical line marks the theoretical model.")
        aria = "Inside diameter by model, in inches,"
        tip = lambda k, x: (f"{NOMES[k][1]}: {fnum(val(x))} in ({fnum(x['D'] * 1000, 0)} mm) · "
                            f"V = {fnum(x['V'], 1)} m/s")
        col = "ID (in)"
        nome_cards = "Results by model"
    else:
        big = max(x["Q"] for x in m.values()) * 86400 >= 1e6
        div, un = (1e6, "10⁶ Nm³/d") if big else (1e3, "10³ Nm³/d")
        ext = "million" if big else "thousand"
        val = lambda x: x["Q"] * 86400 / div
        un_h = un
        rot = "Standard flow rate, theoretical model (Colebrook)"
        sub_u = "across the nine models: %s to %s " + ext + " Nm³/d"
        sub_extra = f"{fmil(teo['Q'] * 86400)} Nm³/d"
        t1 = (f"Flow rate by model ({un})", "The vertical line marks the theoretical model.")
        aria = f"Standard flow rate by model, in {ext} Nm³/d,"
        tip = lambda k, x: (f"{NOMES[k][1]}: {fnum(val(x))} {ext} Nm³/d · "
                            f"V = {fnum(x['V'], 1)} m/s")
        col = f"Q ({un})"
        nome_cards = "Results by model"
    todos = [val(x) for x in m.values()]
    hero = _hero(rot, fnum(val(teo)), un_h,
                 sub_extra + " · " + sub_u % (fnum(min(todos)), fnum(max(todos))))
    tiles = [("Erosion velocity", fnum(Ve, 1), "m/s"),
             ("Max. recommended", fnum(Vm, 1), "m/s"),
             ("Avg. gradient", fnum(r["dPdx"], 1), "kPa/km")]
    cartoes, linhas_tab = [], []
    for chave, x in m.items():
        curto, longo, cat = NOMES[chave]
        st, txt = situacao(x["V"], Vm, Ve)
        f_txt = fsci(x["f"], 3) if x["f_sci"] else "—"
        if modo == "D":
            sec = f"{fnum(x['D'] * 1000, 0)} mm"
            sec_rot = "Diameter"
        else:
            sec = f"{fmil(x['Q'] * 86400)} Nm³/d"
            sec_rot = "Flow rate"
        cartoes.append(
            f'<article class="card"><header><i class="sw c-{cat}"></i>'
            f'<h3>{_esc(longo)}</h3><span class="cat">{CATEGORIAS[cat]}</span>'
            f'<div class="big">{fnum(val(x))}<small> {_esc(un)}</small></div></header>'
            f'<dl class="kv"><div><dt>{sec_rot}</dt><dd>{sec}</dd></div>'
            f'<div><dt>Velocity</dt><dd>{fnum(x["V"], 1)} m/s</dd></div>'
            f'<div><dt>Mach</dt><dd>{fnum(x["Ma"], 3)}</dd></div>'
            f'<div><dt>Reynolds</dt><dd>{fsci(x["Re"])}</dd></div>'
            f'<div><dt>Friction f</dt><dd>{f_txt}</dd></div></dl>'
            f'<p class="chip">{ICONES[st]}<span>{txt}</span></p></article>')
        linhas_tab.append(
            f"<tr><th scope='row'>{_esc(longo)}</th><td>{CATEGORIAS[cat]}</td>"
            f"<td>{fnum(val(x))}</td><td>{f_txt}</td>"
            f"<td>{fsci(x['Re'])}</td><td>{fnum(x['Ma'], 3)}</td><td>{fnum(x['V'], 1)}</td>"
            f"<td>{txt}</td></tr>")
    tabela = (f"<table><thead><tr><th>Model</th><th>Type</th><th>{col}</th><th>f</th>"
              f"<th>Re</th><th>Mach</th><th>V (m/s)</th><th>Status</th></tr></thead>"
              f"<tbody>{''.join(linhas_tab)}</tbody></table>")
    paineis = [(t1[0], t1[1], _grafico_pontos(m, val, aria, tip)),
               ("Outlet velocity (m/s)",
                "Compared with the maximum recommended and the erosion velocity.",
                _grafico_barras(m, Vm, Ve))]
    return _secao(prog, brutos, r, hero, tiles, _legenda_modelos(), paineis,
                  nome_cards, "".join(cartoes), tabela)


def _resultado_loop(prog, brutos, r):
    v = r["entrada"]
    PB, PE, PF = r["P_barg"]["B"], r["P_barg"]["E"], r["P_barg"]["F"]
    PA = v["P1"] - P_ATM
    Q = v["Q9"]
    qc, qd = r["Q_bce"] * 86400, r["Q_bde"] * 86400
    hero = _hero("Pressure at point F (outlet)", fnum(PF, 1), "bar(g)",
                 f"Total drop of {fnum(PA - PF, 1)} bar between A ({fnum(PA, 1)} bar(g)) and F")
    dv, du = (1e6, "10⁶") if Q >= 2e6 else (1e3, "10³")
    tiles = [("Branch C", fnum(qc / dv, 2), f"{du} Nm³/d"),
             ("Branch D", fnum(qd / dv, 2), f"{du} Nm³/d"),
             ("Split C / D", f"{fnum(100 * qc / Q, 0)} / {fnum(100 * qd / Q, 0)}", "%")]
    seg = [("ab", "Segment 1", "A → B", "L1", "D1", Q),
           ("bce", "Segment 2", "B → C → E", "L2", "D2", qc),
           ("bde", "Segment 3", "B → D → E", "L3", "D3", qd),
           ("ef", "Segment 4", "E → F", "L4", "D4", Q)]
    rot_c = {"ab": "A–B", "bce": "B–C–E", "bde": "B–D–E", "ef": "E–F"}
    segs_graf = [(rot_c[k], r["V"][k], r["Vers"][k]) for k, *_ in seg]
    perfil = [("A", 0.0, PA), ("B", v["L1"], PB), ("E", v["L1"] + v["L2"], PE),
              ("F", v["L1"] + v["L2"] + v["L4"], PF)]
    cartoes, tab = [], []
    for k, tit, trajeto, lk, dk, qs in seg:
        V, Ve = r["V"][k], r["Vers"][k]
        st, txt = situacao(V, 0.4 * Ve, Ve)
        cartoes.append(
            f'<article class="card"><header><i class="sw c-teo"></i>'
            f'<h3>{tit}</h3><span class="cat">{trajeto}</span>'
            f'<div class="big">{fnum(v[dk])}<small> in</small></div></header>'
            f'<dl class="kv"><div><dt>Length</dt><dd>{fnum(v[lk], 1)} km</dd></div>'
            f'<div><dt>Flow rate</dt><dd>{fmil(qs)} Nm³/d</dd></div>'
            f'<div><dt>Velocity</dt><dd>{fnum(V, 1)} m/s</dd></div>'
            f'<div><dt>Erosion limit</dt><dd>{fnum(Ve, 1)} m/s</dd></div>'
            f'<div><dt>Gradient</dt><dd>{fnum(r["dPdx"][k], 1)} kPa/km</dd></div>'
            f'<div><dt>Kt</dt><dd>{fsci(r["Kt"][k], 3)}</dd></div></dl>'
            f'<p class="chip">{ICONES[st]}<span>{txt}</span></p></article>')
        tab.append(
            f"<tr><th scope='row'>{tit} ({rot_c[k]})</th><td>{fnum(v[lk], 1)}</td>"
            f"<td>{fnum(v[dk])}</td><td>{fmil(qs)}</td><td>{fnum(V, 1)}</td>"
            f"<td>{fnum(Ve, 1)}</td><td>{fnum(r['dPdx'][k], 1)}</td><td>{txt}</td></tr>")
    tabela = (
        "<table><thead><tr><th>Point</th><th>Elevation (m)</th><th>P (bar g)</th></tr></thead><tbody>"
        + "".join(f"<tr><th scope='row'>{n}</th><td>{fnum(z, 0)}</td><td>{fnum(p, 1)}</td></tr>"
                  for n, z, p in (("A", v["Za"], PA), ("B", v["Zb"], PB),
                                  ("E", v["Ze"], PE), ("F", v["Zf"], PF)))
        + "</tbody></table><br>"
        "<table><thead><tr><th>Segment</th><th>L (km)</th><th>D (in)</th><th>Q (Nm³/d)</th>"
        "<th>V (m/s)</th><th>V erosion</th><th>dP/dx (kPa/km)</th><th>Status</th></tr></thead>"
        f"<tbody>{''.join(tab)}</tbody></table>")
    fech = r["fechamento_Q_bde"]
    if fech == fech:
        tabela += (f"<p class='dica'>Flow balance closure on branch D: {fnum(100 * fech, 3)} %"
                   " (difference between the flow obtained from Q − Q<sub>C</sub> and the direct calculation).</p>")
    legenda = ('<span class="lg"><i class="lm warn"></i>Max. recommended (0.4·Ve)</span>'
               '<span class="lg"><i class="lm crit"></i>Erosion velocity (Ve)</span>')
    paineis = [("Pressure profile (gauge, bar)",
                "Distance measured along branch B–C–E.", _grafico_perfil(perfil)),
               ("Velocity at the end of each segment (m/s)",
                "Marks per segment: amber = max. recommended, red = erosion.",
                _grafico_trechos(segs_graf))]
    return _secao(prog, brutos, r, hero, tiles, legenda, paineis,
                  "Results by segment", "".join(cartoes), tabela)


def _res_diam(prog, brutos, r):
    return _resultado_modelos(prog, brutos, r, "D")


def _res_vaz(prog, brutos, r):
    return _resultado_modelos(prog, brutos, r, "Q")


# ----------------------------------------------------------------------
# Available programs
# ----------------------------------------------------------------------
def _prog(slug, titulo, resumo, entra, sai, botao, grupos, calc, render,
          diagrama=""):
    return dict(slug=slug, diagrama=diagrama, titulo=titulo, resumo=resumo, entra=entra, sai=sai,
                botao=botao, grupos=grupos,
                campos=[c for _, cs in grupos for c in cs], calc=calc,
                render=render)


DIAGRAMA_LOOP = '''<svg viewBox="0 0 340 84" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Looped pipeline schematic: segment 1 from A to B, segments 2 and 3 in parallel from B to E, segment 4 from E to F.">
  <style>.p{stroke:#8FA8B5;stroke-width:2.5;fill:none;} .t{fill:#B0B0B0;font-size:10px;font-family:sans-serif;} .lbl{fill:#E64A19;font-size:11px;font-weight:600;font-family:sans-serif;}</style>
  <line class="p" x1="6" y1="55" x2="120" y2="55"/>
  <path class="p" d="M120,55 C 140,20 190,20 210,55"/>
  <line class="p" x1="120" y1="55" x2="210" y2="55"/>
  <line class="p" x1="210" y1="55" x2="290" y2="55"/>
  <circle cx="6" cy="55" r="3" fill="#E64A19"/>
  <circle cx="120" cy="55" r="3" fill="#E64A19"/>
  <circle cx="210" cy="55" r="3" fill="#E64A19"/>
  <circle cx="290" cy="55" r="3" fill="#8FA8B5"/>
  <path d="M290,49 L302,55 L290,61" fill="#8FA8B5"/>
  <text x="0" y="42" class="lbl">PA</text>
  <text x="2" y="70" class="t">A</text>
  <text x="65" y="42" class="lbl">D1, L1</text>
  <text x="113" y="70" class="t">B</text>
  <text x="140" y="15" class="lbl">D2, L2</text>
  <text x="140" y="70" class="lbl">D3, L3</text>
  <text x="205" y="70" class="t">E</text>
  <text x="235" y="42" class="lbl">D4, L4</text>
  <text x="283" y="70" class="t">F</text>
</svg>'''


PROGRAMAS = {
    p["slug"]: p for p in (
        _prog("diameter", "Gas pipeline diameter",
              "Compressible flow, nine models compared. Enter the data and tap "
              "calculate.",
              "standard flow rate", "inside diameter", "Calculate diameter",
              _grupos_diametro(), _calc_diametro, _res_diam),
        _prog("flow", "Gas pipeline flow rate",
              "Compressible flow, nine models compared. Enter the diameter and "
              "the pressures and tap calculate.",
              "inside diameter", "standard flow rate", "Calculate flow rate",
              _grupos_vazao(), _calc_vazao, _res_vaz),
        _prog("loop", "Looped gas pipeline",
              "Segments A–B, B–C–E / B–D–E in parallel, and E–F. Enter the "
              "network and tap calculate.",
              "network, diameters and total flow rate",
              "pressures, flow rates and velocities",
              "Calculate loop", _grupos_loop(), _calc_loop,
              lambda prog, brutos, r: _resultado_loop(prog, brutos, r),
              diagrama=DIAGRAMA_LOOP),
    )
}


def valores_padrao(prog):
    return {c["nome"]: c["padrao"] for c in prog["campos"]}


CSS = """
:root{--plane:#0b1118;--surface:#141b24;--raised:#1b2531;--line:#263140;--grid:#222d3b;
--axis:#3a4959;--ink:#f2f6fa;--ink2:#c3cdd8;--muted:#8b98a9;
--s1:#3987e5;--s2:#d95926;--s3:#199e70;--good:#0ca30c;--warn:#fab219;--crit:#d03b3b;
--sans:system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif}
*{box-sizing:border-box}
html{color-scheme:dark;-webkit-text-size-adjust:100%}
body{margin:0;background:var(--plane);color:var(--ink);font:16px/1.45 var(--sans)}
main{max-width:640px;margin:0 auto;padding:20px 16px 56px}
h1{font-size:26px;line-height:1.15;margin:8px 0 4px;font-weight:650;letter-spacing:-.01em}
.lead{color:var(--muted);margin:0 0 20px;font-size:15px}
.sr{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0)}
fieldset{border:1px solid var(--line);background:var(--surface);border-radius:14px;
margin:0 0 14px;padding:6px 14px 14px}
legend{padding:0 6px;font-weight:600;font-size:15px;color:var(--ink2)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(148px,1fr));gap:12px 12px}
.f.largo{grid-column:1/-1}
.f{display:flex;flex-direction:column;align-self:start}
.f label{display:flex;flex-direction:column;gap:1px;white-space:nowrap;
font-size:13px;color:var(--ink2);margin:6px 0 5px}
.f label em{font-style:normal;color:var(--muted);font-size:11.5px;white-space:nowrap}
.in{display:flex;align-items:center;background:var(--raised);border:1px solid var(--line);
border-radius:10px}
.in:focus-within{border-color:var(--s1);box-shadow:0 0 0 3px rgba(57,135,229,.28)}
.in input{flex:1;min-width:0;background:none;border:0;color:var(--ink);font:inherit;
font-size:17px;padding:11px 4px 11px 12px;outline:0;font-variant-numeric:tabular-nums}
.u{color:var(--muted);font-size:12.5px;padding:0 12px 0 4px;white-space:nowrap}
.dica{display:block;color:var(--muted);font-size:12px;margin-top:5px}
.msg{display:none;color:var(--ink);font-size:12.5px;margin-top:5px}
.f.bad .in{border-color:var(--crit)}
.f.bad .msg{display:block}
.f.bad .msg::before{content:"";display:inline-block;width:8px;height:8px;border-radius:50%;
background:var(--crit);margin-right:6px}
.acoes{display:flex;gap:10px;flex-wrap:wrap;margin:18px 0 8px}
button,.sec{font:inherit;font-weight:600;border-radius:12px;padding:13px 18px;cursor:pointer;
text-decoration:none;text-align:center}
.prim{flex:1 1 220px;background:var(--s1);color:#04101e;border:0;font-size:17px}
.prim:active{filter:brightness(.92)}
.sec{background:none;color:var(--ink2);border:1px solid var(--line)}
:focus-visible{outline:2px solid var(--s1);outline-offset:2px}
.banner{display:flex;gap:10px;background:var(--surface);border:1px solid var(--crit);
border-radius:12px;padding:12px 14px;margin:0 0 14px}
.aviso{display:flex;gap:10px;background:var(--surface);border:1px solid var(--line);
border-left:3px solid var(--warn);border-radius:12px;padding:12px 14px;margin:0 0 14px;font-size:14px}
.aviso ul{margin:4px 0 0;padding-left:18px;color:var(--ink2)}
.voltar{display:inline-block;color:var(--muted);text-decoration:none;font-size:14px;
padding:6px 0;margin:-4px 0 2px}
.voltar:hover{color:var(--ink2)}
.menu{display:grid;gap:12px;margin:8px 0 0}
.item{display:block;text-decoration:none;color:inherit;background:var(--surface);
border:1px solid var(--line);border-radius:14px;padding:16px 16px 14px}
.item:active{background:var(--raised)}
.item h2{margin:0 0 4px;font-size:19px;font-weight:620}
.item p{margin:0;color:var(--ink2);font-size:14.5px}
.item .es{margin-top:10px;color:var(--muted);font-size:13px}
.diagrama{background:var(--surface);border:1px solid var(--line);border-radius:14px;
padding:10px 12px;margin:6px 0 14px}
.diagrama svg{display:block;width:100%;max-width:520px;height:auto;margin:0 auto}
.lm{display:inline-block;width:3px;height:14px;border-radius:2px;flex:none}
.lm.warn{background:var(--warn)}.lm.crit{background:var(--crit)}
.ln{fill:none;stroke:var(--s1);stroke-width:2;stroke-linejoin:round;stroke-linecap:round}
.ic{width:18px;height:18px;flex:none;fill:none;stroke-width:1.6;stroke-linecap:round;stroke-linejoin:round}
.ic.ok{stroke:var(--good)}.ic.warn{stroke:var(--warn)}.ic.crit{stroke:var(--crit)}
#resultado{margin-top:30px;scroll-margin-top:12px}
.hero{padding:6px 2px 4px}
.hero .rot{margin:0;color:var(--ink2);font-size:15px}
.hero .num{margin:2px 0 0;font-size:clamp(52px,17vw,76px);line-height:1.02;font-weight:650;
letter-spacing:-.03em}
.hero .num span{font-size:.32em;font-weight:500;color:var(--ink2);letter-spacing:0}
.hero .sub{margin:6px 0 0;color:var(--muted);font-size:14.5px}
.tiles{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;margin:18px 0}
.tiles div{background:var(--surface);border:1px solid var(--line);border-radius:12px;
padding:10px 11px;font-size:12.5px;color:var(--ink2);line-height:1.3}
.tiles span{display:block;color:var(--muted);font-size:12px;margin-bottom:3px}
.tiles b{font-size:22px;color:var(--ink);font-weight:650}
.legenda{display:flex;gap:16px;flex-wrap:wrap;margin:6px 2px 10px;font-size:13.5px;color:var(--ink2)}
.lg{display:inline-flex;align-items:center;gap:7px}
.sw{display:inline-block;width:10px;height:10px;border-radius:50%;flex:none}
.sw.c-teo,.dot.c-teo,.bar.c-teo{background:var(--s1);fill:var(--s1)}
.sw.c-emp,.dot.c-emp,.bar.c-emp{background:var(--s2);fill:var(--s2)}
.sw.c-aga,.dot.c-aga,.bar.c-aga{background:var(--s3);fill:var(--s3)}
.painel{margin:0 0 14px;background:var(--surface);border:1px solid var(--line);
border-radius:14px;padding:14px 8px 8px}
.painel figcaption{padding:0 8px 6px;font-weight:600;font-size:15px}
.painel small{display:block;font-weight:400;color:var(--muted);font-size:12.5px;margin-top:2px}
.chart{display:block;width:100%;max-width:480px;margin:0 auto;overflow:visible}
.gl{stroke:var(--grid);stroke-width:1}
.bl{stroke:var(--axis);stroke-width:1}
.ref{stroke:var(--muted);stroke-width:1}
.dot{stroke:var(--surface);stroke-width:2}
.ax{fill:var(--muted);font-size:11px}
.rl{fill:var(--ink2);font-size:12px}
.vl{fill:var(--ink);font-size:11.5px;font-variant-numeric:tabular-nums}
.tl{fill:var(--ink2);font-size:11px}
.th{stroke-width:1.5}.th.warn{stroke:var(--warn)}.th.crit{stroke:var(--crit)}
.hit{fill:transparent;cursor:pointer}
.hit:focus-visible{outline:none;fill:rgba(57,135,229,.12)}
.tit{font-size:18px;margin:24px 2px 10px}
.cards{display:grid;gap:10px}
.card{background:var(--surface);border:1px solid var(--line);border-radius:14px;padding:12px 14px}
.card header{display:grid;grid-template-columns:auto 1fr auto;column-gap:9px;align-items:center}
.card h3{margin:0;font-size:16px;font-weight:600}
.card .cat{grid-column:2;color:var(--muted);font-size:12.5px}
.card .big{grid-column:3;grid-row:1/3;font-size:28px;font-weight:650;letter-spacing:-.02em;text-align:right}
.card .big small{font-size:13px;color:var(--ink2);font-weight:500}
.card .sw{grid-row:1/3}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(92px,1fr));gap:8px;margin:10px 0 8px}
.kv dt{color:var(--muted);font-size:12px}.kv dd{margin:1px 0 0;font-size:14.5px;font-variant-numeric:tabular-nums}
.chip{display:flex;align-items:center;gap:8px;margin:0;font-size:13.5px;color:var(--ink2)}
.tab{margin:14px 0 0;background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:0 14px}
.tab summary{padding:13px 0;cursor:pointer;font-weight:600;font-size:15px}
.rolar{overflow-x:auto;padding-bottom:12px}
table{border-collapse:collapse;font-size:13px;min-width:560px;width:100%}
th,td{padding:7px 8px;text-align:center;border-bottom:1px solid var(--line);white-space:nowrap;
font-variant-numeric:tabular-nums}
th:first-child,td:first-child,th[scope=row]{text-align:left}
thead th{color:var(--muted);font-weight:500}
.kv2{margin:0 0 12px;display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:8px 16px}
.kv2 dt{color:var(--muted);font-size:12.5px}.kv2 dd{margin:0;font-size:14.5px}
.nota{color:var(--muted);font-size:12.5px;margin:22px 2px 0}
#tip{position:fixed;z-index:9;pointer-events:none;background:#e9eef4;color:#0b1118;
border-radius:8px;padding:7px 10px;font-size:12.5px;max-width:260px;display:none;
box-shadow:0 4px 14px rgba(0,0,0,.4)}
@media (max-width:380px){.tiles{grid-template-columns:1fr 1fr}.tiles div:last-child{grid-column:1/-1}}
@media print{
:root{--plane:#fff;--surface:#fff;--raised:#f3f4f6;--line:#d5d9df;--grid:#e5e7eb;--axis:#c3c9d1;
--ink:#0b0b0b;--ink2:#3a4350;--muted:#5b6674}
.noprint,#tip{display:none!important}
body{background:#fff}.card,.painel{break-inside:avoid}
}
"""

JS = r"""
(function(){
  function num(s){
    s=String(s).replace(/\s/g,'');
    if(!s)return NaN;
    var c=s.indexOf(',')>=0,d=s.indexOf('.')>=0;
    if(c&&d){
      if(s.lastIndexOf(',')>s.lastIndexOf('.'))s=s.replace(/\./g,'').replace(',','.');
      else s=s.replace(/,/g,'');
    }else if(c){
      var p=s.replace(/^[+-]/,'').split(',');
      if(p.length>2||(p.length===2&&p[1].length===3&&p[0].length>=1&&p[0].length<=3&&p[0].charAt(0)!=='0'))s=s.replace(/,/g,'');
      else s=s.replace(',','.');
    }
    else if((s.match(/\./g)||[]).length>1)s=s.replace(/\./g,'');
    return /^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/i.test(s)?parseFloat(s):NaN;
  }
  var form=document.getElementById('form');
  if(form){
    var campos=form.querySelectorAll('input[data-sym]');
    var msgDe=function(inp){
      var v=num(inp.value);
      if(isNaN(v))return 'Enter a number.';
      var lo=inp.getAttribute('data-lo'),hi=inp.getAttribute('data-hi');
      var lok=inp.getAttribute('data-lo-ok')==='1',hok=inp.getAttribute('data-hi-ok')==='1';
      if(lo!==''&&(lok?v<+lo:v<=+lo))return 'Out of range: '+inp.getAttribute('data-faixa');
      if(hi!==''&&(hok?v>+hi:v>=+hi))return 'Out of range: '+inp.getAttribute('data-faixa');
      var rel=inp.getAttribute('data-rel');
      if(rel){var p=rel.split(':'),o=form.elements[p[1]],ov=o?num(o.value):NaN;
        if(!isNaN(ov)){
          var os=o.getAttribute('data-sym'),me=inp.getAttribute('data-rel-sym');
          if(p[0]==='<'&&v>=ov)return me+' must be lower than '+os+'.';
          if(p[0]==='<='&&v>ov)return me+' must be lower than or equal to '+os+'.';
        }}
      return '';
    };
    var marcar=function(inp){
      var m=msgDe(inp),box=inp.parentNode.parentNode;
      box.className=box.className.replace(/\s*\bbad\b/g,'')+(m?' bad':'');
      box.querySelector('.msg').textContent=m;
      return !m;
    };
    Array.prototype.forEach.call(campos,function(inp){
      inp.addEventListener('input',function(){marcar(inp);
        Array.prototype.forEach.call(campos,function(o){
          var r=o.getAttribute('data-rel');
          if(r&&r.split(':')[1]===inp.name)marcar(o);});});
    });
    form.addEventListener('submit',function(e){
      var ok=true,first=null;
      Array.prototype.forEach.call(campos,function(inp){
        if(!marcar(inp)){ok=false;first=first||inp;}
      });
      if(!ok){e.preventDefault();first.focus();}
    });
  }
  var tip=document.getElementById('tip'),timer=null;
  function mostrar(el,x,y){
    tip.textContent=el.getAttribute('data-tip');tip.style.display='block';
    var w=tip.offsetWidth,h=tip.offsetHeight;
    var left=Math.max(8,Math.min(x-w/2,window.innerWidth-w-8));
    var top=y-h-16;if(top<8)top=y+20;
    tip.style.left=left+'px';tip.style.top=top+'px';
    clearTimeout(timer);timer=setTimeout(function(){tip.style.display='none';},4500);
  }
  function alvo(e){return e.target.closest?e.target.closest('[data-tip]'):null;}
  document.addEventListener('pointerdown',function(e){var a=alvo(e);if(a)mostrar(a,e.clientX,e.clientY);else tip.style.display='none';});
  document.addEventListener('pointerover',function(e){if(e.pointerType==='mouse'){var a=alvo(e);if(a)mostrar(a,e.clientX,e.clientY);}});
  document.addEventListener('pointerout',function(e){if(e.pointerType==='mouse'&&alvo(e))tip.style.display='none';});
  document.addEventListener('focusin',function(e){var a=alvo(e);if(a){var r=a.getBoundingClientRect();mostrar(a,r.left+r.width/2,r.top);}});
  window.addEventListener('scroll',function(){tip.style.display='none';},{passive:true});
  var imp=document.getElementById('imprimir');
  if(imp)imp.addEventListener('click',function(){window.print();});
  window.addEventListener('beforeprint',function(){
    Array.prototype.forEach.call(document.querySelectorAll('details'),function(d){d.setAttribute('open','');});
  });
  var res=document.getElementById('resultado');
  if(res){
    var reduz=window.matchMedia&&window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    window.addEventListener('load',function(){res.scrollIntoView({behavior:reduz?'auto':'smooth',block:'start'});});
  }
})();
"""


def _cabecalho(titulo):
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="dark">
<meta name="theme-color" content="#0b1118">
<link rel="icon" href="data:,">
<title>{_esc(titulo)}</title>
<style>{CSS}</style>
</head>
<body>
<main>
"""


RODAPE = """
</main>
<div id="tip" role="status"></div>
<script>{JS}</script>
</body>
</html>"""


def render_page(prog, brutos=None, resultado=None, erro="", erros_campo=None):
    brutos = brutos or valores_padrao(prog)
    erros_campo = erros_campo or {}
    banner = ""
    if erro:
        banner = (f'<div class="banner" role="alert">{ICONES["crit"]}'
                  f'<div>{_esc(erro)}</div></div>')
    corpo = prog["render"](prog, brutos, resultado) if resultado else ""
    return (_cabecalho(prog["titulo"])
            + '<a class="voltar noprint" href="/">‹ Menu</a>\n'
            + f'<h1>{_esc(prog["titulo"])}</h1>\n'
            + (f'<div class="diagrama">{prog["diagrama"]}</div>\n' if prog["diagrama"] else "")
            + f'<p class="lead noprint">{_esc(prog["resumo"])}</p>\n'
            + banner + _formulario(prog, brutos, erros_campo) + corpo
            + RODAPE.replace("{JS}", JS))


def render_menu():
    itens = "".join(
        f'<a class="item" href="/{p["slug"]}"><h2>{_esc(p["titulo"])}</h2>'
        f'<p>Input: {_esc(p["entra"])}.</p><p>Output: {_esc(p["sai"])}.</p></a>'
        for p in PROGRAMAS.values())
    return (_cabecalho("Gas pipeline calculations")
            + '<h1>Gas pipeline calculations</h1>\n'
            + '<p class="lead">Choose a program. All calculations run on this device.</p>\n'
            + f'<nav class="menu" aria-label="Programs">{itens}</nav>\n'
            + '<p class="nota">Calculation author: Antonio Ricardo Andrade Bozolla. '
              'Theoretical basis: Stuckenbruck, <i>Escoamento em Dutos</i> (Pipe Flow), Vol. B (PUC-Rio, 2014).</p>'
            + RODAPE.replace("{JS}", JS))


# ----------------------------------------------------------------------
# Server
# ----------------------------------------------------------------------
CSP = ("default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; "
       "img-src data:; form-action 'self'; base-uri 'none'; frame-ancestors 'none'")


class Handler(BaseHTTPRequestHandler):
    server_version = "GasApp/2.0"

    def _enviar(self, codigo, texto):
        dados = texto.encode("utf-8")
        self.send_response(codigo)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(dados)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", CSP)
        self.end_headers()
        self.wfile.write(dados)

    def _host_ok(self):
        host = (self.headers.get("Host") or "").rsplit(":", 1)[0]
        return host in ("127.0.0.1", "localhost")

    def _rota(self):
        """Returns (ok, slug). slug '' = menu."""
        if not self._host_ok():
            self._enviar(403, "Access denied.")
            return False, None
        slug = self.path.split("?")[0].strip("/")
        if slug == "" or slug in PROGRAMAS:
            return True, slug
        self._enviar(404, "Page not found.")
        return False, None

    def do_GET(self):
        ok, slug = self._rota()
        if not ok:
            return
        self._enviar(200, render_page(PROGRAMAS[slug]) if slug else render_menu())

    def do_POST(self):
        from urllib.parse import parse_qs
        ok, slug = self._rota()
        if not ok:
            return
        if not slug:
            return self._enviar(400, "Invalid request.")
        prog = PROGRAMAS[slug]
        try:
            tam = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            tam = 0
        if tam <= 0 or tam > 20000:
            return self._enviar(400, "Invalid request.")
        corpo = self.rfile.read(tam).decode("utf-8", "replace")
        campos = parse_qs(corpo, keep_blank_values=True)
        brutos = {c["nome"]: (campos.get(c["nome"], [""])[0])[:40]
                  for c in prog["campos"]}
        try:
            r, erro, erros = calcular(prog, brutos)
            pagina = render_page(prog, brutos, r, erro, erros)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return self._enviar(500, "Internal error: %s: %s" % (type(e).__name__, e))
        self._enviar(200, pagina)

    def log_message(self, formato, *args):
        print("[%s] %s" % (self.log_date_time_string(), formato % args))

    def handle_one_request(self):
        try:
            super().handle_one_request()
        except Exception:
            import traceback
            traceback.print_exc()


def main():
    servidor, porta = None, PORTA
    for porta in range(PORTA, PORTA + 11):
        try:
            servidor = ThreadingHTTPServer((HOST, porta), Handler)
            servidor.daemon_threads = True
            break
        except OSError:
            servidor = None
    if servidor is None:
        print("Could not open any port between %d and %d." % (PORTA, PORTA + 10))
        return
    url = f"http://{HOST}:{porta}"
    if porta != PORTA:
        print("NOTE: port %d is busy (is another run still active?)." % PORTA)
    print("Server running. Open this address in your browser:")
    print("  " + url)
    print("To quit, stop the program (in Pydroid: stop button in the editor).")
    try:
        webbrowser.open(url)
    except Exception:
        pass
    try:
        servidor.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        servidor.server_close()


if __name__ == "__main__":
    main()

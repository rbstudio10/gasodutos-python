# Gas pipelines – diameter, flow rate and loop (Python)

Python 3 translation of the HP Prime (HP PPL) programs `GasDiameter`, `GasFlowRate` and
`GasPipeLoop`, by Antonio Ricardo Andrade Bozolla, plus a local web interface (dark theme)
that runs on a phone or a PC.

Theoretical basis: STUCKENBRUCK, S. *Escoamento em Dutos* (Pipe Flow), Volume B. PUC-Rio, 2014.

## Files

| File | Purpose |
|---|---|
| `gasutil.py` | Constants, validations, Colebrook, empirical models and output formatting |
| `gas_diameter.py` | Diameter from the standard flow rate (Theoretical, Weymouth, Panhandle A/B, IGT, Mueller, Fritzsche, AGA-A, AGA-B) |
| `gas_flow_rate.py` | Standard flow rate from the diameter (same models) |
| `gas_pipe_loop.py` | Looped gas pipeline (segments A–B, B–C–E, B–D–E, E–F) |
| `web_app.py` | Local web interface (menu + the three programs above) |
| `test_consistencia.py` | Internal consistency tests |

Standard library only (`math`, `http.server`, …). Keep all files in the same folder.

## Quick start (script)

```python
from gas_diameter import gas_diameter, formatar
r = gas_diameter(L9=100, Q9=3.0e6, eps=0.046, lam=0.6, mu=1.1e-5, k=1.3,
                 P19=70, P29=40, Zm=0.9, Tm9=25, n=0.92, Ca=0.95, h=0)
print(formatar(r))
```

Units are the same as in the PPL programs: L in km; Q in Nm³/d; ε in mm; P in bar (absolute);
T in °C (`gas_diameter`, `gas_flow_rate`) or K (`gas_pipe_loop`); D in m (`gas_flow_rate`)
or in (`gas_pipe_loop`).

The function and parameter names keep the original (Portuguese/PPL) naming on purpose.

## Web interface

```
python web_app.py
```

Then open **http://127.0.0.1:8000** in a browser (on Android with Pydroid 3, use Chrome on the same
phone). The page has a menu with three programs:

- **Gas pipeline diameter** – nine models compared, with charts and per-model cards.
- **Gas pipeline flow rate** – same layout, diameter in inches.
- **Looped gas pipeline** – pressure profile, flow split between the loop branches and velocities
  against the erosion limit, with a schematic of the network.

Input fields are validated while you type, results can be printed or saved as PDF, and all
calculations run locally; the server only accepts connections from the device itself (127.0.0.1).
Numbers can be typed as `3000000`, `3,000,000` or `0.046` (a decimal comma such as `1,5` also works).
If port 8000 is busy, the program uses the next free port and prints the address. To quit, stop
the program (in Pydroid, go back to the editor and tap the stop button).

## Diagnostics

`estrito=True` (default) aborts the calculation with `GasPipelineError` on roughness, Reynolds
and maximum-diameter warnings (equivalent to the PPL `BREAK`). `estrito=False` records the warning
in `r["avisos"]` and carries on (this is what the web interface uses).

## Tests

```
python test_consistencia.py
```

## License

MIT – see `LICENSE`.

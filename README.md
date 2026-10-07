# Gas pipelines: diameter, flow rate and loop calculations (Python)

Python 3 translation of the HP Prime (HP PPL) programs `GasDiameter`, `GasFlowRate` and
`GasPipeLoop`, by Antonio Ricardo Andrade Bozolla.

Theoretical basis: STUCKENBRUCK, S. *Escoamento em Dutos*, Volume B. PUC-Rio, 2014.

## Files

| File | Purpose |
|---|---|
| `gasutil.py` | Constants, input validation, Colebrook, empirical models and output formatting |
| `gas_diameter.py` | Diameter from standard flow rate (Theoretical, Weymouth, Panhandle A/B, IGT, Mueller, Fritzsche, AGA-A, AGA-B) |
| `gas_flow_rate.py` | Standard flow rate from diameter (same models) |
| `gas_pipe_loop.py` | Pipeline with a loop (segments AB, B–C–E, B–D–E, EF) |
| `test_consistencia.py` | Internal consistency tests |

Standard library only (`math`). All files must be kept in the same folder.

## Usage

```python
from gas_diameter import gas_diameter, formatar
r = gas_diameter(L9=100, Q9=3.0e6, eps=0.046, lam=0.6, mu=1.1e-5, k=1.3,
                 P19=70, P29=40, Zm=0.9, Tm9=25, n=0.92, Ca=0.95, h=0)
print(formatar(r))
```

Units are the same as in the PPL programs: length in km; flow rate in Nm³/d;
roughness ε in mm; pressures in bar (absolute); temperature in °C
(`gas_diameter`, `gas_flow_rate`) or K (`gas_pipe_loop`); diameter in m
(`gas_flow_rate`) or inches (`gas_pipe_loop`).

To run on Android, open the file you want in Pydroid 3, edit the example values at
the bottom (inside the `if __name__ == "__main__":` block) and tap ▶.

## Diagnostics

With `estrito=True` (default), the calculation stops with a `GasPipelineError` on
roughness, Reynolds number and maximum-diameter warnings (equivalent to `BREAK` in
the PPL code). With `estrito=False`, the warning is stored in `r["avisos"]` and the
calculation continues.

## Tests

```
python test_consistencia.py
```

## License

MIT. See the `LICENSE` file.

## Disclaimer

These are calculation-support tools. Results must be verified by a qualified
engineer before any use in a design.

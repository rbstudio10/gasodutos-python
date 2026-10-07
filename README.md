# Gasodutos – cálculo de diâmetro, vazão e loop (Python)

Tradução para Python 3 dos programas HP Prime (HP PPL) `GasDiameter`, `GasFlowRate` e
`GasPipeLoop`, de Antonio Ricardo Andrade Bozolla.

Base teórica: STUCKENBRUCK, S. *Escoamento em Dutos*, Volume B. PUC-Rio, 2014.

## Arquivos

| Arquivo | Função |
|---|---|
| `gasutil.py` | Constantes, validações, Colebrook, modelos empíricos e formatação |
| `gas_diameter.py` | Diâmetro a partir da vazão padrão (Teórico, Weymouth, Panhandle A/B, IGT, Mueller, Fritzsche, AGA-A, AGA-B) |
| `gas_flow_rate.py` | Vazão padrão a partir do diâmetro (mesmos modelos) |
| `gas_pipe_loop.py` | Gasoduto com loop (trechos AB, B–C–E, B–D–E, EF) |
| `test_consistencia.py` | Testes de consistência interna |

Somente biblioteca padrão (`math`). Todos os arquivos devem ficar na mesma pasta.

## Uso

```python
from gas_diameter import gas_diameter, formatar
r = gas_diameter(L9=100, Q9=3.0e6, eps=0.046, lam=0.6, mu=1.1e-5, k=1.3,
                 P19=70, P29=40, Zm=0.9, Tm9=25, n=0.92, Ca=0.95, h=0)
print(formatar(r))
```

Unidades iguais às do PPL: L em km; Q em Nm³/d; ε em mm; P em bar (absoluta);
T em °C (`gas_diameter`, `gas_flow_rate`) ou K (`gas_pipe_loop`); D em m
(`gas_flow_rate`) ou pol (`gas_pipe_loop`).

## Diagnósticos

`estrito=True` (padrão) interrompe o cálculo com `GasPipelineError` nos avisos de
rugosidade, Reynolds e diâmetro máximo (equivalente ao `BREAK` do PPL).
`estrito=False` registra o aviso em `r["avisos"]` e segue.

## Testes

```
python test_consistencia.py
```

## Licença

Defina a licença antes de publicar (ex.: MIT).

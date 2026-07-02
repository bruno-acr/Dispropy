# dispropy

`dispropy` calcula métricas de desproporcionalidade em farmacovigilância a
partir de tabelas 2x2 armazenadas em um `pandas.DataFrame`.

|                         | Evento de interesse | Outros eventos |
|-------------------------|---------------------|----------------|
| Medicamento de interesse | A                  | B              |
| Outros medicamentos      | C                  | D              |

## Métricas

- **ROR:** `(A * D) / (B * C)`, com intervalo de confiança de 95%.
- **PRR:** `[A / (A + B)] / [C / (C + D)]`, com intervalo de confiança de 95%.
- **IC:** `log2((Obs + 0.5) / (Esp + 0.5))`, onde `Obs = A` e
  `Esp = ((A+B) * (A+C)) / (A+B+C+D)`. Inclui IC025 e IC975.
- **EBGM:** média geométrica bayesiana empírica do método GPS de DuMouchel.
  Uma mistura de duas distribuições Gama é ajustada conjuntamente a todos os
  pares. A saída inclui o peso posterior `qn`, `ebgm`, `eb05` e `eb95`.

O IC025 é calculado por
`IC - 3.3 * (Obs + 0.5)^(-0.5) - 2 * (Obs + 0.5)^(-1.5)`.

## Instalação

```bash
pip install -e .
```

Para desenvolvimento:

```bash
pip install -e ".[dev]"
pytest
```

## Uso

Um tutorial executável cobrindo todas as métricas, validações, flags e
diagnósticos GPS está disponível em
[`examples/dispropy_tutorial.ipynb`](examples/dispropy_tutorial.ipynb).

```python
import pandas as pd
from dispropy import calculate_disproportionality

df = pd.DataFrame({
    "drug": ["Drug A", "Drug A", "Drug A", "Drug A"],
    "event": ["Event W", "Event X", "Event Y", "Event Z"],
    "A": [2, 5, 10, 20],
    "B": [98, 95, 90, 80],
    "C": [20, 25, 20, 30],
    "D": [880, 875, 880, 870],
})

resultado = calculate_disproportionality(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
    metrics=("ror", "prr", "ic", "ebgm"),
    correction=0.5,
    shrinkage=0.5,
    add_signal_flags=True,
)

print(resultado)
print(resultado.attrs["gps_model"])
```

ROR, PRR e IC são calculados por padrão. O EBGM precisa ser solicitado porque
envolve ajuste numérico global e exige pelo menos dois pares com contagem
esperada positiva. Os parâmetros GPS ajustados ficam em
`resultado.attrs["gps_model"]`.

Também é possível chamar `calculate_ror`, `calculate_prr`, `calculate_ic` e
`calculate_ebgm` individualmente. Os nomes das colunas A, B, C e D são livres.

## Sinalização e interpretação

Com `add_signal_flags=True`, os critérios de triagem são:

- `signal_ror`: `ror_lower_95 > 1`;
- `signal_prr`: `prr_lower_95 > 1`;
- `signal_ic`: `ic025 > 0`;
- `signal_ebgm`: `eb05 > 2` e pelo menos três observações (`A >= 3`).

ROR ou PRR acima de 1 e IC acima de 0 sugerem relato maior que o comparador ou
que o esperado. Essas métricas e flags indicam apenas desproporcionalidade
estatística. Elas não demonstram causalidade, e a avaliação clínica e
farmacológica continua necessária.

## Limitações

- Não corrige vieses inerentes à notificação espontânea.
- Não estima incidência nem risco absoluto.
- O ajuste GPS depende da quantidade e da composição dos pares analisados.
- Não substitui avaliação clínica.
- Não implementa estratificação por idade, sexo, país ou período.

## Quando os resultados podem não ser confiáveis

A validação da biblioteca confirma que as colunas existem, são numéricas,
não contêm valores ausentes ou negativos e não representam uma linha totalmente
vazia. Isso não comprova que `A`, `B`, `C` e `D` foram construídos corretamente.
O usuário deve verificar a unidade de contagem, deduplicação, definição de
medicamento e evento, população comparadora e consistência dos totais.

ROR, PRR e IC podem ser estatisticamente instáveis quando as contagens são
pequenas. A correção de continuidade e o shrinkage evitam operações indefinidas
e reduzem parte da instabilidade, mas não criam informação nem substituem a
avaliação da precisão e da relevância clínica.

O EBGM ajusta cinco hiperparâmetros usando todos os pares válidos. A biblioteca
emite `GPSFitWarning` quando há menos de 50 pares válidos ou quando algum
parâmetro termina a menos de 1% de um bound do otimizador. O limite de 50 é um
diagnóstico operacional conservador, equivalente a 10 pares por hiperparâmetro,
e não um corte formal estabelecido na literatura. O GPS foi desenvolvido para
grandes tabelas de frequência, e implementações de referência também verificam
convergência dentro do espaço de parâmetros e estabilidade entre soluções.
Consulte [DuMouchel (1999)](https://doi.org/10.1080/00031305.1999.10474456) e
[Canida e Ihrie (2017)](https://journal.r-project.org/articles/RJ-2017-063/).

Quando um warning aparecer, examine
`resultado.attrs["gps_model"]`. Os campos `parameters_near_bounds` e
`near_bound_parameters` indicam proximidade dos limites; `valid_pair_count`
informa quantos pares sustentaram o ajuste. Não interprete o EBGM isoladamente:
revise a construção da tabela, aumente e diversifique o conjunto de pares
quando possível e realize análise de sensibilidade ou validação estatística
independente antes de usar o resultado para decisão.

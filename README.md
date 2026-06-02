# dispropy

`dispropy` e uma biblioteca Python simples para calculo vetorizado de metricas
de desproporcionalidade em farmacovigilancia a partir de tabelas 2x2.

Ela recebe um `pandas.DataFrame` com colunas correspondentes a `A`, `B`, `C` e
`D` e calcula:

- Reporting Odds Ratio (`ROR`)
- Proportional Reporting Ratio (`PRR`)
- Information Component (`IC`)
- Limites `IC025` e `IC975`

A biblioteca foi feita para bases ja preparadas de farmacovigilancia. Ou seja:
ela calcula as metricas a partir das contagens `A`, `B`, `C` e `D`, mas nao
monta a tabela 2x2 automaticamente a partir de uma base bruta de notificacoes.

## Tabela 2x2

Para cada par medicamento-evento:

```text
                                Evento de interesse      Outros eventos
Medicamento de interesse                 A                     B
Outros medicamentos                       C                     D
```

Onde:

- `A`: notificacoes com o medicamento de interesse e o evento de interesse.
- `B`: notificacoes com o medicamento de interesse e outros eventos.
- `C`: notificacoes com outros medicamentos e o evento de interesse.
- `D`: notificacoes com outros medicamentos e outros eventos.

Os nomes das colunas podem ser livres. A planilha pode ter colunas chamadas
`A`, `B`, `C`, `D`, ou nomes mais descritivos como
`drug_event_count`, `drug_other_events`, `other_drugs_event` e
`other_drugs_other_events`.

## Instalacao

### Pelo GitHub

Depois que o repositório estiver publicado no GitHub, a biblioteca pode ser
instalada em qualquer computador com:

```bash
python3 -m pip install git+https://github.com/SEU-USUARIO/dispropy.git
```

Substitua `SEU-USUARIO` pelo nome da sua conta no GitHub.

### Localmente

Entre na pasta do projeto:

```bash
cd /Users/bruno/Documents/Developer/dispropy
```

Instale a biblioteca em modo editavel:

```bash
python3 -m pip install -e .
```

Para instalar tambem as dependencias de desenvolvimento:

```bash
python3 -m pip install -e ".[dev]"
```

### Futuramente pelo PyPI

Quando a biblioteca estiver publicada no PyPI, a instalacao sera:

```bash
python3 -m pip install dispropy
```

## Como importar

O jeito recomendado e importar como `disp`:

```python
import dispropy as disp
```

Com isso, voce pode usar:

```python
disp.ror(...)
disp.prr(...)
disp.ic(...)
disp.disproportionality(...)
```

Tambem e possivel importar as funcoes completas:

```python
from dispropy import (
    calculate_ror,
    calculate_prr,
    calculate_ic,
    calculate_disproportionality,
)
```

## Exemplo basico

```python
import pandas as pd
import dispropy as disp

df = pd.DataFrame({
    "drug": ["Drug A", "Drug B"],
    "event": ["Event X", "Event Y"],
    "A": [10, 5],
    "B": [90, 45],
    "C": [20, 10],
    "D": [880, 940],
})

resultado = disp.disproportionality(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
    correction=0.5,
    shrinkage=0.5,
    add_signal_flags=True,
)

print(resultado)
```

O resultado mantem as colunas originais e adiciona as metricas calculadas.

## Calcular metricas separadamente

### ROR

```python
resultado_ror = disp.ror(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
)

print(resultado_ror[[
    "ror",
    "log_ror",
    "se_log_ror",
    "ror_lower_95",
    "ror_upper_95",
]])
```

### PRR

```python
resultado_prr = disp.prr(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
)

print(resultado_prr[[
    "prr",
    "log_prr",
    "se_log_prr",
    "prr_lower_95",
    "prr_upper_95",
]])
```

### IC

```python
resultado_ic = disp.ic(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
)

print(resultado_ic[[
    "observed_count",
    "expected_count",
    "ic",
    "ic025",
    "ic975",
]])
```

## Calcular todas as metricas

Use `disp.disproportionality` para calcular tudo de uma vez:

```python
resultado = disp.disproportionality(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
    add_signal_flags=True,
)
```

Por padrao, essa funcao calcula:

- `ror`
- `prr`
- `ic`

E adiciona tambem as colunas auxiliares:

- `log_ror`
- `se_log_ror`
- `ror_lower_95`
- `ror_upper_95`
- `log_prr`
- `se_log_prr`
- `prr_lower_95`
- `prr_upper_95`
- `observed_count`
- `expected_count`
- `ic025`
- `ic975`

Quando `add_signal_flags=True`, tambem adiciona:

- `signal_ror`: `ror_lower_95 > 1`
- `signal_prr`: `prr_lower_95 > 1`
- `signal_ic`: `ic025 > 0`

Essas flags indicam apenas desproporcionalidade estatistica para triagem. Elas
nao indicam causalidade.

## Escolher quais metricas calcular

Voce pode calcular apenas uma ou duas metricas usando o argumento `metrics`.

Somente ROR:

```python
resultado = disp.disproportionality(
    df,
    "A",
    "B",
    "C",
    "D",
    metrics=("ror",),
)
```

ROR e IC:

```python
resultado = disp.disproportionality(
    df,
    "A",
    "B",
    "C",
    "D",
    metrics=("ror", "ic"),
)
```

Valores validos:

- `"ror"`
- `"prr"`
- `"ic"`

## Usar com uma planilha Excel

Se voce tem uma planilha com muitos pares medicamento-evento, carregue a
planilha com pandas:

```python
import pandas as pd
import dispropy as disp

df = pd.read_excel("minha_planilha.xlsx")

resultado = disp.disproportionality(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
    add_signal_flags=True,
)

resultado.to_excel("resultado_disproporcionalidade.xlsx", index=False)
```

Exemplo de planilha de entrada:

```text
drug       event      A    B    C    D
Drug A     Event X    10   90   20   880
Drug B     Event Y    5    45   10   940
Drug C     Event Z    0    12   4    1200
```

## Usar com CSV

```python
import pandas as pd
import dispropy as disp

df = pd.read_csv("minha_base.csv")

resultado = disp.disproportionality(
    df,
    a_col="A",
    b_col="B",
    c_col="C",
    d_col="D",
    add_signal_flags=True,
)

resultado.to_csv("resultado_disproporcionalidade.csv", index=False)
```

## Usar com nomes de colunas diferentes

```python
resultado = disp.disproportionality(
    df,
    a_col="notificacoes_medicamento_evento",
    b_col="notificacoes_medicamento_outros_eventos",
    c_col="notificacoes_outros_medicamentos_evento",
    d_col="notificacoes_outros_medicamentos_outros_eventos",
    add_signal_flags=True,
)
```

## `inplace`

Por padrao, `inplace=False`. Isso significa que a funcao retorna uma copia do
DataFrame original e nao modifica `df`.

```python
resultado = disp.disproportionality(df, "A", "B", "C", "D")
```

Se quiser modificar o DataFrame original:

```python
disp.disproportionality(df, "A", "B", "C", "D", inplace=True)
```

## Parametros principais

### `correction`

Usado em `ROR` e `PRR` para evitar divisao por zero.

Padrao:

```python
correction=0.5
```

Para reproduzir formulas sem correcao em exemplos teoricos:

```python
disp.ror(df, "A", "B", "C", "D", correction=0)
disp.prr(df, "A", "B", "C", "D", correction=0)
```

### `shrinkage`

Usado em `IC`, `IC025` e `IC975`.

Padrao:

```python
shrinkage=0.5
```

## Validacoes

A biblioteca valida automaticamente:

- se `df` e um `pandas.DataFrame`
- se as colunas `A`, `B`, `C` e `D` existem
- se as colunas sao numericas
- se nao ha valores ausentes
- se nao ha valores negativos
- se `A + B + C + D` nao e zero no calculo do `IC`

Exemplos de problemas que geram erro:

- coluna informada nao existe
- coluna com texto em vez de numero
- valores negativos
- valores ausentes
- linha com `A = B = C = D = 0`

## Equacoes

### ROR

```text
A_calc = A + correction
B_calc = B + correction
C_calc = C + correction
D_calc = D + correction

ROR = (A_calc * D_calc) / (B_calc * C_calc)
log_ROR = ln(ROR)
SE_log_ROR = sqrt(1/A_calc + 1/B_calc + 1/C_calc + 1/D_calc)

ROR_lower_95 = exp(log_ROR - 1.96 * SE_log_ROR)
ROR_upper_95 = exp(log_ROR + 1.96 * SE_log_ROR)
```

### PRR

```text
A_calc = A + correction
B_calc = B + correction
C_calc = C + correction
D_calc = D + correction

PRR = [A_calc / (A_calc + B_calc)] / [C_calc / (C_calc + D_calc)]
log_PRR = ln(PRR)

SE_log_PRR = sqrt(
    1/A_calc
    - 1/(A_calc + B_calc)
    + 1/C_calc
    - 1/(C_calc + D_calc)
)

PRR_lower_95 = exp(log_PRR - 1.96 * SE_log_PRR)
PRR_upper_95 = exp(log_PRR + 1.96 * SE_log_PRR)
```

### IC

```text
Obs = A
N = A + B + C + D
Esp = ((A + B) * (A + C)) / N

IC = log2((Obs + shrinkage) / (Esp + shrinkage))
```

### IC025 e IC975

```text
IC025 = IC - 3.3 * (Obs + shrinkage)^(-0.5)
          - 2.0 * (Obs + shrinkage)^(-1.5)

IC975 = IC + 2.4 * (Obs + shrinkage)^(-0.5)
          - 0.5 * (Obs + shrinkage)^(-1.5)
```

## Testar a biblioteca

Na pasta do projeto:

```bash
cd /Users/bruno/Documents/Developer/dispropy
python3 -m pytest
```

Ou, se ainda nao instalou as dependencias de desenvolvimento:

```bash
python3 -m pip install -e ".[dev]"
python3 -m pytest
```

## Interpretacao cautelosa

- `ROR > 1` sugere maior odds de relato.
- `PRR > 1` sugere maior proporcao de relato.
- `IC > 0` sugere observado maior que esperado.
- `ROR_lower_95 > 1` pode ser usado como criterio estatistico de triagem.
- `PRR_lower_95 > 1` pode ser usado como criterio estatistico de triagem.
- `IC025 > 0` pode ser usado como criterio estatistico de triagem.
- Nenhuma dessas metricas prova causalidade.
- A avaliacao clinica e farmacologica continua necessaria.

## Limitacoes

- A biblioteca nao corrige vieses de notificacao espontanea.
- A biblioteca nao estima incidencia ou risco absoluto.
- A biblioteca nao substitui avaliacao clinica.
- A biblioteca nao implementa EBGM nesta versao.
- A biblioteca nao implementa graficos nesta versao.
- A biblioteca nao implementa interface web nesta versao.
- A biblioteca nao implementa leitura interna de arquivos nesta versao.
- A biblioteca nao implementa estratificacao por idade, sexo, pais ou periodo
  nesta primeira versao.

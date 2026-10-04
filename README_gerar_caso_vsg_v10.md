# Gerador de casos VSG v10 — `gerar_caso_vsg_v10.py`

## Visão geral

Este utilitário transforma requisitos de engenharia em um arquivo YAML pronto para `vsg_2a_ordem_degrau_carga_v10.py`. Ele:

- interpreta potências, tensões e frequências com unidades;
- converte potência ativa e fator de potência em cargas RL;
- dimensiona `Lf`, `Cf` e `Rf` por regras de ripple, reativo e X/R;
- calcula a inércia virtual a partir de um limite de RoCoF;
- estima o VDC mínimo e o ponto de operação ilhado;
- em modo rede, verifica ressonância LCL e usa a análise modal completa do simulador para dimensionar `D_w`;
- prepara eventos de fase, frequência, falta e sincronização;
- configura limitador de corrente e anti-windup;
- na v10, inclui opcionalmente o modelo Thévenin da bateria e o barramento CC dinâmico;
- gera uma matriz de dez casos de aceitação.

O gerador **não executa a simulação completa**. Sua saída principal é um YAML e um resumo de projeto no console.

## Relação com o simulador

Mantenha estes arquivos na mesma pasta:

```text
gerar_caso_vsg_v10.py
vsg_2a_ordem_degrau_carga_v10.py
```

O simulador é carregado pelo gerador quando necessário para análise modal em modo rede e para `--gerar-exemplos`. Se não estiver na mesma pasta, essas operações falham.

## Requisitos

- Python 3;
- NumPy;
- o simulador v10 na mesma pasta;
- SciPy para a análise modal acionada pelo simulador em casos de rede;
- PyYAML para executar posteriormente o YAML no simulador;
- Matplotlib apenas se a execução posterior do simulador gerar PNGs.

Instalação recomendada:

```bash
python -m venv .venv
source .venv/bin/activate       # Linux/macOS
# .venv\Scripts\Activate.ps1  # Windows PowerShell
python -m pip install numpy scipy matplotlib pyyaml
```

## Início rápido

### Gerar o caso padrão

```bash
python gerar_caso_vsg_v10.py
```

O arquivo padrão é `meu_caso.yaml`.

### Gerar e simular

```bash
python gerar_caso_vsg_v10.py -o casos/meu_caso.yaml --prefixo resultados/meu_caso
python vsg_2a_ordem_degrau_carga_v10.py --config casos/meu_caso.yaml
```

### Consultar todas as opções

```bash
python gerar_caso_vsg_v10.py --help
```

## Unidades aceitas

| Grandeza | Formatos aceitos | Observações |
|---|---|---|
| Potência ativa | `W`, `kW`, `MW`, `GW`, `pu` | `pu` usa `Sn` como base. Sem unidade, assume W. |
| Potência aparente | `VA`, `kVA`, `MVA`, `GVA` | Sem unidade, assume VA. |
| Tensão | `V`, `kV`, `pu` | `pu` só é aceito quando há uma base disponível, como em `--v-inversor`. |
| Frequência | `Hz`, `kHz` | Sem unidade, assume Hz. |
| Listas | valores separados por vírgula | Usado por `--ocv-soc-pu`, `--ocv-v`, `--scr-faixa` e `--h-faixa`. |

Ponto ou vírgula podem ser usados como separador decimal em valores simples. Em listas, a vírgula é o separador entre elementos; use ponto decimal nos itens.

## Regras de projeto

### Bases

```text
Zb = VLL² / Sn
Lb = Zb / ω0
Cb = 1 / (ω0·Zb)
În = √2·Sn / (√3·VLL)
```

### Cargas RL

A partir de `P` e fator de potência indutivo:

```text
Q = P·tan(arccos(fp))
R = VLL²·P / (P² + Q²)
L = [VLL²·Q / (P² + Q²)] / ω0
```

A carga inicial usa `--carga`/`--fp`; o bloco conectado no degrau usa `--degrau`/`--fp-degrau`.

### Filtro

```text
Vdc(auto) = 1,3·√2·VLL
ΔI        = ripple·În
Lf        = Vdc / (6·fsw·ΔI)
Cf        = qc_max·Cb
Xf        = ω0·Lf
Rf        = Xf / xr
fr,LC     = 1 / (2π√(Lf·Cf))
Rd        = 1 / (3·2π·fr·Cf)
```

O gerador avisa quando `fr` fica fora de `10·f0 < fr < fsw/2`, quando `Xf > 0,15 pu` e quando o VDC informado fica abaixo do mínimo estimado para SVPWM.

> `--vdc` participa do **projeto do filtro e da verificação de margem**. `--vdc-inicial-v` pertence ao **modelo dinâmico do lado CC**. São parâmetros diferentes.

### Inércia

```text
Hmin = ΔPmax_pu·f0 / (2·RoCoFmax)
H    = Hmin·margem_h, arredondado para cima em passos de 0,1 s
Hmax = (pmax_pcs − P1_pu)·f0 / (2·RoCoF_rede)
```

`--h-faixa min,max` limita o H adotado e produz aviso se o requisito de RoCoF não puder ser atendido.

### Rede e washout

Para cada SCR da faixa:

1. calcula `Rg` e `Lg` com base em `Zb/SCR` e `X/R`;
2. verifica a frequência de ressonância LCL;
3. lineariza o modelo completo do simulador no equilíbrio;
4. identifica o modo eletromecânico dominante;
5. se `--dw` não foi fornecido, procura `D_w` entre 0 e 2000 e refina por bisseção até atender aproximadamente `ζ_alvo + 0,005` no SCR mais forte;
6. relata o amortecimento com e sem washout em todos os SCRs considerados.

### Eventos de frequência

Para rampas, o resumo estima:

```text
Heff       = H + Dw·Tw/2
Tr         = |Δf| / RoCoF
ΔPinercial = (RoCoF/f0)·[2H + Dw·Tw·(1 − exp(−Tr/Tw))]
ΔPdroop    = Dp·|Δf|/f0
Ppico      = P1_pu + ΔPinercial + ΔPdroop
```

É emitido aviso quando `Ppico > pmax_pcs` ou quando o washout domina a inércia efetiva.

### Limitador de corrente

Com `--imax-pu > 0`, o gerador:

- usa `i_on = 0,98·imax` se `--i-on-pu` não for informado;
- usa `rv_max = 3 pu` se `--rv-max-pu` for omitido;
- faz uma pré-estimativa conservadora da corrente e da resistência virtual necessária;
- avisa se `rv_max` parecer insuficiente ou se a resposta inercial prevista for recortada.

A aceitação final do limitador depende da simulação EMT; a pré-verificação do gerador não substitui os resultados do simulador.

## Exemplos

### Caso ilhado dimensionado automaticamente

```bash
python gerar_caso_vsg_v10.py \
  --sn "100 kVA" --vll "380 V" \
  --carga "50 kW" --fp 0.958 \
  --degrau "30 kW" --fp-degrau 0.958 \
  --evento carga \
  --prefixo resultados/ilhado \
  -o casos/ilhado.yaml
```

### Rede com faixa de SCR e amortecimento automático

```bash
python gerar_caso_vsg_v10.py \
  --modo rede --scr 5 --scr-faixa 2,20 --xr-rede 10 \
  --zeta-alvo 0.30 --tw 1 \
  --evento carga \
  --prefixo resultados/rede \
  -o casos/rede.yaml
```

### Rampa de frequência

```bash
python gerar_caso_vsg_v10.py \
  --modo rede --evento freq_rampa \
  --df-rede -0.5 --rocof-rede 0.5 \
  --pmax-pcs 1.1 \
  -o casos/freq_rampa.yaml
```

### Falta trifásica e limitação de corrente

```bash
python gerar_caso_vsg_v10.py \
  --modo rede --evento falta_3f \
  --t-step 1 --t-clear 1.15 --t-end 3 --vg-falta 0 \
  --imax-pu 1.2 --rv-max-pu 3 --k-aw 20 \
  -o casos/falta.yaml
```

### Fechamento ativo

```bash
python gerar_caso_vsg_v10.py \
  --modo rede --evento fechamento \
  --disjuntor-inicial aberto --estrategia-sync ativo \
  --delta-g0-graus 20 --df-g0-hz 0.12 \
  --t-step 0.5 --t-end 8 \
  -o casos/fechamento_ativo.yaml
```

### Bateria Thévenin v10

```bash
python gerar_caso_vsg_v10.py \
  --modelo-cc thevenin \
  --c-dc-f 0.10 --vdc-inicial-v auto \
  --capacidade-ah 100 --soc-inicial 0.60 --r0-ohm 0.05 \
  --ocv-soc-pu 0,0.25,0.5,0.75,1 \
  --ocv-v 660,680,700,720,740 \
  --vdc-min-oper-v 500 --vdc-max-oper-v 800 \
  --m-max-pu 1.0 \
  --prefixo resultados/bess \
  -o casos/bess.yaml
```

## Opções completas

### Geral

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--gerar-exemplos PASTA` | `—` | Gera os dez casos de aceitação v10 e encerra. |
| `-o, --saida ARQUIVO` | `meu_caso.yaml` | Arquivo YAML de saída. |
| `--prefixo PREFIXO` | `vsg` | Prefixo que será gravado na seção `saida` do caso. |

### Sistema

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--sn VALOR` | `100 kVA` | Potência aparente nominal; aceita VA, kVA, MVA e GVA. |
| `--vll VALOR` | `380 V` | Tensão de linha RMS; aceita V e kV. |
| `--f0 HZ` | `60` | Frequência nominal. |
| `--modo MODO` | `ilhado` | `ilhado` ou `rede`. |

### Rede, amortecimento e evento

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--scr VALOR` | `5` | SCR do caso simulado. |
| `--scr-faixa MIN,MAX` | `null` | Extremos de SCR usados nas verificações e no projeto de `D_w`. |
| `--xr-rede VALOR` | `10` | Relação X/R da rede. |
| `--tw S` | `1` | Constante de tempo do washout. |
| `--zeta-alvo VALOR` | `0.3` | Amortecimento modal mínimo desejado. |
| `--dw VALOR` | `auto` | Fixa `D_w`; quando ausente, o gerador o dimensiona em rede. |
| `--evento NOME` | `carga` | Evento: `nenhum`, `carga`, `pref`, `fase`, `freq_degrau`, `freq_rampa`, `falta_3f`, `fechamento` ou `sincronizacao`. |
| `--df-rede HZ` | `0` | Desvio final/associado da frequência da rede. |
| `--d-fase GRAUS` | `5` | Salto angular; também serve à interface legada `sincronizacao`. |
| `--vg-falta PU` | `0` | Multiplicador da tensão durante a falta. |
| `--t-clear S` | `null` | Instante de eliminação da falta. |

### Cargas

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--carga VALOR` | `50 kW` | Potência ativa inicial; aceita W/kW/MW/GW/pu. |
| `--fp VALOR` | `0.958` | Fator de potência indutivo da carga inicial. |
| `--degrau VALOR` | `30 kW` | Potência ativa do segundo bloco. |
| `--fp-degrau VALOR` | `igual a --fp` | Fator de potência indutivo do segundo bloco. |

### Filtro

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--vdc VALOR` | `auto` | VDC usado no projeto do filtro; `auto` ou V/kV. |
| `--fsw VALOR` | `8 kHz` | Frequência de chaveamento; Hz/kHz. |
| `--ripple FRAÇÃO` | `0.15` | Ripple de corrente usado para dimensionar `Lf`. |
| `--qc-max PU` | `0.05` | Limite de reativo do capacitor; define `Cf`. |
| `--xr VALOR` | `40` | Relação X/R do indutor do filtro; define `Rf`. |
| `--margem-vdc FRAÇÃO` | `0.10` | Margem usada na verificação de VDC mínimo. |

### Inércia e capacidade do PCS

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--rocof-max HZ/S` | `2` | RoCoF máximo para calcular `H_min`. |
| `--dpmax VALOR` | `igual ao degrau` | Maior desequilíbrio de potência; W/kW/MW/GW/pu. |
| `--margem-h FATOR` | `1` | Margem aplicada a `H_min`. |
| `--h-faixa MIN,MAX` | `null` | Faixa ajustável de H do PCS. |
| `--pmax-pcs PU` | `1.1` | Limite de potência usado para calcular `H_max` e emitir avisos. |
| `--rocof-rede HZ/S` | `igual a --rocof-max` | RoCoF da rede para `H_max` e para `freq_rampa`. |

### Controle e limitador

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--v-inversor VALOR` | `auto` | Referência do inversor; `auto`, V, kV ou pu. |
| `--mp PU` | `0.05` | Droop P-f. |
| `--nq PU` | `0.05` | Droop Q-V. |
| `--fc HZ` | `10` | Corte da medição de P/Q. |
| `--imax-pu PU` | `0` | Limite de corrente; zero desabilita. |
| `--i-on-pu PU` | `0,98·Imax` | Limiar de início da impedância virtual. |
| `--rv-max-pu PU` | `2 sem limitador; 3 com limitador` | Teto de resistência virtual; automático quando omitido. |
| `--xv-rv VALOR` | `0` | Relação Xv/Rv. |
| `--k-aw 1/S` | `20` | Ganho de anti-windup. |

### Lado CC e bateria v10

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--modelo-cc {ideal,thevenin}` | `ideal` | Seleciona o lado CC. |
| `--c-dc-f F` | `null` | Capacitância do barramento CC. |
| `--vdc-inicial-v V|auto` | `auto` | VDC inicial do modelo dinâmico. |
| `--permitir-desequilibrio-inicial` | `false` | Aceita VDC manual fora do equilíbrio. |
| `--vdc-min-oper-v V` | `null` | Limite diagnóstico inferior. |
| `--vdc-max-oper-v V` | `null` | Limite diagnóstico superior. |
| `--m-max-pu PU` | `1` | Índice máximo de modulação. |
| `--vdc-piso-numerico-v V|auto` | `auto` | Piso terminal de VDC. |
| `--capacidade-ah AH` | `null` | Capacidade do pack. |
| `--soc-inicial PU` | `null` | SOC inicial. |
| `--r0-ohm OHM` | `null` | Resistência série Thévenin. |
| `--ocv-soc-pu LISTA` | `null` | Nós de SOC separados por vírgula. |
| `--ocv-v LISTA` | `null` | OCVs correspondentes, em volts. |
| `--soc-min-oper-pu PU` | `0.10` | Limite operacional inferior de SOC. |
| `--soc-max-oper-pu PU` | `0.90` | Limite operacional superior de SOC. |
| `--i-desc-max-a A` | `0` | Limite diagnóstico de descarga; zero desabilita. |
| `--i-carga-max-a A` | `0` | Limite diagnóstico de carga; zero desabilita. |

### Sincronismo e disjuntor

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--delta-g0-graus GRAUS` | `0` | Ângulo inicial da rede para `fechamento`. |
| `--df-g0-hz HZ` | `0` | Desvio inicial de frequência para `fechamento`. |
| `--disjuntor-inicial {aberto,fechado}` | `automático` | Aberto por padrão em `fechamento`; fechado nos demais eventos. |
| `--estrategia-sync {forcado,passivo,ativo}` | `passivo` | Estratégia do supervisor. |
| `--dv-sync-max-pu PU` | `0.05` | Janela de ΔV. |
| `--df-sync-max-hz HZ` | `0.10` | Janela de Δf. |
| `--dtheta-sync-max-graus GRAUS` | `5` | Janela angular atual e prevista. |
| `--t-sync-hold-s S` | `0.10` | Tempo de qualificação. |
| `--t-fechamento-s S` | `0.06` | Tempo mecânico. |
| `--t-sync-timeout-s S` | `10` | Timeout da manobra. |
| `--dt-rele-s S` | `0.001` | Passo fixo do relé 25. |
| `--antecipar-fechamento {true,false}` | `true` | Habilita antecipação angular. |
| `--vmin-medicao-pu PU` | `0.20` | Tensão mínima para medição válida. |
| `--kp-theta-hz-rad` | `0.60` | Ganho P angular. |
| `--ki-theta-hz-rad-s` | `0.20` | Ganho I angular. |
| `--df-sync-lim-hz` | `0.50` | Limite do bias de frequência. |
| `--kp-v` | `0.80` | Ganho P de tensão. |
| `--ki-v-s` | `0.30` | Ganho I de tensão. |
| `--de-sync-lim-pu` | `0.10` | Limite do bias de tensão. |
| `--t-release-sync-s` | `0.20` | Retirada dos biases após o contato. |
| `--pll-sync-bw-hz` | `5` | Banda do PLL dedicado. |
| `--pre-sync {0,1}` | `1` | Compatibilidade com `sincronizacao`. |
| `--kp-sync-p` | `4` | Ganho P angular legado. |
| `--ki-sync-p` | `2` | Ganho I angular legado. |
| `--p-sync-max-pu` | `0.50` | Limite de viés ativo legado. |
| `--kp-sync-v` | `0.50` | Ganho P de tensão legado. |
| `--ki-sync-v` | `1` | Ganho I de tensão legado. |
| `--e-sync-max-pu` | `0.20` | Limite de correção de tensão legado. |
| `--sync-dv-max-pu` | `0.05` | Janela ΔV legada. |
| `--sync-df-max-hz` | `0.10` | Janela Δf legada. |
| `--sync-delta-max-deg` | `10` | Janela angular prevista legada. |
| `--sync-hold-s` | `0.10` | Qualificação legada. |
| `--breaker-delay-s` | `0.05` | Tempo mecânico legado. |

### Simulação

| Opção | Padrão | Finalidade |
|---|---:|---|
| `--t-step S` | `1` | Instante do evento. |
| `--t-end S` | `4` | Tempo total. |
| `--dt-out S` | `5e-5` | Passo de amostragem da saída. |

## Regras e validações relevantes

- `--modo` deve ser `ilhado` ou `rede`;
- `fase`, `freq_degrau`, `freq_rampa`, `falta_3f`, `fechamento` e `sincronizacao` exigem rede;
- `fechamento` exige disjuntor inicialmente aberto;
- potências de carga devem ser positivas e fatores de potência devem estar em `(0,1]`;
- `--t-end` deve ser maior que `--t-step`;
- em `falta_3f`, `t_step < t_clear < t_end` e `vg_falta >= 0`;
- eventos de frequência exigem `--df-rede != 0` e RoCoF positivo;
- com Thévenin, são obrigatórios `c_dc_f`, `capacidade_ah`, `soc_inicial`, `r0_ohm`, `ocv_soc_pu` e `ocv_v`;
- ganhos e limites de sincronização são verificados quanto ao sinal;
- a estratégia ativa exige ao menos um canal de PI útil no simulador;
- mensagens iniciadas por `(!)` são advertências de engenharia; o YAML ainda pode ser gerado.

## Saída YAML

O YAML contém:

```text
sistema
rede
filtro
carga
controle
sincronismo
[lado_cc e bateria, quando a interface v10 foi explicitamente usada]
evento
simulacao
saida
```

No modo rede, `pref` e `qref` são gravados como `null`. O simulador então resolve um equilíbrio com troca inicial nula com a rede, no qual a carga local é suprida pelo VSG.

O cabeçalho do arquivo registra:

- o comando usado para gerar o caso;
- o resumo de bases, cargas, filtro e inércia;
- resultados modais e de rede quando aplicáveis;
- avisos de engenharia.

Strings são gravadas entre aspas simples para evitar que caminhos Windows com barra invertida sejam interpretados como sequências de escape.

## Opt-in do lado CC v10

Para preservar regressões antigas, as seções `lado_cc` e `bateria` só são materializadas quando alguma opção v10 do lado CC aparece explicitamente na linha de comando. Isso inclui `--modelo-cc`, mesmo quando o valor informado é `ideal`.

Consequências práticas:

- chamada sem opções CC: YAML compatível com a interface v5–v9; o simulador usa o caminho ideal legado;
- chamada com `--modelo-cc ideal`: YAML explicita a interface v10, mas sem dinâmica de bateria;
- chamada com `--modelo-cc thevenin`: YAML inclui todos os parâmetros do barramento CC e da bateria.

## Dez exemplos de aceitação

```bash
python gerar_caso_vsg_v10.py --gerar-exemplos exemplos_v10
```

Arquivos gerados:

1. `caso_cc_ideal_regressao_v10.yaml` — regressão do caminho ideal;
2. `caso_bess_degrau_carga_v10.yaml` — BESS ilhado com degrau;
3. `caso_bess_cdc_baixo_v10.yaml` — sensibilidade a Cdc baixo;
4. `caso_bess_cdc_alto_v10.yaml` — sensibilidade a Cdc alto;
5. `caso_bess_regeneracao_v10.yaml` — potência regenerativa em rede;
6. `caso_bess_saturacao_modulacao_v10.yaml` — saturação de modulação;
7. `caso_bess_falta_3f_v10.yaml` — falta trifásica;
8. `caso_bess_freq_rampa_v10.yaml` — rampa de frequência;
9. `caso_bess_fechamento_ativo_v10.yaml` — fechamento com pré-sincronização ativa;
10. `caso_bess_inviavel_v10.yaml` — caso deliberadamente inviável para validação.

A opção apenas grava a matriz; não executa os dez casos.

## Fluxo recomendado

1. Defina bases, carga e limites do PCS.
2. Gere o YAML e leia todos os avisos.
3. Revise `Lf`, `Cf`, `Rf`, ressonância, VDC mínimo, `H`, `Hmax` e `D_w`.
4. Em v10, confirme que a curva OCV é do pack completo e que `R0`, `Cdc` e `QAh` usam as unidades esperadas.
5. Execute o simulador com o YAML.
6. Confira o JSON, principalmente `modal`, `limitador`, `sincronizacao` e `lado_cc`.
7. Faça varreduras de SCR, SOC, Cdc, RoCoF, falta e tolerâncias de sincronismo antes de usar o caso como evidência de projeto.

## Solução de problemas

| Mensagem/sintoma | Causa provável | Correção |
|---|---|---|
| Simulador não encontrado | Os dois `.py` não estão na mesma pasta. | Reúna os arquivos sem alterar o nome `vsg_2a_ordem_degrau_carga_v10.py`. |
| Unidade inválida | Sufixo não suportado ou espaço/formato incorreto. | Use as unidades da tabela; coloque valores com espaço entre aspas no shell. |
| Não foi possível atingir ζ com `D_w <= 2000` | Meta modal incompatível com o caso. | Reveja `zeta-alvo`, filtro, SCR, X/R e `Tw`; ou fixe `--dw` conscientemente. |
| Ressonância fora da faixa | Combinação de `fsw`, ripple, `qc_max`, filtro e SCR. | Redimensione o filtro e considere amortecimento físico. |
| Resposta inercial limitada | `H`, `D_w·T_w`, RoCoF ou Δf exigem mais que o PCS. | Ajuste o requisito ou aumente a capacidade; valide no simulador. |
| Parâmetros Thévenin ausentes | `--modelo-cc thevenin` sem todos os campos obrigatórios. | Informe Cdc, QAh, SOC, R0 e as duas listas OCV. |
| O YAML não contém `lado_cc` | Nenhuma opção CC v10 foi passada. | Inclua `--modelo-cc ideal` ou `--modelo-cc thevenin` explicitamente. |

## Limitações

- As regras de filtro e VDC são critérios de pré-dimensionamento, não substituem projeto de hardware, perdas, térmica, EMC ou coordenação de proteção.
- O cálculo modal depende do modelo e do equilíbrio implementados no simulador; mudanças no simulador podem alterar o `D_w` calculado.
- O fator de potência de carga é tratado como indutivo.
- O gerador não estima automaticamente `Cdc`, capacidade Ah, `R0` ou curva OCV; esses dados precisam vir do projeto do BESS.
- A pré-estimativa de corrente é conservadora e não substitui a simulação temporal.
- O YAML não inclui opções que o gerador não expõe, como `d_pref`, `max_step` e `rtol`; ajuste-as depois no YAML ou na linha de comando do simulador.

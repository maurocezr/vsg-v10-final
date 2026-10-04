# Simulador VSG v10 — `vsg_2a_ordem_degrau_carga_v10.py`

## Visão geral

Este programa executa uma simulação EMT de **modelo médio** de um inversor trifásico baseado em recursos (IBR) com controle VSM/VSG de segunda ordem. O mesmo núcleo atende a:

- operação **ilhada**, alimentando uma carga RL local;
- operação **conectada à rede**, representada por um barramento infinito atrás de uma impedância definida por SCR e X/R;
- eventos de carga, referência de potência, fase, frequência, falta trifásica e fechamento/sincronização;
- limitação de corrente por impedância virtual adaptativa e anti-windup;
- disjuntor, relé 25 e pré-sincronizador;
- na v10, um lado CC opcional com bateria Thévenin, capacitor de barramento CC, SOC e saturação de modulação SVPWM.

O modo `modelo_cc: ideal` preserva o caminho numérico legado. O modo `modelo_cc: thevenin` adiciona os estados dinâmicos `VDC` e `SOC` e publica balanços de potência e energia do lado CC.

## Topologia e modelo

### Circuito por fase

```text
 e(VSG) -- [Rf + Lf] --+-- vc (PCC) -- [Rg + Lg] -- vg (rede)
                       |
                       +-- Cf
                       +-- [RL1 + LL1]
                       +-- S -- [RL2 + LL2]   (degrau de carga)
```

O ramo `Rg + Lg` só participa quando há rede e o disjuntor está fechado. No evento `fechamento`, a simulação começa em ilha, mantém o ramo da rede aberto e o conecta no instante de contato.

### Controle VSG

As equações centrais em pu são:

```text
2H · dΔω/dt = Pref − Pf − Dp·Δω − Dw·(Δω − z)
dz/dt       = (Δω − z) / Tw
dδv/dt      = ω0·Δω
E            = Eref − nq·(Qf − Qref)
dPf/dt       = ωc·(P − Pf)
dQf/dt       = ωc·(Q − Qf)
Dp           = 1/mp
S            = 1,5·Vc·conj(If) / Sn
```

A integração é feita em referencial síncrono fixo em `ω0`. Em rede, o filtro e a impedância da rede formam um LCL. A análise modal usa um Jacobiano numérico no equilíbrio inicial.

### Lado CC v10

No modo Thévenin:

```text
Ibat       = (OCV(SOC) − VDC) / R0
Iconv      = Pconv / VDC
Cdc·dVDC/dt = Ibat − Iconv
dSOC/dt    = −Ibat / (3600·QAh)
```

Convenções:

- `Ibat > 0`: descarga da bateria;
- `Pconv > 0`: fluxo de potência CC → CA;
- a OCV é interpolada linearmente e **não é extrapolada**;
- a saturação SVPWM é radial, com `Emax = m_max·VDC/(√2·VLL_base)`;
- a integração termina se `VDC` atingir o piso numérico ou se o SOC alcançar o limite do domínio da curva OCV.

Quando `vdc_inicial_v: auto`, o programa procura o ponto de alta tensão que equilibra bateria, conversor e limite de modulação. Um valor manual fora do equilíbrio só é aceito com `permitir_desequilibrio_inicial: true`.

## Requisitos

- Python 3;
- NumPy;
- SciPy;
- Matplotlib, apenas para gerar gráficos;
- PyYAML, apenas para ler `.yaml`/`.yml`;
- Python 3.11+ para leitura nativa de `.toml` por `tomllib`.

Instalação recomendada:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib pyyaml
```

## Início rápido

### 1. Executar o caso padrão

```bash
python vsg_2a_ordem_degrau_carga_v10.py
```

O padrão é ilhado, com evento de degrau de carga em `t = 1 s`, término em `4 s` e prefixo `vsg`.

### 2. Gerar um modelo de configuração

```bash
python vsg_2a_ordem_degrau_carga_v10.py --gerar-config modelo_v10.yaml
```

### 3. Executar um arquivo de caso

```bash
python vsg_2a_ordem_degrau_carga_v10.py --config meu_caso.yaml
```

### 4. Sobrescrever valores do arquivo pela linha de comando

```bash
python vsg_2a_ordem_degrau_carga_v10.py \
  --config meu_caso.yaml \
  --t-end 8 \
  --prefixo resultados/ensaio_01 \
  --sem-graficos
```

A precedência é:

```text
valores padrão < arquivo de configuração < linha de comando
```

## Arquivos de configuração

São aceitos `.yaml`, `.yml`, `.json` e `.toml`. YAML requer PyYAML. O conteúdo pode ser organizado por seções ou fornecido como pares planos; as seções são achatadas internamente. Uma chave repetida após o achatamento ou uma chave desconhecida causa erro.

Exemplo mínimo:

```yaml
sistema:
  modo: 'rede'

rede:
  scr: 5
  xr_rede: 10

controle:
  H: 5
  dw: 100

evento:
  evento: 'freq_rampa'
  df_g: -0.5
  rocof_g: 0.5

simulacao:
  t_step: 1
  t_end: 4

saida:
  prefixo: 'resultados/freq_rampa'
```

`null`, `none`, `auto` e string vazia são aceitos para parâmetros anuláveis. Listas como `ocv_soc_pu` e `ocv_v` podem ser listas reais no arquivo ou valores separados por vírgula na linha de comando.

## Eventos

| Evento | Modo | Efeito principal | Requisitos relevantes |
|---|---|---|---|
| `nenhum` | ilhado ou rede | Mantém o ponto inicial durante todo o intervalo. | Nenhum requisito adicional. |
| `carga` | ilhado ou rede | Conecta a segunda carga RL em `t_step`. | `p2 > 0`; `t_end > t_step`. |
| `pref` | ilhado ou rede | Aplica `d_pref` à referência ativa em `t_step`. | `t_end > t_step`. |
| `fase` | rede | Aplica salto `d_fase` ao ângulo da rede. | `modo: rede`. |
| `freq_degrau` | rede | Aplica degrau `df_g` à frequência da rede. | `df_g != 0`. |
| `freq_rampa` | rede | Aplica rampa até `df_g`, com módulo de RoCoF `rocof_g`. | `df_g != 0`; `rocof_g > 0`. |
| `falta_3f` | rede | Multiplica a tensão de rede por `vg_falta` entre `t_step` e `t_clear`. | `t_step < t_clear < t_end`; `vg_falta >= 0`. |
| `sincronizacao` | rede | Interface de compatibilidade do pré-sincronizador anterior. | Usa `pre_sync`, `kp_sync_*`, `sync_*` e `breaker_delay_s`. |
| `fechamento` | rede | Manobra completa com disjuntor aberto, supervisor, relé 25 e contato. | `disjuntor_inicial: aberto`; estratégia `forcado`, `passivo` ou `ativo`. |

### Estratégias do evento `fechamento`

- `forcado`: emite comando no início da manobra sem exigir aprovação do relé 25;
- `passivo`: aguarda que ΔV, Δf e ângulo atual/previsto permaneçam nas janelas pelo tempo de qualificação;
- `ativo`: usa PI de ângulo/frequência e PI de tensão para reduzir os erros antes de aplicar o mesmo sync-check.

A antecipação angular considera o tempo mecânico de fechamento. Depois do contato, os biases ativos são retirados linearmente em `t_release_sync_s`.

## Exemplos de execução

### Caso ilhado com degrau de carga

```bash
python vsg_2a_ordem_degrau_carga_v10.py \
  --modo ilhado --evento carga \
  --p1 50000 --q1 15000 --p2 30000 --q2 10000 \
  --t-step 1 --t-end 4 \
  --prefixo resultados/ilhado
```

### Rede com rampa de frequência

```bash
python vsg_2a_ordem_degrau_carga_v10.py \
  --modo rede --scr 5 --xr-rede 10 \
  --evento freq_rampa --df-g -0.5 --rocof-g 0.5 \
  --dw 100 --tw 1 \
  --prefixo resultados/freq_rampa
```

### Falta trifásica com limitador de corrente

```bash
python vsg_2a_ordem_degrau_carga_v10.py \
  --modo rede --evento falta_3f \
  --t-step 1 --t-clear 1.15 --t-end 3 \
  --vg-falta 0 \
  --imax-pu 1.2 --rv-max-pu 3 --k-aw 20 \
  --prefixo resultados/falta
```

### Fechamento ativo

```bash
python vsg_2a_ordem_degrau_carga_v10.py \
  --modo rede --evento fechamento \
  --disjuntor-inicial aberto --estrategia-sync ativo \
  --delta-g0-graus 20 --df-g0-hz 0.12 \
  --t-step 0.5 --t-end 8 \
  --prefixo resultados/fechamento
```

### Bateria Thévenin e barramento CC dinâmico

```bash
python vsg_2a_ordem_degrau_carga_v10.py \
  --modo ilhado --evento carga \
  --modelo-cc thevenin \
  --c-dc-f 0.10 --vdc-inicial-v auto --m-max-pu 1.0 \
  --capacidade-ah 100 --soc-inicial 0.60 --r0-ohm 0.05 \
  --ocv-soc-pu 0,0.25,0.5,0.75,1 \
  --ocv-v 660,680,700,720,740 \
  --vdc-min-oper-v 500 --vdc-max-oper-v 800 \
  --prefixo resultados/bess
```

## Parâmetros completos

As chaves YAML usam sublinhado; as opções de linha de comando usam hífen. A opção `--H` mantém o `H` maiúsculo.

### `sistema`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `sn` / `--sn` | `VA` | `100000.0` | Potência nominal do IBR (base) |
| `vll` / `--vll` | `V` | `380.0` | Tensão de linha nominal (rms) |
| `f0` / `--f0` | `Hz` | `60.0` | Frequência nominal |
| `modo` / `--modo` | `-` | `ilhado` | Modo de operação: ilhado \| rede |

### `rede`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `scr` / `--scr` | `-` | `5.0` | Relação de curto-circuito da rede (base S_n) |
| `xr_rede` / `--xr-rede` | `-` | `10.0` | Relação X/R da impedância da rede |
| `vg` / `--vg` | `pu` | `1.0` | Tensão da rede |

### `filtro`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `rf` / `--rf` | `pu` | `0.01` | Resistência série do filtro |
| `xf` / `--xf` | `pu` | `0.15` | Reatância série do filtro (L_f) |
| `bc` / `--bc` | `pu` | `0.05` | Susceptância do capacitor (C_f) |

### `carga`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `p1` / `--p1` | `W` | `50000.0` | Carga inicial - potência ativa trifásica |
| `q1` / `--q1` | `var` | `15000.0` | Carga inicial - potência reativa trifásica |
| `p2` / `--p2` | `W` | `30000.0` | Degrau de carga - potência ativa trifásica |
| `q2` / `--q2` | `var` | `10000.0` | Degrau de carga - potência reativa trifásica |

### `controle`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `pref` / `--pref` | `pu` | `null / auto` | Setpoint de P (null = automático) |
| `qref` / `--qref` | `pu` | `null / auto` | Setpoint de Q (null = automático) |
| `eref` / `--eref` | `pu` | `1.0` | Referência da tensão interna (E_ref) |
| `H` / `--H` | `s` | `5.0` | Constante de inércia virtual |
| `mp` / `--mp` | `pu` | `0.05` | Droop P-f (D_p = 1/m_p) |
| `nq` / `--nq` | `pu` | `0.05` | Droop Q-V |
| `fc` / `--fc` | `Hz` | `10.0` | Corte do filtro de medição de P e Q |
| `dw` / `--dw` | `pu` | `0.0` | Ganho de amortecimento washout D_w |
| `tw` / `--tw` | `s` | `1.0` | Constante de tempo do washout T_w |
| `imax_pu` / `--imax-pu` | `pu` | `0.0` | Limite da envoltória de I_f (0 desabilita) |
| `i_on_pu` / `--i-on-pu` | `pu` | `null / auto` | Início da atuação (auto = 0,98 I_max) |
| `rv_max_pu` / `--rv-max-pu` | `pu` | `2.0` | Teto da resistência virtual |
| `xv_rv` / `--xv-rv` | `-` | `0.0` | Relação X_v/R_v |
| `k_aw` / `--k-aw` | `1/s` | `20.0` | Ganho do anti-windup angular |

### `sincronismo`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `delta_g0_graus` / `--delta-g0-graus` | `graus` | `0.0` | Ângulo inicial da rede na manobra |
| `df_g0_hz` / `--df-g0-hz` | `Hz` | `0.0` | Desvio inicial de frequência da rede |
| `disjuntor_inicial` / `--disjuntor-inicial` | `-` | `fechado` | Estado inicial: aberto \| fechado |
| `estrategia_sync` / `--estrategia-sync` | `-` | `passivo` | forcado \| passivo \| ativo |
| `dv_sync_max_pu` / `--dv-sync-max-pu` | `pu` | `0.05` | Janela máxima de diferença de tensão |
| `df_sync_max_hz` / `--df-sync-max-hz` | `Hz` | `0.1` | Janela máxima de diferença de frequência |
| `dtheta_sync_max_graus` / `--dtheta-sync-max-graus` | `graus` | `5.0` | Janela angular atual e prevista |
| `t_sync_hold_s` / `--t-sync-hold-s` | `s` | `0.1` | Tempo contínuo dentro das janelas |
| `t_fechamento_s` / `--t-fechamento-s` | `s` | `0.06` | Tempo mecânico até o contato |
| `t_sync_timeout_s` / `--t-sync-timeout-s` | `s` | `10.0` | Prazo máximo da manobra |
| `dt_rele_s` / `--dt-rele-s` | `s` | `0.001` | Passo fixo do relé 25 |
| `antecipar_fechamento` / `--antecipar-fechamento` | `true|false` | `true` | Antecipa o erro angular no contato |
| `vmin_medicao_pu` / `--vmin-medicao-pu` | `pu` | `0.2` | Tensão mínima para medição válida |
| `kp_theta_hz_rad` / `--kp-theta-hz-rad` | `Hz/rad` | `0.6` | Ganho P do PI angular |
| `ki_theta_hz_rad_s` / `--ki-theta-hz-rad-s` | `Hz/(rad.s)` | `0.2` | Ganho I do PI angular |
| `df_sync_lim_hz` / `--df-sync-lim-hz` | `Hz` | `0.5` | Limite do bias de frequência |
| `kp_v` / `--kp-v` | `pu/pu` | `0.8` | Ganho P do PI de tensão |
| `ki_v_s` / `--ki-v-s` | `1/s` | `0.3` | Ganho I do PI de tensão |
| `de_sync_lim_pu` / `--de-sync-lim-pu` | `pu` | `0.1` | Limite do bias de tensão |
| `t_release_sync_s` / `--t-release-sync-s` | `s` | `0.2` | Retirada linear dos biases |
| `pll_sync_bw_hz` / `--pll-sync-bw-hz` | `Hz` | `5.0` | Banda do PLL dedicado |
| `pre_sync` / `--pre-sync` | `0|1` | `1.0` | Compatibilidade: habilita pré-sync antigo |
| `kp_sync_p` / `--kp-sync-p` | `pu/rad` | `4.0` | Ganho P do sincronizador angular |
| `ki_sync_p` / `--ki-sync-p` | `pu/(rad.s)` | `2.0` | Ganho I do sincronizador angular |
| `p_sync_max_pu` / `--p-sync-max-pu` | `pu` | `0.5` | Limite do viés de potência do sincronizador |
| `kp_sync_v` / `--kp-sync-v` | `pu/pu` | `0.5` | Ganho P do sincronizador de tensão |
| `ki_sync_v` / `--ki-sync-v` | `1/s` | `1.0` | Ganho I do sincronizador de tensão |
| `e_sync_max_pu` / `--e-sync-max-pu` | `pu` | `0.2` | Limite da correção de tensão do sincronizador |
| `sync_dv_max_pu` / `--sync-dv-max-pu` | `pu` | `0.05` | Janela máxima de diferença de tensão |
| `sync_df_max_hz` / `--sync-df-max-hz` | `Hz` | `0.1` | Janela máxima de diferença de frequência |
| `sync_delta_max_deg` / `--sync-delta-max-deg` | `graus` | `10.0` | Janela máxima de ângulo previsto |
| `sync_hold_s` / `--sync-hold-s` | `s` | `0.1` | Tempo contínuo exigido dentro das janelas |
| `breaker_delay_s` / `--breaker-delay-s` | `s` | `0.05` | Tempo mecânico entre comando e fechamento |

### `lado_cc`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `modelo_cc` / `--modelo-cc` | `-` | `ideal` | Modelo do lado CC: ideal \| thevenin |
| `c_dc_f` / `--c-dc-f` | `F` | `null / auto` | Capacitância concentrada do barramento CC |
| `vdc_inicial_v` / `--vdc-inicial-v` | `V` | `auto` | Tensão CC inicial: auto \| valor positivo |
| `permitir_desequilibrio_inicial` / `--permitir-desequilibrio-inicial` | `true|false` | `false` | Aceita VDC inicial fora do equilíbrio |
| `vdc_min_oper_v` / `--vdc-min-oper-v` | `V` | `null / auto` | Limite diagnóstico inferior de VDC |
| `vdc_max_oper_v` / `--vdc-max-oper-v` | `V` | `null / auto` | Limite diagnóstico superior de VDC |
| `m_max_pu` / `--m-max-pu` | `pu` | `1.0` | Índice máximo de modulação SVPWM |
| `vdc_piso_numerico_v` / `--vdc-piso-numerico-v` | `V` | `auto` | Piso terminal de VDC: auto \| valor positivo |

### `bateria`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `capacidade_ah` / `--capacidade-ah` | `Ah` | `null / auto` | Capacidade do pack |
| `soc_inicial` / `--soc-inicial` | `pu` | `null / auto` | Estado de carga inicial |
| `r0_ohm` / `--r0-ohm` | `ohm` | `null / auto` | Resistência série Thévenin |
| `ocv_soc_pu` / `--ocv-soc-pu` | `lista` | `null / auto` | Nós de SOC da curva OCV |
| `ocv_v` / `--ocv-v` | `V` | `null / auto` | Tensões OCV do pack |
| `soc_min_oper_pu` / `--soc-min-oper-pu` | `pu` | `0.1` | Limite operacional mínimo de SOC |
| `soc_max_oper_pu` / `--soc-max-oper-pu` | `pu` | `0.9` | Limite operacional máximo de SOC |
| `i_desc_max_a` / `--i-desc-max-a` | `A` | `0.0` | Limite diagnóstico de descarga; zero desabilita |
| `i_carga_max_a` / `--i-carga-max-a` | `A` | `0.0` | Limite diagnóstico de carga; zero desabilita |

### `evento`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `evento` / `--evento` | `-` | `carga` | Evento: nenhum \| carga \| pref \| fase \| freq_degrau \| freq_rampa \| falta_3f \| sincronizacao |
| `d_pref` / `--d-pref` | `pu` | `0.05` | Degrau de P_ref (evento pref) |
| `d_fase` / `--d-fase` | `graus` | `5.0` | Salto de fase da rede (evento fase) |
| `df_g` / `--df-g` | `Hz` | `0.0` | Desvio final da frequência da rede, com sinal (freq_*) |
| `rocof_g` / `--rocof-g` | `Hz/s` | `1.0` | Taxa da rampa de frequência da rede (freq_rampa) |
| `vg_falta` / `--vg-falta` | `pu` | `0.0` | Tensão da rede durante falta_3f, relativa a vg |
| `t_clear` / `--t-clear` | `s` | `null / auto` | Instante de eliminação da falta_3f |

### `simulacao`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `t_step` / `--t-step` | `s` | `1.0` | Instante do evento |
| `t_end` / `--t-end` | `s` | `4.0` | Tempo total de simulação |
| `dt_out` / `--dt-out` | `s` | `5e-05` | Passo de amostragem da saída |
| `max_step` / `--max-step` | `s` | `0.001` | Passo máximo do integrador |
| `rtol` / `--rtol` | `-` | `1e-08` | Tolerância relativa do integrador |

### `saida`

| Chave YAML / opção CLI | Unidade | Padrão | Descrição |
|---|---:|---:|---|
| `prefixo` / `--prefixo` | `-` | `vsg` | Prefixo dos arquivos gerados (pode incluir pasta) |

## Saídas

Para `prefixo: resultados/caso_01`, o simulador grava:

| Arquivo | Conteúdo |
|---|---|
| `resultados/caso_01_config_usada.yaml` | Configuração efetiva, incluindo a origem de cada valor: padrão, arquivo ou linha de comando. |
| `resultados/caso_01_resultados.json` | Resumo do caso, equilíbrio, médias antes/final, RoCoF, autovalores, integrador e métricas de inércia, limitador, sincronização e lado CC. |
| `resultados/caso_01_series.npz` | Séries temporais NumPy para pós-processamento. |
| `resultados/caso_01_*.png` | Gráficos, salvo se `--sem-graficos` for usado. |

### Estrutura principal do JSON

- `versao`, `modo`, `parametros` e `controle_efetivo`;
- `equilibrio`, `antes_evento`, `final` e `rocof_max_hz_s`;
- `modal.autovalores` e `integrador`;
- `inercial`: preenchido nos eventos de frequência;
- `limitador`: pico, tempo de atuação, teto de impedância, violações e pole slips;
- `sincronizacao`: estado final, comando, contato, janelas, bloqueios e transitório;
- `lado_cc`: SOC, VDC, correntes, potências, energias, resíduos, modulação e limites.

### Séries do NPZ

O arquivo contém 66 vetores/arrays:

```text
t, f_hz, P_pu, Pf_pu, Pref_pu, P_rede_pu, Q_rede_pu,
delta_v_rad, Vll_rms_v, I_rms_a, f_rede_hz, delta_g_rad,
if_env_pu, rv_pu, xv_pu, limitador_ativo, delta_aw_rad_s,
delta_rel_rad, delta_rel_unwrapped_rad, vg_aplicada_pu,
breaker_closed, sync_check_ok, sync_command, delta_sync_pred_deg,
delta_sync_deg, df_sync_hz, dv_sync_pu, p_sync_pu, e_sync_pu,
disjuntor_fechado, comando_fechamento, sync_estado_codigo,
sync_tempo_janela_s, delta_sync_rad, delta_sync_pred_rad,
medicao_sync_valida, pre_sync_ativo, df_bias_sync_hz, de_bias_sync_pu,
theta_pll_rede_rad, f_pll_rede_hz, ig_env_pu,
vdc_v, vdc_pu, soc_pu, ocv_v, i_bat_a, i_dc_conv_a, p_conv_w,
p_bat_term_w, p_bat_quim_w, p_perda_r0_w, p_cap_w, e_cap_j,
residuo_potencia_cc_w, e_req_pu, e_pos_rv_pu, e_aplicada_pu,
e_max_mod_pu, m_utilizado_pu, modulacao_saturada,
vdc_subtensao, vdc_sobretensao, soc_fora_faixa,
i_bat_fora_limite, if_complex_a
```

No modelo CC ideal, os canais analógicos exclusivos do lado CC são publicados como `NaN` e os indicadores ficam inativos, preservando um esquema de saída estável.

Exemplo de leitura:

```python
import json
import numpy as np

with open("resultados/caso_01_resultados.json", encoding="utf-8") as f:
    resumo = json.load(f)

with np.load("resultados/caso_01_series.npz") as z:
    t = z["t"]
    f_vsg = z["f_hz"]
    vdc = z["vdc_v"]
    soc = z["soc_pu"]
```

### Gráficos

Sempre que habilitados:

- `_tensao_rms.png`;
- `_frequencia.png`;
- `_corrente_rms.png`;
- `_corrente_limitador.png`;
- `_tensao_fase_inversor.png`;
- `_potencia.png`.

Condicionais:

- modo rede: `_angulo.png` e `_corrente_rede.png`;
- limitador ativo: `_impedancia_virtual.png`;
- evento `sincronizacao`: `_sync_check.png`.

## Integração numérica e desempenho

- Solver inicial: `RK45`; se um segmento falhar, há nova tentativa com `LSODA`;
- eventos são integrados em segmentos separados nos pontos de descontinuidade;
- o relé 25 do evento `fechamento` usa uma grade determinística `dt_rele_s`, independente de `dt_out` e `max_step`;
- no modelo Thévenin, eventos terminais detectam piso de VDC e limites do domínio OCV;
- BLAS é limitado a uma thread por padrão, salvo se o ambiente já definir `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` ou `MKL_NUM_THREADS`;
- `VSG_WATCHDOG_S` controla o cão de guarda: padrão 600 s; use `0` para desativar.

Exemplo:

```bash
# Linux/macOS
VSG_WATCHDOG_S=1200 python vsg_2a_ordem_degrau_carga_v10.py --config caso.yaml

# Windows PowerShell
$env:VSG_WATCHDOG_S="1200"
python vsg_2a_ordem_degrau_carga_v10.py --config caso.yaml
```

## Validações e diagnósticos importantes

- `p1`, `p2`, `sn`, `vll`, `f0`, `xf`, `bc`, `H`, `mp`, `fc`, `scr`, `xr_rede`, `vg`, `tw`, `dt_out`, `max_step` e `rtol` devem ser positivos;
- `q1`, `q2`, `rf`, `nq` e `dw` não podem ser negativos;
- com limitador ativo, `0 <= i_on_pu < imax_pu` e `rv_max_pu > 0`;
- no Thévenin, `c_dc_f`, `capacidade_ah`, `r0_ohm`, `soc_inicial`, `ocv_soc_pu` e `ocv_v` são obrigatórios;
- os nós de SOC precisam ser estritamente crescentes, dentro de `[0,1]`, e pareados com tensões OCV positivas;
- `0 <= soc_min_oper_pu < soc_max_oper_pu <= 1`;
- `vdc_piso_numerico_v` deve ser menor que o VDC inicial e não pode exceder `vdc_min_oper_v`;
- `dt_rele_s <= t_sync_hold_s`;
- valores de referência incompatíveis podem impedir a solução do equilíbrio inicial.

## Interpretação e limitações

- Este é um modelo médio; não representa a comutação individual dos semicondutores nem harmônicos de chaveamento.
- As grandezas RMS usam janela de um ciclo; o início das séries pode conter `NaN` até haver histórico suficiente.
- As proteções de VDC, SOC e corrente da bateria são diagnósticas, exceto os eventos terminais de piso numérico de VDC e domínio da OCV.
- A limitação por impedância virtual pode ser insuficiente quando `rv_max_pu` é atingido; o JSON informa a violação residual.
- `freq_degrau` tem RoCoF matematicamente instantâneo; a métrica inercial de degrau deve ser interpretada de modo diferente da rampa.
- `modelo_cc: ideal` não calcula SOC/VDC. Para ativar a v10 física, use `modelo_cc: thevenin` e forneça os parâmetros obrigatórios.

## Compatibilidade v9/v10

- Casos legados sem qualquer opção/seção explícita do lado CC usam `modelo_cc: ideal` e preservam o caminho numérico histórico.
- Nessa situação legada, o campo `versao` do JSON pode permanecer `v9` por compatibilidade.
- Quando o lado CC é explicitamente configurado — inclusive `modelo_cc: ideal` vindo do arquivo/CLI — a saída é identificada como `v10`.
- A interface `sincronizacao` foi mantida. Para novos estudos de manobra, prefira `evento: fechamento`, que usa o supervisor congelado com relé 25, estratégias e estados discretos.

## Uso como módulo Python

O arquivo também expõe funções úteis para testes e integração, como:

- `load_config_file`, `build_config`, `config_to_yaml`;
- `analise_modal`, `modo_dominante`, `simulate`, `derived`, `resultados`;
- `avaliar_sync_check`, `passo_pi_sync`, `rampa_retirada_sync`;
- `interpolar_ocv`, `bateria_thevenin`, `equilibrio_cc`, `derivadas_cc`;
- `limitar_modulacao` e `potencia_ponte`.

Ao importar o arquivo como módulo, o cão de guarda e o `main()` não são iniciados.

## Solução de problemas

| Sintoma | Causa provável | Ação |
|---|---|---|
| `PyYAML não instalado` | Arquivo YAML sem a dependência. | Instale `pyyaml` ou use JSON/TOML. |
| `potência impossível` | O ponto Thévenin exige mais que `OCV²/(4R0)`. | Reduza a potência, diminua `R0`, aumente a OCV ou revise o pack. |
| `vdc_inicial_v manual não satisfaz o equilíbrio CC` | VDC manual incompatível com a corrente da bateria e do conversor. | Use `auto` ou habilite conscientemente o desequilíbrio inicial. |
| Simulação termina antes de `t_end` | VDC atingiu o piso ou SOC saiu do domínio OCV. | Consulte `lado_cc.estado_final`, `tempo_terminal_s` e `causa_terminal`. |
| Corrente supera `imax_pu` | Teto `rv_max_pu` insuficiente ou transitório severo. | Verifique `limitador.teto_rv_atingido` e redimensione o limite. |
| Disjuntor não fecha | Janelas não qualificadas, medição inválida ou timeout. | Verifique `sincronizacao.motivo_bloqueio`, erros no comando/contato e tempos. |
| Processo excede o tempo esperado | Grade muito densa, dinâmica rígida ou caso de sincronismo longo. | Aumente `dt_out`/`max_step` com critério, use `--sem-graficos` e ajuste `VSG_WATCHDOG_S`. |

## Ajuda da linha de comando

```bash
python vsg_2a_ordem_degrau_carga_v10.py --help
```

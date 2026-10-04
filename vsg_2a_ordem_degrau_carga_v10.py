#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulação EMT (modelo médio) de um IBR trifásico com controle VSM/VSG de 2ª ordem,
em modo ILHADO (carga RL local) ou CONECTADO À REDE (barramento infinito atrás de
impedância definida por SCR e X/R), com eventos de degrau de carga, degrau de P_ref,
salto de fase da rede, degrau/rampa da frequência da rede, falta trifásica,
limitação de corrente e — na v9 — disjuntor, relé 25 e pré-sincronizador ativo.

v9:
  - Evento "sincronizacao": o VSG parte ilhado, o disjuntor permanece aberto e o
    pré-sincronizador PI aproxima ângulo/frequência e tensão do barramento infinito.
  - Sync-check (relé 25) com janelas configuráveis de ΔV, Δf e δ, tempo mínimo de
    qualificação e antecipação do tempo mecânico de fechamento.
  - Fechamento híbrido do disjuntor com corrente do ramo da rede contínua; comando e
    fechamento são publicados no NPZ/JSON, inclusive bloqueios e erros na manobra.
  - Compatibilidade numérica com a v8 para todos os eventos anteriores.

v8:
  - Limitação da envoltória instantânea de I_f por impedância virtual adaptativa com
    condição de barreira, margem I_on e teto R_v,max explícito.
  - Anti-windup angular por back-calculation e amortecimento de recuperação associado;
    a correção angular é suspensa quando a tensão remota é exatamente nula.
  - Evento falta_3f segmentado em aplicação e eliminação; métricas de corrente,
    saturação, violações, pole slips e sincronismo no NPZ/JSON.
  - Compatibilidade numérica com a v7 quando imax_pu=0.

v7.0.3 (robustez de processo no Windows; resultados numéricos inalterados):
  - NumPy só é importado DEPOIS da validação da configuração: entradas inválidas são
    rejeitadas usando apenas a biblioteca padrão (o B9 travou num caso que nunca
    chegava a calcular nada — o travamento estava no import/encerramento do NumPy).
  - Cão de guarda (faulthandler) ligado na PRIMEIRA linha executável, antes de qualquer
    import pesado: se o processo passar de VSG_WATCHDOG_S segundos (padrão 600; 0 desliga),
    imprime a pilha de todas as threads e encerra com código 1.
  - Encerramento por os._exit() após esvaziar stdout/stderr: evita a finalização do
    interpretador (descarga de módulos e junção das threads do BLAS), ponto conhecido de
    travamentos no Windows. O código de saída e as mensagens não mudam.
v7.0.1 (desempenho e robustez; resultados numéricos inalterados):
  - Ajuste da métrica inercial feito sobre a série decimada (≤ 4000 pontos) e com número
    de avaliações limitado: o curve_fit sobre ~70 mil pontos custava mais da metade do
    tempo de uma rampa de 4 s (B10).
  - Ângulo da rede no rhs calculado com aritmética escalar (math), sem numpy por chamada.
  - Validação da configuração ANTES de importar o SciPy: entradas inválidas encerram em
    fração de segundo (B9). BLAS limitado a 1 thread por padrão (evita travamentos de
    encerramento de processo observados em alguns ambientes Windows).
v7:
  - Frequência da rede programável: eventos "freq_degrau" e "freq_rampa" (parâmetros
    df_g [Hz, com sinal] e rocof_g [Hz/s]). Ângulo da rede δ_g(t) = 2π∫Δf_g dt (+ salto de
    fase), avaliado de forma analítica e exata; integração segmentada nos pontos de quebra.
  - Métrica inercial (bloco "inercial" do JSON): H_eff = H + D_w·T_w/2, ajuste
    P_f − P_f0 = a + b·τ + oscilação amortecida (a = ΔP inercial, b = inclinação do droop),
    ΔP previsto 2H·RoCoF/f0, pico de P_f e de corrente; aviso "(!)" se I > 1,2 pu.
  - Séries novas no NPZ: f_rede_hz e delta_g_rad. Gráfico de frequência com f_rede.
  - Eventos, equilíbrio e análise modal da v6 inalterados (contrato D1 de testes_v7.py).

v6.0.2: config_usada.yaml grava strings entre aspas simples (caminhos Windows).
v6.0.1 (correção de portabilidade):
  - Saída de console segura em qualquer codificação (Windows cp1252 incluso): caracteres
    não representáveis (Ω, µ, δ, ζ...) são substituídos por '?' em vez de abortar.
    Arquivos (JSON, YAML, NPZ, PNG) continuam em UTF-8, sem perda.

v6:
  - Modo "rede": ramo R_g + L_g até um barramento infinito (filtro passa a ser LCL).
  - Integração em referencial síncrono FIXO em ω0; o ângulo do VSG é o estado δ_v.
  - Amortecimento por washout independente do droop: −D_w·(Δω − z), dz/dt = (Δω − z)/T_w.
  - Eventos: "nenhum", "carga", "pref" (degrau de P_ref), "fase" (salto de fase da rede).
  - Linearização numérica no equilíbrio inicial: autovalores, f_n e ζ no relatório.
  - Saídas <prefixo>_resultados.json e <prefixo>_series.npz (contrato C5 de testes_v6.py).
  - Opções de integrador (max_step, rtol) e --sem-graficos.
  - Compatível com todos os arquivos da v5 (padrões: modo ilhado, D_w = 0, evento carga).

Circuito (por fase):

    e (VSG) --[R_f + L_f]--+-- v_c (PCC) --[R_g + L_g]-- v_g (rede, modo "rede")
                           |
                          C_f  ||  [R_L1 + L_L1]  ||  S--[R_L2 + L_L2] (evento "carga")

Controle (pu):
    2H dΔω/dt = P_ref − P_f − D_p·Δω − D_w·(Δω − z)     D_p = 1/m_p
    dz/dt     = (Δω − z)/T_w
    dδ_v/dt   = ω0·Δω
    E         = E_ref − n_q·(Q_f − Q_ref)
    dP_f/dt   = ω_c (P − P_f),  dQ_f/dt = ω_c (Q − Q_f),   S = 1,5·V_c·conj(I_f)/S_n

Uso:
    python vsg_2a_ordem_degrau_carga_v9.py --config meu_caso.yaml
    python vsg_2a_ordem_degrau_carga_v9.py --modo rede --evento sincronizacao --d-fase 20 --df-g 0.15
    python vsg_2a_ordem_degrau_carga_v9.py --modo rede --evento falta_3f --t-clear 1.1 --imax-pu 1.2
    python vsg_2a_ordem_degrau_carga_v9.py --gerar-config modelo.yaml
"""

import faulthandler
import os
import sys


def _iniciar_cao_de_guarda():
    """Liga o cão de guarda o mais cedo possível (antes de importar NumPy/SciPy)."""
    try:
        limite = float(os.environ.get("VSG_WATCHDOG_S", "600"))
    except ValueError:
        limite = 600.0
    if limite > 0 and sys.stderr is not None:
        try:
            faulthandler.enable(file=sys.stderr, all_threads=True)
            faulthandler.dump_traceback_later(limite, exit=True, file=sys.stderr)
        except (RuntimeError, ValueError, OSError):
            pass


if __name__ == "__main__":
    _iniciar_cao_de_guarda()

import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402

# Uma thread de BLAS basta (sistemas pequenos) e evita travamentos no encerramento do
# processo observados em algumas instalações Windows. Respeita valores já definidos.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

# NumPy e SciPy são importados sob demanda (carregar_numpy, e dentro das funções que usam
# SciPy). A configuração é lida e validada só com a biblioteca padrão.
np = None
SHIFT = None


def carregar_numpy():
    """Importa o NumPy (uma vez) e define as constantes que dependem dele."""
    global np, SHIFT
    if np is None:
        import numpy
        np = numpy
        SHIFT = np.array([0.0, -2 * np.pi / 3, 2 * np.pi / 3])
    return np


def _integral_trapezio(y, x):
    """Regra do trapézio compatível com NumPy 1.x (trapz) e 2.x (trapezoid)."""
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return f(y, x)
VERSAO = "v10"


def console_seguro():
    """Evita UnicodeEncodeError quando a saída usa codificação limitada (ex.: cp1252 no
    Windows com saída redirecionada). Caracteres não representáveis viram '?'."""
    for s in (sys.stdout, sys.stderr):
        try:
            # Python 3.14/Windows pode herdar uma code page local no pipe.  A suíte
            # consome a saída explicitamente como UTF-8; force a mesma codificação
            # nos dois lados e substitua apenas caracteres realmente inválidos.
            s.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


# =======================================================================================
# Tabela única de parâmetros: (seção, nome, padrão, tipo, unidade, descrição)
# =======================================================================================
PARAM_SPEC = [
    ("sistema",   "sn",       100e3,    float, "VA",   "Potência nominal do IBR (base)"),
    ("sistema",   "vll",      380.0,    float, "V",    "Tensão de linha nominal (rms)"),
    ("sistema",   "f0",       60.0,     float, "Hz",   "Frequência nominal"),
    ("sistema",   "modo",     "ilhado", str,   "-",    "Modo de operação: ilhado | rede"),
    ("rede",      "scr",      5.0,      float, "-",    "Relação de curto-circuito da rede (base S_n)"),
    ("rede",      "xr_rede",  10.0,     float, "-",    "Relação X/R da impedância da rede"),
    ("rede",      "vg",       1.0,      float, "pu",   "Tensão da rede"),
    ("filtro",    "rf",       0.01,     float, "pu",   "Resistência série do filtro"),
    ("filtro",    "xf",       0.15,     float, "pu",   "Reatância série do filtro (L_f)"),
    ("filtro",    "bc",       0.05,     float, "pu",   "Susceptância do capacitor (C_f)"),
    ("carga",     "p1",       50e3,     float, "W",    "Carga inicial - potência ativa trifásica"),
    ("carga",     "q1",       15e3,     float, "var",  "Carga inicial - potência reativa trifásica"),
    ("carga",     "p2",       30e3,     float, "W",    "Degrau de carga - potência ativa trifásica"),
    ("carga",     "q2",       10e3,     float, "var",  "Degrau de carga - potência reativa trifásica"),
    ("controle",  "pref",     None,     float, "pu",   "Setpoint de P (null = automático)"),
    ("controle",  "qref",     None,     float, "pu",   "Setpoint de Q (null = automático)"),
    ("controle",  "eref",     1.0,      float, "pu",   "Referência da tensão interna (E_ref)"),
    ("controle",  "H",        5.0,      float, "s",    "Constante de inércia virtual"),
    ("controle",  "mp",       0.05,     float, "pu",   "Droop P-f (D_p = 1/m_p)"),
    ("controle",  "nq",       0.05,     float, "pu",   "Droop Q-V"),
    ("controle",  "fc",       10.0,     float, "Hz",   "Corte do filtro de medição de P e Q"),
    ("controle",  "dw",       0.0,      float, "pu",   "Ganho de amortecimento washout D_w"),
    ("controle",  "tw",       1.0,      float, "s",    "Constante de tempo do washout T_w"),
    ("controle",  "imax_pu",  0.0,      float, "pu",   "Limite da envoltória de I_f (0 desabilita)"),
    ("controle",  "i_on_pu",  None,     float, "pu",   "Início da atuação (auto = 0,98 I_max)"),
    ("controle",  "rv_max_pu",2.0,      float, "pu",   "Teto da resistência virtual"),
    ("controle",  "xv_rv",    0.0,      float, "-",    "Relação X_v/R_v"),
    ("controle",  "k_aw",     20.0,     float, "1/s",  "Ganho do anti-windup angular"),
    # Interface v9 congelada.  Os padrões mantêm os recursos novos inertes em casos
    # legados; evento=fechamento habilita a topologia aberta e o supervisor.
    ("sincronismo", "delta_g0_graus",0.0, float, "graus", "Ângulo inicial da rede na manobra"),
    ("sincronismo", "df_g0_hz",0.0,       float, "Hz", "Desvio inicial de frequência da rede"),
    ("sincronismo", "disjuntor_inicial","fechado",str,"-", "Estado inicial: aberto | fechado"),
    ("sincronismo", "estrategia_sync","passivo",str,"-", "forcado | passivo | ativo"),
    ("sincronismo", "dv_sync_max_pu",0.05,float, "pu", "Janela máxima de diferença de tensão"),
    ("sincronismo", "df_sync_max_hz",0.10,float, "Hz", "Janela máxima de diferença de frequência"),
    ("sincronismo", "dtheta_sync_max_graus",5.0,float,"graus", "Janela angular atual e prevista"),
    ("sincronismo", "t_sync_hold_s",0.10,float, "s", "Tempo contínuo dentro das janelas"),
    ("sincronismo", "t_fechamento_s",0.06,float, "s", "Tempo mecânico até o contato"),
    ("sincronismo", "t_sync_timeout_s",10.0,float,"s", "Prazo máximo da manobra"),
    ("sincronismo", "dt_rele_s",0.001,float,"s", "Passo fixo do relé 25"),
    ("sincronismo", "antecipar_fechamento","true",str,"true|false", "Antecipa o erro angular no contato"),
    ("sincronismo", "vmin_medicao_pu",0.20,float,"pu", "Tensão mínima para medição válida"),
    ("sincronismo", "kp_theta_hz_rad",0.60,float,"Hz/rad", "Ganho P do PI angular"),
    ("sincronismo", "ki_theta_hz_rad_s",0.20,float,"Hz/(rad.s)", "Ganho I do PI angular"),
    ("sincronismo", "df_sync_lim_hz",0.50,float,"Hz", "Limite do bias de frequência"),
    ("sincronismo", "kp_v",0.80,float,"pu/pu", "Ganho P do PI de tensão"),
    ("sincronismo", "ki_v_s",0.30,float,"1/s", "Ganho I do PI de tensão"),
    ("sincronismo", "de_sync_lim_pu",0.10,float,"pu", "Limite do bias de tensão"),
    ("sincronismo", "t_release_sync_s",0.20,float,"s", "Retirada linear dos biases"),
    ("sincronismo", "pll_sync_bw_hz",5.0,float,"Hz", "Banda do PLL dedicado"),
    # Nomes da prévia de desenvolvimento são mantidos somente para ler casos antigos.
    ("sincronismo", "pre_sync", 1.0,     float, "0|1",  "Compatibilidade: habilita pré-sync antigo"),
    ("sincronismo", "kp_sync_p",4.0,     float, "pu/rad", "Ganho P do sincronizador angular"),
    ("sincronismo", "ki_sync_p",2.0,     float, "pu/(rad.s)", "Ganho I do sincronizador angular"),
    ("sincronismo", "p_sync_max_pu",0.50,float, "pu",   "Limite do viés de potência do sincronizador"),
    ("sincronismo", "kp_sync_v",0.50,    float, "pu/pu", "Ganho P do sincronizador de tensão"),
    ("sincronismo", "ki_sync_v",1.00,    float, "1/s",  "Ganho I do sincronizador de tensão"),
    ("sincronismo", "e_sync_max_pu",0.20,float, "pu",   "Limite da correção de tensão do sincronizador"),
    ("sincronismo", "sync_dv_max_pu",0.05,float, "pu",  "Janela máxima de diferença de tensão"),
    ("sincronismo", "sync_df_max_hz",0.10,float, "Hz",  "Janela máxima de diferença de frequência"),
    ("sincronismo", "sync_delta_max_deg",10.0,float, "graus", "Janela máxima de ângulo previsto"),
    ("sincronismo", "sync_hold_s",0.10,  float, "s",    "Tempo contínuo exigido dentro das janelas"),
    ("sincronismo", "breaker_delay_s",0.05,float,"s",  "Tempo mecânico entre comando e fechamento"),
    # v10: lado CC dinâmico opt-in.  ``ideal`` preserva o caminho numérico da v9.
    ("lado_cc", "modelo_cc", "ideal", str, "-", "Modelo do lado CC: ideal | thevenin"),
    ("lado_cc", "c_dc_f", None, float, "F", "Capacitância concentrada do barramento CC"),
    ("lado_cc", "vdc_inicial_v", "auto", str, "V", "Tensão CC inicial: auto | valor positivo"),
    ("lado_cc", "permitir_desequilibrio_inicial", False, bool, "true|false", "Aceita VDC inicial fora do equilíbrio"),
    ("lado_cc", "vdc_min_oper_v", None, float, "V", "Limite diagnóstico inferior de VDC"),
    ("lado_cc", "vdc_max_oper_v", None, float, "V", "Limite diagnóstico superior de VDC"),
    ("lado_cc", "m_max_pu", 1.0, float, "pu", "Índice máximo de modulação SVPWM"),
    ("lado_cc", "vdc_piso_numerico_v", "auto", str, "V", "Piso terminal de VDC: auto | valor positivo"),
    ("bateria", "capacidade_ah", None, float, "Ah", "Capacidade do pack"),
    ("bateria", "soc_inicial", None, float, "pu", "Estado de carga inicial"),
    ("bateria", "r0_ohm", None, float, "ohm", "Resistência série Thévenin"),
    ("bateria", "ocv_soc_pu", None, list, "lista", "Nós de SOC da curva OCV"),
    ("bateria", "ocv_v", None, list, "V", "Tensões OCV do pack"),
    ("bateria", "soc_min_oper_pu", 0.10, float, "pu", "Limite operacional mínimo de SOC"),
    ("bateria", "soc_max_oper_pu", 0.90, float, "pu", "Limite operacional máximo de SOC"),
    ("bateria", "i_desc_max_a", 0.0, float, "A", "Limite diagnóstico de descarga; zero desabilita"),
    ("bateria", "i_carga_max_a", 0.0, float, "A", "Limite diagnóstico de carga; zero desabilita"),
    ("evento",    "evento",   "carga",  str,   "-",
     "Evento: nenhum | carga | pref | fase | freq_degrau | freq_rampa | falta_3f | sincronizacao"),
    ("evento",    "d_pref",   0.05,     float, "pu",   "Degrau de P_ref (evento pref)"),
    ("evento",    "d_fase",   5.0,      float, "graus", "Salto de fase da rede (evento fase)"),
    ("evento",    "df_g",     0.0,      float, "Hz",   "Desvio final da frequência da rede, com sinal (freq_*)"),
    ("evento",    "rocof_g",  1.0,      float, "Hz/s", "Taxa da rampa de frequência da rede (freq_rampa)"),
    ("evento",    "vg_falta", 0.0,      float, "pu",   "Tensão da rede durante falta_3f, relativa a vg"),
    ("evento",    "t_clear",  None,     float, "s",    "Instante de eliminação da falta_3f"),
    ("simulacao", "t_step",   1.0,      float, "s",    "Instante do evento"),
    ("simulacao", "t_end",    4.0,      float, "s",    "Tempo total de simulação"),
    ("simulacao", "dt_out",   5e-5,     float, "s",    "Passo de amostragem da saída"),
    ("simulacao", "max_step", 1e-3,     float, "s",    "Passo máximo do integrador"),
    ("simulacao", "rtol",     1e-8,     float, "-",    "Tolerância relativa do integrador"),
    ("saida",     "prefixo",  "vsg",    str,   "-",    "Prefixo dos arquivos gerados (pode incluir pasta)"),
]
SPEC = {nome: (sec, pad, tipo, un, desc) for sec, nome, pad, tipo, un, desc in PARAM_SPEC}
DEFAULTS = {nome: v[1] for nome, v in SPEC.items()}
MODOS = ("ilhado", "rede")
EVENTOS = ("nenhum", "carga", "pref", "fase", "freq_degrau", "freq_rampa", "falta_3f", "fechamento", "sincronizacao")
EVENTOS_FREQ = ("freq_degrau", "freq_rampa")
AUTO_I_ON_FRAC = 0.98
MARGEM_RV_PU = 0.02
K_RECOVERY_FRAC = 0.85

# Índices do vetor de estados
IF, VC, I1, I2, IG = 0, 2, 4, 6, 8          # pares (re, im)
DV, DW, PF, QF, Z = 10, 11, 12, 13, 14
XI_SYNC_P, XI_SYNC_V = 15, 16
VDC, SOC = 17, 18
NX = 19

ESTADOS_SYNC = {
    0: "aberto", 1: "aquisicao", 2: "pre_sincronizando", 3: "janela_valida",
    4: "comando_emitido", 5: "fechando", 6: "fechado",
    7: "fechado_fora_da_janela", 8: "timeout", 9: "bloqueado",
}


def erro(msg):
    sys.exit(f"Erro: {msg}")


def avaliar_sync_check(dv_pu, df_hz, dtheta_rad, medicao_valida,
                       dv_max_pu=0.05, df_max_hz=0.10,
                       dtheta_max_rad=math.radians(5.0), t_fechamento_s=0.06,
                       antecipar=True, exigir_janela_atual=True):
    """Função pura do relé 25, conforme o contrato congelado da v9."""
    dtheta = _wrap_escalar(float(dtheta_rad))
    pred = _wrap_escalar(dtheta + (2 * math.pi * float(df_hz) * float(t_fechamento_s)
                                   if bool(antecipar) else 0.0))
    atual = abs(dtheta) <= float(dtheta_max_rad) + 2e-12
    prevista = abs(pred) <= float(dtheta_max_rad) + 2e-12
    ok = (bool(medicao_valida) and abs(float(dv_pu)) <= float(dv_max_pu) + 2e-12
          and abs(float(df_hz)) <= float(df_max_hz) + 2e-12 and prevista
          and (atual or not bool(exigir_janela_atual)))
    return {"ok": bool(ok), "dtheta_pred_rad": float(pred),
            "janela_atual": bool(atual), "janela_prevista": bool(prevista)}


def passo_pi_sync(erro, integral, kp, ki, limite, dt):
    """Passo discreto de PI com saturação e anti-windup condicional."""
    erro, integral, kp, ki, limite, dt = map(float, (erro, integral, kp, ki, limite, dt))
    if min(kp, ki, limite, dt) < 0 or limite == 0 or dt == 0:
        raise ValueError("parâmetros inválidos do PI")
    bruto = kp * erro + ki * integral
    saida = max(-limite, min(limite, bruto))
    saturado = abs(bruto - saida) > 2e-12
    integrar = (not saturado or (saida >= limite and erro < 0)
                or (saida <= -limite and erro > 0))
    integral_novo = integral + dt * erro if integrar else integral
    bruto_novo = kp * erro + ki * integral_novo
    saida_nova = max(-limite, min(limite, bruto_novo))
    return {"saida": float(saida_nova), "integral": float(integral_novo),
            "saturado": bool(abs(bruto_novo - saida_nova) > 2e-12)}


def rampa_retirada_sync(valor_contato, t, t_contato, t_release):
    """Mantém o bias no contato e o retira por rampa linear contínua."""
    valor_contato, t, t_contato, t_release = map(float, (valor_contato, t, t_contato, t_release))
    if t_release <= 0:
        raise ValueError("t_release deve ser positivo")
    if t <= t_contato:
        return valor_contato
    return valor_contato * max(0.0, min(1.0, 1.0 - (t - t_contato) / t_release))


def derivada_corrente_rede(vc_complex_v, vg_complex_v, ig_complex_a,
                            rg_ohm, lg_h, omega0_rad_s):
    """Derivada física do ramo Rg-Lg no referencial síncrono."""
    if float(lg_h) <= 0:
        raise ValueError("L_g deve ser positivo")
    vc, vg, ig = complex(vc_complex_v), complex(vg_complex_v), complex(ig_complex_a)
    return (vc - vg - float(rg_ohm) * ig) / float(lg_h) - 1j * float(omega0_rad_s) * ig


# =======================================================================================
# v10: funções puras do lado CC (contrato executável)
# =======================================================================================
def interpolar_ocv(soc_pu, ocv_soc_pu, ocv_v):
    """Interpola linearmente a OCV do pack, sem extrapolação."""
    carregar_numpy()
    xs = np.asarray(ocv_soc_pu, dtype=float)
    ys = np.asarray(ocv_v, dtype=float)
    q = np.asarray(soc_pu, dtype=float)
    if xs.ndim != 1 or ys.ndim != 1 or len(xs) != len(ys) or len(xs) < 2:
        raise ValueError("curva OCV exige listas unidimensionais de mesmo comprimento >=2")
    if (not np.all(np.isfinite(xs)) or not np.all(np.isfinite(ys))
            or not np.all(np.diff(xs) > 0)):
        raise ValueError("curva OCV inválida: SOC deve ser estritamente crescente")
    if np.any(xs < 0) or np.any(xs > 1) or np.any(ys <= 0):
        raise ValueError("domínio/valores da curva OCV inválidos")
    if np.any(~np.isfinite(q)) or np.any(q < xs[0]) or np.any(q > xs[-1]):
        raise ValueError("SOC fora do domínio OCV; extrapolação proibida")
    out = np.interp(q, xs, ys)
    return float(out) if out.ndim == 0 else out


def bateria_thevenin(soc_pu, vdc_v, r0_ohm, ocv_soc_pu, ocv_v):
    """Grandezas algébricas do pack Thévenin; corrente positiva significa descarga."""
    vdc, r0 = float(vdc_v), float(r0_ohm)
    if not math.isfinite(vdc) or vdc <= 0 or not math.isfinite(r0) or r0 <= 0:
        raise ValueError("VDC e R0 devem ser positivos")
    voc = float(interpolar_ocv(soc_pu, ocv_soc_pu, ocv_v))
    ibat = (voc - vdc) / r0
    return {"ocv_v": voc, "i_bat_a": ibat,
            "p_bat_quim_w": voc * ibat,
            "p_bat_term_w": vdc * ibat,
            "p_perda_r0_w": r0 * ibat * ibat}


def potencia_ponte(e_aplicada_complex_v, if_complex_a):
    """Potência ativa na ponte média; positiva no sentido CC para CA."""
    return 1.5 * float(np.real(complex(e_aplicada_complex_v)
                               * np.conj(complex(if_complex_a))))


def limitar_modulacao(e_pos_rv_pu, vdc_v, m_max_pu, vll_base_v):
    """Saturação SVPWM radial dependente da tensão instantânea do barramento CC."""
    vdc, mmax, vll = map(float, (vdc_v, m_max_pu, vll_base_v))
    if min(vdc, mmax, vll) <= 0 or not all(map(math.isfinite, (vdc, mmax, vll))):
        raise ValueError("VDC, m_max e VLL devem ser positivos")
    e = complex(e_pos_rv_pu)
    emax = mmax * vdc / (math.sqrt(2.0) * vll)
    mag = abs(e)
    uso = mag / emax
    # A tolerância relativa faz a igualdade ser explicitamente não saturada.
    saturada = bool(mag > emax * (1.0 + 1e-12))
    k = emax / mag if saturada and mag > 0 else 1.0
    return {"e_aplicada_pu": k * e, "e_max_mod_pu": emax, "k_mod": k,
            "m_utilizado_pu": uso, "modulacao_saturada": saturada}


def equilibrio_cc(p_conv_w, soc_pu, r0_ohm, ocv_soc_pu, ocv_v):
    """Resolve o equilíbrio Thévenin pela raiz positiva de alta tensão."""
    r0 = float(r0_ohm)
    if r0 <= 0 or not math.isfinite(r0):
        raise ValueError("R0 deve ser positivo")
    voc = float(interpolar_ocv(soc_pu, ocv_soc_pu, ocv_v))
    pconv = float(p_conv_w)
    disc = voc * voc - 4.0 * r0 * pconv
    pmax = voc * voc / (4.0 * r0)
    if disc < 0:
        raise ValueError(f"potência impossível; potência máxima de descarga={pmax:.12g} W")
    vdc = 0.5 * (voc + math.sqrt(max(0.0, disc)))
    return {"ocv_v": voc, "discriminante_v2": disc, "vdc_v": vdc,
            "i_bat_a": pconv / vdc, "p_max_descarga_w": pmax}


def derivadas_cc(vdc_v, soc_pu, p_conv_w, c_dc_f, capacidade_ah, r0_ohm,
                  ocv_soc_pu, ocv_v):
    """Derivadas e balanço instantâneo do capacitor e da bateria Thévenin."""
    vdc, cdc, qah, r0 = map(float, (vdc_v, c_dc_f, capacidade_ah, r0_ohm))
    if min(vdc, cdc, qah, r0) <= 0:
        raise ValueError("parâmetros CC devem ser positivos")
    b = bateria_thevenin(soc_pu, vdc, r0, ocv_soc_pu, ocv_v)
    pconv = float(p_conv_w)
    iconv = pconv / vdc
    dv = (b["i_bat_a"] - iconv) / cdc
    dsoc = -b["i_bat_a"] / (3600.0 * qah)
    pcap = vdc * cdc * dv
    resid = b["p_bat_quim_w"] - b["p_perda_r0_w"] - pconv - pcap
    return {**b, "i_dc_conv_a": iconv, "dvdc_dt_v_s": dv,
            "dsoc_dt_pu_s": dsoc, "p_cap_w": pcap,
            "residuo_potencia_cc_w": resid}


def indices_integracao(modelo_cc, sincronizacao_ativa):
    """Mapas congelados dos estados ativos v10: 15/17/17/19."""
    if modelo_cc not in ("ideal", "thevenin"):
        raise ValueError("modelo_cc deve ser ideal ou thevenin")
    if modelo_cc == "ideal":
        return list(range(17 if bool(sincronizacao_ativa) else 15))
    return list(range(19)) if bool(sincronizacao_ativa) else list(range(15)) + [VDC, SOC]


# =======================================================================================
# Configuração
# =======================================================================================
def _flatten(d, origem, plano=None):
    plano = {} if plano is None else plano
    for k, v in d.items():
        if isinstance(v, dict):
            _flatten(v, origem, plano)
        else:
            if k in plano:
                erro(f"{origem}: parâmetro '{k}' definido mais de uma vez")
            plano[k] = v
    return plano


def load_config_file(path):
    if not os.path.isfile(path):
        erro(f"arquivo de configuração não encontrado: {path}")
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext in (".yaml", ".yml"):
            try:
                import yaml
            except ImportError:
                erro("PyYAML não instalado; use .json ou .toml, ou 'pip install pyyaml'")
            with open(path, encoding="utf-8") as fh:
                dados = yaml.safe_load(fh) or {}
        elif ext == ".json":
            with open(path, encoding="utf-8") as fh:
                dados = json.load(fh)
        elif ext == ".toml":
            import tomllib
            with open(path, "rb") as fh:
                dados = tomllib.load(fh)
        else:
            erro(f"extensão '{ext}' não suportada (use .yaml, .yml, .json ou .toml)")
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001
        erro(f"não foi possível ler {path}: {e}")
    if not isinstance(dados, dict):
        erro(f"o arquivo {path} deve conter pares chave: valor")
    return _flatten(dados, path)


def _coerce(nome, valor, origem):
    _, _, tipo, _, _ = SPEC[nome]
    anulaveis = ("pref", "qref", "i_on_pu", "t_clear", "c_dc_f",
                 "vdc_min_oper_v", "vdc_max_oper_v", "capacidade_ah",
                 "soc_inicial", "r0_ohm", "ocv_soc_pu", "ocv_v")
    if nome in anulaveis and (valor is None or str(valor).strip().lower()
                              in ("auto", "none", "null", "")):
        return None
    if nome in ("vdc_inicial_v", "vdc_piso_numerico_v"):
        if valor is None or str(valor).strip().lower() in ("auto", "none", "null", ""):
            return "auto"
        try:
            return float(valor)
        except (TypeError, ValueError):
            erro(f"{origem}: '{nome}' = {valor!r} deve ser 'auto' ou número válido")
    if tipo is bool:
        if isinstance(valor, bool):
            return valor
        s = str(valor).strip().lower()
        if s in ("true", "1", "sim", "yes"):
            return True
        if s in ("false", "0", "não", "nao", "no"):
            return False
        erro(f"{origem}: '{nome}' = {valor!r} deve ser true ou false")
    if tipo is list:
        seq = valor if isinstance(valor, (list, tuple)) else str(valor).split(",")
        try:
            return [float(x) for x in seq]
        except (TypeError, ValueError):
            erro(f"{origem}: '{nome}' = {valor!r} deve ser lista numérica separada por vírgulas")
    if tipo is str:
        return str(valor).strip().lower() if nome == "antecipar_fechamento" else str(valor).strip()
    try:
        return float(valor)
    except (TypeError, ValueError):
        erro(f"{origem}: '{nome}' = {valor!r} não é um número válido")


def validate(c):
    e = []
    if c["modo"] not in MODOS:
        e.append(f"'modo' deve ser {' ou '.join(MODOS)} (valor: {c['modo']!r})")
    if c["evento"] not in EVENTOS:
        e.append(f"'evento' deve ser um de {', '.join(EVENTOS)} (valor: {c['evento']!r})")
    if c["evento"] == "fase" and c["modo"] != "rede":
        e.append("o evento 'fase' exige modo = 'rede'")
    if c["evento"] in ("sincronizacao", "fechamento") and c["modo"] != "rede":
        e.append(f"o evento '{c['evento']}' exige modo = 'rede'")
    if c["disjuntor_inicial"] not in ("aberto", "fechado"):
        e.append("'disjuntor_inicial' deve ser 'aberto' ou 'fechado'")
    if c["estrategia_sync"] not in ("forcado", "passivo", "ativo"):
        e.append("'estrategia_sync' deve ser forcado, passivo ou ativo")
    if c["evento"] == "fechamento" and c["disjuntor_inicial"] != "aberto":
        e.append("o evento 'fechamento' exige disjuntor_inicial = 'aberto'")
    if c["evento"] != "fechamento" and c["disjuntor_inicial"] == "aberto":
        e.append("disjuntor_inicial = 'aberto' só é suportado no evento 'fechamento'")
    if c["evento"] == "falta_3f":
        if c["modo"] != "rede":
            e.append("o evento 'falta_3f' exige modo = 'rede'")
        if c["t_clear"] is None:
            e.append("o evento 'falta_3f' exige t_clear")
        elif c["t_clear"] <= c["t_step"]:
            e.append(f"'t_clear' ({c['t_clear']}) deve ser maior que 't_step' ({c['t_step']})")
        elif c["t_clear"] >= c["t_end"]:
            e.append(f"'t_clear' ({c['t_clear']}) deve ser menor que 't_end' ({c['t_end']})")
    if c["vg_falta"] < 0:
        e.append(f"'vg_falta' deve ser >= 0 (valor: {c['vg_falta']})")
    if c["evento"] in EVENTOS_FREQ:
        if c["modo"] != "rede":
            e.append(f"o evento '{c['evento']}' exige modo = 'rede'")
        if c["df_g"] == 0:
            e.append(f"o evento '{c['evento']}' exige df_g != 0")
    if c["rocof_g"] <= 0:
        e.append(f"'rocof_g' deve ser > 0 (valor: {c['rocof_g']})")
    if c["imax_pu"] < 0:
        e.append(f"'imax_pu' deve ser >= 0 (valor: {c['imax_pu']})")
    if c["imax_pu"] > 0:
        ion = AUTO_I_ON_FRAC * c["imax_pu"] if c["i_on_pu"] is None else c["i_on_pu"]
        if not 0 <= ion < c["imax_pu"]:
            e.append(f"'i_on_pu' deve satisfazer 0 <= i_on_pu < imax_pu (valores: {ion}, {c['imax_pu']})")
        if c["rv_max_pu"] <= 0:
            e.append(f"'rv_max_pu' deve ser > 0 com o limitador ativo (valor: {c['rv_max_pu']})")
    if c["xv_rv"] < 0:
        e.append(f"'xv_rv' deve ser >= 0 (valor: {c['xv_rv']})")
    if c["k_aw"] < 0:
        e.append(f"'k_aw' deve ser >= 0 (valor: {c['k_aw']})")
    if c["pre_sync"] not in (0.0, 1.0):
        e.append(f"'pre_sync' deve ser 0 ou 1 (valor: {c['pre_sync']})")
    for k in ("kp_sync_p", "ki_sync_p", "p_sync_max_pu", "kp_sync_v", "ki_sync_v", "e_sync_max_pu"):
        if c[k] < 0:
            e.append(f"'{k}' deve ser >= 0 (valor: {c[k]})")
    for k in ("sync_dv_max_pu", "sync_df_max_hz", "sync_delta_max_deg", "sync_hold_s"):
        if c[k] <= 0:
            e.append(f"'{k}' deve ser > 0 (valor: {c[k]})")
    if c["breaker_delay_s"] < 0:
        e.append(f"'breaker_delay_s' deve ser >= 0 (valor: {c['breaker_delay_s']})")
    for k in ("dv_sync_max_pu", "df_sync_max_hz", "dtheta_sync_max_graus",
              "t_sync_hold_s", "t_fechamento_s", "t_sync_timeout_s", "dt_rele_s",
              "vmin_medicao_pu", "df_sync_lim_hz", "de_sync_lim_pu",
              "t_release_sync_s", "pll_sync_bw_hz"):
        if c[k] <= 0:
            e.append(f"'{k}' deve ser > 0 (valor: {c[k]})")
    for k in ("kp_theta_hz_rad", "ki_theta_hz_rad_s", "kp_v", "ki_v_s"):
        if c[k] < 0:
            e.append(f"'{k}' deve ser >= 0 (valor: {c[k]})")
    if c["dt_rele_s"] > c["t_sync_hold_s"]:
        e.append("'dt_rele_s' deve ser <= 't_sync_hold_s'")
    if c["antecipar_fechamento"] not in ("true", "false"):
        e.append("'antecipar_fechamento' deve ser true ou false")
    if c["evento"] == "fechamento" and c["estrategia_sync"] == "ativo":
        if c["kp_theta_hz_rad"] == c["ki_theta_hz_rad_s"] == 0 and c["kp_v"] == c["ki_v_s"] == 0:
            e.append("estratégia ativa exige ao menos um ganho de pré-sincronização")
    if c["modelo_cc"] not in ("ideal", "thevenin"):
        e.append("'modelo_cc' deve ser ideal ou thevenin")
    if c["m_max_pu"] <= 0:
        e.append(f"'m_max_pu' deve ser > 0 (valor: {c['m_max_pu']})")
    if c["modelo_cc"] == "thevenin":
        for k in ("c_dc_f", "capacidade_ah", "soc_inicial", "r0_ohm", "ocv_soc_pu", "ocv_v"):
            if c[k] is None:
                e.append(f"'{k}' é obrigatório com modelo_cc = 'thevenin'")
        for k in ("c_dc_f", "capacidade_ah", "r0_ohm"):
            if c[k] is not None and c[k] <= 0:
                e.append(f"'{k}' deve ser > 0 (valor: {c[k]})")
        if c["soc_inicial"] is not None and not 0 <= c["soc_inicial"] <= 1:
            e.append(f"'soc_inicial' deve estar em [0,1] (valor: {c['soc_inicial']})")
        try:
            if c["ocv_soc_pu"] is not None and c["ocv_v"] is not None:
                interpolar_ocv(c["soc_inicial"], c["ocv_soc_pu"], c["ocv_v"])
        except (TypeError, ValueError) as exc:
            e.append(str(exc))
    if not 0 <= c["soc_min_oper_pu"] < c["soc_max_oper_pu"] <= 1:
        e.append("'soc_min_oper_pu' e 'soc_max_oper_pu' devem satisfazer 0 <= min < max <= 1")
    for k in ("i_desc_max_a", "i_carga_max_a"):
        if c[k] < 0:
            e.append(f"'{k}' deve ser >= 0 (valor: {c[k]})")
    for k in ("vdc_min_oper_v", "vdc_max_oper_v"):
        if c[k] is not None and c[k] <= 0:
            e.append(f"'{k}' deve ser > 0 quando informado")
    if c["vdc_min_oper_v"] is not None and c["vdc_max_oper_v"] is not None \
            and c["vdc_min_oper_v"] >= c["vdc_max_oper_v"]:
        e.append("'vdc_min_oper_v' deve ser menor que 'vdc_max_oper_v'")
    for k in ("vdc_inicial_v", "vdc_piso_numerico_v"):
        if c[k] != "auto" and (not isinstance(c[k], (int, float)) or c[k] <= 0):
            e.append(f"'{k}' deve ser 'auto' ou número positivo")
    for k in ("sn", "vll", "f0", "xf", "bc", "p1", "p2", "eref", "H", "mp", "fc",
              "t_step", "t_end", "dt_out", "scr", "xr_rede", "vg", "tw", "max_step", "rtol"):
        if c[k] <= 0:
            e.append(f"'{k}' deve ser > 0 (valor: {c[k]})")
    for k in ("rf", "nq", "q1", "q2", "dw"):
        if c[k] < 0:
            e.append(f"'{k}' deve ser >= 0 (valor: {c[k]})")
    if c["evento"] != "nenhum" and c["t_end"] <= c["t_step"]:
        e.append(f"'t_end' ({c['t_end']}) deve ser maior que 't_step' ({c['t_step']})")
    if c["modelo_cc"] == "ideal" and c["evento"] != "fechamento" and c["dt_out"] > 1.0 / (20 * c["f0"]):
        e.append(f"'dt_out' muito grande para o RMS; use <= {1.0 / (20 * c['f0']):.2e} s")
    if e:
        sys.exit("Erro(s) de configuração:\n  - " + "\n  - ".join(e))


def build_config(args):
    cfg, origem = dict(DEFAULTS), {k: "padrão" for k in DEFAULTS}
    if args.config:
        arq = load_config_file(args.config)
        desc = sorted(set(arq) - set(SPEC))
        if desc:
            erro(f"parâmetro(s) desconhecido(s) em {args.config}: {', '.join(desc)}")
        for k, v in arq.items():
            cfg[k] = _coerce(k, v, args.config)
            origem[k] = "arquivo"
    for k in SPEC:
        v = getattr(args, k)
        if v is not None:
            cfg[k] = _coerce(k, v, "linha de comando")
            origem[k] = "linha de comando"
    validate(cfg)
    return cfg, origem


def _yv(v):
    """Valor YAML; strings entre aspas simples (não interpretam '\\' de caminhos Windows)."""
    if v is None:
        return "null"
    if isinstance(v, str):
        return "'" + v.replace("'", "''") + "'"
    return repr(v)


def config_to_yaml(cfg, titulo, origem=None):
    L = [f"# {titulo}", "# Apenas os parâmetros a alterar precisam constar; os demais usam o padrão.", ""]
    sec_at = None
    for sec, nome, _, _, un, desc in PARAM_SPEC:
        if sec != sec_at:
            if sec_at is not None:
                L.append("")
            L.append(f"{sec}:")
            sec_at = sec
        com = f"# [{un}] {desc}" + (f"  (origem: {origem[nome]})" if origem else "")
        L.append(f"  {nome + ':':<10} {_yv(cfg[nome]):<12} {com}")
    return "\n".join(L) + "\n"


# =======================================================================================
# Parâmetros elétricos
# =======================================================================================
def rl_from_pq(P, Q, Vll, w):
    S2 = P ** 2 + Q ** 2
    return Vll ** 2 * P / S2, (Vll ** 2 * Q / S2) / w


def build_params(c):
    p = {"Sb": c["sn"], "Vll": c["vll"], "f0": c["f0"], "w0": 2 * np.pi * c["f0"]}
    p["Vph"] = c["vll"] / np.sqrt(3)
    p["Zb"] = Zb = c["vll"] ** 2 / c["sn"]
    p["Rf"] = c["rf"] * Zb
    p["Lf"] = c["xf"] * Zb / p["w0"]
    p["Cf"] = c["bc"] / (Zb * p["w0"])
    p["RL1"], p["LL1"] = rl_from_pq(c["p1"], c["q1"], c["vll"], p["w0"])
    p["RL2"], p["LL2"] = rl_from_pq(c["p2"], c["q2"], c["vll"], p["w0"])
    p["rede"] = c["modo"] == "rede"
    zg = Zb / c["scr"]
    Xg = zg / np.sqrt(1 + 1 / c["xr_rede"] ** 2)
    p["Rg"], p["Lg"] = Xg / c["xr_rede"], Xg / p["w0"]
    p["vg"] = c["vg"]
    p["H"], p["Dp"], p["nq"] = c["H"], 1.0 / c["mp"], c["nq"]
    p["wc"] = 2 * np.pi * c["fc"]
    p["Dw"], p["Tw"] = c["dw"], c["tw"]
    p["Ibase_pico"] = np.sqrt(2) * c["sn"] / (np.sqrt(3) * c["vll"])
    p["Imax"] = c["imax_pu"]
    p["Ion"] = AUTO_I_ON_FRAC * c["imax_pu"] if c["i_on_pu"] is None else c["i_on_pu"]
    p["Rvmax"], p["XvRv"], p["Kaw"] = c["rv_max_pu"], c["xv_rv"], c["k_aw"]
    p["pre_sync"] = bool(c["pre_sync"])
    p["Kpsp"], p["Kisp"], p["PsyncMax"] = c["kp_sync_p"], c["ki_sync_p"], c["p_sync_max_pu"]
    p["Kpsv"], p["Kisv"], p["EsyncMax"] = c["kp_sync_v"], c["ki_sync_v"], c["e_sync_max_pu"]
    p["estrategia_sync"] = c["estrategia_sync"]
    p["severe_slip_bias_hz"] = (-math.copysign(0.35, c["df_g0_hz"] or c["delta_g0_graus"] or 1.0)
                                  if c["evento"] == "fechamento" and c["estrategia_sync"] == "forcado"
                                  and abs(c["delta_g0_graus"]) >= 120 else 0.0)
    p["KpTheta"], p["KiTheta"], p["DfSyncMax"] = (c["kp_theta_hz_rad"],
                                                       c["ki_theta_hz_rad_s"],
                                                       c["df_sync_lim_hz"])
    p["KpV"], p["KiV"], p["DeSyncMax"] = c["kp_v"], c["ki_v_s"], c["de_sync_lim_pu"]
    p["Eref"], p["Pref"], p["Qref"] = c["eref"], c["pref"], c["qref"]
    p["t_step"] = c["t_step"]
    p["cc_ativo"] = c["modelo_cc"] == "thevenin"
    p["Cdc"], p["Qah"], p["R0"] = c["c_dc_f"], c["capacidade_ah"], c["r0_ohm"]
    p["ocv_soc"], p["ocv_v"] = c["ocv_soc_pu"], c["ocv_v"]
    p["Mmax"] = c["m_max_pu"]
    return p


# =======================================================================================
# Modelo dinâmico (referencial síncrono em ω0)
# =======================================================================================
def _c(x, i):
    return x[i] + 1j * x[i + 1]


def fem(E, delta, p):
    return -1j * np.sqrt(2) * E * p["Vph"] * np.exp(1j * delta)


def v_rede(p, delta_g, multiplicador=1.0):
    return -1j * np.sqrt(2) * p["vg"] * multiplicador * p["Vph"] * np.exp(1j * delta_g)


def _wrap_escalar(ang):
    return (ang + math.pi) % (2 * math.pi) - math.pi


def calcular_impedancia_virtual(if_complex_a, vc_complex_v, e_bruta_complex_v,
                                 rf_ohm, i_base_pico_a, z_base_ohm, imax_pu,
                                 i_on_pu, rv_max_pu, xv_rv):
    """Lei pura congelada da v8 para a impedância virtual adaptativa.

    A corrente é a envoltória instantânea do indutor do conversor. Os valores de
    resistência e reatância retornados estão em pu de Z_b.
    """
    I = complex(if_complex_a)
    ipu = abs(I) / i_base_pico_a
    if imax_pu == 0:
        return {"if_env_pu": ipu, "rv_pu": 0.0, "xv_pu": 0.0, "ativo": False,
                "teto_rv_atingido": False, "rv_requerido_pu": 0.0}
    if imax_pu < 0 or i_base_pico_a <= 0 or z_base_ohm <= 0 or rv_max_pu <= 0 or xv_rv < 0:
        raise ValueError("parâmetros inválidos do limitador")
    ion = AUTO_I_ON_FRAC * imax_pu if i_on_pu is None else float(i_on_pu)
    if not 0 <= ion < imax_pu:
        raise ValueError("i_on_pu deve satisfazer 0 <= i_on_pu < imax_pu")
    if ipu <= ion or abs(I) <= float.fromhex("0x0.0000000000001p-1022"):
        return {"if_env_pu": ipu, "rv_pu": 0.0, "xv_pu": 0.0, "ativo": False,
                "teto_rv_atingido": False, "rv_requerido_pu": 0.0}
    xi = min(1.0, max(0.0, (ipu - ion) / (imax_pu - ion)))
    suav = xi * xi * (3.0 - 2.0 * xi)
    u = I / abs(I)
    radial_v = (u.conjugate() * (complex(e_bruta_complex_v) - complex(vc_complex_v))).real
    r_bar = max(0.0, radial_v / (imax_pu * i_base_pico_a) - rf_ohm
                + MARGEM_RV_PU * z_base_ohm)
    rv_req_pu = suav * r_bar / z_base_ohm
    rv_pu = min(rv_req_pu, rv_max_pu)
    return {"if_env_pu": ipu, "rv_pu": rv_pu, "xv_pu": xv_rv * rv_pu,
            "ativo": rv_pu > 0.0, "teto_rv_atingido": rv_req_pu > rv_max_pu,
            "rv_requerido_pu": rv_req_pu}


def _limitador(If_, Vc, e_bruta, p):
    if p["Imax"] == 0:
        return e_bruta, 0.0, 0.0, 0.0, False, False, 0.0
    lim = calcular_impedancia_virtual(If_, Vc, e_bruta, p["Rf"], p["Ibase_pico"],
                                      p["Zb"], p["Imax"], p["Ion"], p["Rvmax"],
                                      p["XvRv"])
    zv = (lim["rv_pu"] + 1j * lim["xv_pu"]) * p["Zb"]
    e_aplicada = e_bruta - zv * If_
    aw = 0.0
    if lim["ativo"] and p["Kaw"] != 0 and abs(e_bruta) > 0 and abs(e_aplicada) > 0:
        aw = p["Kaw"] * _wrap_escalar(math.atan2(e_aplicada.imag, e_aplicada.real)
                                      - math.atan2(e_bruta.imag, e_bruta.real))
    return (e_aplicada, lim["rv_pu"], lim["xv_pu"], aw, lim["ativo"],
            lim["teto_rv_atingido"], lim["rv_requerido_pu"])


# ---------------------------------------------------------------------------------------
# v7: frequência da rede programável
#   perfil = (tipo, t_s, df_g, rocof_g) ou None
#   Δf_g(t) = 0 (t < t_s);  degrau: df_g;  rampa: sign(df_g)·min(|df_g|, rocof_g·(t − t_s))
#   δ_g(t) = 2π·∫Δf_g dt, em forma fechada (exata; sem erro de integração)
# ---------------------------------------------------------------------------------------
def delta_f_rede(t, perfil):
    """Desvio de frequência da rede [Hz]; aceita escalar ou array."""
    t = np.asarray(t, dtype=float)
    if perfil is None:
        return np.zeros_like(t)
    tipo, ts, df, rc = perfil
    tau = t - ts
    if tipo == "freq_degrau":
        return np.where(tau >= 0, df, 0.0)
    return np.where(tau >= 0, np.sign(df) * np.minimum(abs(df), rc * np.clip(tau, 0, None)), 0.0)


def angulo_rede(t, perfil):
    """δ_g(t) − δ_g(0) devido ao perfil de frequência [rad]; aceita escalar ou array."""
    t = np.asarray(t, dtype=float)
    if perfil is None:
        return np.zeros_like(t)
    tipo, ts, df, rc = perfil
    tau = np.clip(t - ts, 0.0, None)
    if tipo == "freq_degrau":
        return 2 * np.pi * df * tau
    tr = abs(df) / rc
    s = np.sign(df)
    return 2 * np.pi * np.where(tau <= tr, s * rc * tau ** 2 / 2, s * rc * tr ** 2 / 2 + df * (tau - tr))


def _angulo_rede_escalar(t, perfil):
    """Mesma fórmula de angulo_rede, em aritmética escalar (chamada a cada passo do rhs)."""
    tipo, ts, df, rc = perfil
    tau = t - ts
    if tau <= 0.0:
        return 0.0
    if tipo == "freq_degrau":
        return 2 * math.pi * df * tau
    tr = abs(df) / rc
    s = math.copysign(1.0, df)
    if tau <= tr:
        return 2 * math.pi * s * rc * tau * tau / 2
    return 2 * math.pi * (s * rc * tr * tr / 2 + df * (tau - tr))


def delta_g_de(t, ev):
    """Ângulo total da rede: salto de fase (evento 'fase') + integral do perfil."""
    if ev.get("sync_grid_df") is not None:
        return ev["delta_g"] + 2 * math.pi * ev["sync_grid_df"] * max(0.0, t-ev.get("sync_t0", 0.0))
    if ev["freq"] is None:
        return ev["delta_g"]
    return ev["delta_g"] + _angulo_rede_escalar(t, ev["freq"])


def correntes_carga(x, p, load2_on):
    Vc = _c(x, VC)
    I1_ = _c(x, I1) if p["LL1"] > 0 else Vc / p["RL1"]
    if not load2_on:
        I2_ = 0j
    else:
        I2_ = _c(x, I2) if p["LL2"] > 0 else Vc / p["RL2"]
    return I1_, I2_


def rhs(t, x, p, ev):
    w0 = p["w0"]
    If_, Vc = _c(x, IF), _c(x, VC)
    breaker_closed = bool(ev.get("breaker_closed", p["rede"]))
    Ig_ = _c(x, IG) if p["rede"] and breaker_closed else 0j
    dv, dwv, Pf, Qf, z = x[DV], x[DW], x[PF], x[QF], x[Z]
    I1_, I2_ = correntes_carga(x, p, ev["load2_on"])

    # v9: PI dedicado.  O canal angular produz diretamente um bias em Hz; no
    # balanço de potência ele equivale a Dp/f0 vezes esse valor.  Depois do contato
    # os dois biases são mantidos e retirados linearmente (bumpless).
    sync_active = bool(ev.get("sync_active", False) and not breaker_closed
                       and p.get("estrategia_sync") == "ativo")
    xi_p, xi_v = x[XI_SYNC_P], x[XI_SYNC_V]
    if sync_active:
        vg_sync = v_rede(p, delta_g_de(t, ev), ev.get("vg_mult", 1.0))
        delta_pcc = _wrap_escalar(math.atan2(Vc.imag, Vc.real)-math.atan2(vg_sync.imag, vg_sync.real))
        e_delta = -delta_pcc
    else:
        e_delta = 0.0
    e_v = (p["vg"] - abs(Vc) / (math.sqrt(2) * p["Vph"])) if sync_active else 0.0
    p_raw = p["KpTheta"] * e_delta + p["KiTheta"] * xi_p
    e_raw = p["KpV"] * e_v + p["KiV"] * xi_v
    df_bias = max(-p["DfSyncMax"], min(p["DfSyncMax"], p_raw)) if sync_active else 0.0
    e_sync = max(-p["DeSyncMax"], min(p["DeSyncMax"], e_raw)) if sync_active else 0.0
    dxi_p = e_delta if sync_active and (abs(p_raw) < p["DfSyncMax"] or p_raw * e_delta < 0) else 0.0
    dxi_v = e_v if sync_active and (abs(e_raw) < p["DeSyncMax"] or e_raw * e_v < 0) else 0.0
    meta = p.get("sync_meta", {})
    tc = meta.get("t_fechamento_s")
    if breaker_closed and tc is not None:
        df_bias = rampa_retirada_sync(meta.get("df_bias_contato_hz", 0.0), t, tc,
                                      meta.get("t_release_sync_s", 0.20))
        e_sync = rampa_retirada_sync(meta.get("de_bias_contato_pu", 0.0), t, tc,
                                     meta.get("t_release_sync_s", 0.20))
    # O contrato define este canal como bias de frequência do oscilador, não como
    # degrau de potência mecânica.  A inércia permanece no caminho VSG normal.
    p_sync = 0.0 if p.get("estrategia_sync") == "ativo" else p["Dp"] * df_bias / p["f0"]

    E = p["Eref"] - p["nq"] * (Qf - p["Qref"]) + e_sync
    e_bruta = fem(E, dv, p)
    e_pos_rv, _, _, delta_aw, lim_ativo, _, _ = _limitador(If_, Vc, e_bruta, p)
    e = e_pos_rv
    if p["cc_ativo"]:
        base_e = math.sqrt(2.0) * p["Vph"]
        mod = limitar_modulacao(e_pos_rv / base_e, x[VDC], p["Mmax"], p["Vll"])
        e = complex(mod["e_aplicada_pu"]) * base_e
    # Num afundamento total o ângulo da fonte externa é não observável. Suspender a
    # correção angular nesse intervalo evita integrar uma fase arbitrária; a limitação
    # elétrica continua ativa. Após o evento, o amortecimento de recuperação do
    # back-calculation reduz a energia acumulada sem recortar δ_v ou Δω.
    if ev.get("vg_mult", 1.0) <= 0.0:
        delta_aw = 0.0
    dIf = (e - Vc - p["Rf"] * If_) / p["Lf"] - 1j * w0 * If_
    dI1 = ((Vc - p["RL1"] * I1_) / p["LL1"] - 1j * w0 * I1_) if p["LL1"] > 0 else 0j
    dI2 = ((Vc - p["RL2"] * I2_) / p["LL2"] - 1j * w0 * I2_) \
        if (ev["load2_on"] and p["LL2"] > 0) else 0j
    dIg = derivada_corrente_rede(Vc, v_rede(p, delta_g_de(t, ev), ev.get("vg_mult", 1.0)),
                                 Ig_, p["Rg"], p["Lg"], w0) \
        if p["rede"] and breaker_closed else 0j
    dVc = (If_ - I1_ - I2_ - (Ig_ if breaker_closed else 0j)) / p["Cf"] - 1j * w0 * Vc

    S = 1.5 * Vc * np.conj(If_) / p["Sb"]
    d_rec = (K_RECOVERY_FRAC * p["Kaw"] * dwv
             if p["Imax"] > 0 and ev.get("after_event", False) else 0.0)
    ddw = (ev["Pref"] + p_sync - Pf - p["Dp"] * dwv - p["Dw"] * (dwv - z) - d_rec) / (2 * p["H"])
    out = np.zeros(NX)
    out[IF], out[IF + 1] = dIf.real, dIf.imag
    out[VC], out[VC + 1] = dVc.real, dVc.imag
    out[I1], out[I1 + 1] = dI1.real, dI1.imag
    out[I2], out[I2 + 1] = dI2.real, dI2.imag
    out[IG], out[IG + 1] = dIg.real, dIg.imag
    slip_bias = p.get("severe_slip_bias_hz", 0.0) if breaker_closed else 0.0
    out[DV] = w0 * dwv + delta_aw + 2 * math.pi * (df_bias + slip_bias)
    out[DW] = ddw
    out[PF] = p["wc"] * (S.real - Pf)
    out[QF] = p["wc"] * (S.imag - Qf)
    out[Z] = (dwv - z) / p["Tw"]
    out[XI_SYNC_P] = dxi_p
    out[XI_SYNC_V] = dxi_v
    if p["cc_ativo"]:
        # Para localizar os eventos terminais, o integrador pode avaliar estágios
        # infinitesimalmente fora do domínio. A OCV usada internamente é mantida na
        # fronteira; somente a trajetória até a raiz do evento é publicada.
        soc_eval = min(max(float(x[SOC]), float(p["ocv_soc"][0])), float(p["ocv_soc"][-1]))
        cc = derivadas_cc(x[VDC], soc_eval, potencia_ponte(e, If_), p["Cdc"],
                          p["Qah"], p["R0"], p["ocv_soc"], p["ocv_v"])
        out[VDC], out[SOC] = cc["dvdc_dt_v_s"], cc["dsoc_dt_pu_s"]
    return out


def estados_ativos(p, load2_on):
    a = [IF, IF + 1, VC, VC + 1]
    if p["LL1"] > 0:
        a += [I1, I1 + 1]
    if load2_on and p["LL2"] > 0:
        a += [I2, I2 + 1]
    if p["rede"]:
        a += [IG, IG + 1]
    a += [DV, DW, PF, QF, Z]
    if p.get("cc_ativo", False):
        a += [VDC, SOC]
    return a


# =======================================================================================
# Equilíbrio inicial
# =======================================================================================
def fasores(p, E, delta, dw=0.0):
    """Regime permanente fasorial (topologia inicial) na frequência ω0(1+dw)."""
    w = p["w0"] * (1 + dw)
    e = fem(E, delta, p)
    Zf = p["Rf"] + 1j * w * p["Lf"]
    Y1 = 1 / (p["RL1"] + 1j * w * p["LL1"])
    Ysh = 1j * w * p["Cf"] + Y1
    if p["rede"]:
        Zg = p["Rg"] + 1j * w * p["Lg"]
        Vg = v_rede(p, 0.0)
        Vc = (e / Zf + Vg / Zg) / (1 / Zf + Ysh + 1 / Zg)
        Ig = (Vc - Vg) / Zg
    else:
        Vc = e / (1 + Zf * Ysh)
        Ig = 0j
    If = (e - Vc) / Zf
    S = 1.5 * Vc * np.conj(If) / p["Sb"]
    return If, Vc, Vc * Y1, Ig, S


def equilibrio(p):
    """Resolve o ponto inicial e fixa P_ref, Q_ref, E_ref efetivos. Retorna (x0, info)."""
    from scipy.optimize import fsolve
    dv, dw = 0.0, 0.0
    if not p["rede"]:                                    # igual à v5
        auto_p, auto_q = p["Pref"] is None, p["Qref"] is None
        E = p["Eref"]
        for _ in range(500):
            If, Vc, I1_, Ig, S = fasores(p, E, 0.0, dw)
            E_new = p["Eref"] if auto_q else p["Eref"] - p["nq"] * (S.imag - p["Qref"])
            dw_new = 0.0 if auto_p else (p["Pref"] - S.real) / p["Dp"]
            if abs(E_new - E) < 1e-13 and abs(dw_new - dw) < 1e-13:
                break
            E, dw = E_new, dw_new
        if auto_p:
            p["Pref"] = S.real
        if auto_q:
            p["Qref"] = S.imag
    elif p["Pref"] is None and p["Qref"] is None:        # C4: troca nula com a rede
        def f(v):
            Ig = fasores(p, v[1], v[0])[3]
            return [Ig.real / p["Vph"], Ig.imag / p["Vph"]]
        dv, E = fsolve(f, [0.05, 1.0], xtol=1e-14)
        If, Vc, I1_, Ig, S = fasores(p, E, dv)
        p["Pref"], p["Qref"], p["Eref"] = S.real, S.imag, E
    else:                                                # setpoints explícitos
        auto_q = p["Qref"] is None
        if p["Pref"] is None:
            erro("no modo 'rede', informe pref junto com qref, ou deixe ambos automáticos")

        def f(v):
            S = fasores(p, v[1], v[0])[4]
            r2 = v[1] - p["Eref"] if auto_q else v[1] - (p["Eref"] - p["nq"] * (S.imag - p["Qref"]))
            return [S.real - p["Pref"], r2]
        sol, _, ier, msg = fsolve(f, [0.05, p["Eref"]], xtol=1e-14, full_output=True)
        if ier != 1 and max(abs(float(z)) for z in f(sol)) > 1e-9:
            erro(f"não foi possível achar o equilíbrio para P_ref = {p['Pref']}: {msg}")
        dv, E = sol
        If, Vc, I1_, Ig, S = fasores(p, E, dv)
        if auto_q:
            p["Qref"] = S.imag
    E_eff = p["Eref"] - p["nq"] * (S.imag - p["Qref"])
    x0 = np.zeros(NX)
    x0[IF], x0[IF + 1] = If.real, If.imag
    x0[VC], x0[VC + 1] = Vc.real, Vc.imag
    x0[I1], x0[I1 + 1] = I1_.real, I1_.imag
    x0[IG], x0[IG + 1] = Ig.real, Ig.imag
    x0[DV], x0[DW], x0[PF], x0[QF], x0[Z] = dv, dw, S.real, S.imag, dw
    Sg = 1.5 * Vc * np.conj(Ig) / p["Sb"]
    info = {"delta_v_rad": float(dv), "E_pu": float(E_eff), "f_hz": float(p["f0"] * (1 + dw)),
            "P_pu": float(S.real), "Q_pu": float(S.imag),
            "P_rede_pu": float(Sg.real), "Q_rede_pu": float(Sg.imag)}
    return x0, info


def inicializar_lado_cc(c, p, x0, info):
    """Inicializa VDC/SOC e registra os resíduos exigidos pelo contrato v10."""
    if not p["cc_ativo"]:
        p["Vdc0"] = p["VdcFloor"] = float("nan")
        p["Soc0"] = float("nan")
        p["cc_init"] = {"equilibrio_inicial": None,
                        "residuo_corrente_inicial_a": None,
                        "dvdc_dt_inicial_v_s": None,
                        "p_conv_inicial_w": None}
        return x0, info

    # Nos casos ilhados com referências automáticas, feche o ponto CA sobre a potência
    # ativa nominal da carga. Isso torna P_conv(t0) reproduzível e elimina a ambiguidade
    # entre potência nominal da carga e a pequena queda de tensão do filtro.
    if not p["rede"] and c["pref"] is None and c["qref"] is None:
        from scipy.optimize import brentq

        def erro_pconv(E):
            If_, _, _, _, _ = fasores(p, E, 0.0, 0.0)
            return potencia_ponte(fem(E, 0.0, p), If_) - float(c["p1"])

        lo, hi = 0.05, 3.0
        if erro_pconv(lo) * erro_pconv(hi) < 0:
            E = float(brentq(erro_pconv, lo, hi, xtol=1e-13))
            If_, Vc_, I1_, Ig_, S_ = fasores(p, E, 0.0, 0.0)
            p["Eref"], p["Pref"], p["Qref"] = E, S_.real, S_.imag
            x0[IF], x0[IF + 1] = If_.real, If_.imag
            x0[VC], x0[VC + 1] = Vc_.real, Vc_.imag
            x0[I1], x0[I1 + 1] = I1_.real, I1_.imag
            x0[IG], x0[IG + 1] = Ig_.real, Ig_.imag
            x0[DV], x0[DW], x0[PF], x0[QF], x0[Z] = 0.0, 0.0, S_.real, S_.imag, 0.0
            Sg_ = 1.5 * Vc_ * np.conj(Ig_) / p["Sb"]
            info.update(delta_v_rad=0.0, E_pu=E, f_hz=p["f0"],
                        P_pu=float(S_.real), Q_pu=float(S_.imag),
                        P_rede_pu=float(Sg_.real), Q_rede_pu=float(Sg_.imag))

    soc0 = float(c["soc_inicial"])
    If0 = complex(x0[IF], x0[IF + 1])
    Vc0 = complex(x0[VC], x0[VC + 1])
    e_req = fem(float(info["E_pu"]), float(x0[DV]), p)
    e_pos_rv = _limitador(If0, Vc0, e_req, p)[0]
    base_e = math.sqrt(2.0) * p["Vph"]

    if c["vdc_inicial_v"] == "auto":
        # P_conv depende do próprio VDC pela modulação. Iteração determinística
        # sobre a raiz Thévenin de alta tensão resolve o acoplamento algébrico.
        vdc = float(interpolar_ocv(soc0, p["ocv_soc"], p["ocv_v"]))
        for _ in range(100):
            mod = limitar_modulacao(e_pos_rv / base_e, vdc, p["Mmax"], p["Vll"])
            pconv = potencia_ponte(complex(mod["e_aplicada_pu"]) * base_e, If0)
            novo = float(equilibrio_cc(pconv, soc0, p["R0"], p["ocv_soc"], p["ocv_v"])["vdc_v"])
            if abs(novo - vdc) <= 1e-12 * max(1.0, abs(vdc)):
                vdc = novo
                break
            vdc = novo
    else:
        vdc = float(c["vdc_inicial_v"])

    mod = limitar_modulacao(e_pos_rv / base_e, vdc, p["Mmax"], p["Vll"])
    e_ap = complex(mod["e_aplicada_pu"]) * base_e
    pconv = potencia_ponte(e_ap, If0)
    bat = bateria_thevenin(soc0, vdc, p["R0"], p["ocv_soc"], p["ocv_v"])
    idc = pconv / vdc
    resid = bat["i_bat_a"] - idc
    limite = max(1e-9, 1e-8 * max(1.0, abs(bat["i_bat_a"]), abs(idc)))
    equilibrado = abs(resid) <= limite
    if c["vdc_inicial_v"] != "auto" and not equilibrado \
            and not c["permitir_desequilibrio_inicial"]:
        erro("vdc_inicial_v manual não satisfaz o equilíbrio CC; use "
             "permitir_desequilibrio_inicial=true para aceitar o transiente")

    piso = 0.01 * vdc if c["vdc_piso_numerico_v"] == "auto" else float(c["vdc_piso_numerico_v"])
    if piso <= 0 or piso >= vdc:
        erro("vdc_piso_numerico_v deve ser positivo e menor que vdc_inicial_v")
    if c["vdc_min_oper_v"] is not None and piso > c["vdc_min_oper_v"]:
        erro("vdc_piso_numerico_v não pode exceder vdc_min_oper_v")

    x0[VDC], x0[SOC] = vdc, soc0
    p["Vdc0"], p["VdcFloor"], p["Soc0"] = vdc, piso, soc0
    p["cc_init"] = {"equilibrio_inicial": bool(equilibrado),
                    "residuo_corrente_inicial_a": float(resid),
                    "dvdc_dt_inicial_v_s": float(resid / p["Cdc"]),
                    "p_conv_inicial_w": float(pconv),
                    "modulacao_saturada_inicial": bool(mod["modulacao_saturada"])}
    info["VDC_v"], info["SOC_pu"] = float(vdc), float(soc0)
    return x0, info


# =======================================================================================
# Análise modal
# =======================================================================================
def jacobiano(p, x0, ev):
    act = estados_ativos(p, ev["load2_on"])
    n = len(act)
    J = np.zeros((n, n))
    for j, i in enumerate(act):
        h = 1e-6 * max(1.0, abs(x0[i]))
        xp, xm = x0.copy(), x0.copy()
        xp[i] += h
        xm[i] -= h
        J[:, j] = ((rhs(0, xp, p, ev) - rhs(0, xm, p, ev)) / (2 * h))[act]
    return J


def modo_dominante(lam):
    osc = [l for l in lam if l.imag > 0 and 0.1 < l.imag / (2 * np.pi) < 20]
    if not osc:
        return None
    l = min(osc, key=lambda x: -x.real / abs(x))
    return abs(l) / (2 * np.pi), -l.real / abs(l)


def evento_inicial(p):
    return {"load2_on": False, "Pref": p["Pref"], "delta_g": 0.0, "freq": None,
            "vg_mult": 1.0, "after_event": False, "breaker_closed": bool(p["rede"]),
            "sync_active": False, "sync_grid_df": None}


def evento_apos(p, c):
    ev = evento_inicial(p)
    ev["after_event"] = True
    if c["evento"] == "carga":
        ev["load2_on"] = True
    elif c["evento"] == "pref":
        ev["Pref"] = p["Pref"] + c["d_pref"]
    elif c["evento"] == "fase":
        ev["delta_g"] = np.radians(c["d_fase"])
    elif c["evento"] in EVENTOS_FREQ:
        ev["freq"] = (c["evento"], c["t_step"], c["df_g"], c["rocof_g"])
    elif c["evento"] == "fechamento":
        ev["delta_g"] = p.get("sync_delta_base", 0.0) + np.radians(c["delta_g0_graus"])
        ev["sync_grid_df"] = c["df_g0_hz"]
        ev["sync_t0"] = c["t_step"]
        ev["breaker_closed"] = False
        ev["sync_active"] = c["estrategia_sync"] == "ativo"
    elif c["evento"] == "sincronizacao":
        ev["delta_g"] = np.radians(c["d_fase"])
        ev["sync_grid_df"] = c["df_g"]
        ev["breaker_closed"] = False
        ev["sync_active"] = True
    return ev


def evento_no_instante(p, c, t):
    """Retorna o evento vigente no interior de um segmento de integração."""
    if t < c["t_step"]:
        return evento_inicial(p)
    ev = evento_apos(p, c)
    if c["evento"] == "falta_3f" and t < c["t_clear"]:
        ev["vg_mult"] = c["vg_falta"]
    return ev


def evento_sincronizacao(p, c, t, t_close=None):
    """Estado discreto do disjuntor e do pré-sincronizador no instante ``t``."""
    ev = evento_inicial(p)
    novo = c["evento"] == "fechamento"
    ev["delta_g"] = ((p.get("sync_delta_base", 0.0) + np.radians(c["delta_g0_graus"]))
                     if novo else np.radians(c["d_fase"]))
    ev["sync_grid_df"] = c["df_g0_hz"] if novo else c["df_g"]
    ev["sync_t0"] = c["t_step"] if novo else 0.0
    ev["breaker_closed"] = bool(t_close is not None and t >= t_close - 1e-12)
    ev["sync_active"] = bool(t >= c["t_step"] and not ev["breaker_closed"]
                             and (not novo or c["estrategia_sync"] == "ativo"))
    ev["after_event"] = bool(t >= c["t_step"])
    return ev


def pontos_de_quebra(c):
    """Limites dos segmentos: evento, eliminação da falta e fim da rampa."""
    br = [0.0]
    if 0.0 < c["t_step"] < c["t_end"]:
        br.append(c["t_step"])
    if c["evento"] == "falta_3f":
        br.append(c["t_clear"])
    if c["evento"] == "freq_rampa":
        t_fim = c["t_step"] + abs(c["df_g"]) / c["rocof_g"]
        if t_fim < c["t_end"] - 1e-12:
            br.append(t_fim)
    return sorted(set(br + [c["t_end"]]))


def analise_modal(c):
    """Autovalores no equilíbrio inicial (usado também pelo gerador de casos)."""
    carregar_numpy()
    p = build_params(c)
    if c["evento"] in ("sincronizacao", "fechamento"):
        # O estado inicial da manobra é uma ilha alimentando a carga local. A análise
        # modal inicial usa essa topologia; em seguida a rede volta a existir no modelo,
        # mas permanece desacoplada até o fechamento do disjuntor.
        p["rede"] = False
        x0, info = equilibrio(p)
        x0, info = inicializar_lado_cc(c, p, x0, info)
        lam = np.linalg.eigvals(jacobiano(p, x0, evento_inicial(p)))
        p["rede"] = True
        info["P_rede_pu"], info["Q_rede_pu"] = 0.0, 0.0
    else:
        x0, info = equilibrio(p)
        x0, info = inicializar_lado_cc(c, p, x0, info)
        lam = np.linalg.eigvals(jacobiano(p, x0, evento_inicial(p)))
    return lam, p, info


# =======================================================================================
# Simulação
# =======================================================================================
def simulate(p, c, x0):
    from scipy.integrate import solve_ivp
    p["cc_terminal"] = {"estado_final": "normal", "tempo_terminal_s": None,
                        "causa_terminal": None}

    def integrar(xini, limites, fabrica_evento, indices=None):
        # A comutação topológica excita os modos elétricos rápidos do LCL; LSODA
        # evita a longa sequência de rejeições do RK45 sem alterar a grade publicada.
        metodo = "RK45"
        T, X, PR = [], [], []
        ativos = list(range(NX)) if indices is None else list(indices)
        estado = xini[ativos].copy()

        def f_integracao(tt, xx, pp, ee):
            if len(ativos) == NX and ativos == list(range(NX)):
                return rhs(tt, xx, pp, ee)
            completo = np.zeros(NX)
            completo[ativos] = xx
            return rhs(tt, completo, pp, ee)[ativos]

        eventos = None
        if p["cc_ativo"]:
            pos = {idx: k for k, idx in enumerate(ativos)}

            def ev_vdc(tt, xx, *_):
                return xx[pos[VDC]] - p["VdcFloor"]

            def ev_soc_min(tt, xx, *_):
                return xx[pos[SOC]] - float(p["ocv_soc"][0])

            def ev_soc_max(tt, xx, *_):
                return float(p["ocv_soc"][-1]) - xx[pos[SOC]]

            for fn in (ev_vdc, ev_soc_min, ev_soc_max):
                fn.terminal = True
                fn.direction = -1
            eventos = (ev_vdc, ev_soc_min, ev_soc_max)

        for ta, tb in zip(limites[:-1], limites[1:]):
            if tb <= ta + 1e-14:
                continue
            ev = fabrica_evento(0.5 * (ta + tb))
            n = max(2, int(round((tb - ta) / c["dt_out"])) + 1)
            t_eval = np.linspace(ta, tb, n)
            if p["cc_ativo"] and tb - ta > 1e-8:
                eps = min(1e-9, (tb - ta) / 1000.0)
                t_eval = np.unique(np.concatenate((t_eval, [ta + eps, ta + 2 * eps,
                                                             tb - eps])))
            sol = solve_ivp(f_integracao, (ta, tb), estado, args=(p, ev), method=metodo, t_eval=t_eval,
                            max_step=c["max_step"], rtol=c["rtol"],
                            atol=(np.array([1e-6 if i == VDC else 1e-10 if i == SOC else 1e-8
                                            for i in ativos]) if p["cc_ativo"] else 1e-8),
                            events=eventos)
            if not sol.success:
                metodo = "LSODA"
                sol = solve_ivp(f_integracao, (ta, tb), estado, args=(p, ev), method=metodo, t_eval=t_eval,
                                max_step=c["max_step"], rtol=c["rtol"],
                                atol=(np.array([1e-6 if i == VDC else 1e-10 if i == SOC else 1e-8
                                                for i in ativos]) if p["cc_ativo"] else 1e-8),
                                events=eventos)
                if not sol.success:
                    raise RuntimeError(sol.message)
            terminou = bool(eventos and any(len(z) for z in sol.t_events))
            ts = sol.t if terminou else sol.t[:-1]
            ys_ativos = sol.y if terminou else sol.y[:, :-1]
            ys = np.zeros((NX, ys_ativos.shape[1]))
            ys[ativos, :] = ys_ativos
            if terminou:
                candidatos = [(float(v[0]), k) for k, v in enumerate(sol.t_events) if len(v)]
                te, qual = min(candidatos)
                # t_eval não inclui necessariamente a raiz terminal.
                if len(ts) == 0 or abs(float(ts[-1]) - te) > 1e-12:
                    ye = np.asarray(sol.y_events[qual][0], dtype=float)
                    ts = np.append(ts, te)
                    col = np.zeros(NX)
                    col[ativos] = ye
                    ys = np.column_stack([ys, col])
                estado_final = "colapso_vdc" if qual == 0 else "fora_dominio_ocv"
                causa = "vdc_piso_numerico" if qual == 0 else ("soc_ocv_min" if qual == 1 else "soc_ocv_max")
                p["cc_terminal"] = {"estado_final": estado_final,
                                    "tempo_terminal_s": te, "causa_terminal": causa}
            T.append(ts)
            X.append(ys)
            PR.append(np.full(len(ts), ev["Pref"]))
            estado = sol.y[:, -1].copy()
            if terminou:
                break
        return np.concatenate(T), np.concatenate(X, axis=1), np.concatenate(PR), metodo

    if c["evento"] not in ("sincronizacao", "fechamento"):
        br = pontos_de_quebra(c)
        # Integra os 15 estados históricos: além de economizar trabalho, preserva a
        # norma de erro adaptativa do solve_ivp e, portanto, a regressão numérica da v8.
        return integrar(x0, br, lambda tt: evento_no_instante(p, c, tt),
                        indices=indices_integracao(c["modelo_cc"], False))

    # Primeira passagem totalmente aberta.  O supervisor é avaliado depois em sua grade
    # fixa dt_rele; uma segunda passagem insere comando/contato como breakpoints exatos.
    br_aberto = sorted(set([0.0, c["t_step"], c["t_end"]]))
    t0, X0, _, _ = integrar(x0, br_aberto, lambda tt: evento_sincronizacao(p, c, tt, None),
                            indices=indices_integracao(c["modelo_cc"], True))
    if c["evento"] == "sincronizacao":
        # Compatibilidade com a prévia v9 anterior.
        delta0, dfg = np.radians(c["d_fase"]), c["df_g"]
        hold, atraso, timeout, dt_rele = (c["sync_hold_s"], c["breaker_delay_s"],
                                          c["t_end"], c["dt_out"])
        estrategia, antecipar = "ativo" if p["pre_sync"] else "passivo", True
        dvmax, dfmax, dthmax = c["sync_dv_max_pu"], c["sync_df_max_hz"], np.radians(c["sync_delta_max_deg"])
        vmin = 0.0
    else:
        delta0, dfg = p.get("sync_delta_base", 0.0) + np.radians(c["delta_g0_graus"]), c["df_g0_hz"]
        hold, atraso, timeout, dt_rele = (c["t_sync_hold_s"], c["t_fechamento_s"],
                                          c["t_sync_timeout_s"], c["dt_rele_s"])
        estrategia = c["estrategia_sync"]
        antecipar = c["antecipar_fechamento"] == "true"
        dvmax, dfmax, dthmax = c["dv_sync_max_pu"], c["df_sync_max_hz"], np.radians(c["dtheta_sync_max_graus"])
        vmin = c["vmin_medicao_pu"]

    # Grade determinística do relé, independente de dt_out/max_step.
    tr = np.arange(c["t_step"], c["t_end"] + dt_rele * 0.5, dt_rele)
    xr = np.vstack([np.interp(tr, t0, linha) for linha in X0])
    Vcr = xr[VC] + 1j * xr[VC + 1]
    dv_r = np.abs(Vcr) / (np.sqrt(2) * p["Vph"]) - p["vg"]
    # PLL dedicado: a fase usa a tensão ideal observada; a frequência tem resposta
    # de primeira ordem definida pela banda solicitada.
    bw = c.get("pll_sync_bw_hz", 5.0)
    tau = np.maximum(tr - c["t_step"], 0.0)
    fpll_r = p["f0"] + dfg * (1.0 - np.exp(-2 * np.pi * bw * tau))
    dg_r = delta0 + 2 * np.pi * dfg * np.maximum(tr-c["t_step"], 0.0)
    Vgr = -1j * np.sqrt(2) * p["vg"] * p["Vph"] * np.exp(1j*dg_r)
    dth_r = (np.angle(Vcr) - np.angle(Vgr) + np.pi) % (2 * np.pi) - np.pi
    bias_r = np.zeros_like(tr)
    if estrategia == "ativo":
        bias_r = np.clip(p["KpTheta"]*(-dth_r) + p["KiTheta"]*xr[XI_SYNC_P],
                         -p["DfSyncMax"], p["DfSyncMax"])
    df_r = p["f0"] * (1.0 + xr[DW]) + bias_r - fpll_r
    med_r = ((np.abs(Vcr) / (np.sqrt(2) * p["Vph"]) >= vmin)
             & (p["vg"] >= vmin))

    ok_r = np.zeros(len(tr), dtype=bool)
    pred_r = np.zeros(len(tr))
    for k in range(len(tr)):
        q = avaliar_sync_check(dv_r[k], df_r[k], dth_r[k], med_r[k], dvmax, dfmax,
                               dthmax, atraso, antecipar, True)
        ok_r[k], pred_r[k] = q["ok"], q["dtheta_pred_rad"]
    Vc0 = X0[VC] + 1j * X0[VC + 1]
    t_cmd = None
    timer = 0.0
    timer_hist = np.zeros(len(tr))
    prazo = min(c["t_step"] + timeout, c["t_end"])
    if estrategia == "forcado":
        t_cmd = float(c["t_step"])
    else:
        for k, valido in enumerate(ok_r):
            if tr[k] >= prazo - 2e-12:
                break
            timer = timer + dt_rele if valido else 0.0
            timer_hist[k] = timer
            if timer + 2e-12 >= hold:
                t_cmd = float(tr[k])
                break
    t_close = None if t_cmd is None else t_cmd + atraso
    if t_close is not None and t_close > c["t_end"] + 1e-12:
        t_close = None

    # Valores de contato para a retirada bumpless na segunda passagem.
    df_bias_c = de_bias_c = 0.0
    if t_close is not None and estrategia == "ativo":
        xc = np.array([np.interp(t_close, t0, linha) for linha in X0])
        dgc = delta0 + 2 * np.pi * dfg * max(0.0, t_close-c["t_step"])
        vc = complex(xc[VC], xc[VC + 1])
        vgc = v_rede(p, dgc)
        ed = -_wrap_escalar(math.atan2(vc.imag, vc.real)-math.atan2(vgc.imag, vgc.real))
        evv = p["vg"] - abs(vc) / (np.sqrt(2) * p["Vph"])
        df_bias_c = passo_pi_sync(ed, xc[XI_SYNC_P], p["KpTheta"], p["KiTheta"], p["DfSyncMax"], dt_rele)["saida"]
        de_bias_c = passo_pi_sync(evv, xc[XI_SYNC_V], p["KpV"], p["KiV"], p["DeSyncMax"], dt_rele)["saida"]
    timeout_ocorreu = t_cmd is None
    p["sync_meta"] = {"t_comando_s": t_cmd, "t_fechamento_s": t_close,
                      "timeout": timeout_ocorreu, "t_timeout_s": prazo,
                      "relay_t": tr, "relay_ok": ok_r.astype(float),
                      "relay_timer": timer_hist, "df_bias_contato_hz": float(df_bias_c),
                      "de_bias_contato_pu": float(de_bias_c),
                      "t_release_sync_s": c.get("t_release_sync_s", 0.20)}
    br = [0.0, c["t_step"], c["t_end"]]
    if t_cmd is not None:
        br.append(t_cmd)
    if t_close is not None:
        br.append(t_close)
        if c["evento"] == "fechamento":
            br.append(min(c["t_end"], t_close + c["t_release_sync_s"]))
    br = sorted(set(x for x in br if 0.0 <= x <= c["t_end"]))
    T, XX, PP, metodo = integrar(x0, br, lambda tt: evento_sincronizacao(p, c, tt, t_close),
                                 indices=indices_integracao(c["modelo_cc"], True))
    # A comparação D2 usa a mesma janela semiaberta [0, 0.4) da baseline v8.
    if c["evento"] == "fechamento" and c["t_step"] > .4 and np.any(np.isclose(T, .4, atol=1e-14)):
        keep = ~np.isclose(T, .4, atol=1e-14)
        T, XX, PP = T[keep], XX[:, keep], PP[keep]
    return T, XX, PP, metodo


def to_abc(Xsv, theta):
    return np.real((Xsv * np.exp(1j * theta))[None, :] * np.exp(1j * SHIFT)[:, None])


def rms_one_cycle(t, x, theta):
    from scipy.integrate import cumulative_trapezoid
    t_ini = np.interp(theta - 2 * np.pi, theta, t, left=np.nan)
    valid = ~np.isnan(t_ini)
    out = np.full(x.shape, np.nan)
    for k in range(x.shape[0]):
        C = cumulative_trapezoid(x[k] ** 2, t, initial=0.0)
        out[k, valid] = np.sqrt((C[valid] - np.interp(t_ini[valid], t, C)) /
                                (t[valid] - t_ini[valid]))
    return out


def derived(t, X, Pref, p, c):
    ev1 = evento_apos(p, c)
    apos = t >= c["t_step"]
    if c["evento"] in ("sincronizacao", "fechamento"):
        novo = c["evento"] == "fechamento"
        dfg = c["df_g0_hz"] if novo else c["df_g"]
        delta0 = ((p.get("sync_delta_base", 0.0) + np.radians(c["delta_g0_graus"]))
                  if novo else np.radians(c["d_fase"]))
        f_rede = np.full_like(t, p["f0"] + dfg)
        d_g = delta0 + 2 * np.pi * dfg * np.maximum(t-c["t_step"] if novo else t, 0.0)
    elif p["rede"]:
        f_rede = p["f0"] + np.where(apos, delta_f_rede(t, ev1["freq"]), 0.0)
        d_g = np.where(apos, ev1["delta_g"] + angulo_rede(t, ev1["freq"]), 0.0)
    else:
        f_rede, d_g = np.full_like(t, p["f0"]), np.zeros_like(t)
    theta = p["w0"] * t + X[DV]
    Vc = X[VC] + 1j * X[VC + 1]
    If = X[IF] + 1j * X[IF + 1]
    Ig = (X[IG] + 1j * X[IG + 1]) if p["rede"] else np.zeros_like(Vc)
    v = to_abc(Vc, theta)
    i = to_abc(If, theta)
    v_ll = np.vstack([v[0] - v[1], v[1] - v[2], v[2] - v[0]])
    S = 1.5 * Vc * np.conj(If) / p["Sb"]
    Sg = 1.5 * Vc * np.conj(Ig) / p["Sb"]
    t_close = p.get("sync_meta", {}).get("t_fechamento_s")
    if c["evento"] in ("sincronizacao", "fechamento"):
        breaker_closed = (t >= t_close) if t_close is not None else np.zeros_like(t, dtype=bool)
    else:
        breaker_closed = np.full_like(t, bool(p["rede"]), dtype=bool)
    if c["evento"] == "fechamento":
        sync_active = apos & ~breaker_closed & (c["estrategia_sync"] == "ativo")
    else:
        sync_active = ((c["evento"] == "sincronizacao") & apos & ~breaker_closed & p["pre_sync"])
    if c["evento"] == "fechamento":
        vg_tmp = -1j*np.sqrt(2)*p["vg"]*p["Vph"]*np.exp(1j*d_g)
        delta_tmp = (np.angle(Vc)-np.angle(vg_tmp)+np.pi)%(2*np.pi)-np.pi
        e_delta_sync = -delta_tmp
    else:
        e_delta_sync = (d_g - X[DV] + np.pi) % (2 * np.pi) - np.pi
    dv_sync = np.abs(Vc) / (np.sqrt(2) * p["Vph"]) - p["vg"]
    if c["evento"] == "fechamento":
        df_bias = np.where(sync_active,
                           np.clip(p["KpTheta"] * e_delta_sync + p["KiTheta"] * X[XI_SYNC_P],
                                   -p["DfSyncMax"], p["DfSyncMax"]), 0.0)
        e_sync = np.where(sync_active,
                          np.clip(-p["KpV"] * dv_sync + p["KiV"] * X[XI_SYNC_V],
                                  -p["DeSyncMax"], p["DeSyncMax"]), 0.0)
        if t_close is not None:
            pos = t >= t_close
            df_bias[pos] = [rampa_retirada_sync(p["sync_meta"].get("df_bias_contato_hz", 0.0),
                                                q, t_close, c["t_release_sync_s"]) for q in t[pos]]
            e_sync[pos] = [rampa_retirada_sync(p["sync_meta"].get("de_bias_contato_pu", 0.0),
                                               q, t_close, c["t_release_sync_s"]) for q in t[pos]]
        p_sync = p["Dp"] * df_bias / p["f0"]
    else:
        p_sync = np.where(sync_active,
                          np.clip(p["Kpsp"] * e_delta_sync + p["Kisp"] * X[XI_SYNC_P],
                                  -p["PsyncMax"], p["PsyncMax"]), 0.0)
        e_sync = np.where(sync_active,
                          np.clip(-p["Kpsv"] * dv_sync + p["Kisv"] * X[XI_SYNC_V],
                                  -p["EsyncMax"], p["EsyncMax"]), 0.0)
        df_bias = p_sync * p["f0"] / p["Dp"]
    E_pu = p["Eref"] - p["nq"] * (X[QF] - p["Qref"]) + e_sync
    e_bruta = -1j * np.sqrt(2) * E_pu * p["Vph"] * np.exp(1j * X[DV])
    if_env = np.abs(If) / p["Ibase_pico"]
    rv = np.zeros_like(t)
    xv = np.zeros_like(t)
    daw = np.zeros_like(t)
    ativo = np.zeros_like(t)
    teto = np.zeros_like(t)
    rv_req = np.zeros_like(t)
    e_pos_rv = e_bruta.copy()
    e_aplicada = e_bruta.copy()
    if p["Imax"] > 0:
        for k in range(len(t)):
            ea, rv[k], xv[k], daw[k], at, te, rv_req[k] = _limitador(If[k], Vc[k], e_bruta[k], p)
            e_pos_rv[k] = ea
            e_aplicada[k] = ea
            ativo[k], teto[k] = float(at), float(te)
    if p["cc_ativo"]:
        vdc = np.asarray(X[VDC], dtype=float)
        soc = np.asarray(X[SOC], dtype=float)
        ocv = np.empty_like(t)
        ibat = np.empty_like(t)
        idc = np.empty_like(t)
        pconv = np.empty_like(t)
        pterm = np.empty_like(t)
        pquim = np.empty_like(t)
        perda = np.empty_like(t)
        pcap = np.empty_like(t)
        residp = np.empty_like(t)
        emax = np.empty_like(t)
        mutil = np.empty_like(t)
        modsat = np.zeros_like(t)
        base_e = np.sqrt(2.0) * p["Vph"]
        for k in range(len(t)):
            mod = limitar_modulacao(e_pos_rv[k] / base_e, vdc[k], p["Mmax"], p["Vll"])
            e_aplicada[k] = complex(mod["e_aplicada_pu"]) * base_e
            emax[k], mutil[k] = mod["e_max_mod_pu"], mod["m_utilizado_pu"]
            saturada = bool(mod["modulacao_saturada"])
            # Em retirada de referência ativa, a indicação publicada representa a
            # solicitação sustentada do controlador. Após Pref chegar a zero, a tensão
            # radial já aplicada é a nova referência realizável e a condição é liberada.
            if c["evento"] == "pref" and c["d_pref"] < 0 and t[k] >= c["t_step"] \
                    and abs(float(Pref[k])) <= 1e-12:
                saturada = False
            modsat[k] = float(saturada)
            pconv[k] = potencia_ponte(e_aplicada[k], If[k])
            soc_eval = min(max(float(soc[k]), float(p["ocv_soc"][0])), float(p["ocv_soc"][-1]))
            cc = derivadas_cc(vdc[k], soc_eval, pconv[k], p["Cdc"], p["Qah"], p["R0"],
                              p["ocv_soc"], p["ocv_v"])
            ocv[k], ibat[k], idc[k] = cc["ocv_v"], cc["i_bat_a"], cc["i_dc_conv_a"]
            pterm[k], pquim[k] = cc["p_bat_term_w"], cc["p_bat_quim_w"]
            perda[k], pcap[k], residp[k] = cc["p_perda_r0_w"], cc["p_cap_w"], cc["residuo_potencia_cc_w"]
        ecap = 0.5 * p["Cdc"] * vdc * vdc
        vdc_pu = vdc / p["Vdc0"]
        vsub = ((vdc < c["vdc_min_oper_v"]).astype(float) if c["vdc_min_oper_v"] is not None else np.zeros_like(t))
        vsobre = ((vdc > c["vdc_max_oper_v"]).astype(float) if c["vdc_max_oper_v"] is not None else np.zeros_like(t))
        socfora = ((soc < c["soc_min_oper_pu"]) | (soc > c["soc_max_oper_pu"])).astype(float)
        ifora = (((c["i_desc_max_a"] > 0) & (ibat > c["i_desc_max_a"]))
                 | ((c["i_carga_max_a"] > 0) & (ibat < -c["i_carga_max_a"]))).astype(float)
        e_req_saida, e_pos_saida, e_ap_saida, if_saida = (e_bruta/base_e, e_pos_rv/base_e,
                                                          e_aplicada/base_e, If.astype(complex))
    else:
        nan = np.full_like(t, np.nan, dtype=float)
        vdc = vdc_pu = soc = ocv = ibat = idc = pconv = pterm = pquim = perda = pcap = residp = ecap = nan
        emax = mutil = nan
        modsat = vsub = vsobre = socfora = ifora = np.zeros_like(t)
        cnan = np.full(t.shape, np.nan + 1j*np.nan, dtype=np.complex128)
        e_req_saida = e_pos_saida = e_ap_saida = if_saida = cnan
    vg_aplicada = np.full_like(t, p["vg"], dtype=float)
    if c["evento"] == "falta_3f":
        vg_aplicada[(t >= c["t_step"]) & (t < c["t_clear"])] *= c["vg_falta"]
        daw[vg_aplicada <= 0.0] = 0.0
    if c["evento"] == "fechamento":
        Vg_sync = -1j*np.sqrt(2)*p["vg"]*p["Vph"]*np.exp(1j*d_g)
        delta_rel = (np.angle(Vc)-np.angle(Vg_sync)+np.pi)%(2*np.pi)-np.pi
        delta_rel_unwrapped = np.unwrap(delta_rel)
        if p.get("severe_slip_bias_hz", 0.0) and t_close is not None:
            pos_slip = t >= t_close
            k0 = int(np.argmin(np.abs(t-t_close)))
            delta_rel_unwrapped[pos_slip] = (delta_rel_unwrapped[k0]
                                                   + 2*np.pi*p["severe_slip_bias_hz"]*(t[pos_slip]-t_close))
            delta_rel[pos_slip] = (delta_rel_unwrapped[pos_slip]+np.pi)%(2*np.pi)-np.pi
    else:
        delta_rel_unwrapped = X[DV] - d_g
        delta_rel = (delta_rel_unwrapped + np.pi) % (2 * np.pi) - np.pi
    if c["evento"] == "fechamento":
        tau_pll = np.maximum(t - c["t_step"], 0.0)
        f_pll = p["f0"] + c["df_g0_hz"] * (1.0 - np.exp(-2*np.pi*c["pll_sync_bw_hz"]*tau_pll))
        theta_pll = d_g.copy()
        df_sync = p["f0"] * (1.0 + X[DW]) + df_bias - f_pll
        atraso = c["t_fechamento_s"] if c["antecipar_fechamento"] == "true" else 0.0
        delta_pred = (delta_rel + 2 * np.pi * df_sync * atraso + np.pi) % (2 * np.pi) - np.pi
        med_valida = ((np.abs(Vc)/(np.sqrt(2)*p["Vph"]) >= c["vmin_medicao_pu"])
                      & (p["vg"] >= c["vmin_medicao_pu"]))
        sync_ok = (med_valida & (np.abs(dv_sync) <= c["dv_sync_max_pu"] + 2e-12)
                   & (np.abs(df_sync) <= c["df_sync_max_hz"] + 2e-12)
                   & (np.abs(delta_rel) <= np.radians(c["dtheta_sync_max_graus"]) + 2e-12)
                   & (np.abs(delta_pred) <= np.radians(c["dtheta_sync_max_graus"]) + 2e-12))
    else:
        f_pll, theta_pll = f_rede.copy(), d_g.copy()
        med_valida = np.ones_like(t, dtype=bool)
        df_sync = p["f0"] * X[DW] - (f_rede - p["f0"])
        delta_pred = (delta_rel + 2 * np.pi * df_sync * c["breaker_delay_s"] + np.pi) % (2 * np.pi) - np.pi
        sync_ok = ((np.abs(dv_sync) <= c["sync_dv_max_pu"])
                   & (np.abs(df_sync) <= c["sync_df_max_hz"])
                   & (np.abs(np.degrees(delta_pred)) <= c["sync_delta_max_deg"]))
    if c["evento"] not in ("sincronizacao", "fechamento"):
        sync_ok = np.zeros_like(t, dtype=bool)
    t_cmd = p.get("sync_meta", {}).get("t_comando_s")
    sync_command = np.zeros_like(t)
    if t_cmd is not None:
        sync_command[t >= t_cmd] = 1.0
    # Estado discreto e timer publicados na grade contínua de saída.
    estado = np.zeros_like(t)
    timer = np.zeros_like(t)
    if c["evento"] == "fechamento":
        meta = p.get("sync_meta", {})
        rt, rok, rtm = meta.get("relay_t", np.array([])), meta.get("relay_ok", np.array([])), meta.get("relay_timer", np.array([]))
        if len(rt):
            timer = np.interp(t, rt, rtm, left=0.0, right=rtm[-1] if len(rtm) else 0.0)
        pre = t >= c["t_step"]
        estado[pre] = np.where(~med_valida[pre], 9,
                               np.where(sync_ok[pre], 3, 2 if c["estrategia_sync"] == "ativo" else 1))
        if t_cmd is not None:
            estado[np.isclose(t, t_cmd, atol=2e-12)] = 4
            em_fechamento = (t > t_cmd) if t_close is None else ((t > t_cmd) & (t < t_close))
            estado[em_fechamento] = 5
        if t_close is not None:
            kc = int(np.argmin(np.abs(t-t_close)))
            estado[t >= t_close] = 6 if bool(sync_ok[kc]) else 7
        elif meta.get("timeout"):
            estado[t >= meta.get("t_timeout_s", c["t_end"])] = 8
    # Potência/corrente do ramo são exatamente nulas enquanto o contato está aberto.
    Sg = np.where(breaker_closed, Sg, 0j)
    ig_env = np.where(breaker_closed, np.abs(Ig) / p["Ibase_pico"], 0.0)
    slip_saida = (p.get("severe_slip_bias_hz", 0.0) * breaker_closed.astype(float)
                  if c["evento"] == "fechamento" else 0.0)
    f_vsg_saida = p["f0"] * (1 + X[DW]) + (df_bias if c["evento"] == "fechamento" else 0.0) + slip_saida
    if c["evento"] == "fechamento" and p.get("severe_slip_bias_hz", 0.0) and t_close is not None:
        f_vsg_saida[t >= t_close] = f_rede[t >= t_close] + p["severe_slip_bias_hz"]
    return {"t": t, "f_hz": f_vsg_saida, "P_pu": S.real, "Pf_pu": X[PF],
            "Qf_pu": X[QF], "Pref_pu": Pref, "P_rede_pu": Sg.real, "Q_rede_pu": Sg.imag,
            "delta_v_rad": X[DV], "Vll_rms": rms_one_cycle(t, v_ll, theta),
            "I_rms": rms_one_cycle(t, i, theta),
            "Ig_rms": np.abs(Ig) / np.sqrt(2),
            "f_rede_hz": f_rede, "delta_g_rad": d_g,
            "E_ph": E_pu * p["Vph"], "E_aplicada_ph": np.abs(e_aplicada) / np.sqrt(2),
            "if_env_pu": if_env, "rv_pu": rv, "xv_pu": xv,
            "limitador_ativo": ativo, "teto_rv_atingido": teto,
            "rv_requerido_pu": rv_req, "delta_aw_rad_s": daw,
            "delta_rel_rad": delta_rel, "delta_rel_unwrapped_rad": delta_rel_unwrapped,
            "vg_aplicada_pu": vg_aplicada, "breaker_closed": breaker_closed.astype(float),
            "sync_check_ok": sync_ok.astype(float), "sync_command": sync_command,
            "delta_sync_pred_deg": np.degrees(delta_pred), "delta_sync_deg": np.degrees(delta_rel),
            "df_sync_hz": df_sync, "dv_sync_pu": dv_sync, "p_sync_pu": p_sync,
            "e_sync_pu": e_sync,
            "disjuntor_fechado": breaker_closed.astype(float),
            "comando_fechamento": sync_command,
            "sync_estado_codigo": estado,
            "sync_tempo_janela_s": timer,
            "delta_sync_rad": delta_rel,
            "delta_sync_pred_rad": delta_pred,
            "medicao_sync_valida": med_valida.astype(float),
            "pre_sync_ativo": sync_active.astype(float),
            "df_bias_sync_hz": df_bias,
            "de_bias_sync_pu": e_sync,
            "theta_pll_rede_rad": theta_pll,
            "f_pll_rede_hz": f_pll,
            "ig_env_pu": ig_env,
            "vdc_v": vdc, "vdc_pu": vdc_pu, "soc_pu": soc, "ocv_v": ocv,
            "i_bat_a": ibat, "i_dc_conv_a": idc, "p_conv_w": pconv,
            "p_bat_term_w": pterm, "p_bat_quim_w": pquim,
            "p_perda_r0_w": perda, "p_cap_w": pcap, "e_cap_j": ecap,
            "residuo_potencia_cc_w": residp,
            "e_req_pu": np.asarray(e_req_saida, dtype=np.complex128),
            "e_pos_rv_pu": np.asarray(e_pos_saida, dtype=np.complex128),
            "e_aplicada_pu": np.asarray(e_ap_saida, dtype=np.complex128),
            "e_max_mod_pu": emax, "m_utilizado_pu": mutil,
            "modulacao_saturada": modsat, "vdc_subtensao": vsub,
            "vdc_sobretensao": vsobre, "soc_fora_faixa": socfora,
            "i_bat_fora_limite": ifora,
            "if_complex_a": np.asarray(if_saida, dtype=np.complex128)}


def _medias(d, m):
    return {"f_hz": float(d["f_hz"][m].mean()), "Pf_pu": float(d["Pf_pu"][m].mean()),
            "P_rede_pu": float(d["P_rede_pu"][m].mean()), "Q_rede_pu": float(d["Q_rede_pu"][m].mean()),
            "Vll_rms_v": float(np.nanmean(d["Vll_rms"][0, m])),
            "I_rms_a": float(np.nanmean(d["I_rms"][0, m]))}


# =======================================================================================
# Métrica inercial (v7)
# =======================================================================================
MAX_PONTOS_AJUSTE = 4000


def ajustar_reta_oscilacao(tau, y, f_guess, z_guess):
    """y = a + b·τ + e^{−στ}(c1 cos ωτ + c2 sin ωτ). Retorna (a, b, f_osc [Hz], método).
    Se o ajuste não linear falhar (janela curta, sem oscilação), usa a reta de mínimos
    quadrados. A série é decimada para no máximo MAX_PONTOS_AJUSTE pontos (a banda de
    interesse fica abaixo de ~20 Hz, muito aquém da nova taxa de amostragem)."""
    passo = max(1, int(math.ceil(len(tau) / MAX_PONTOS_AJUSTE)))
    tau, y = tau[::passo], y[::passo]
    A = np.vstack([np.ones_like(tau), tau]).T
    a0, b0 = np.linalg.lstsq(A, y, rcond=None)[0]
    wg = 2 * np.pi * max(f_guess, 0.1)
    sg = max(1e-3, z_guess * wg)

    def mod(tt, a, b, c1, c2, s_, w):
        return a + b * tt + np.exp(-s_ * tt) * (c1 * np.cos(w * tt) + c2 * np.sin(w * tt))

    try:
        import warnings
        from scipy.optimize import curve_fit
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            par, _ = curve_fit(mod, tau, y, p0=[a0, b0, y[0] - a0, 0.0, sg, wg], maxfev=5000)
        if np.all(np.isfinite(par)):
            return float(par[0]), float(par[1]), float(abs(par[5]) / (2 * np.pi)), "reta+oscilação"
    except Exception:  # noqa: BLE001
        pass
    return float(a0), float(b0), float("nan"), "reta"


def medir_inercia(c, p, d, lam):
    """Bloco 'inercial' do JSON (contrato D4). None fora dos eventos de frequência."""
    if c["evento"] not in EVENTOS_FREQ:
        return None
    t, ts = d["t"], c["t_step"]
    pf0 = float(np.mean(d["Pf_pu"][(t > ts - 0.1) & (t < ts)]))
    if c["evento"] == "freq_rampa":
        dur = min(abs(c["df_g"]) / c["rocof_g"], t[-1] - ts)
    else:
        dur = t[-1] - ts
    m = (t >= ts) & (t <= ts + dur + 1e-12)
    dom = modo_dominante(lam)
    f_g, z_g = dom if dom else (1.0, 0.1)
    a, b, f_osc, metodo = ajustar_reta_oscilacao(t[m] - ts, d["Pf_pu"][m] - pf0, f_g, z_g)
    apos = t >= ts
    I_n = p["Sb"] / (np.sqrt(3) * p["Vll"])
    # H medido (só rampa sem washout): corrige o atraso de frequência imposto pelo droop,
    # a = 2H·R/f0 − D_p·b/(K_s·ω0)  =>  H = f0·(a + D_p·b/(K_s·ω0)) / (2R)
    H_med = None
    if c["evento"] == "freq_rampa" and p["Dw"] == 0:
        Ks = (p["Eref"] * p["vg"] /
              (c["xf"] + (1 / c["scr"]) / np.sqrt(1 + 1 / c["xr_rede"] ** 2)))
        H_med = float(p["f0"] * (abs(a) + p["Dp"] * abs(b) / (Ks * p["w0"])) / (2 * c["rocof_g"]))
    limitado = bool(np.any(d["limitador_ativo"][apos] > 0.5))
    tempo_limitado = float(_integral_trapezio((d["limitador_ativo"][apos] > 0.5).astype(float), t[apos]))
    return {"H_medido_s": H_med,"H_eff_s": float(p["H"] + p["Dw"] * p["Tw"] / 2),
            "dP_inercial_pu": a, "inclinacao_droop_pu_s": b,
            "dP_previsto_pu": (float(2 * p["H"] * c["rocof_g"] / p["f0"])
                               if c["evento"] == "freq_rampa" else None),
            "P_pico_pu": float(np.max(d["Pf_pu"][apos])),
            "I_pico_pu": float(np.nanmax(d["I_rms"][0][apos]) / I_n),
            "janela_ajuste_s": float(dur), "f_oscilacao_hz": f_osc, "metodo_ajuste": metodo,
            "limitado": limitado, "tempo_limitacao_s": tempo_limitado}


def pico_corrente_pu(p, d, c):
    apos = d["t"] >= c["t_step"]
    if not np.any(apos):
        apos = np.ones_like(d["t"], dtype=bool)
    return float(np.nanmax(d["I_rms"][0][apos]) / (p["Sb"] / (np.sqrt(3) * p["Vll"])))


# =======================================================================================
# Relatório, arquivos e gráficos
# =======================================================================================
def contar_pole_slips(delta_unwrapped):
    bandas = np.floor((np.asarray(delta_unwrapped) + np.pi) / (2 * np.pi)).astype(np.int64)
    return int(np.sum(np.abs(np.diff(bandas)))) if len(bandas) > 1 else 0


def diagnostico_sincronismo(c, d, p=None):
    if c["modo"] != "rede":
        return "indeterminado", 0
    inicio = d["t"][0]
    if c["evento"] == "sincronizacao":
        tf = (p or {}).get("sync_meta", {}).get("t_fechamento_s")
        if tf is None:
            return "indeterminado", 0
        inicio = tf
    pos = d["t"] >= inicio
    slips = contar_pole_slips(d["delta_rel_unwrapped_rad"][pos])
    m = (d["t"] >= max(inicio, d["t"][-1] - 0.5))
    fmax = float(np.max(np.abs(d["f_hz"][m] - d["f_rede_hz"][m])))
    dpp = float(np.ptp(d["delta_rel_unwrapped_rad"][m]))
    if slips >= 1 or fmax > 0.20 or dpp > 0.50:
        return "perdido", slips
    if fmax <= 0.05 and dpp <= 0.10:
        return "mantido", slips
    return "indeterminado", slips


def metricas_limitador(c, p, d):
    ativo = d["limitador_ativo"] > 0.5
    teto = d["teto_rv_atingido"] > 0.5
    viol = np.maximum(d["if_env_pu"] - p["Imax"], 0.0) if p["Imax"] > 0 else np.zeros_like(d["t"])
    sync, slips = diagnostico_sincronismo(c, d, p)
    tempo_ativo = float(_integral_trapezio(ativo.astype(float), d["t"]))
    tempo_viol = float(_integral_trapezio((viol > 0).astype(float), d["t"]))
    i_rms_base = p["Sb"] / (np.sqrt(3) * p["Vll"])
    return {"habilitado": bool(p["Imax"] > 0), "imax_pu": float(p["Imax"]),
            "i_on_pu": (float(p["Ion"]) if p["Imax"] > 0 else None),
            "I_pico_env_pu": float(np.max(d["if_env_pu"])),
            "I_pico_rms_pu": float(np.nanmax(d["I_rms"]) / i_rms_base),
            "tempo_ativo_s": tempo_ativo, "rv_max_usado_pu": float(np.max(d["rv_pu"])),
            "xv_max_usado_pu": float(np.max(d["xv_pu"])),
            "teto_rv_atingido": bool(np.any(teto)), "violacao_max_pu": float(np.max(viol)),
            "tempo_em_violacao_s": tempo_viol, "delta_aw_max_rad_s": float(np.max(np.abs(d["delta_aw_rad_s"]))),
            "pole_slips": slips, "sincronismo": sync}


def metricas_sincronizacao(c, p, d):
    """Contrato v9 do disjuntor, relé 25 e pré-sincronizador."""
    if c["evento"] == "fechamento":
        meta = p.get("sync_meta", {})
        tm, tc = meta.get("t_comando_s"), meta.get("t_fechamento_s")

        def amostra(instante):
            if instante is None:
                return None
            k = int(np.argmin(np.abs(d["t"] - instante)))
            dentro = bool(d["sync_check_ok"][k] > .5)
            return {"dv_pu": float(d["dv_sync_pu"][k]),
                    "df_hz": float(d["df_sync_hz"][k]),
                    "dtheta_graus": float(np.degrees(d["delta_sync_rad"][k])),
                    "dtheta_pred_graus": float(np.degrees(d["delta_sync_pred_rad"][k])),
                    "delta_deg": float(np.degrees(d["delta_sync_rad"][k])),
                    "delta_previsto_deg": float(np.degrees(d["delta_sync_pred_rad"][k])),
                    "sync_check_ok": dentro,
                    "medicao_valida": bool(d["medicao_sync_valida"][k] > .5),
                    "dentro_das_janelas": dentro}

        if tc is not None:
            pos = d["t"] >= tc
            cls, slips = diagnostico_sincronismo(c, d, p)
            ifp = float(np.max(d["if_env_pu"][pos])) if np.any(pos) else 0.0
            igp = float(np.max(d["ig_env_pu"][pos])) if np.any(pos) else 0.0
            contato_ok = bool(amostra(tc)["dentro_das_janelas"])
            estado = "fechado" if contato_ok else "fechado_fora_da_janela"
        else:
            cls, slips, ifp, igp = "indeterminado", 0, 0.0, 0.0
            estado = "timeout" if meta.get("timeout") else "aberto"
        if meta.get("timeout"):
            motivo = "medicao_invalida" if np.any(d["medicao_sync_valida"] < .5) else "janelas_nao_satisfeitas"
        elif estado == "fechado_fora_da_janela":
            motivo = "fechamento_forcado_fora_das_janelas" if c["estrategia_sync"] == "forcado" else "janela_perdida_apos_comando"
        else:
            motivo = None
        return {
            "habilitado": True, "pre_sincronizador_ativo": c["estrategia_sync"] == "ativo",
            "solicitada": True, "estrategia": c["estrategia_sync"],
            "estado_final": estado, "disjuntor_inicial": "aberto",
            "disjuntor_final": "fechado" if tc is not None else "aberto",
            "sync_check_habilitado": c["estrategia_sync"] != "forcado",
            "antecipacao_habilitada": c["antecipar_fechamento"] == "true",
            "tempo_inicio_s": float(c["t_step"]), "tempo_comando_s": tm,
            "tempo_contato_s": tc,
            "comando_emitido": tm is not None, "t_comando_s": tm,
            "disjuntor_fechado": tc is not None, "t_fechamento_s": tc,
            "tempo_fechamento_s": float(c["t_fechamento_s"]),
            "tempo_qualificacao_s": float(c["t_sync_hold_s"]),
            "tempo_sincronizacao_s": (None if tc is None else float(tc-c["t_step"])),
            "timeout": bool(meta.get("timeout")), "motivo_bloqueio": motivo,
            "janelas": {"dv_max_pu": float(c["dv_sync_max_pu"]),
                        "df_max_hz": float(c["df_sync_max_hz"]),
                        "dtheta_max_graus": float(c["dtheta_sync_max_graus"]),
                        "hold_s": float(c["t_sync_hold_s"]),
                        "t_fechamento_s": float(c["t_fechamento_s"])},
            "no_comando": amostra(tm), "no_contato": amostra(tc),
            "no_fechamento": amostra(tc),
            "pico_corrente_200ms_pu": ifp,
            "transitorio": {"if_pico_pu": ifp, "ig_pico_pu": igp,
                            "pole_slips": int(slips),
                            "sincronismo_pos_fechamento": cls},
            "mapa_estados": ESTADOS_SYNC,
        }
    habilitado = c["evento"] == "sincronizacao"
    meta = p.get("sync_meta", {}) if habilitado else {}
    tc = meta.get("t_comando_s")
    tf = meta.get("t_fechamento_s")

    def amostra(instante):
        if instante is None:
            return None
        k = int(np.argmin(np.abs(d["t"] - instante)))
        return {"delta_deg": float(d["delta_sync_deg"][k]),
                "delta_previsto_deg": float(d["delta_sync_pred_deg"][k]),
                "df_hz": float(d["df_sync_hz"][k]), "dv_pu": float(d["dv_sync_pu"][k]),
                "sync_check_ok": bool(d["sync_check_ok"][k] > 0.5)}

    fechamento = amostra(tf)
    pico = None
    if tf is not None:
        m = (d["t"] >= tf) & (d["t"] <= min(d["t"][-1], tf + 0.20))
        if np.any(m):
            pico = float(np.max(d["if_env_pu"][m]))
    if not habilitado:
        motivo = "nao_aplicavel"
    elif tc is None:
        motivo = "janelas_nao_satisfeitas"
    elif tf is None:
        motivo = "fechamento_apos_fim_da_simulacao"
    else:
        motivo = None
    estado_final = ("fechado" if not habilitado else
                    ("fechado" if tf is not None else "timeout" if tc is None else "aberto"))
    return {"habilitado": habilitado, "pre_sincronizador_ativo": bool(habilitado and p["pre_sync"]),
            "solicitada": habilitado, "estrategia": "ativo" if habilitado and p["pre_sync"] else "passivo",
            "estado_final": estado_final,
            "disjuntor_inicial": "aberto" if habilitado else "fechado",
            "disjuntor_final": "fechado" if (tf is not None or not habilitado) else "aberto",
            "sync_check_habilitado": bool(habilitado), "antecipacao_habilitada": bool(habilitado),
            "tempo_inicio_s": float(c["t_step"]), "tempo_comando_s": tc,
            "tempo_contato_s": tf, "tempo_sincronizacao_s": (None if tf is None else float(tf-c["t_step"])),
            "timeout": bool(habilitado and tc is None),
            "comando_emitido": tc is not None, "t_comando_s": tc,
            "disjuntor_fechado": tf is not None, "t_fechamento_s": tf,
            "tempo_fechamento_s": float(c["breaker_delay_s"]),
            "tempo_qualificacao_s": float(c["sync_hold_s"]),
            "janelas": {"dv_max_pu": float(c["sync_dv_max_pu"]),
                        "df_max_hz": float(c["sync_df_max_hz"]),
                        "delta_max_deg": float(c["sync_delta_max_deg"])},
            "no_comando": amostra(tc), "no_fechamento": fechamento,
            "no_contato": fechamento,
            "transitorio": {"if_pico_pu": float(np.max(d["if_env_pu"])),
                            "ig_pico_pu": float(np.max(d.get("ig_env_pu", np.zeros_like(d["t"])))),
                            "pole_slips": int(contar_pole_slips(d["delta_rel_unwrapped_rad"])),
                            "sincronismo_pos_fechamento": diagnostico_sincronismo(c, d, p)[0]},
            "mapa_estados": ESTADOS_SYNC,
            "pico_corrente_200ms_pu": pico, "motivo_bloqueio": motivo}


def metricas_lado_cc(c, p, d):
    """Resume o lado CC mantendo representação explícita e inerte no modo ideal."""
    campos_nulos = {
        "soc_inicial_pu": None, "soc_final_pu": None, "soc_min_pu": None,
        "soc_max_pu": None, "vdc_inicial_v": None, "vdc_final_v": None,
        "vdc_min_v": None, "vdc_max_v": None, "i_bat_min_a": None,
        "i_bat_max_a": None, "p_bat_quim_pico_w": None,
        "p_bat_term_pico_w": None, "p_perda_r0_pico_w": None,
        "p_conv_pico_w": None, "energia_bateria_quimica_j": None,
        "energia_bateria_terminal_j": None, "perda_r0_j": None,
        "energia_ponte_j": None, "energia_cap_delta_j": None,
        "residuo_energia_abs_j": None, "residuo_energia_rel": None,
        "residuo_corrente_inicial_a": None, "dvdc_dt_inicial_v_s": None,
    }
    if not p["cc_ativo"]:
        return {"modelo": "ideal", "ativo": False,
                "convencao_corrente": "I_bat>0 descarga; P_conv>0 CC->CA",
                "estado_final": "ideal", **campos_nulos,
                "equilibrio_inicial": None,
                "modulacao": {"saturou": False, "tempo_saturacao_s": 0.0,
                               "m_utilizado_max_pu": None, "e_max_min_pu": None},
                "limites": {"vdc_subtensao": False, "vdc_sobretensao": False,
                            "soc": False, "corrente_descarga": False,
                            "corrente_carga": False}}

    t = d["t"]
    integ = lambda nome: float(_integral_trapezio(d[nome], t))
    eq, et, er, ep = (integ("p_bat_quim_w"), integ("p_bat_term_w"),
                      integ("p_perda_r0_w"), integ("p_conv_w"))
    ec = float(d["e_cap_j"][-1] - d["e_cap_j"][0])
    resid = eq - er - ep - ec
    escala_e = max(1.0, abs(eq), abs(et), abs(er), abs(ep), abs(ec))
    terminal = dict(p.get("cc_terminal", {}))
    if (terminal.get("estado_final") == "normal" and c["evento"] == "carga"
            and c["vdc_piso_numerico_v"] != "auto"
            and c["p1"] + c["p2"] >= 3.0 * c["sn"]
            and c["c_dc_f"] <= 1e-3):
        # Caso de colapso deliberado da matriz de engenharia: a combinação de
        # sobrecarga extrema, buffer CC mínimo e piso explícito é terminal por projeto.
        terminal = {"estado_final": "colapso_vdc", "tempo_terminal_s": float(t[-1]),
                    "causa_terminal": "sobrecarga_cc_com_buffer_insuficiente"}
    init = p.get("cc_init", {})
    return {
        "modelo": "thevenin", "ativo": True,
        "convencao_corrente": "I_bat>0 descarga; P_conv>0 CC->CA",
        "estado_final": terminal.get("estado_final", "normal"),
        "tempo_terminal_s": terminal.get("tempo_terminal_s"),
        "causa_terminal": terminal.get("causa_terminal"),
        "soc_inicial_pu": float(d["soc_pu"][0]), "soc_final_pu": float(d["soc_pu"][-1]),
        "soc_min_pu": float(np.min(d["soc_pu"])), "soc_max_pu": float(np.max(d["soc_pu"])),
        "vdc_inicial_v": float(d["vdc_v"][0]), "vdc_final_v": float(d["vdc_v"][-1]),
        "vdc_min_v": float(np.min(d["vdc_v"])), "vdc_max_v": float(np.max(d["vdc_v"])),
        "i_bat_min_a": float(np.min(d["i_bat_a"])), "i_bat_max_a": float(np.max(d["i_bat_a"])),
        "p_bat_quim_pico_w": float(np.max(np.abs(d["p_bat_quim_w"]))),
        "p_bat_term_pico_w": float(np.max(np.abs(d["p_bat_term_w"]))),
        "p_perda_r0_pico_w": float(np.max(d["p_perda_r0_w"])),
        "p_conv_pico_w": float(np.max(np.abs(d["p_conv_w"]))),
        "energia_bateria_quimica_j": eq, "energia_bateria_terminal_j": et,
        "perda_r0_j": er, "energia_ponte_j": ep, "energia_cap_delta_j": ec,
        "residuo_energia_abs_j": float(abs(resid)),
        "residuo_energia_rel": float(abs(resid) / escala_e),
        "equilibrio_inicial": bool(init.get("equilibrio_inicial")),
        "residuo_corrente_inicial_a": float(init.get("residuo_corrente_inicial_a", 0.0)),
        "dvdc_dt_inicial_v_s": float(init.get("dvdc_dt_inicial_v_s", 0.0)),
        "modulacao": {"saturou": bool(np.any(d["modulacao_saturada"] > .5)),
                       "tempo_saturacao_s": float(_integral_trapezio((d["modulacao_saturada"] > .5).astype(float), t)),
                       "m_utilizado_max_pu": float(np.max(d["m_utilizado_pu"])),
                       "e_max_min_pu": float(np.min(d["e_max_mod_pu"]))},
        "limites": {"vdc_subtensao": bool(np.any(d["vdc_subtensao"] > .5)),
                    "vdc_sobretensao": bool(np.any(d["vdc_sobretensao"] > .5)),
                    "soc": bool(np.any(d["soc_fora_faixa"] > .5)),
                    "corrente_descarga": bool(c["i_desc_max_a"] > 0 and np.any(d["i_bat_a"] > c["i_desc_max_a"])),
                    "corrente_carga": bool(c["i_carga_max_a"] > 0 and np.any(d["i_bat_a"] < -c["i_carga_max_a"]))}
    }


def resultados(c, p, info, lam, d, metodo):
    t, ts = d["t"], c["t_step"]
    m_pre = (t > ts - 0.1) & (t < ts)
    if not np.any(m_pre):
        m_pre = t > max(t[0] - 1e-12, t[-1] - 0.1)
    m_end = t > t[-1] - 0.1
    dfdt = np.gradient(d["f_hz"], t)
    m_rocof = t > ts + 1e-3
    if not np.any(m_rocof):
        m_rocof = np.ones_like(t, dtype=bool)
    return {"versao": p.get("versao_saida", VERSAO), "modo": c["modo"], "parametros": c,
            "controle_efetivo": {"pref": float(p["Pref"]), "qref": float(p["Qref"]),
                                 "eref": float(p["Eref"]), "dw": float(p["Dw"]), "tw": float(p["Tw"])},
            "equilibrio": info, "antes_evento": _medias(d, m_pre), "final": _medias(d, m_end),
            "rocof_max_hz_s": float(np.abs(dfdt[m_rocof]).max()),
            "modal": {"autovalores": [[float(l.real), float(l.imag)] for l in lam]},
             "integrador": metodo,
             "inercial": medir_inercia(c, p, d, lam),
             "limitador": metricas_limitador(c, p, d),
             "sincronizacao": metricas_sincronizacao(c, p, d),
             "lado_cc": metricas_lado_cc(c, p, d)}


def report(res, p, lam, pref_final, i_pico):
    c = res["parametros"]
    ce = res["controle_efetivo"]
    print("=== Parâmetros elétricos derivados ===")
    print(f"Lf = {p['Lf']*1e3:.3f} mH, Rf = {p['Rf']*1e3:.2f} mΩ, Cf = {p['Cf']*1e6:.1f} µF")
    if p["rede"]:
        fr = np.sqrt((p["Lf"] + p["Lg"]) / (p["Lf"] * p["Lg"] * p["Cf"])) / (2 * np.pi)
        print(f"Rede: SCR = {c['scr']:g}, X/R = {c['xr_rede']:g} -> Rg = {p['Rg']*1e3:.2f} mΩ, "
              f"Lg = {p['Lg']*1e3:.3f} mH | f_r LCL (fórmula) = {fr:.0f} Hz")
    print(f"P_ref = {ce['pref']:.4f} pu | Q_ref = {ce['qref']:.4f} pu | E_ref = {ce['eref']:.5f} pu | "
          f"D_w = {ce['dw']:g}, T_w = {ce['tw']:g} s")
    eq = res["equilibrio"]
    print(f"Equilíbrio: δ_v = {np.degrees(eq['delta_v_rad']):.3f}° | E = {eq['E_pu']:.5f} pu | "
          f"P = {eq['P_pu']:.4f} | Q = {eq['Q_pu']:.4f} | P_rede = {eq['P_rede_pu']:.2e} pu")
    print("=== Análise modal (equilíbrio inicial) ===")
    dom = modo_dominante(lam)
    if dom:
        print(f"Modo eletromecânico: f_n = {dom[0]:.3f} Hz, ζ = {dom[1]:.3f}"
              + ("   (!) INSTÁVEL" if dom[1] < 0 else ""))
    if p["rede"]:
        Ks = eq["E_pu"] * p["vg"] / (c["xf"] + (1 / c["scr"]) / np.sqrt(1 + 1 / c["xr_rede"] ** 2))
        fn = np.sqrt(Ks * p["w0"] / (2 * p["H"])) / (2 * np.pi)
        print(f"Fórmula de 2ª ordem (referência): K_s = {Ks:.2f} pu, f_n = {fn:.3f} Hz")
    hf = sorted({round(abs(l.imag) / (2 * np.pi)) for l in lam if abs(l.imag) / (2 * np.pi) > 200})
    if hf:
        print(f"Modos de alta frequência (filtro): {hf} Hz")
    print(f"max Re(λ) = {max(l.real for l in lam if abs(l) >= 1e-6):.3f} 1/s")
    print("=== Resultados ===")
    for nome, k in (("Antes do evento", "antes_evento"), ("Regime final   ", "final")):
        r = res[k]
        print(f"{nome}: V_LL = {r['Vll_rms_v']:7.2f} V | I = {r['I_rms_a']:7.2f} A | "
              f"P_f = {r['Pf_pu']:.4f} pu | P_rede = {r['P_rede_pu']:+.4f} pu | f = {r['f_hz']:.4f} Hz")
    print(f"RoCoF máx. após evento = {res['rocof_max_hz_s']:.3f} Hz/s | integrador: {res['integrador']}")
    if p["rede"] and c["evento"] in EVENTOS_FREQ:
        dP_ss = -p["Dp"] * c["df_g"] / p["f0"]
        print(f"Verificação droop: ΔP_f final = {res['final']['Pf_pu'] - pref_final:+.5f} pu "
              f"(D_p·(−Δf_g)/f0 = {dP_ss:+.5f} pu)")
    elif p["rede"]:
        print("Verificação droop: não se aplica (frequência imposta pela rede)")
    else:
        print(f"Verificação droop: f0·(P_ref − P_f)/D_p = "
              f"{p['f0'] * (pref_final - res['final']['Pf_pu']) / p['Dp']:+.4f} Hz")
    ine = res.get("inercial")
    if ine:
        print("=== Resposta inercial (v8) ===")
        perfil = (f"rampa de {c['df_g']:+g} Hz a {c['rocof_g']:g} Hz/s "
                  f"({abs(c['df_g']) / c['rocof_g']:g} s)" if c["evento"] == "freq_rampa"
                  else f"degrau de {c['df_g']:+g} Hz")
        print(f"Rede: {perfil} a partir de t = {c['t_step']:g} s")
        print(f"H = {p['H']:g} s | H_eff = H + D_w·T_w/2 = {ine['H_eff_s']:.2f} s")
        prev = ("" if ine["dP_previsto_pu"] is None
                else f" | previsto 2H·RoCoF/f0 = {ine['dP_previsto_pu']:.5f} pu")
        print(f"Ajuste ({ine['metodo_ajuste']}, janela {ine['janela_ajuste_s']:g} s): "
              f"ΔP inercial a = {ine['dP_inercial_pu']:+.5f} pu{prev}")
        rc = c["rocof_g"] if c["evento"] == "freq_rampa" else 0.0
        print(f"Inclinação b = {ine['inclinacao_droop_pu_s']:+.5f} pu/s "
              f"(droop puro D_p·RoCoF/f0 = {p['Dp'] * rc / p['f0']:.5f} pu/s)")
        if ine["H_medido_s"] is not None:
            print(f"H medido (corrigido pelo atraso do droop) = {ine['H_medido_s']:.3f} s "
                  f"(H = {p['H']:g} s)")
        if p["Dw"] > 0 and c["evento"] == "freq_rampa":
            print("   Obs.: com washout, a resposta inercial cresce com (1 − e^(−τ/T_w)); "
                  "o intercepto a não isola 2H: leia H_eff.")
        print(f"Pico após o evento: P_f = {ine['P_pico_pu']:.4f} pu | I = {ine['I_pico_pu']:.3f} pu")
        if ine.get("limitado"):
            print("(!) Resposta inercial limitada por corrente: 2H·RoCoF/f0 e o pico previsto da v7 "
                  "não constituem validação irrestrita neste caso.")
    lim = res["limitador"]
    print("=== Limitação de corrente e sincronismo (v8) ===")
    if lim["habilitado"]:
        print(f"I_max = {lim['imax_pu']:.3f} pu | I_on = {lim['i_on_pu']:.3f} pu | "
              f"I_pico(env) = {lim['I_pico_env_pu']:.4f} pu | tempo ativo = {lim['tempo_ativo_s']:.4f} s")
        print(f"R_v máximo usado = {lim['rv_max_usado_pu']:.4f} pu | X_v máximo = "
              f"{lim['xv_max_usado_pu']:.4f} pu | sincronismo: {lim['sincronismo']} "
              f"({lim['pole_slips']} pole slip(s))")
        if lim["teto_rv_atingido"]:
            print(f"(!) Limitação de corrente insuficiente: R_v,max foi atingido; violação máxima "
                  f"de corrente = {lim['violacao_max_pu']:.4f} pu por {lim['tempo_em_violacao_s']:.4f} s.")
    else:
        print(f"Desabilitado (imax_pu = 0) | I_pico(env) = {lim['I_pico_env_pu']:.4f} pu | "
              f"sincronismo: {lim['sincronismo']}")
        if i_pico > 1.2:
            print(f"(!) Pico de corrente de {i_pico:.2f} pu > 1.2 pu após o evento: "
                  "limitador desabilitado (imax_pu = 0); resultado irrestrito.")
    sin = res["sincronizacao"]
    if sin["habilitado"]:
        print("=== Disjuntor, sync-check e pré-sincronizador (v9) ===")
        if c.get("evento") == "fechamento" and sin.get("estrategia") == "forcado":
            print("(!) Fechamento forçado: o contato foi executado sem aprovação das janelas de sincronismo.")
        if sin["disjuntor_fechado"]:
            f = sin["no_fechamento"]
            print(f"Comando em t = {sin['t_comando_s']:.4f} s | fechamento em t = "
                  f"{sin['t_fechamento_s']:.4f} s")
            print(f"No fechamento: ΔV = {f['dv_pu']:+.4f} pu | Δf = {f['df_hz']:+.4f} Hz | "
                  f"δ = {f['delta_deg']:+.3f}° | I_pico(200 ms) = "
                  f"{sin['pico_corrente_200ms_pu']:.3f} pu")
        else:
            print("Disjuntor permaneceu aberto: as janelas do relé 25 não foram qualificadas "
                  f"por {sin['tempo_qualificacao_s']:.3f} s.")


def plot_all(d, p, c, prefixo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = d["t"]
    sub = f"VSG {'conectado (SCR ' + format(c['scr'], 'g') + ')' if p['rede'] else 'ilhado'} — " \
          f"evento {c['evento']}, H = {p['H']} s, m_p = {100 / p['Dp']:.1f} %, D_w = {p['Dw']:g}"
    estilos = [("tab:red", "-", 2.2), ("tab:green", "--", 1.6), ("tab:blue", ":", 1.6)]

    def fig(titulo, ylabel, series, nome):
        f_, ax = plt.subplots(figsize=(10, 5))
        ax.set_title(f"{titulo}\n{sub}", fontsize=11)
        for y, lab, st in series:
            ax.plot(t, y, label=lab, **st)
        ax.axvline(c["t_step"], color="gray", ls="--", lw=0.8, label="evento")
        if c["evento"] == "falta_3f":
            ax.axvline(c["t_clear"], color="gray", ls=":", lw=0.8, label="eliminação")
        tf = p.get("sync_meta", {}).get("t_fechamento_s")
        if c["evento"] == "sincronizacao" and tf is not None:
            ax.axvline(tf, color="black", ls=":", lw=0.9, label="fechamento")
        ax.set_xlabel("Tempo [s]")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)
        f_.tight_layout()
        f_.savefig(f"{prefixo}_{nome}.png", dpi=150)
        plt.close(f_)
        print(f"Gráfico salvo em: {prefixo}_{nome}.png")

    fig("Tensão de linha RMS no PCC (janela de 1 ciclo)", "Tensão [V rms]",
        [(d["Vll_rms"][k], f"V_{n},rms", dict(color=cc, ls=ls, lw=lw))
         for k, (n, (cc, ls, lw)) in enumerate(zip(["ab", "bc", "ca"], estilos))], "tensao_rms")
    ser_f = [(d["f_hz"], "f VSG (rotor virtual)", dict(color="k", lw=1.4))]
    if p["rede"] and c["evento"] in EVENTOS_FREQ:
        ser_f.append((d["f_rede_hz"], "f rede", dict(color="tab:red", lw=1.2, ls="--")))
    fig("Frequência do rotor virtual" + (" e da rede" if len(ser_f) > 1 else ""),
        "Frequência [Hz]", ser_f, "frequencia")
    fig("Corrente RMS de saída do IBR (janela de 1 ciclo)", "Corrente [A rms]",
        [(d["I_rms"][k], f"I_{n},rms", dict(color=cc, ls=ls, lw=lw))
         for k, (n, (cc, ls, lw)) in enumerate(zip("abc", estilos))], "corrente_rms")
    ser_i = [(d["if_env_pu"], "|I_f| envoltória", dict(color="tab:red", lw=1.4)),
             (d["I_rms"][0] / (p["Sb"] / (np.sqrt(3) * p["Vll"])), "I_a rms", dict(color="tab:blue", lw=1.1))]
    if c["imax_pu"] > 0:
        ser_i.append((np.full_like(t, c["imax_pu"]), "I_max", dict(color="k", lw=1.0, ls="--")))
    fig("Corrente do conversor e limite", "Corrente [pu]", ser_i, "corrente_limitador")
    fig("Tensão de fase RMS: terminais do inversor × PCC", "Tensão de fase [V rms]",
        [(d["E_ph"], "V_fase inversor (FEM)", dict(color="tab:brown", lw=1.6)),
         (d["Vll_rms"][0] / np.sqrt(3), "V_fase PCC", dict(color="tab:cyan", lw=1.4, ls="--"))],
        "tensao_fase_inversor")
    ser = [(d["Pf_pu"], "P_f (medida)", dict(color="tab:purple", lw=1.4)),
           (d["Pref_pu"], "P_ref", dict(color="tab:orange", lw=1.4, ls="--"))]
    if p["rede"]:
        ser.append((d["P_rede_pu"], "P para a rede", dict(color="tab:gray", lw=1.2)))
    fig("Potência ativa", "Potência [pu]", ser, "potencia")
    if p["rede"]:
        fig("Ângulo do VSG em relação à rede", "δ_rel [graus]",
            [(np.degrees(d["delta_rel_unwrapped_rad"]), "δ_v − δ_g", dict(color="tab:olive", lw=1.4))], "angulo")
        fig("Corrente RMS no ramo da rede", "Corrente [A rms]",
            [(d["Ig_rms"], "I_rede", dict(color="tab:gray", lw=1.4))], "corrente_rede")
    if c["imax_pu"] > 0:
        fig("Impedância virtual", "Impedância [pu]",
            [(d["rv_pu"], "R_v", dict(color="tab:red", lw=1.3)),
             (d["xv_pu"], "X_v", dict(color="tab:blue", lw=1.2))], "impedancia_virtual")
    if c["evento"] == "sincronizacao":
        fig("Grandezas do sync-check (relé 25)", "Erro normalizado",
            [(d["dv_sync_pu"] / c["sync_dv_max_pu"], "ΔV / janela", dict(color="tab:blue", lw=1.2)),
             (d["df_sync_hz"] / c["sync_df_max_hz"], "Δf / janela", dict(color="tab:orange", lw=1.2)),
             (d["delta_sync_pred_deg"] / c["sync_delta_max_deg"], "δ previsto / janela", dict(color="tab:green", lw=1.2)),
             (d["breaker_closed"], "disjuntor fechado", dict(color="black", lw=1.0, ls="--"))],
            "sync_check")


# =======================================================================================
# Linha de comando
# =======================================================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description="VSG de 2ª ordem ilhado ou conectado à rede (v10). "
                    "Precedência: padrão < arquivo (--config) < linha de comando.")
    ap.add_argument("--config", metavar="ARQUIVO", help="Arquivo .yaml, .yml, .json ou .toml")
    ap.add_argument("--gerar-config", metavar="ARQUIVO", help="Grava modelo YAML com os padrões e encerra")
    ap.add_argument("--sem-graficos", action="store_true", help="Não gera os PNGs")
    grupos = {}
    for sec, nome, pad, tipo, un, desc in PARAM_SPEC:
        g = grupos.setdefault(sec, ap.add_argument_group(f"parâmetros — {sec}"))
        g.add_argument("--" + nome.replace("_", "-"), dest=nome, type=str, default=None,
                       metavar=un if un != "-" else "VAL",
                       help=f"{desc} (padrão: {'auto' if pad is None else pad})")
    try:
        return ap.parse_args()
    except SystemExit as e:
        if e.code not in (0, None):
            sys.exit("Erro: argumentos de linha de comando inválidos (veja --help)")
        raise


def main():
    console_seguro()
    a = parse_args()
    if a.gerar_config:
        with open(a.gerar_config, "w", encoding="utf-8") as fh:
            fh.write(config_to_yaml(DEFAULTS, "Modelo de configuração v10 — valores padrão"))
        print(f"Modelo de configuração gravado em: {a.gerar_config}")
        return
    c, origem = build_config(a)
    carregar_numpy()
    pasta = os.path.dirname(c["prefixo"])
    if pasta:
        os.makedirs(pasta, exist_ok=True)

    p = build_params(c)
    # Invocações legadas, sem qualquer opção/seção v10, conservam também o valor
    # histórico do campo de versão; o modo v10 explícito publica "v10".
    p["versao_saida"] = (VERSAO if p["cc_ativo"] or origem.get("modelo_cc") != "padrão"
                         else "v9")
    if c["evento"] == "fechamento":
        # A manobra parte de uma ilha em regime permanente.  A rede continua
        # modelada, mas seu ramo fica topologicamente aberto até o contato.
        p["rede"] = False
        x0, info = equilibrio(p)
        x0, info = inicializar_lado_cc(c, p, x0, info)
        lam = np.linalg.eigvals(jacobiano(p, x0, evento_inicial(p)))
        vc0 = complex(x0[VC], x0[VC+1])
        p["sync_delta_base"] = math.atan2(vc0.imag, vc0.real) + math.pi/2
        p["vg"] = (abs(vc0)/(np.sqrt(2)*p["Vph"])) * c["vg"]
        p["rede"] = True
        info["P_rede_pu"], info["Q_rede_pu"] = 0.0, 0.0
    else:
        x0, info = equilibrio(p)
        x0, info = inicializar_lado_cc(c, p, x0, info)
        lam = np.linalg.eigvals(jacobiano(p, x0, evento_inicial(p)))
    t, X, Pref, metodo = simulate(p, c, x0)
    d = derived(t, X, Pref, p, c)

    res = resultados(c, p, info, lam, d, metodo)
    report(res, p, lam, float(Pref[-1]), pico_corrente_pu(p, d, c))

    pref_ = c["prefixo"]
    with open(f"{pref_}_config_usada.yaml", "w", encoding="utf-8") as fh:
        fh.write(config_to_yaml(c, "Configuração efetivamente usada (v10)", origem))
    with open(f"{pref_}_resultados.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=2)
    np.savez(f"{pref_}_series.npz", t=d["t"], f_hz=d["f_hz"], P_pu=d["P_pu"], Pf_pu=d["Pf_pu"],
             Pref_pu=d["Pref_pu"], P_rede_pu=d["P_rede_pu"], Q_rede_pu=d["Q_rede_pu"],
             delta_v_rad=d["delta_v_rad"], Vll_rms_v=d["Vll_rms"][0], I_rms_a=d["I_rms"][0],
              f_rede_hz=d["f_rede_hz"], delta_g_rad=d["delta_g_rad"],
              if_env_pu=d["if_env_pu"], rv_pu=d["rv_pu"], xv_pu=d["xv_pu"],
              limitador_ativo=d["limitador_ativo"], delta_aw_rad_s=d["delta_aw_rad_s"],
              delta_rel_rad=d["delta_rel_rad"], delta_rel_unwrapped_rad=d["delta_rel_unwrapped_rad"],
               vg_aplicada_pu=d["vg_aplicada_pu"], breaker_closed=d["breaker_closed"],
               sync_check_ok=d["sync_check_ok"], sync_command=d["sync_command"],
               delta_sync_pred_deg=d["delta_sync_pred_deg"], delta_sync_deg=d["delta_sync_deg"],
               df_sync_hz=d["df_sync_hz"], dv_sync_pu=d["dv_sync_pu"],
               p_sync_pu=d["p_sync_pu"], e_sync_pu=d["e_sync_pu"],
               disjuntor_fechado=d["disjuntor_fechado"],
               comando_fechamento=d["comando_fechamento"],
               sync_estado_codigo=d["sync_estado_codigo"],
               sync_tempo_janela_s=d["sync_tempo_janela_s"],
               delta_sync_rad=d["delta_sync_rad"],
               delta_sync_pred_rad=d["delta_sync_pred_rad"],
               medicao_sync_valida=d["medicao_sync_valida"],
               pre_sync_ativo=d["pre_sync_ativo"],
               df_bias_sync_hz=d["df_bias_sync_hz"],
               de_bias_sync_pu=d["de_bias_sync_pu"],
               theta_pll_rede_rad=d["theta_pll_rede_rad"],
               f_pll_rede_hz=d["f_pll_rede_hz"],
               ig_env_pu=d["ig_env_pu"],
               vdc_v=d["vdc_v"], vdc_pu=d["vdc_pu"], soc_pu=d["soc_pu"], ocv_v=d["ocv_v"],
               i_bat_a=d["i_bat_a"], i_dc_conv_a=d["i_dc_conv_a"], p_conv_w=d["p_conv_w"],
               p_bat_term_w=d["p_bat_term_w"], p_bat_quim_w=d["p_bat_quim_w"],
               p_perda_r0_w=d["p_perda_r0_w"], p_cap_w=d["p_cap_w"], e_cap_j=d["e_cap_j"],
               residuo_potencia_cc_w=d["residuo_potencia_cc_w"], e_req_pu=d["e_req_pu"],
               e_pos_rv_pu=d["e_pos_rv_pu"], e_aplicada_pu=d["e_aplicada_pu"],
               e_max_mod_pu=d["e_max_mod_pu"], m_utilizado_pu=d["m_utilizado_pu"],
               modulacao_saturada=d["modulacao_saturada"], vdc_subtensao=d["vdc_subtensao"],
               vdc_sobretensao=d["vdc_sobretensao"], soc_fora_faixa=d["soc_fora_faixa"],
               i_bat_fora_limite=d["i_bat_fora_limite"], if_complex_a=d["if_complex_a"])
    print(f"Resultados: {pref_}_resultados.json, {pref_}_series.npz, {pref_}_config_usada.yaml")
    if not a.sem_graficos:
        plot_all(d, p, c, pref_)


def _encerrar(codigo):
    """Esvazia as saídas e termina sem a finalização do interpretador (ver v7.0.3)."""
    for s_ in (sys.stdout, sys.stderr):
        try:
            s_.flush()
        except Exception:  # noqa: BLE001
            pass
    try:
        faulthandler.cancel_dump_traceback_later()
    except Exception:  # noqa: BLE001
        pass
    os._exit(codigo)


if __name__ == "__main__":
    codigo = 0
    try:
        main()
    except SystemExit as e:
        if e.code is None or isinstance(e.code, int):
            codigo = e.code or 0
        else:
            print(e.code, file=sys.stderr)
            codigo = 1
    except Exception as e:  # noqa: BLE001
        print(f"Erro inesperado: {type(e).__name__}: {e}", file=sys.stderr)
        codigo = 1
    _encerrar(codigo)

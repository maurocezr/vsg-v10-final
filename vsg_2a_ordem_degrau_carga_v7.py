#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulação EMT (modelo médio) de um IBR trifásico com controle VSM/VSG de 2ª ordem,
em modo ILHADO (carga RL local) ou CONECTADO À REDE (barramento infinito atrás de
impedância definida por SCR e X/R), com eventos de degrau de carga, degrau de P_ref,
salto de fase da rede e — desde a v7 — degrau ou rampa da FREQUÊNCIA da rede.

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
    python vsg_2a_ordem_degrau_carga_v7.py --config meu_caso.yaml
    python vsg_2a_ordem_degrau_carga_v7.py --modo rede --scr 5 --dw 100 --evento pref
    python vsg_2a_ordem_degrau_carga_v7.py --modo rede --scr 2 --evento freq_rampa --df-g -0.5 --rocof-g 0.5
    python vsg_2a_ordem_degrau_carga_v7.py --gerar-config modelo.yaml
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
VERSAO = "v7"


def console_seguro():
    """Evita UnicodeEncodeError quando a saída usa codificação limitada (ex.: cp1252 no
    Windows com saída redirecionada). Caracteres não representáveis viram '?'."""
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(errors="replace")
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
    ("evento",    "evento",   "carga",  str,   "-",
     "Evento: nenhum | carga | pref | fase | freq_degrau | freq_rampa"),
    ("evento",    "d_pref",   0.05,     float, "pu",   "Degrau de P_ref (evento pref)"),
    ("evento",    "d_fase",   5.0,      float, "graus", "Salto de fase da rede (evento fase)"),
    ("evento",    "df_g",     0.0,      float, "Hz",   "Desvio final da frequência da rede, com sinal (freq_*)"),
    ("evento",    "rocof_g",  1.0,      float, "Hz/s", "Taxa da rampa de frequência da rede (freq_rampa)"),
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
EVENTOS = ("nenhum", "carga", "pref", "fase", "freq_degrau", "freq_rampa")
EVENTOS_FREQ = ("freq_degrau", "freq_rampa")
LIMITE_I_PU = 1.2

# Índices do vetor de estados
IF, VC, I1, I2, IG = 0, 2, 4, 6, 8          # pares (re, im)
DV, DW, PF, QF, Z = 10, 11, 12, 13, 14
NX = 15


def erro(msg):
    sys.exit(f"Erro: {msg}")


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
    if nome in ("pref", "qref") and (valor is None or str(valor).strip().lower()
                                     in ("auto", "none", "null", "")):
        return None
    if tipo is str:
        return str(valor).strip()
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
    if c["evento"] in EVENTOS_FREQ:
        if c["modo"] != "rede":
            e.append(f"o evento '{c['evento']}' exige modo = 'rede'")
        if c["df_g"] == 0:
            e.append(f"o evento '{c['evento']}' exige df_g != 0")
    if c["rocof_g"] <= 0:
        e.append(f"'rocof_g' deve ser > 0 (valor: {c['rocof_g']})")
    for k in ("sn", "vll", "f0", "xf", "bc", "p1", "p2", "eref", "H", "mp", "fc",
              "t_step", "t_end", "dt_out", "scr", "xr_rede", "vg", "tw", "max_step", "rtol"):
        if c[k] <= 0:
            e.append(f"'{k}' deve ser > 0 (valor: {c[k]})")
    for k in ("rf", "nq", "q1", "q2", "dw"):
        if c[k] < 0:
            e.append(f"'{k}' deve ser >= 0 (valor: {c[k]})")
    if c["t_end"] <= c["t_step"]:
        e.append(f"'t_end' ({c['t_end']}) deve ser maior que 't_step' ({c['t_step']})")
    if c["dt_out"] > 1.0 / (20 * c["f0"]):
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
    p["Eref"], p["Pref"], p["Qref"] = c["eref"], c["pref"], c["qref"]
    p["t_step"] = c["t_step"]
    return p


# =======================================================================================
# Modelo dinâmico (referencial síncrono em ω0)
# =======================================================================================
def _c(x, i):
    return x[i] + 1j * x[i + 1]


def fem(E, delta, p):
    return -1j * np.sqrt(2) * E * p["Vph"] * np.exp(1j * delta)


def v_rede(p, delta_g):
    return -1j * np.sqrt(2) * p["vg"] * p["Vph"] * np.exp(1j * delta_g)


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
    Ig_ = _c(x, IG) if p["rede"] else 0j
    dv, dwv, Pf, Qf, z = x[DV], x[DW], x[PF], x[QF], x[Z]
    I1_, I2_ = correntes_carga(x, p, ev["load2_on"])

    E = p["Eref"] - p["nq"] * (Qf - p["Qref"])
    e = fem(E, dv, p)
    dIf = (e - Vc - p["Rf"] * If_) / p["Lf"] - 1j * w0 * If_
    dI1 = ((Vc - p["RL1"] * I1_) / p["LL1"] - 1j * w0 * I1_) if p["LL1"] > 0 else 0j
    dI2 = ((Vc - p["RL2"] * I2_) / p["LL2"] - 1j * w0 * I2_) \
        if (ev["load2_on"] and p["LL2"] > 0) else 0j
    dIg = ((Vc - p["Rg"] * Ig_ - v_rede(p, delta_g_de(t, ev))) / p["Lg"] - 1j * w0 * Ig_) \
        if p["rede"] else 0j
    dVc = (If_ - I1_ - I2_ - Ig_) / p["Cf"] - 1j * w0 * Vc

    S = 1.5 * Vc * np.conj(If_) / p["Sb"]
    ddw = (ev["Pref"] - Pf - p["Dp"] * dwv - p["Dw"] * (dwv - z)) / (2 * p["H"])
    out = np.empty(NX)
    out[IF], out[IF + 1] = dIf.real, dIf.imag
    out[VC], out[VC + 1] = dVc.real, dVc.imag
    out[I1], out[I1 + 1] = dI1.real, dI1.imag
    out[I2], out[I2 + 1] = dI2.real, dI2.imag
    out[IG], out[IG + 1] = dIg.real, dIg.imag
    out[DV] = w0 * dwv
    out[DW] = ddw
    out[PF] = p["wc"] * (S.real - Pf)
    out[QF] = p["wc"] * (S.imag - Qf)
    out[Z] = (dwv - z) / p["Tw"]
    return out


def estados_ativos(p, load2_on):
    a = [IF, IF + 1, VC, VC + 1]
    if p["LL1"] > 0:
        a += [I1, I1 + 1]
    if load2_on and p["LL2"] > 0:
        a += [I2, I2 + 1]
    if p["rede"]:
        a += [IG, IG + 1]
    return a + [DV, DW, PF, QF, Z]


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
        if ier != 1:
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
    return {"load2_on": False, "Pref": p["Pref"], "delta_g": 0.0, "freq": None}


def evento_apos(p, c):
    ev = evento_inicial(p)
    if c["evento"] == "carga":
        ev["load2_on"] = True
    elif c["evento"] == "pref":
        ev["Pref"] = p["Pref"] + c["d_pref"]
    elif c["evento"] == "fase":
        ev["delta_g"] = np.radians(c["d_fase"])
    elif c["evento"] in EVENTOS_FREQ:
        ev["freq"] = (c["evento"], c["t_step"], c["df_g"], c["rocof_g"])
    return ev


def pontos_de_quebra(c):
    """Limites dos segmentos de integração: 0, t_step, fim da rampa (se houver), t_end."""
    br = [0.0, c["t_step"]]
    if c["evento"] == "freq_rampa":
        t_fim = c["t_step"] + abs(c["df_g"]) / c["rocof_g"]
        if t_fim < c["t_end"] - 1e-12:
            br.append(t_fim)
    return br + [c["t_end"]]


def analise_modal(c):
    """Autovalores no equilíbrio inicial (usado também pelo gerador de casos)."""
    carregar_numpy()
    p = build_params(c)
    x0, info = equilibrio(p)
    lam = np.linalg.eigvals(jacobiano(p, x0, evento_inicial(p)))
    return lam, p, info


# =======================================================================================
# Simulação
# =======================================================================================
def simulate(p, c, x0):
    from scipy.integrate import solve_ivp
    ev0, ev1 = evento_inicial(p), evento_apos(p, c)
    metodo = "RK45"
    T, X, PR = [], [], []
    br = pontos_de_quebra(c)
    segs = [(br[0], br[1], ev0)] + [(a, b, ev1) for a, b in zip(br[1:-1], br[2:])]
    for ta, tb, ev in segs:
        t_eval = np.linspace(ta, tb, int(round((tb - ta) / c["dt_out"])) + 1)
        sol = solve_ivp(rhs, (ta, tb), x0, args=(p, ev), method=metodo, t_eval=t_eval,
                        max_step=c["max_step"], rtol=c["rtol"], atol=1e-8)
        if not sol.success:
            metodo = "LSODA"
            sol = solve_ivp(rhs, (ta, tb), x0, args=(p, ev), method=metodo, t_eval=t_eval,
                            max_step=c["max_step"], rtol=c["rtol"], atol=1e-8)
            if not sol.success:
                raise RuntimeError(sol.message)
        T.append(sol.t[:-1])
        X.append(sol.y[:, :-1])
        PR.append(np.full(len(sol.t) - 1, ev["Pref"]))
        x0 = sol.y[:, -1].copy()
    return np.concatenate(T), np.concatenate(X, axis=1), np.concatenate(PR), metodo


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
    if p["rede"]:
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
    return {"t": t, "f_hz": p["f0"] * (1 + X[DW]), "P_pu": S.real, "Pf_pu": X[PF],
            "Qf_pu": X[QF], "Pref_pu": Pref, "P_rede_pu": Sg.real, "Q_rede_pu": Sg.imag,
            "delta_v_rad": X[DV], "Vll_rms": rms_one_cycle(t, v_ll, theta),
            "I_rms": rms_one_cycle(t, i, theta),
            "Ig_rms": np.abs(Ig) / np.sqrt(2),
            "f_rede_hz": f_rede, "delta_g_rad": d_g,
            "E_ph": (p["Eref"] - p["nq"] * (X[QF] - p["Qref"])) * p["Vph"]}


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
    return {"H_medido_s": H_med,"H_eff_s": float(p["H"] + p["Dw"] * p["Tw"] / 2),
            "dP_inercial_pu": a, "inclinacao_droop_pu_s": b,
            "dP_previsto_pu": (float(2 * p["H"] * c["rocof_g"] / p["f0"])
                               if c["evento"] == "freq_rampa" else None),
            "P_pico_pu": float(np.max(d["Pf_pu"][apos])),
            "I_pico_pu": float(np.nanmax(d["I_rms"][0][apos]) / I_n),
            "janela_ajuste_s": float(dur), "f_oscilacao_hz": f_osc, "metodo_ajuste": metodo}


def pico_corrente_pu(p, d, c):
    apos = d["t"] >= c["t_step"]
    return float(np.nanmax(d["I_rms"][0][apos]) / (p["Sb"] / (np.sqrt(3) * p["Vll"])))


# =======================================================================================
# Relatório, arquivos e gráficos
# =======================================================================================
def resultados(c, p, info, lam, d, metodo):
    t, ts = d["t"], c["t_step"]
    m_pre = (t > ts - 0.1) & (t < ts)
    m_end = t > t[-1] - 0.1
    dfdt = np.gradient(d["f_hz"], t)
    return {"versao": VERSAO, "modo": c["modo"], "parametros": c,
            "controle_efetivo": {"pref": float(p["Pref"]), "qref": float(p["Qref"]),
                                 "eref": float(p["Eref"]), "dw": float(p["Dw"]), "tw": float(p["Tw"])},
            "equilibrio": info, "antes_evento": _medias(d, m_pre), "final": _medias(d, m_end),
            "rocof_max_hz_s": float(np.abs(dfdt[t > ts + 1e-3]).max()),
            "modal": {"autovalores": [[float(l.real), float(l.imag)] for l in lam]},
            "integrador": metodo,
            "inercial": medir_inercia(c, p, d, lam)}


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
        print("=== Resposta inercial (v7) ===")
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
    if i_pico > LIMITE_I_PU:
        print(f"(!) Pico de corrente de {i_pico:.2f} pu > {LIMITE_I_PU} pu após o evento: "
              f"o modelo não tem limitação de corrente (v8); resultado otimista.")


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
        fig("Ângulo do VSG em relação à rede", "δ_v [graus]",
            [(np.degrees(d["delta_v_rad"]), "δ_v", dict(color="tab:olive", lw=1.4))], "angulo")
        fig("Corrente RMS no ramo da rede", "Corrente [A rms]",
            [(d["Ig_rms"], "I_rede", dict(color="tab:gray", lw=1.4))], "corrente_rede")


# =======================================================================================
# Linha de comando
# =======================================================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description="VSG de 2ª ordem ilhado ou conectado à rede (v7). "
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
            fh.write(config_to_yaml(DEFAULTS, "Modelo de configuração v7 — valores padrão"))
        print(f"Modelo de configuração gravado em: {a.gerar_config}")
        return
    c, origem = build_config(a)
    carregar_numpy()
    pasta = os.path.dirname(c["prefixo"])
    if pasta:
        os.makedirs(pasta, exist_ok=True)

    p = build_params(c)
    x0, info = equilibrio(p)
    lam = np.linalg.eigvals(jacobiano(p, x0, evento_inicial(p)))
    t, X, Pref, metodo = simulate(p, c, x0)
    d = derived(t, X, Pref, p, c)

    res = resultados(c, p, info, lam, d, metodo)
    report(res, p, lam, float(Pref[-1]), pico_corrente_pu(p, d, c))

    pref_ = c["prefixo"]
    with open(f"{pref_}_config_usada.yaml", "w", encoding="utf-8") as fh:
        fh.write(config_to_yaml(c, "Configuração efetivamente usada (v7)", origem))
    with open(f"{pref_}_resultados.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=2)
    np.savez(f"{pref_}_series.npz", t=d["t"], f_hz=d["f_hz"], P_pu=d["P_pu"], Pf_pu=d["Pf_pu"],
             Pref_pu=d["Pref_pu"], P_rede_pu=d["P_rede_pu"], Q_rede_pu=d["Q_rede_pu"],
             delta_v_rad=d["delta_v_rad"], Vll_rms_v=d["Vll_rms"][0], I_rms_a=d["I_rms"][0],
             f_rede_hz=d["f_rede_hz"], delta_g_rad=d["delta_g_rad"])
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

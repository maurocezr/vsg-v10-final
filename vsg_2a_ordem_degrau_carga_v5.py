#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulação EMT (modelo médio) de um IBR trifásico com controle VSM/VSG de 2ª ordem
alimentando uma carga linear RL trifásica em modo ilhado, com degrau de carga.

v5.0.2 (portabilidade Windows; sem alteração de cálculo):
  - Saída de console segura em qualquer codificação (cp1252 incluso).
  - config_usada.yaml grava strings entre aspas simples (caminhos com '\').

v5:
  - Novos parâmetros E_ref (eref) e Q_ref (qref; null = casado com a carga inicial).
  - Cargas puramente resistivas (Q = 0) passam a ser aceitas (corrente algébrica).
  - Relatório e gráfico da tensão de fase nos terminais do inversor (FEM interna).
  - Compatível com os arquivos gerados por gerar_caso_vsg.py.

v4:
  - Parâmetros podem ser lidos de um arquivo de configuração (YAML, JSON ou TOML).
  - Todo parâmetro tem valor padrão; o arquivo só precisa conter o que se quer alterar.
  - Precedência: valores padrão  <  arquivo (--config)  <  linha de comando.
  - --gerar-config cria um modelo de arquivo YAML comentado com todos os padrões.
  - A configuração efetivamente usada é salva junto dos gráficos (<prefixo>_config_usada.yaml).

v3 (mantido):
  - Cada gráfico é salvo em um PNG separado (nenhuma janela é aberta).
  - Tensão e corrente em valor RMS (janela deslizante de um ciclo, sincronizada com θ).

Topologia (por fase, estrela aterrada):

    e_abc (VSG) --[R_f + L_f]--+-- v_abc (PCC) --+----------------+
                               |                 |                |
                              C_f           [R_L1 + L_L1]   S --[R_L2 + L_L2]
                               |                 |                |   (fecha em t_step)

Controle VSG (p.u.):
    2H dΔω/dt = P_ref - P_f - D_p Δω          (swing; D_p = 1/m_p -> droop P-f)
    dθ/dt     = ω0 (1 + Δω)
    E         = E_ref - n_q (Q_f - Q_ref)     (droop Q-V)
    P_f, Q_f  : potências medidas no PCC, filtro passa-baixa de 1ª ordem (ω_c)

Uso:
    python vsg_2a_ordem_degrau_carga_v5.py                          # só padrões
    python vsg_2a_ordem_degrau_carga_v5.py --gerar-config caso.yaml # cria modelo
    python vsg_2a_ordem_degrau_carga_v5.py --config caso.yaml       # lê o arquivo
    python vsg_2a_ordem_degrau_carga_v5.py --config caso.yaml --H 2 # arquivo + ajuste
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")                      # backend sem interface gráfica
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

SHIFT = np.array([0.0, -2 * np.pi / 3, 2 * np.pi / 3])


def console_seguro():
    """Evita UnicodeEncodeError em consoles/pipes com codificação limitada (Windows cp1252)."""
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


# =======================================================================================
# Definição central dos parâmetros: (seção, nome, padrão, tipo, unidade, descrição)
# =======================================================================================
PARAM_SPEC = [
    ("sistema",   "sn",      100e3,  float, "VA",   "Potência nominal do IBR (base)"),
    ("sistema",   "vll",     380.0,  float, "V",    "Tensão de linha nominal (rms)"),
    ("sistema",   "f0",      60.0,   float, "Hz",   "Frequência nominal"),
    ("filtro",    "rf",      0.01,   float, "pu",   "Resistência série do filtro"),
    ("filtro",    "xf",      0.15,   float, "pu",   "Reatância série do filtro (L_f)"),
    ("filtro",    "bc",      0.05,   float, "pu",   "Susceptância do capacitor (C_f)"),
    ("carga",     "p1",      50e3,   float, "W",    "Carga inicial - potência ativa trifásica"),
    ("carga",     "q1",      15e3,   float, "var",  "Carga inicial - potência reativa trifásica"),
    ("carga",     "p2",      30e3,   float, "W",    "Degrau de carga - potência ativa trifásica"),
    ("carga",     "q2",      10e3,   float, "var",  "Degrau de carga - potência reativa trifásica"),
    ("controle",  "pref",    None,   float, "pu",   "Setpoint de potência ativa (null = casado com a carga inicial)"),
    ("controle",  "qref",    None,   float, "pu",   "Setpoint de potência reativa (null = casado com a carga inicial)"),
    ("controle",  "eref",    1.0,    float, "pu",   "Referência da tensão interna do inversor (E_ref)"),
    ("controle",  "H",       5.0,    float, "s",    "Constante de inércia virtual"),
    ("controle",  "mp",      0.05,   float, "pu",   "Droop P-f (D_p = 1/m_p)"),
    ("controle",  "nq",      0.05,   float, "pu",   "Droop Q-V"),
    ("controle",  "fc",      10.0,   float, "Hz",   "Frequência de corte do filtro de medição de P e Q"),
    ("simulacao", "t_step",  1.0,    float, "s",    "Instante do degrau de carga"),
    ("simulacao", "t_end",   4.0,    float, "s",    "Tempo total de simulação"),
    ("simulacao", "dt_out",  5e-5,   float, "s",    "Passo de amostragem da saída"),
    ("saida",     "prefixo", "vsg",  str,   "-",    "Prefixo dos arquivos gerados (pode incluir pasta)"),
]
SPEC = {nome: (sec, pad, tipo, un, desc) for sec, nome, pad, tipo, un, desc in PARAM_SPEC}
DEFAULTS = {nome: v[1] for nome, v in SPEC.items()}


# =======================================================================================
# Leitura, mescla e validação da configuração
# =======================================================================================
def _flatten(d, origem):
    plano = {}
    for k, v in d.items():
        if isinstance(v, dict):
            plano.update(_flatten(v, origem))
        else:
            if k in plano:
                raise ValueError(f"{origem}: parâmetro '{k}' definido mais de uma vez")
            plano[k] = v
    return plano


def load_config_file(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError:
            sys.exit("Erro: PyYAML não instalado. Use um arquivo .json ou .toml, "
                     "ou instale com 'pip install pyyaml'.")
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
        sys.exit(f"Erro: extensão '{ext}' não suportada (use .yaml, .yml, .json ou .toml).")
    if not isinstance(dados, dict):
        sys.exit(f"Erro: o arquivo {path} deve conter pares chave: valor.")
    return _flatten(dados, path)


def _coerce(nome, valor, origem):
    _, _, tipo, _, _ = SPEC[nome]
    if nome in ("pref", "qref") and (valor is None or str(valor).strip().lower() in ("auto", "none", "null", "")):
        return None
    try:
        return tipo(valor)
    except (TypeError, ValueError):
        sys.exit(f"Erro: {origem}: '{nome}' = {valor!r} não é um valor {tipo.__name__} válido")


def validate(cfg):
    erros = []
    positivos = ["sn", "vll", "f0", "xf", "bc", "p1", "p2", "eref",
                 "H", "mp", "fc", "t_step", "t_end", "dt_out"]
    for k in positivos:
        if cfg[k] <= 0:
            erros.append(f"'{k}' deve ser > 0 (valor: {cfg[k]})")
    for k in ("rf", "nq", "q1", "q2"):
        if cfg[k] < 0:
            erros.append(f"'{k}' deve ser >= 0 (valor: {cfg[k]})")
    if cfg["t_end"] <= cfg["t_step"]:
        erros.append(f"'t_end' ({cfg['t_end']}) deve ser maior que 't_step' ({cfg['t_step']})")
    if cfg["dt_out"] > 1.0 / (20 * cfg["f0"]):
        erros.append(f"'dt_out' ({cfg['dt_out']}) muito grande para calcular o RMS; "
                     f"use <= {1.0 / (20 * cfg['f0']):.2e} s")
    if erros:
        sys.exit("Erro(s) de configuração:\n  - " + "\n  - ".join(erros))


def build_config(args):
    cfg = dict(DEFAULTS)
    origem = {k: "padrão" for k in cfg}
    if args.config:
        arquivo = load_config_file(args.config)
        desconhecidos = sorted(set(arquivo) - set(SPEC))
        if desconhecidos:
            sys.exit(f"Erro: parâmetro(s) desconhecido(s) em {args.config}: "
                     f"{', '.join(desconhecidos)}\nVálidos: {', '.join(SPEC)}")
        for k, v in arquivo.items():
            cfg[k] = _coerce(k, v, args.config)
            origem[k] = "arquivo"
    for k in SPEC:
        v = getattr(args, k)
        if v is not None:
            cfg[k] = _coerce(k, v, "linha de comando")
            origem[k] = "linha de comando"
    validate(cfg)
    return cfg, origem


def _yaml_value(v):
    if v is None:
        return "null"
    if isinstance(v, str):
        return "'" + v.replace("'", "''") + "'"      # aspas simples: '\' literal
    return repr(v)


def config_to_yaml(cfg, titulo, origem=None):
    linhas = [f"# {titulo}",
              "# Formato: YAML. Apenas os parâmetros que se deseja alterar precisam constar;",
              "# os demais assumem os valores padrão do script.",
              ""]
    secao_atual = None
    for sec, nome, _, _, un, desc in PARAM_SPEC:
        if sec != secao_atual:
            if secao_atual is not None:
                linhas.append("")
            linhas.append(f"{sec}:")
            secao_atual = sec
        coment = f"# [{un}] {desc}"
        if origem is not None:
            coment += f"  (origem: {origem[nome]})"
        linhas.append(f"  {nome + ':':<9} {_yaml_value(cfg[nome]):<12} {coment}")
    return "\n".join(linhas) + "\n"


def print_config(cfg, origem):
    print("=== Parâmetros da simulação ===")
    for sec, nome, _, _, un, _ in PARAM_SPEC:
        val = "auto (casado com a carga inicial)" if (nome in ("pref", "qref") and cfg[nome] is None) \
            else cfg[nome]
        marca = "" if origem[nome] == "padrão" else f"   <- {origem[nome]}"
        print(f"  {sec:<10} {nome:<8} = {val} {'' if un == '-' else un}{marca}")


# =======================================================================================
# Parâmetros elétricos derivados
# =======================================================================================
def rl_from_pq(P, Q, Vll, w):
    S2 = P ** 2 + Q ** 2
    return Vll ** 2 * P / S2, (Vll ** 2 * Q / S2) / w


def build_params(c):
    p = {"Sb": c["sn"], "Vll": c["vll"], "f0": c["f0"], "w0": 2 * np.pi * c["f0"]}
    p["Vph"] = c["vll"] / np.sqrt(3)
    p["Zb"] = c["vll"] ** 2 / c["sn"]
    p["Rf"] = c["rf"] * p["Zb"]
    p["Lf"] = c["xf"] * p["Zb"] / p["w0"]
    p["Cf"] = c["bc"] / (p["Zb"] * p["w0"])
    p["RL1"], p["LL1"] = rl_from_pq(c["p1"], c["q1"], c["vll"], p["w0"])
    p["RL2"], p["LL2"] = rl_from_pq(c["p2"], c["q2"], c["vll"], p["w0"])
    p["H"] = c["H"]
    p["Dp"] = 1.0 / c["mp"]
    p["nq"] = c["nq"]
    p["wc"] = 2 * np.pi * c["fc"]
    p["Eref"] = c["eref"]
    p["Qref"] = c["qref"]
    p["t_step"] = c["t_step"]
    p["Pref"] = c["pref"]
    return p


# =======================================================================================
# Modelo dinâmico
# =======================================================================================
def emf_sv(E, p):
    return -1j * np.sqrt(2) * E * p["Vph"]


def rhs(t, x, p, load2_on):
    If, Vc = x[0] + 1j * x[1], x[2] + 1j * x[3]
    I1, I2 = x[4] + 1j * x[5], x[6] + 1j * x[7]
    dw, Pf, Qf = x[9], x[10], x[11]
    w = p["w0"] * (1 + dw)
    E = p["Eref"] - p["nq"] * (Qf - p["Qref"])
    e = emf_sv(E, p)
    dIf = (e - Vc - p["Rf"] * If) / p["Lf"] - 1j * w * If
    if p["LL1"] > 0:
        dI1 = (Vc - p["RL1"] * I1) / p["LL1"] - 1j * w * I1
    else:
        I1, dI1 = Vc / p["RL1"], 0j
    if not load2_on:
        I2, dI2 = 0j, 0j
    elif p["LL2"] > 0:
        dI2 = (Vc - p["RL2"] * I2) / p["LL2"] - 1j * w * I2
    else:
        I2, dI2 = Vc / p["RL2"], 0j
    dVc = (If - I1 - I2) / p["Cf"] - 1j * w * Vc
    S = 1.5 * Vc * np.conj(If) / p["Sb"]
    dPf = p["wc"] * (S.real - Pf)
    dQf = p["wc"] * (S.imag - Qf)
    ddw = (p["Pref"] - Pf - p["Dp"] * dw) / (2 * p["H"])
    return np.array([dIf.real, dIf.imag, dVc.real, dVc.imag, dI1.real, dI1.imag,
                     dI2.real, dI2.imag, w, ddw, dPf, dQf])


def to_abc(Xsv, theta):
    return np.real((Xsv * np.exp(1j * theta))[None, :] * np.exp(1j * SHIFT)[:, None])


def solve_phasors(E, dw, p):
    w = p["w0"] * (1 + dw)
    Zf = p["Rf"] + 1j * w * p["Lf"]
    Z1 = p["RL1"] + 1j * w * p["LL1"]
    Ypar = 1j * w * p["Cf"] + 1 / Z1
    e = emf_sv(E, p)
    Vc = e / (1 + Zf * Ypar)
    If = (e - Vc) / Zf
    S = 1.5 * Vc * np.conj(If) / p["Sb"]
    return If, Vc, Vc / Z1, S


def steady_state_init(p):
    auto_p = p["Pref"] is None
    auto_q = p["Qref"] is None
    E, dw = p["Eref"], 0.0
    for _ in range(500):
        If, Vc, I1, S = solve_phasors(E, dw, p)
        E_new = p["Eref"] if auto_q else p["Eref"] - p["nq"] * (S.imag - p["Qref"])
        dw_new = 0.0 if auto_p else (p["Pref"] - S.real) / p["Dp"]
        if abs(E_new - E) < 1e-13 and abs(dw_new - dw) < 1e-13:
            break
        E, dw = E_new, dw_new
    if auto_p:
        p["Pref"] = S.real
    if auto_q:
        p["Qref"] = S.imag
    return np.array([If.real, If.imag, Vc.real, Vc.imag, I1.real, I1.imag, 0.0, 0.0,
                     0.0, dw, S.real, S.imag])


def simulate(p, t_end, dt_out):
    x0 = steady_state_init(p)
    T, X = [], []
    for ta, tb, on in [(0.0, p["t_step"], False), (p["t_step"], t_end, True)]:
        t_eval = np.linspace(ta, tb, int(round((tb - ta) / dt_out)) + 1)
        sol = solve_ivp(rhs, (ta, tb), x0, args=(p, on), method="RK45", t_eval=t_eval,
                        max_step=1e-3, rtol=1e-8, atol=1e-8)
        if not sol.success:
            raise RuntimeError(sol.message)
        T.append(sol.t[:-1])
        X.append(sol.y[:, :-1])
        x0 = sol.y[:, -1]
    return np.concatenate(T), np.concatenate(X, axis=1)


def rms_one_cycle(t, x, theta):
    t_ini = np.interp(theta - 2 * np.pi, theta, t, left=np.nan)
    valid = ~np.isnan(t_ini)
    out = np.full(x.shape, np.nan)
    for k in range(x.shape[0]):
        C = cumulative_trapezoid(x[k] ** 2, t, initial=0.0)
        C_ini = np.interp(t_ini[valid], t, C)
        out[k, valid] = np.sqrt((C[valid] - C_ini) / (t[valid] - t_ini[valid]))
    return out


def derived_quantities(t, X, p):
    theta = X[8]
    v = to_abc(X[2] + 1j * X[3], theta)
    i = to_abc(X[0] + 1j * X[1], theta)
    v_ll = np.vstack([v[0] - v[1], v[1] - v[2], v[2] - v[0]])
    return {"t": t,
            "Vll_rms": rms_one_cycle(t, v_ll, theta),
            "I_rms": rms_one_cycle(t, i, theta),
            "f": p["f0"] * (1 + X[9]),
            "Pf": X[10],
            "Qf": X[11],
            "E_ph": (p["Eref"] - p["nq"] * (X[11] - p["Qref"])) * p["Vph"]}


def report(d, p):
    t, f = d["t"], d["f"]
    ts = p["t_step"]
    m_pre = (t > ts - 0.1) & (t < ts)
    m_end = t > t[-1] - 0.1
    dfdt = np.gradient(f, t)
    print("=== Parâmetros elétricos derivados ===")
    print(f"Lf = {p['Lf']*1e3:.3f} mH, Rf = {p['Rf']*1e3:.2f} mΩ, Cf = {p['Cf']*1e6:.1f} µF")
    print(f"Carga 1: R = {p['RL1']:.3f} Ω, L = {p['LL1']*1e3:.3f} mH")
    print(f"Carga 2: R = {p['RL2']:.3f} Ω, L = {p['LL2']*1e3:.3f} mH (entra em t = {ts} s)")
    print(f"P_ref = {p['Pref']:.4f} pu ({p['Pref'] * p['Sb'] / 1e3:.1f} kW) | "
          f"Q_ref = {p['Qref']:.4f} pu ({p['Qref'] * p['Sb'] / 1e3:.1f} kvar) | "
          f"E_ref = {p['Eref']:.4f} pu")
    print("=== Resultados ===")
    for nome, m in [("Antes do degrau", m_pre), ("Regime final   ", m_end)]:
        print(f"{nome}: V_LL,rms = {np.nanmean(d['Vll_rms'][0, m]):7.2f} V | "
              f"I_rms = {np.nanmean(d['I_rms'][0, m]):7.2f} A | "
              f"P_f = {d['Pf'][m].mean():.4f} pu | f = {f[m].mean():.4f} Hz | "
              f"V_fase,inv = {d['E_ph'][m].mean():.2f} V")
    print(f"RoCoF máx. após degrau = {np.abs(dfdt[t > ts + 1e-3]).max():.3f} Hz/s")
    print(f"Verificação droop (regime final): f0·(P_ref − P_f)/D_p = "
          f"{p['f0'] * (p['Pref'] - d['Pf'][m_end].mean()) / p['Dp']:+.4f} Hz")


PHASE_STYLE = [("tab:red", "-", 2.2), ("tab:green", "--", 1.6), ("tab:blue", ":", 1.6)]


def _new_fig(p, titulo):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_title(f"{titulo}\nVSG 2ª ordem ilhado — P_ref = {p['Pref']:.3f} pu, "
                 f"H = {p['H']} s, m_p = {100 / p['Dp']:.1f} %, n_q = {100 * p['nq']:.1f} %",
                 fontsize=11)
    ax.axvline(p["t_step"], color="gray", ls="--", lw=0.8, label="degrau de carga")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Tempo [s]")
    return fig, ax


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Gráfico salvo em: {path}")


def plot_all(d, p, prefixo):
    t = d["t"]
    fig, ax = _new_fig(p, "Tensão de linha RMS no PCC (janela de 1 ciclo)")
    for k, nome in enumerate(["ab", "bc", "ca"]):
        c, ls, lw = PHASE_STYLE[k]
        ax.plot(t, d["Vll_rms"][k], color=c, ls=ls, lw=lw, label=f"V_{nome},rms")
    ax.set_ylabel("Tensão [V rms]")
    ax.legend(loc="best", fontsize=9)
    _save(fig, f"{prefixo}_tensao_rms.png")

    fig, ax = _new_fig(p, "Frequência do rotor virtual")
    ax.plot(t, d["f"], color="k", lw=1.4, label="f")
    ax.set_ylabel("Frequência [Hz]")
    ax.legend(loc="best", fontsize=9)
    _save(fig, f"{prefixo}_frequencia.png")

    fig, ax = _new_fig(p, "Corrente RMS de saída do IBR (janela de 1 ciclo)")
    for k, nome in enumerate("abc"):
        c, ls, lw = PHASE_STYLE[k]
        ax.plot(t, d["I_rms"][k], color=c, ls=ls, lw=lw, label=f"I_{nome},rms")
    ax.set_ylabel("Corrente [A rms]")
    ax.legend(loc="best", fontsize=9)
    _save(fig, f"{prefixo}_corrente_rms.png")

    fig, ax = _new_fig(p, "Tensão de fase RMS: terminais do inversor × PCC")
    ax.plot(t, d["E_ph"], color="tab:brown", lw=1.6, label="V_fase inversor (FEM interna)")
    ax.plot(t, d["Vll_rms"][0] / np.sqrt(3), color="tab:cyan", lw=1.4, ls="--", label="V_fase PCC")
    ax.set_ylabel("Tensão de fase [V rms]")
    ax.legend(loc="best", fontsize=9)
    _save(fig, f"{prefixo}_tensao_fase_inversor.png")

    fig, ax = _new_fig(p, "Potência ativa medida (P_f) e setpoint (P_ref)")
    ax.plot(t, d["Pf"], color="tab:purple", lw=1.4, label="P_f (medida)")
    ax.axhline(p["Pref"], color="tab:orange", ls="--", lw=1.4, label="P_ref")
    ax.set_ylabel("Potência [pu]")
    ax.legend(loc="best", fontsize=9)
    _save(fig, f"{prefixo}_potencia.png")


def parse_args():
    ap = argparse.ArgumentParser(
        description="VSG de 2ª ordem + carga RL trifásica com degrau. "
                    "Precedência: padrão < arquivo (--config) < linha de comando.")
    ap.add_argument("--config", metavar="ARQUIVO",
                    help="Arquivo de parâmetros (.yaml, .yml, .json ou .toml)")
    ap.add_argument("--gerar-config", metavar="ARQUIVO",
                    help="Grava um modelo YAML comentado com os valores padrão e encerra")
    grupos = {}
    for sec, nome, pad, tipo, un, desc in PARAM_SPEC:
        g = grupos.setdefault(sec, ap.add_argument_group(f"parâmetros — {sec}"))
        flag = "--" + nome.replace("_", "-")
        g.add_argument(flag, dest=nome, type=str if nome in ("pref", "qref") else tipo, default=None,
                       metavar=un if un != "-" else "TXT",
                       help=f"{desc} (padrão: {'auto' if pad is None else pad})")
    return ap.parse_args()


def main():
    console_seguro()
    args = parse_args()
    if args.gerar_config:
        with open(args.gerar_config, "w", encoding="utf-8") as fh:
            fh.write(config_to_yaml(DEFAULTS, "Modelo de configuração — valores padrão"))
        print(f"Modelo de configuração gravado em: {args.gerar_config}")
        return
    cfg, origem = build_config(args)
    print_config(cfg, origem)
    pasta = os.path.dirname(cfg["prefixo"])
    if pasta:
        os.makedirs(pasta, exist_ok=True)
    cfg_path = f"{cfg['prefixo']}_config_usada.yaml"
    with open(cfg_path, "w", encoding="utf-8") as fh:
        fh.write(config_to_yaml(cfg, "Configuração efetivamente usada nesta simulação", origem))
    print(f"Configuração usada salva em: {cfg_path}")
    p = build_params(cfg)
    t, X = simulate(p, cfg["t_end"], cfg["dt_out"])
    d = derived_quantities(t, X, p)
    report(d, p)
    plot_all(d, p, cfg["prefixo"])


if __name__ == "__main__":
    main()

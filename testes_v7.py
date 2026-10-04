#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
testes_v7.py — Suíte de aceitação CONGELADA da v7 do simulador VSG
                (frequência da rede programável: degrau e rampa).
=====================================================================================
REGRA DE CONGELAMENTO
=====================================================================================
Este arquivo foi escrito ANTES da implementação da v7 e define, de forma fixa:
  (1) o CONTRATO de interface que a v7 deve cumprir (seção CONTRATO abaixo);
  (2) os testes de aceitação R1–R3 e B1–B10 e suas tolerâncias;
  (3) um ORÁCULO independente: simulação NÃO LINEAR escrita aqui, num referencial que
      gira com a REDE (e não fixo em ω0, como a implementação), sem reutilizar código
      da implementação. Como a formulação é outra, erros de convenção não se cancelam.
Durante a implementação da v7 este arquivo NÃO pode ser alterado. A integridade é
verificada pelo hash SHA-256 publicado em testes_v7.sha256 (e registrado fora do
ambiente pelo usuário). Qualquer mudança exige uma nova versão (testes_v7_1.py) com
justificativa escrita — nunca uma edição silenciosa deste arquivo.

Diferenças em relação ao tasks_v7.md, decididas na elaboração desta suíte (com base no
oráculo) e por isso já fazem parte do contrato congelado:
  * B3/B4: o intercepto da reta de P_f numa rampa NÃO é exatamente 2H·|R|/f0. Em regime
    de rampa, o ângulo precisa crescer para a potência crescer; com isso o VSG fica
    atrasado em frequência de b/(K_s·ω0) e o droop desconta D_p·b/(K_s·ω0):
        a = 2H·|R|/f0 − D_p·b/(K_s·ω0),   b = D_p·|R|/f0
    O "H medido" é, portanto, H_med = f0·(a + D_p·b/(K_s·ω0)) / (2|R|). No caso do B3
    (SCR 2) a correção vale 7 % do intercepto — maior que a tolerância de 5 %.
  * B4 usa rampa de −0,25 Hz/s por 4 s (janela de 1 s é curta demais para H = 9 s,
    cujo modo eletromecânico tem 0,92 Hz).
  * B2 não verifica a frequência da oscilação: com washout (ζ ≈ 0,5) ela não é
    mensurável com a precisão exigida; verifica o regime e a comparação com o oráculo.
  * B7 compara a simulação com o oráculo (2 %) e a previsão do gerador com o oráculo (5 %).

=====================================================================================
CONTRATO DE INTERFACE DA v7 (vsg_2a_ordem_degrau_carga_v7.py)
=====================================================================================
D1. Compatibilidade: todo parâmetro e arquivo aceito pela v6 continua aceito, com o mesmo
    significado e os MESMOS resultados (R1–R3). Todo o contrato C1–C7 de testes_v6.py
    continua valendo, exceto o campo "versao", que passa a ser "v7".
D2. Novos parâmetros (arquivo e linha de comando, "_" <-> "-"), seção "evento":
      evento   passa a aceitar também "freq_degrau" e "freq_rampa"
      df_g     desvio final da frequência da rede [Hz], com sinal     (padrão 0)
      rocof_g  módulo da taxa da rampa [Hz/s], > 0                    (padrão 1.0)
    Eventos de frequência exigem modo="rede" e df_g ≠ 0 (senão: erro de configuração).
    rocof_g ≤ 0 é sempre erro de configuração.
D3. Perfil da rede (t_s = t_step):
      Δf_g(t) = 0                                            para t < t_s
      freq_degrau: Δf_g(t) = df_g                            para t ≥ t_s
      freq_rampa:  Δf_g(t) = sign(df_g)·min(|df_g|, rocof_g·(t − t_s))  para t ≥ t_s
    Tensão da rede no referencial fixo em ω0:
      V_g = −j·√2·v_g·V_ph·exp(j·δ_g),   dδ_g/dt = 2π·Δf_g(t),   δ_g(0) = 0
    O evento "fase" continua somando d_fase a δ_g em t_step. Os pontos de quebra
    (t_step e, na rampa, t_step + |df_g|/rocof_g) devem ser limites de segmento de
    integração. A análise modal continua no equilíbrio inicial (autovalores da v6).
D4. Saídas adicionais:
    <prefixo>_series.npz ganha  f_rede_hz = f0 + Δf_g(t)  e  delta_g_rad = δ_g(t)
      (no modo ilhado: f_rede_hz = f0 e delta_g_rad = 0).
    <prefixo>_resultados.json: "versao": "v7" e, nos eventos de frequência, o bloco
      "inercial": {"H_eff_s", "dP_inercial_pu", "inclinacao_droop_pu_s",
                   "dP_previsto_pu", "P_pico_pu", "I_pico_pu"}
      - H_eff_s = H + D_w·T_w/2
      - dP_inercial_pu, inclinacao_droop_pu_s: a e b do ajuste
            P_f − P_f0 = a + b·τ + e^{−στ}(c1·cos ωτ + c2·sin ωτ),  τ = t − t_step,
        feito durante a rampa (freq_rampa) ou em [t_step, t_end] (freq_degrau);
        valores com o sinal físico (positivos quando a frequência da rede cai)
      - dP_previsto_pu = 2H·rocof_g/f0 (freq_rampa) ou null (freq_degrau)
      - P_pico_pu = max(P_f) para t ≥ t_step;  I_pico_pu = max(I_rms_a)/I_n,
        I_n = S_n/(√3·V_LL)
      Nos demais eventos, "inercial" pode ser null ou ausente.
    Quando I_pico_pu > 1,2, o terminal mostra uma linha iniciada por "(!)" contendo
    a palavra "corrente".
D5. Gerador (gerar_caso_vsg.py v7): --evento aceita freq_degrau e freq_rampa;
    --df-rede [Hz] grava df_g; a taxa da rampa (rocof_g) é --rocof-rede (ou --rocof-max
    se ausente). Nos eventos de frequência imprime, em linhas próprias:
        "H_eff = <x> s"            (H + D_w·T_w/2)
        "P_pico previsto = <y> pu" (P1 + (|R|/f0)·[2H + D_w·T_w·(1 − e^{−T_r/T_w})]
                                    + D_p·|Δf|/f0,  T_r = |Δf|/|R|; no degrau, T_r → 0
                                    e o termo inercial é omitido)
    e emite um aviso "(!)" contendo "H_eff" quando y > --pmax-pcs.
    O YAML gerado roda no simulador sem edição.

=====================================================================================
USO
=====================================================================================
  python testes_v7.py                         # todos os testes
  python testes_v7.py --so R1,B1,B3           # subconjunto (desenvolvimento incremental)
  python testes_v7.py --autoteste             # valida o oráculo e os critérios (sem v7)
  python testes_v7.py --script v7.py --gerador gerar.py --v6 v6.py --v5 v5.py
Tempo típico da suíte completa: 10 a 20 minutos.
Saídas: tabela no terminal, relatorio_testes_v7.md e relatorio_testes_v7.json.
Código de saída 0 somente se TODOS os testes selecionados forem aprovados e o arquivo
estiver ÍNTEGRO.
"""
import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit, fsolve

# =======================================================================================
# TOLERÂNCIAS (congeladas)
# =======================================================================================
TOL = {
    "R1_rel": 1e-3,            # 0,1 % em V, I e RoCoF
    "R1_f_abs": 1e-3,          # Hz
    "R2_pu": 1e-6,             # pu   (v7 × v6, séries)
    "R2_hz": 1e-6,             # Hz
    "R2_rad": 1e-6,            # rad
    "R3_rel": 1e-5,            # autovalores v7 × v6 (relativo a max(1, |λ|))
    "R3_oraculo": 2e-3,        # f [Hz] e ζ contra ORACULO_V6
    "B1_hz": 1e-9,             # f_rede × perfil analítico
    "B1_rad": 1e-6,            # δ_g × integral analítica
    "B2_f_abs": 1e-4,          # Hz
    "B2_P_rel": 5e-3,          # 0,5 %
    "B3_rel": 0.05,            # a, b e H medido: 5 %
    "B3_json_rel": 0.05,       # a e b relatados no JSON × ajuste da suíte
    "B5_rel": 0.05,            # ΔP_f em τ = 1, 2, 3 s × fórmula quase estática
    "B6_frac": 0.02,           # |P_f impl − P_f oráculo| ≤ 2 % de max|ΔP_f oráculo|
    "B7_orac": 0.02,           # pico simulado × oráculo
    "B7_prev": 0.05,           # pico previsto pelo gerador × oráculo
    "B8_rel": 1e-3,            # 0,1 %
    "B10_razao": 1.5,          # tempo v7 (rampa) ≤ 1,5 × v6 (pref), mesmo caso
}

# Caso-base (mesmo de testes_v6.py)
BASE = {
    "sn": 100000.0, "vll": 380.0, "f0": 60.0,
    "rf": 0.002947423967, "xf": 0.1178969587, "bc": 0.05,
    "p1": 50000.0, "q1": 14967.025, "p2": 30000.0, "q2": 8980.215002,
    "H": 4.5, "mp": 0.05, "nq": 0.05, "fc": 10.0,
    "xr_rede": 10.0, "vg": 1.0,
}
DW_REF, TW_REF = 103.0, 1.0          # D_w do caso_rede.yaml da v6
PMAX_PCS = 1.1

YAML_V5 = """sistema:
  sn:       100000
  vll:      380
  f0:       60
filtro:
  rf:       0.002947423967
  xf:       0.1178969587
  bc:       0.05
carga:
  p1:       50000
  q1:       14967.025
  p2:       30000
  q2:       8980.215002
controle:
  pref:     0.5
  qref:     0.09967025004
  eref:     1.014920847
  H:        4.5
  mp:       0.05
  nq:       0.05
  fc:       10
simulacao:
  t_step:   1
  t_end:    4
  dt_out:   0.0001
saida:
  prefixo:  "vsg"
"""
REF_V5 = {"antes": {"Vll_rms_v": 380.00, "I_rms_a": 77.46, "f_hz": 60.0000},
          "final": {"Vll_rms_v": 373.54, "I_rms_a": 122.98, "f_hz": 59.1770},
          "rocof_max_hz_s": 1.608}

# Modo eletromecânico publicado na v6 (testes_v6.py): (scr, dw) -> (f_n [Hz], ζ)
ORACULO_V6 = {
    (20, 0): (2.455, -0.050), (5, 0): (1.803, 0.005), (2, 0): (1.302, 0.067),
    (20, 100): (2.387, 0.295), (5, 100): (1.713, 0.497), (2, 100): (1.159, 0.804),
}

# Valores do oráculo v7 calculados na elaboração deste arquivo (verificados em --autoteste)
ORACULO_V7 = {
    "B3_a": 0.06924, "B3_b": 0.16730,                  # SCR 2, R = −0,5 Hz/s, 1 s
    "B5_dP": (0.15486, 0.23141, 0.27904),              # τ = 1, 2, 3 s
    "B2_P_final": 0.06674, "B2_f_final": 59.79999,     # degrau −0,2 Hz, 8 s
    "B7_pico": 1.2943,                                 # P_f máx., rampa −0,5 Hz/s até −0,5 Hz
}

# =======================================================================================
# ORÁCULO — simulação não linear independente, referencial que gira com a REDE
#   x' = x·e^{−jδ_g}:   termos cruzados −j(ω0 + 2πΔf_g)·x',  V_g' constante,
#   ângulo relativo φ = δ_v − δ_g com dφ/dt = ω0·Δω − 2πΔf_g
# =======================================================================================
def perfil(t, ev):
    tipo, ts, df, rc = ev
    if tipo not in ("freq_degrau", "freq_rampa") or t < ts:
        return 0.0
    if tipo == "freq_degrau":
        return df
    return math.copysign(min(abs(df), rc * (t - ts)), df)


def integral_perfil(t, ev):
    """∫ Δf_g dt de 0 a t (analítico) [Hz·s]."""
    tipo, ts, df, rc = ev
    if tipo not in ("freq_degrau", "freq_rampa"):
        return np.zeros_like(t)
    tau = np.clip(t - ts, 0.0, None)
    if tipo == "freq_degrau":
        return df * tau
    tr = abs(df) / rc
    s = np.sign(df)
    return np.where(tau <= tr, s * rc * tau ** 2 / 2, s * rc * tr ** 2 / 2 + df * (tau - tr))


def _eletrica(c, scr):
    Sn, V, f0 = c["sn"], c["vll"], c["f0"]
    w0 = 2 * np.pi * f0
    Zb = V ** 2 / Sn
    d = dict(Sn=Sn, w0=w0, f0=f0, Vph=V / np.sqrt(3),
             Rf=c["rf"] * Zb, Lf=c["xf"] * Zb / w0, Cf=c["bc"] / (Zb * w0))
    S2 = c["p1"] ** 2 + c["q1"] ** 2
    d["R1"], d["L1"] = V ** 2 * c["p1"] / S2, V ** 2 * c["q1"] / S2 / w0
    xr = c["xr_rede"]
    Xg = (Zb / scr) / np.sqrt(1 + 1 / xr ** 2)
    d["Rg"], d["Lg"], d["vg"] = Xg / xr, Xg / w0, c["vg"]
    return d


class Oraculo:
    def __init__(self, c, scr, dw=0.0, tw=1.0):
        self.c, self.scr, self.dw, self.tw = c, scr, dw, tw
        d = self.d = _eletrica(c, scr)
        w0, Vph = d["w0"], d["Vph"]
        self.Vg = -1j * np.sqrt(2) * d["vg"] * Vph
        phi0, E0 = fsolve(lambda v: self._ig_nulo(v), [0.05, 1.0], xtol=1e-14)
        If, Vc, I1, Ig = self._fasores(phi0, E0)
        S = 1.5 * Vc * np.conj(If) / d["Sn"]
        self.Pref, self.Qref, self.Eref = S.real, S.imag, E0
        self.x0 = np.array([If.real, If.imag, Vc.real, Vc.imag, I1.real, I1.imag,
                            Ig.real, Ig.imag, phi0, 0.0, S.real, S.imag, 0.0])
        self.Ks = E0 * d["vg"] / (c["xf"] + (1 / scr) / np.sqrt(1 + 1 / c["xr_rede"] ** 2))
        self.w0 = w0

    def _fasores(self, phi, E):
        d = self.d
        w0 = d["w0"]
        e = -1j * np.sqrt(2) * E * d["Vph"] * np.exp(1j * phi)
        Zf = d["Rf"] + 1j * w0 * d["Lf"]
        Zg = d["Rg"] + 1j * w0 * d["Lg"]
        Y1 = 1 / (d["R1"] + 1j * w0 * d["L1"])
        Vc = (e / Zf + self.Vg / Zg) / (1 / Zf + 1j * w0 * d["Cf"] + Y1 + 1 / Zg)
        return (e - Vc) / Zf, Vc, Vc * Y1, (Vc - self.Vg) / Zg

    def _ig_nulo(self, v):
        Ig = self._fasores(v[0], v[1])[3]
        return [Ig.real / self.d["Vph"], Ig.imag / self.d["Vph"]]

    def F(self, t, x, ev):
        d, c = self.d, self.c
        w0, Vph = d["w0"], d["Vph"]
        dfg = perfil(t, ev)
        wg = w0 + 2 * np.pi * dfg
        If_, Vc_ = x[0] + 1j * x[1], x[2] + 1j * x[3]
        I1_, Ig_ = x[4] + 1j * x[5], x[6] + 1j * x[7]
        phi, dwv, Pf, Qf, z = x[8:13]
        E = self.Eref - c["nq"] * (Qf - self.Qref)
        e = -1j * np.sqrt(2) * E * Vph * np.exp(1j * phi)
        dIf = (e - Vc_ - d["Rf"] * If_) / d["Lf"] - 1j * wg * If_
        dI1 = (Vc_ - d["R1"] * I1_) / d["L1"] - 1j * wg * I1_
        dIg = (Vc_ - d["Rg"] * Ig_ - self.Vg) / d["Lg"] - 1j * wg * Ig_
        dVc = (If_ - I1_ - Ig_) / d["Cf"] - 1j * wg * Vc_
        S = 1.5 * Vc_ * np.conj(If_) / d["Sn"]
        ddw = (self.Pref - Pf - dwv / c["mp"] - self.dw * (dwv - z)) / (2 * c["H"])
        wc = 2 * np.pi * c["fc"]
        return np.array([dIf.real, dIf.imag, dVc.real, dVc.imag, dI1.real, dI1.imag,
                         dIg.real, dIg.imag, w0 * dwv - 2 * np.pi * dfg, ddw,
                         wc * (S.real - Pf), wc * (S.imag - Qf), (dwv - z) / self.tw])

    def autovalores(self):
        ev = ("nenhum", 0.0, 0.0, 1.0)
        n = len(self.x0)
        J = np.zeros((n, n))
        for i in range(n):
            h = 1e-6 * max(1.0, abs(self.x0[i]))
            e = np.zeros(n)
            e[i] = h
            J[:, i] = (self.F(0, self.x0 + e, ev) - self.F(0, self.x0 - e, ev)) / (2 * h)
        return np.linalg.eigvals(J), np.abs(self.F(0, self.x0, ev)).max()

    def simular(self, ev, t_end, dt_out=1e-4):
        """Retorna t, P_f [pu], f_VSG [Hz]."""
        tipo, ts, df, rc = ev
        br = [0.0]
        if tipo in ("freq_degrau", "freq_rampa"):
            br.append(ts)
        if tipo == "freq_rampa":
            br.append(ts + abs(df) / rc)
        br = [b for b in br if b < t_end] + [t_end]
        x0, T, X = self.x0.copy(), [], []
        for a, b in zip(br[:-1], br[1:]):
            te = np.linspace(a, b, int(round((b - a) / dt_out)) + 1)
            s = solve_ivp(self.F, (a, b), x0, args=(ev,), t_eval=te, method="RK45",
                          rtol=1e-10, atol=1e-10, max_step=2e-4)
            if not s.success:
                raise RuntimeError(f"oráculo: {s.message}")
            T.append(s.t[:-1])
            X.append(s.y[:, :-1])
            x0 = s.y[:, -1]
        T.append(np.array([t_end]))
        X.append(x0[:, None])
        t, X = np.concatenate(T), np.concatenate(X, axis=1)
        return t, X[10], self.c["f0"] * (1 + X[9])


def pico_previsto(P1, H, dw, tw, R, df, f0=60.0, Dp=20.0):
    """Fórmula do contrato D5 (quase estática)."""
    tr = abs(df) / abs(R)
    return P1 + (abs(R) / f0) * (2 * H + dw * tw * (1 - math.exp(-tr / tw))) + Dp * abs(df) / f0


# =======================================================================================
# Ferramentas de análise
# =======================================================================================
def autovalores(js):
    lam = np.array([complex(r, i) for r, i in js["modal"]["autovalores"]])
    return lam[np.abs(lam) >= 1e-6]


def dominante(lam):
    osc = [l for l in lam if l.imag > 0 and 0.1 < l.imag / (2 * np.pi) < 20]
    if not osc:
        raise AssertionError("nenhum modo oscilatório entre 0,1 e 20 Hz")
    l = min(osc, key=lambda x: -x.real / abs(x))
    return abs(l) / (2 * np.pi), -l.real / abs(l)


def ajustar_reta_oscilacao(tau, y, f_guess, z_guess):
    """y = a + b·τ + e^{−στ}(c1 cos ωτ + c2 sin ωτ). Retorna (a, b, f_osc [Hz])."""
    A = np.vstack([np.ones_like(tau), tau]).T
    a0, b0 = np.linalg.lstsq(A, y, rcond=None)[0]
    wg = 2 * np.pi * f_guess
    sg = max(1e-3, z_guess * wg)

    def mod(tt, a, b, c1, c2, s, w):
        return a + b * tt + np.exp(-s * tt) * (c1 * np.cos(w * tt) + c2 * np.sin(w * tt))

    p, _ = curve_fit(mod, tau, y, p0=[a0, b0, y[0] - a0, 0.0, sg, wg], maxfev=40000)
    return p[0], p[1], abs(p[5]) / (2 * np.pi)


def h_medido(a, b, Ks, f0, R, Dp=20.0):
    w0 = 2 * np.pi * f0
    return f0 * (a + Dp * b / (Ks * w0)) / (2 * abs(R))


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-15)


def checar(cond, msg):
    if not cond:
        raise AssertionError(msg)


# =======================================================================================
# Execução do simulador / gerador
# =======================================================================================
class Ambiente:
    def __init__(self, script, gerador, v6, v5, workdir):
        self.script, self.gerador, self.v6, self.v5, self.dir = script, gerador, v6, v5, workdir
        self.base_json = os.path.join(workdir, "base.json")
        with open(self.base_json, "w", encoding="utf-8") as fh:
            json.dump(BASE, fh)
        self.cont = 0
        self.stdout = ""

    def _pref(self, nome):
        self.cont += 1
        return os.path.join(self.dir, f"{self.cont:03d}_{nome}")

    def sim(self, nome, args, config=None, graficos=False, script=None, esperar_erro=False,
            ler=True):
        pref = self._pref(nome)
        cmd = [sys.executable, script or self.script, "--config", config or self.base_json,
               *[str(a) for a in args], "--prefixo", pref]
        if not graficos:
            cmd.append("--sem-graficos")
        t0 = time.perf_counter()
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1800,
                            encoding="utf-8", errors="replace")
        dt = time.perf_counter() - t0
        self.stdout = pr.stdout or ""
        if esperar_erro:
            return pr
        if pr.returncode != 0:
            raise AssertionError(f"simulador falhou (código {pr.returncode}):\n"
                                 f"{(pr.stderr or pr.stdout)[-1500:]}")
        if not ler:
            return None, None, dt
        with open(pref + "_resultados.json", encoding="utf-8") as fh:
            js = json.load(fh)
        npz = np.load(pref + "_series.npz")
        ser = {k: np.asarray(npz[k]) for k in npz.files}
        return js, ser, dt

    def gerar(self, nome, args):
        saida = self._pref(nome) + ".yaml"
        cmd = [sys.executable, self.gerador, *[str(a) for a in args], "-o", saida,
               "--prefixo", self._pref(nome + "_sim")]
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=900,
                            encoding="utf-8", errors="replace")
        return pr, saida


def rede(scr, dw=0.0, tw=1.0):
    return ["--modo", "rede", "--scr", scr, "--dw", dw, "--tw", tw]


def curto(extra=()):
    return ["--evento", "nenhum", "--t-step", 0.2, "--t-end", 0.3, "--dt-out", 1e-4, *extra]


def rampa(df, rc, ts, te, dt=1e-4):
    return ["--evento", "freq_rampa", "--df-g", df, "--rocof-g", rc,
            "--t-step", ts, "--t-end", te, "--dt-out", dt]


def cfg(H=None):
    c = dict(BASE)
    if H is not None:
        c["H"] = H
    return c


def orac_cache(cache, chave, fn):
    k = "orac_" + chave
    if k not in cache:
        cache[k] = fn()
    return cache[k]


# ---------------------------------------------------------------------------------------
# Casos padronizados
# ---------------------------------------------------------------------------------------
B3_EV = ("freq_rampa", 0.5, -0.5, 0.5)                      # SCR 2, D_w = 0, 1 s de rampa
B3_ARGS = [*rede(2), *rampa(-0.5, 0.5, 0.5, 1.6)]
B5_EV = ("freq_rampa", 0.5, -0.3, 0.1)                      # SCR 5, D_w = 103, 3 s de rampa
B5_ARGS = [*rede(5, DW_REF, TW_REF), *rampa(-0.3, 0.1, 0.5, 4.0)]
B2_EV = ("freq_degrau", 0.5, -0.2, 1.0)
B2_ARGS = [*rede(5, DW_REF, TW_REF), "--evento", "freq_degrau", "--df-g", -0.2,
           "--t-step", 0.5, "--t-end", 8.5, "--dt-out", 2e-4]


def medir_rampa(t, pf, ts, dur, pf0, f_guess, z_guess):
    m = (t >= ts) & (t <= ts + dur + 1e-12)
    return ajustar_reta_oscilacao(t[m] - ts, pf[m] - pf0, f_guess, z_guess)


def pf0_de(t, pf, ts):
    return float(np.mean(pf[(t > ts - 0.1) & (t < ts)]))


# =======================================================================================
# TESTES — Regressão
# =======================================================================================
def R1(amb, cache):
    """Regressão do modo ilhado: meu_caso.yaml da v5 reproduz os valores publicados."""
    y = os.path.join(amb.dir, "caso_v5.yaml")
    with open(y, "w", encoding="utf-8") as fh:
        fh.write(YAML_V5)
    js, s, _ = amb.sim("R1", [], config=y)
    checar(js.get("versao") == "v7", f"campo 'versao' = {js.get('versao')!r} (esperado 'v7')")
    checar(js.get("modo") == "ilhado", "modo padrão deveria ser 'ilhado'")
    out = []
    for bloco, chave in (("antes", "antes_evento"), ("final", "final")):
        for k, ref in REF_V5[bloco].items():
            v = js[chave][k]
            ok = abs(v - ref) <= TOL["R1_f_abs"] if k == "f_hz" else rel(v, ref) <= TOL["R1_rel"]
            out.append(f"{bloco}.{k}={v:.4f}")
            checar(ok, f"{bloco}.{k} = {v:.5f} difere da v5 ({ref})")
    r = js["rocof_max_hz_s"]
    checar(rel(r, REF_V5["rocof_max_hz_s"]) <= TOL["R1_rel"], f"RoCoF {r:.4f} ≠ 1,608")
    for k in ("f_rede_hz", "delta_g_rad"):
        checar(k in s, f"série '{k}' ausente do NPZ (contrato D4)")
    checar(np.all(s["f_rede_hz"] == BASE["f0"]) and np.all(s["delta_g_rad"] == 0),
           "no modo ilhado f_rede_hz deve ser f0 e delta_g_rad deve ser 0")
    return "; ".join(out) + f"; rocof={r:.4f}"


R2_CASOS = [
    ("pref_scr2", [*rede(2), "--evento", "pref", "--t-step", 0.5, "--t-end", 3.0, "--dt-out", 1e-4]),
    ("pref_scr5_dw", [*rede(5, DW_REF), "--evento", "pref", "--t-step", 0.5, "--t-end", 3.0,
                      "--dt-out", 1e-4]),
    ("carga_scr5_dw", [*rede(5, 100), "--evento", "carga", "--t-step", 0.5, "--t-end", 3.0,
                       "--dt-out", 2e-4]),
    ("fase_scr20_dw", [*rede(20, 100), "--evento", "fase", "--d-fase", 5, "--t-step", 0.5,
                       "--t-end", 2.0, "--dt-out", 5e-5]),
    ("carga_ilha_dw", ["--evento", "carga", "--dw", 100, "--tw", 0.1, "--t-step", 0.5,
                       "--t-end", 3.0, "--dt-out", 1e-4]),
]


def R2(amb, cache):
    """Eventos da v6 (pref, carga, fase; ilha e rede; com e sem washout): séries v7 = v6."""
    checar(os.path.isfile(amb.v6), f"script da v6 não encontrado: {amb.v6}")
    out = []
    for nome, args in R2_CASOS:
        j6, s6, _ = amb.sim(f"R2_{nome}_v6", args, script=amb.v6)
        j7, s7, _ = amb.sim(f"R2_{nome}_v7", args)
        checar(len(s6["t"]) == len(s7["t"]) and np.allclose(s6["t"], s7["t"], atol=1e-12),
               f"{nome}: vetores de tempo diferentes")
        pior = 0.0
        for k, tol in (("Pf_pu", "R2_pu"), ("P_pu", "R2_pu"), ("P_rede_pu", "R2_pu"),
                       ("Q_rede_pu", "R2_pu"), ("f_hz", "R2_hz"), ("delta_v_rad", "R2_rad")):
            d = float(np.nanmax(np.abs(s6[k] - s7[k])))
            pior = max(pior, d / TOL[tol])
            checar(d <= TOL[tol], f"{nome}: {k} difere da v6 em {d:.2e} (tol {TOL[tol]:.0e})")
        for bloco in ("antes_evento", "final"):
            for k, v in j6[bloco].items():
                checar(abs(v - j7[bloco][k]) <= max(1e-6, 1e-7 * abs(v)),
                       f"{nome}: {bloco}.{k} v7 {j7[bloco][k]} ≠ v6 {v}")
        out.append(f"{nome} ok ({pior:.2f}×tol)")
    return "; ".join(out)


def R3(amb, cache):
    """Análise modal: autovalores v7 = v6 e modo eletromecânico = valores publicados."""
    checar(os.path.isfile(amb.v6), f"script da v6 não encontrado: {amb.v6}")
    out = []
    for (scr, dw), (fe, ze) in ORACULO_V6.items():
        j6, _, _ = amb.sim(f"R3_{scr}_{dw}_v6", [*rede(scr, dw), *curto()], script=amb.v6)
        j7, _, _ = amb.sim(f"R3_{scr}_{dw}_v7", [*rede(scr, dw), *curto()])
        l6, l7 = autovalores(j6), autovalores(j7)
        for l in l6:
            dmin = np.min(np.abs(l7 - l)) / max(1.0, abs(l))
            checar(dmin <= TOL["R3_rel"], f"SCR {scr}, D_w {dw}: autovalor {l:.4f} da v6 ausente na v7")
        checar(len(l7) == len(l6), f"SCR {scr}, D_w {dw}: {len(l7)} autovalores na v7 × {len(l6)} na v6")
        f, z = dominante(l7)
        checar(abs(f - fe) <= TOL["R3_oraculo"] and abs(z - ze) <= TOL["R3_oraculo"],
               f"SCR {scr}, D_w {dw}: f={f:.4f}, ζ={z:.4f} (publicado {fe}, {ze})")
        out.append(f"SCR{scr}/Dw{dw}: {f:.3f} Hz, ζ {z:.3f}")
    return "; ".join(out)


# =======================================================================================
# TESTES — Frequência da rede
# =======================================================================================
def B1(amb, cache):
    """Perfil imposto: rampa de −0,5 Hz/s até −0,5 Hz (SCR 5); f_rede e δ_g analíticos."""
    ev = ("freq_rampa", 0.5, -0.5, 0.5)
    js, s, _ = amb.sim("B1", [*rede(5, DW_REF), *rampa(-0.5, 0.5, 0.5, 2.5)])
    t = s["t"]
    ref = np.array([BASE["f0"] + perfil(tt, ev) for tt in t])
    longe = (np.abs(t - 0.5) > 1e-9) & (np.abs(t - 1.5) > 1e-9)
    ef = float(np.max(np.abs(s["f_rede_hz"][longe] - ref[longe])))
    checar(ef <= TOL["B1_hz"], f"f_rede difere do perfil em {ef:.2e} Hz")
    dg_ref = 2 * np.pi * integral_perfil(t, ev)
    ed = float(np.max(np.abs(s["delta_g_rad"] - dg_ref)))
    checar(ed <= TOL["B1_rad"], f"δ_g difere da integral analítica em {ed:.2e} rad")
    salto = float(np.max(np.abs(np.diff(s["delta_g_rad"]))))
    checar(salto < 2 * np.pi * 0.5 * np.median(np.diff(t)) * 1.01 + 1e-12,
           f"δ_g descontínuo (maior passo {salto:.2e} rad)")
    k = np.searchsorted(t, 1.5 + 1e-9)
    checar(abs(s["f_rede_hz"][k] - 59.5) < 1e-9, "a rampa não termina em t_step + |df_g|/rocof_g = 1,5 s")
    checar(js.get("versao") == "v7", "campo 'versao' diferente de 'v7'")
    return f"erro f_rede {ef:.1e} Hz; erro δ_g {ed:.1e} rad; fim da rampa em 1,5 s"


def B2(amb, cache):
    """Degrau de frequência −0,2 Hz (SCR 5, D_w = 103): regime final e aderência ao oráculo."""
    js, s, _ = amb.sim("B2", B2_ARGS)
    fin = js["final"]
    alvo_f = BASE["f0"] - 0.2
    alvo_dP = (1 / BASE["mp"]) * 0.2 / BASE["f0"]
    dP = fin["Pf_pu"] - js["controle_efetivo"]["pref"]
    checar(abs(fin["f_hz"] - alvo_f) <= TOL["B2_f_abs"], f"f final {fin['f_hz']:.6f} Hz (alvo {alvo_f})")
    checar(rel(dP, alvo_dP) <= TOL["B2_P_rel"], f"ΔP_f final {dP:.5f} pu (alvo {alvo_dP:.5f})")
    o = orac_cache(cache, "B2", lambda: Oraculo(cfg(), 5, DW_REF, TW_REF).simular(B2_EV, 8.5, 2e-4))
    to, pfo, _ = o
    esc = np.max(np.abs(pfo - pfo[0]))
    m = s["t"] >= 0.6
    err = np.max(np.abs(s["Pf_pu"][m] - np.interp(s["t"][m], to, pfo)))
    checar(err <= TOL["B6_frac"] * esc, f"P_f difere do oráculo em {err:.4f} pu (> 2 % de {esc:.3f})")
    ine = js.get("inercial") or {}
    checar("H_eff_s" in ine and ine.get("dP_previsto_pu", 0) is None,
           "bloco 'inercial' incompleto (H_eff_s; dP_previsto_pu = null no degrau)")
    return f"f={fin['f_hz']:.5f} Hz, ΔP_f={dP:.5f} (alvo {alvo_dP:.5f}); desvio do oráculo {err:.1e} pu"


# =======================================================================================
# TESTES — Resposta inercial
# =======================================================================================
def _b3_sim(amb, cache, extra=(), nome="B3"):
    js, s, _ = amb.sim(nome, [*B3_ARGS, *extra])
    fi, zi = dominante(autovalores(js))
    ts = 0.5
    pf0 = pf0_de(s["t"], s["Pf_pu"], ts)
    a, b, fo = medir_rampa(s["t"], s["Pf_pu"], ts, 1.0, pf0, fi, zi)
    return js, s, a, b, fo


def B3(amb, cache):
    """ΔP inercial puro: rampa −0,5 Hz/s por 1 s (SCR 2, D_w = 0)."""
    js, s, a, b, fo = _b3_sim(amb, cache)
    cache["B3"] = (js, s, a, b)
    H, R, f0, Dp = BASE["H"], 0.5, BASE["f0"], 1 / BASE["mp"]
    Ks = Oraculo(cfg(), 2).Ks
    b_teo = Dp * R / f0
    a_teo = 2 * H * R / f0 - Dp * b_teo / (Ks * 2 * np.pi * f0)
    Hm = h_medido(a, b, Ks, f0, R)
    checar(rel(b, b_teo) <= TOL["B3_rel"], f"b = {b:.5f} pu/s (teoria {b_teo:.5f})")
    checar(rel(a, a_teo) <= TOL["B3_rel"], f"a = {a:.5f} pu (teoria corrigida {a_teo:.5f})")
    checar(rel(Hm, H) <= TOL["B3_rel"], f"H medido = {Hm:.3f} s (H = {H})")
    ine = js.get("inercial") or {}
    for k in ("H_eff_s", "dP_inercial_pu", "inclinacao_droop_pu_s", "dP_previsto_pu",
              "P_pico_pu", "I_pico_pu"):
        checar(k in ine, f"bloco 'inercial' sem o campo '{k}' (contrato D4)")
    checar(abs(ine["H_eff_s"] - H) < 1e-9, f"H_eff_s = {ine['H_eff_s']} (esperado {H} com D_w = 0)")
    checar(rel(ine["dP_previsto_pu"], 2 * H * R / f0) < 1e-9, "dP_previsto_pu ≠ 2H·rocof_g/f0")
    checar(rel(ine["dP_inercial_pu"], a) <= TOL["B3_json_rel"],
           f"dP_inercial_pu relatado {ine['dP_inercial_pu']:.5f} × ajuste da suíte {a:.5f}")
    checar(rel(ine["inclinacao_droop_pu_s"], b) <= TOL["B3_json_rel"],
           f"inclinacao_droop_pu_s relatada {ine['inclinacao_droop_pu_s']:.5f} × ajuste {b:.5f}")
    ts = 0.5
    pk = float(np.max(s["Pf_pu"][s["t"] >= ts]))
    checar(abs(ine["P_pico_pu"] - pk) <= 1e-4, f"P_pico_pu {ine['P_pico_pu']:.5f} ≠ max(P_f) {pk:.5f}")
    return (f"a={a:.5f} (teo {a_teo:.5f}; 2H|R|/f0 = {2 * H * R / f0:.5f}), b={b:.5f} (teo {b_teo:.5f}), "
            f"H medido {Hm:.3f} s")


def B4(amb, cache):
    """Proporcionalidade com H: rampa −0,25 Hz/s por 4 s, SCR 2, H = 2; 4,5; 9 s."""
    Ks = Oraculo(cfg(), 2).Ks
    out = []
    for H in (2.0, 4.5, 9.0):
        js, s, _ = amb.sim(f"B4_H{H}", [*rede(2), "--H", H, *rampa(-1.0, 0.25, 0.5, 4.6)])
        fi, zi = dominante(autovalores(js))
        pf0 = pf0_de(s["t"], s["Pf_pu"], 0.5)
        a, b, _ = medir_rampa(s["t"], s["Pf_pu"], 0.5, 4.0, pf0, fi, zi)
        Hm = h_medido(a, b, Ks, BASE["f0"], 0.25)
        out.append(f"H={H}: medido {Hm:.3f} s")
        checar(rel(Hm, H) <= TOL["B3_rel"], f"H = {H}: H medido {Hm:.3f} s")
    return "; ".join(out)


def B5(amb, cache):
    """Washout como inércia: rampa −0,1 Hz/s por 3 s (SCR 5, D_w = 103, T_w = 1 s)."""
    js, s, _ = amb.sim("B5", B5_ARGS)
    cache["B5"] = (js, s)
    ts, R, f0, H = 0.5, 0.1, BASE["f0"], BASE["H"]
    dw, tw = js["controle_efetivo"]["dw"], js["controle_efetivo"]["tw"]
    pf0 = pf0_de(s["t"], s["Pf_pu"], ts)
    out = []
    for tau in (1.0, 2.0, 3.0):
        med = float(np.interp(ts + tau, s["t"], s["Pf_pu"])) - pf0
        teo = (R / f0) * (2 * H + dw * tw * (1 - math.exp(-tau / tw))) + (1 / BASE["mp"]) * R * tau / f0
        out.append(f"τ={tau:g}: {med:.4f}/{teo:.4f}")
        checar(rel(med, teo) <= TOL["B5_rel"], f"ΔP_f(τ={tau}) = {med:.4f} pu (fórmula {teo:.4f})")
    ine = js.get("inercial") or {}
    checar(abs(ine.get("H_eff_s", -1) - (H + dw * tw / 2)) < 1e-6,
           f"H_eff_s = {ine.get('H_eff_s')} (esperado {H + dw * tw / 2:.2f})")
    return "; ".join(out) + f"  [simulação/fórmula]; H_eff = {ine['H_eff_s']:.1f} s"


def B6(amb, cache):
    """Aderência ao oráculo independente (referencial da rede): séries P_f de B3 e B5."""
    if "B3" not in cache:
        B3(amb, cache)
    if "B5" not in cache:
        B5(amb, cache)
    casos = [("B3", cache["B3"][1], lambda: Oraculo(cfg(), 2).simular(B3_EV, 1.6)),
             ("B5", cache["B5"][1], lambda: Oraculo(cfg(), 5, DW_REF, TW_REF).simular(B5_EV, 4.0))]
    out = []
    for nome, s, fn in casos:
        to, pfo, fvo = orac_cache(cache, nome, fn)
        esc = np.max(np.abs(pfo - pfo[0]))
        m = s["t"] >= 0.6
        err = float(np.max(np.abs(s["Pf_pu"][m] - np.interp(s["t"][m], to, pfo))))
        errf = float(np.max(np.abs(s["f_hz"][m] - np.interp(s["t"][m], to, fvo))))
        out.append(f"{nome}: ΔP_f máx {err:.1e} pu ({100 * err / esc:.3f} %), Δf máx {errf:.1e} Hz")
        checar(err <= TOL["B6_frac"] * esc, f"{nome}: P_f difere do oráculo em {err:.4f} pu")
    return "; ".join(out)


# =======================================================================================
# TESTES — Critério de H_max e engenharia
# =======================================================================================
_RX_HEFF = re.compile(r"H_eff\s*=\s*([-+0-9.eE]+)\s*s")
_RX_PICO = re.compile(r"P_pico previsto\s*=\s*([-+0-9.eE]+)\s*pu")


def B7(amb, cache):
    """Gerador: H_eff, pico previsto × oráculo × simulação; aviso '(!)' com H_eff."""
    args = ["--modo", "rede", "--scr", 5, "--scr-faixa", "2,20", "--evento", "freq_rampa",
            "--df-rede", -0.5, "--rocof-rede", 0.5, "--t-step", 0.5, "--t-end", 3.0,
            "--dt-out", 2e-4]
    pr, y = amb.gerar("B7", args)
    checar(pr.returncode == 0, f"gerador falhou:\n{(pr.stderr or pr.stdout)[-1200:]}")
    mh, mp_ = _RX_HEFF.search(pr.stdout), _RX_PICO.search(pr.stdout)
    checar(mh and mp_, "o gerador deveria imprimir 'H_eff = x s' e 'P_pico previsto = y pu'")
    heff_g, pico_g = float(mh.group(1)), float(mp_.group(1))
    avisos = [l for l in pr.stdout.splitlines() if l.strip().startswith("(!)")]
    checar(any("H_eff" in l for l in avisos), f"pico previsto {pico_g:.3f} > {PMAX_PCS}: faltou o aviso '(!)' com H_eff")
    js, s, _ = amb.sim("B7_sim", [], config=y)
    dw, tw = js["controle_efetivo"]["dw"], js["controle_efetivo"]["tw"]
    H = js["parametros"]["H"]
    checar(abs(heff_g - (H + dw * tw / 2)) <= 0.05 + 1e-3 * heff_g,
           f"H_eff impresso {heff_g} ≠ H + D_w·T_w/2 = {H + dw * tw / 2:.2f}")
    c = dict(BASE, H=H, rf=js["parametros"]["rf"], xf=js["parametros"]["xf"], bc=js["parametros"]["bc"])
    to, pfo, _ = Oraculo(c, 5, dw, tw).simular(("freq_rampa", 0.5, -0.5, 0.5), 3.0, 2e-4)
    pico_o = float(np.max(pfo))
    pico_s = float(np.max(s["Pf_pu"][s["t"] >= 0.5]))
    checar(rel(pico_s, pico_o) <= TOL["B7_orac"], f"pico simulado {pico_s:.4f} × oráculo {pico_o:.4f}")
    checar(rel(pico_g, pico_o) <= TOL["B7_prev"], f"pico previsto {pico_g:.4f} × oráculo {pico_o:.4f}")
    ine = js.get("inercial") or {}
    In = BASE["sn"] / (math.sqrt(3) * BASE["vll"])
    ipk = float(np.nanmax(s["I_rms_a"][s["t"] >= 0.5])) / In
    checar(abs(ine.get("I_pico_pu", -1) - ipk) <= 0.01 * ipk, f"I_pico_pu {ine.get('I_pico_pu')} ≠ {ipk:.3f}")
    if ipk > 1.2:
        checar(any(l.strip().startswith("(!)") and "corrente" in l.lower() for l in amb.stdout.splitlines()),
               f"I_pico = {ipk:.2f} pu > 1,2 pu: o simulador deveria emitir aviso '(!)' de corrente")
    # controle negativo: rampa branda não deve gerar o aviso de H_eff
    pr2, _ = amb.gerar("B7_brando", ["--modo", "rede", "--scr", 5, "--scr-faixa", "2,20",
                                     "--evento", "freq_rampa", "--df-rede", -0.1, "--rocof-rede", 0.1])
    checar(pr2.returncode == 0, "gerador falhou no caso brando")
    checar(not any(l.strip().startswith("(!)") and "H_eff" in l for l in pr2.stdout.splitlines()),
           "aviso de H_eff emitido num caso cujo pico previsto é < p_max,pcs")
    return (f"H_eff {heff_g:.1f} s; pico previsto {pico_g:.4f} / oráculo {pico_o:.4f} / simulado "
            f"{pico_s:.4f} pu; I_pico {ipk:.2f} pu")


def B8(amb, cache):
    """Convergência: B3 e B5 com max_step e rtol dez vezes menores."""
    if "B3" not in cache:
        B3(amb, cache)
    if "B5" not in cache:
        B5(amb, cache)
    _, _, a1, b1 = cache["B3"]
    _, _, a2, b2, _ = _b3_sim(amb, cache, ["--max-step", 1e-4, "--rtol", 1e-9], "B8_B3")
    checar(rel(a2, a1) <= TOL["B8_rel"] and rel(b2, b1) <= TOL["B8_rel"],
           f"B3: a {a1:.6f}/{a2:.6f}, b {b1:.6f}/{b2:.6f}")
    _, s1 = cache["B5"]
    _, s2, _ = amb.sim("B8_B5", [*B5_ARGS, "--max-step", 1e-4, "--rtol", 1e-9])
    pf01, pf02 = pf0_de(s1["t"], s1["Pf_pu"], 0.5), pf0_de(s2["t"], s2["Pf_pu"], 0.5)
    out = [f"a {a1:.6f}/{a2:.6f}", f"b {b1:.6f}/{b2:.6f}"]
    for tau in (1.0, 2.0, 3.0):
        v1 = np.interp(0.5 + tau, s1["t"], s1["Pf_pu"]) - pf01
        v2 = np.interp(0.5 + tau, s2["t"], s2["Pf_pu"]) - pf02
        checar(rel(v2, v1) <= TOL["B8_rel"], f"B5 τ={tau}: {v1:.6f}/{v2:.6f}")
        out.append(f"ΔP(τ={tau:g}) {v1:.5f}/{v2:.5f}")
    return ", ".join(out)


def B9(amb, cache):
    """Robustez: entradas inválidas dos eventos de frequência geram erro limpo."""
    amb.sim("B9_controle", [*rede(5), *rampa(-0.1, 0.1, 0.2, 0.4)])
    casos = {
        "rocof_g=0": [*rede(5), *rampa(-0.1, 0, 0.2, 0.4)],
        "rocof_g<0": [*rede(5), *rampa(-0.1, -1, 0.2, 0.4)],
        "df_g=0 rampa": [*rede(5), *rampa(0, 0.1, 0.2, 0.4)],
        "df_g=0 degrau": [*rede(5), "--evento", "freq_degrau", "--df-g", 0],
        "rampa em ilha": ["--modo", "ilhado", *rampa(-0.1, 0.1, 0.2, 0.4)],
        "degrau em ilha": ["--modo", "ilhado", "--evento", "freq_degrau", "--df-g", -0.1],
        "evento inválido": [*rede(5), "--evento", "freq_xyz"],
    }
    for nome, args in casos.items():
        pr = amb.sim("B9", args, esperar_erro=True)
        txt = (pr.stdout or "") + (pr.stderr or "")
        checar(pr.returncode != 0, f"{nome}: deveria falhar (código 0)")
        checar("Traceback" not in txt, f"{nome}: traceback exposto")
        checar("Erro" in txt, f"{nome}: mensagem sem 'Erro'")
    return f"{len(casos)} casos inválidos rejeitados corretamente"


def B10(amb, cache):
    """Desempenho: rampa de 4 s na v7 ≤ 1,5 × degrau de P_ref na v6 (SCR 5, D_w = 103, com gráficos)."""
    checar(os.path.isfile(amb.v6), f"script da v6 não encontrado: {amb.v6}")
    comum = [*rede(5, DW_REF), "--t-step", 0.5, "--t-end", 4.0, "--dt-out", 5e-5]
    _, _, t6 = amb.sim("B10_v6", [*comum, "--evento", "pref"], script=amb.v6, graficos=True, ler=False)
    _, _, t7 = amb.sim("B10_v7", [*comum, "--evento", "freq_rampa", "--df-g", -0.3, "--rocof-g", 0.1],
                       graficos=True, ler=False)
    r = t7 / t6
    checar(r <= TOL["B10_razao"], f"v7 {t7:.1f} s × v6 {t6:.1f} s (razão {r:.2f})")
    return f"v6 {t6:.1f} s, v7 {t7:.1f} s, razão {r:.2f}"


TESTES = [("R1", R1), ("R2", R2), ("R3", R3), ("B1", B1), ("B2", B2), ("B3", B3), ("B4", B4),
          ("B5", B5), ("B6", B6), ("B7", B7), ("B8", B8), ("B9", B9), ("B10", B10)]


# =======================================================================================
# Autoteste: valida o oráculo e mostra que os critérios são atingíveis (não usa a v7)
# =======================================================================================
def autoteste():
    ok = True

    def item(cond, msg):
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'ok' if cond else 'FALHA'}] {msg}")

    # 1. Oráculo no referencial da rede reproduz os modos publicados da v6
    for (scr, dw), (fe, ze) in ORACULO_V6.items():
        lam, res = Oraculo(cfg(), scr, dw).autovalores()
        f, z = dominante(lam)
        item(abs(f - fe) < 2e-3 and abs(z - ze) < 2e-3 and res < 1e-6,
             f"modo SCR {scr:>2}, D_w {dw:>3}: f={f:.3f} ζ={z:+.3f} (v6: {fe}, {ze}); resíduo {res:.0e}")
    # 2. Equilíbrio: sem evento, nada se move
    t, pf, fv = Oraculo(cfg(), 5).simular(("nenhum", 0, 0, 1), 0.3)
    item(np.max(np.abs(pf - pf[0])) < 1e-8 and np.max(np.abs(fv - 60)) < 1e-9, "equilíbrio estacionário")
    # 3. Ajuste de reta + senoide amortecida em sinal sintético
    tau = np.arange(0, 1.0, 1e-4)
    y = 0.07 + 0.16 * tau + np.exp(-0.55 * tau) * (-0.05 * np.cos(8.2 * tau) + 0.02 * np.sin(8.2 * tau))
    a, b, fo = ajustar_reta_oscilacao(tau, y, 1.3, 0.07)
    item(abs(a - 0.07) < 1e-6 and abs(b - 0.16) < 1e-6, f"ajuste sintético: a={a:.6f}, b={b:.6f}, f={fo:.4f} Hz")
    # 4. Integral analítica do perfil × integração numérica
    ev = ("freq_rampa", 0.5, -0.5, 0.5)
    tt = np.linspace(0, 2.5, 25001)
    num = np.concatenate([[0], np.cumsum(np.diff(tt) * 0.5 * (np.array([perfil(x, ev) for x in tt[1:]]) +
                                                              np.array([perfil(x, ev) for x in tt[:-1]])))])
    item(np.max(np.abs(num - integral_perfil(tt, ev))) < 1e-9, "integral analítica do perfil")
    # 5. Os critérios de B2–B7 aplicados ao próprio oráculo (o alvo é atingível)
    o2 = Oraculo(cfg(), 2)
    t, pf, _ = o2.simular(B3_EV, 1.6)
    fi, zi = dominante(o2.autovalores()[0])
    a, b, _ = medir_rampa(t, pf, 0.5, 1.0, pf0_de(t, pf, 0.5), fi, zi)
    b_teo = 20 * 0.5 / 60
    a_teo = 0.075 - 20 * b_teo / (o2.Ks * 2 * np.pi * 60)
    Hm = h_medido(a, b, o2.Ks, 60, 0.5)
    item(rel(a, a_teo) < TOL["B3_rel"] and rel(b, b_teo) < TOL["B3_rel"] and rel(Hm, 4.5) < TOL["B3_rel"],
         f"B3: a={a:.5f} (teo {a_teo:.5f}), b={b:.5f}, H medido {Hm:.3f} s")
    item(abs(a - ORACULO_V7["B3_a"]) < 2e-4 and abs(b - ORACULO_V7["B3_b"]) < 2e-4, "B3 = valores congelados")
    for H in (2.0, 9.0):
        oh = Oraculo(cfg(H), 2)
        t, pf, _ = oh.simular(("freq_rampa", 0.5, -1.0, 0.25), 4.6)
        fi, zi = dominante(oh.autovalores()[0])
        a, b, _ = medir_rampa(t, pf, 0.5, 4.0, pf0_de(t, pf, 0.5), fi, zi)
        Hm = h_medido(a, b, oh.Ks, 60, 0.25)
        item(rel(Hm, H) < TOL["B3_rel"], f"B4: H = {H} -> H medido {Hm:.3f} s")
    o5 = Oraculo(cfg(), 5, DW_REF, TW_REF)
    t, pf, _ = o5.simular(B5_EV, 4.0)
    for k, tau in enumerate((1.0, 2.0, 3.0)):
        med = np.interp(0.5 + tau, t, pf) - pf[0]
        teo = (0.1 / 60) * (9 + DW_REF * (1 - math.exp(-tau))) + 20 * 0.1 * tau / 60
        item(rel(med, teo) < TOL["B5_rel"] and abs(med - ORACULO_V7["B5_dP"][k]) < 2e-4,
             f"B5: ΔP_f(τ={tau:g}) = {med:.5f} (fórmula {teo:.5f})")
    t, pf, fv = o5.simular(B2_EV, 8.5, 2e-4)
    m = t > t[-1] - 0.1
    dP, ff = pf[m].mean() - pf[0], fv[m].mean()
    item(rel(dP, 20 * 0.2 / 60) < TOL["B2_P_rel"] and abs(ff - 59.8) < TOL["B2_f_abs"],
         f"B2: ΔP_f final {dP:.5f} pu, f final {ff:.5f} Hz")
    t, pf, _ = o5.simular(("freq_rampa", 0.5, -0.5, 0.5), 3.0, 2e-4)
    pico = float(np.max(pf))
    prev = pico_previsto(pf[0], 4.5, DW_REF, TW_REF, 0.5, -0.5)
    item(rel(prev, pico) < TOL["B7_prev"] and abs(pico - ORACULO_V7["B7_pico"]) < 1e-3,
         f"B7: pico oráculo {pico:.4f} pu, previsto {prev:.4f} pu (> {PMAX_PCS}: aviso esperado)")
    prev_b = pico_previsto(0.5, 4.5, DW_REF, TW_REF, 0.1, -0.1)
    item(prev_b < PMAX_PCS, f"B7 controle negativo: previsto {prev_b:.3f} pu < {PMAX_PCS}")
    return ok


# =======================================================================================
def sha256(caminho):
    h = hashlib.sha256()
    with open(caminho, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def verificar_integridade():
    eu = os.path.abspath(__file__)
    atual = sha256(eu)
    arq = os.path.join(os.path.dirname(eu), "testes_v7.sha256")
    if not os.path.isfile(arq):
        return atual, "SEM REFERÊNCIA (testes_v7.sha256 ausente)"
    esperado = open(arq, encoding="utf-8").read().split()[0].strip().lower()
    return atual, "ÍNTEGRO" if esperado == atual else "ALTERADO — resultados sem validade"


def main():
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description="Suíte de aceitação congelada da v7.")
    aqui = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--script", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v7.py"))
    ap.add_argument("--gerador", default=os.path.join(aqui, "gerar_caso_vsg.py"))
    ap.add_argument("--v6", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v6.py"))
    ap.add_argument("--v5", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v5.py"))
    ap.add_argument("--so", default=None, help="Lista de testes, ex.: R1,B1,B3")
    ap.add_argument("--autoteste", action="store_true", help="Valida oráculo e critérios e encerra")
    ap.add_argument("--manter", action="store_true", help="Mantém a pasta de trabalho")
    a = ap.parse_args()

    hash_atual, estado = verificar_integridade()
    print(f"testes_v7.py  SHA-256 {hash_atual}  [{estado}]")
    if a.autoteste:
        print("Autoteste do oráculo e dos critérios:")
        ok = autoteste()
        print("AUTOTESTE", "APROVADO" if ok else "REPROVADO")
        sys.exit(0 if ok else 1)

    sel = [x.strip() for x in a.so.split(",")] if a.so else [n for n, _ in TESTES]
    desconhecidos = set(sel) - {n for n, _ in TESTES}
    if desconhecidos:
        sys.exit(f"Erro: testes desconhecidos: {', '.join(sorted(desconhecidos))}")
    if not os.path.isfile(a.script):
        sys.exit(f"Erro: implementação não encontrada: {a.script}")

    work = tempfile.mkdtemp(prefix="testes_v7_")
    amb = Ambiente(a.script, a.gerador, a.v6, a.v5, work)
    cache, resultados = {}, []
    t_ini = time.perf_counter()
    for nome, fn in TESTES:
        if nome not in sel:
            continue
        t0 = time.perf_counter()
        try:
            det, st = fn(amb, cache), "APROVADO"
        except AssertionError as e:
            det, st = str(e), "REPROVADO"
        except Exception as e:  # noqa: BLE001
            det, st = f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=2)}", "ERRO"
        dt = time.perf_counter() - t0
        resultados.append({"teste": nome, "status": st, "tempo_s": round(dt, 1),
                           "descricao": fn.__doc__.strip(), "detalhe": det})
        print(f"[{st:^9}] {nome:<4} ({dt:5.1f} s) {fn.__doc__.strip()}\n            {det}")

    n_ok = sum(r["status"] == "APROVADO" for r in resultados)
    total = time.perf_counter() - t_ini
    print(f"\nResultado: {n_ok}/{len(resultados)} aprovados em {total:.0f} s  |  integridade: {estado}")
    rel_json = {"sha256_testes": hash_atual, "integridade": estado, "script": a.script,
                "sha256_script": sha256(a.script), "aprovados": n_ok, "total": len(resultados),
                "resultados": resultados}
    with open("relatorio_testes_v7.json", "w", encoding="utf-8") as fh:
        json.dump(rel_json, fh, ensure_ascii=False, indent=2)
    with open("relatorio_testes_v7.md", "w", encoding="utf-8") as fh:
        fh.write(f"# Relatório de aceitação — v7\n\n- testes_v7.py SHA-256: `{hash_atual}` ({estado})\n"
                 f"- implementação: `{a.script}` SHA-256 `{rel_json['sha256_script']}`\n"
                 f"- resultado: **{n_ok}/{len(resultados)} aprovados** em {total:.0f} s\n\n"
                 "| Teste | Status | Tempo (s) | Descrição | Detalhe |\n|---|---|---|---|---|\n")
        for r in resultados:
            d = r["detalhe"].replace("\n", " ").replace("|", "/")[:400]
            fh.write(f"| {r['teste']} | {r['status']} | {r['tempo_s']} | {r['descricao']} | {d} |\n")
    if not a.manter:
        shutil.rmtree(work, ignore_errors=True)
    else:
        print(f"Pasta de trabalho: {work}")
    sys.exit(0 if n_ok == len(resultados) and estado == "ÍNTEGRO" else 1)


if __name__ == "__main__":
    main()

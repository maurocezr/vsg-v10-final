#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
testes_v8.py — suíte de aceitação CONGELADA da v8 do simulador VSG.

A v8 limita a envoltória instantânea da corrente I_f por impedância virtual
adaptativa e acrescenta anti-windup angular, falta trifásica equilibrada e
métricas de sincronismo. Este arquivo é o contrato executável da T0 e deve ser
congelado antes da implementação de vsg_2a_ordem_degrau_carga_v8.py.

REGRA DE CONGELAMENTO
=====================
O SHA-256 publicado em testes_v8.sha256 identifica esta revisão. Qualquer
mudança de fórmula, caso ou tolerância exige uma nova revisão da suíte e uma
justificativa escrita; referências nunca devem ser ajustadas silenciosamente
para acomodar a implementação.
REVISÃO testes_v8_1 (única alteração de contrato em relação a testes_v8.py)
==========================================================================
Origem: testes_v8.py SHA-256 4ff2ab51c92f791175316b058484abd8f9d3348a4781e5db8c3c73b58b8aefbc.
Teste afetado: C1. Fórmulas, oráculos, demais casos e TODAS as tolerâncias inalterados.
Defeito: C1 executava a v7 com --t-end 0.5 sem informar --t-step; o padrão é 1,0 s.
A v7 exige t_end > t_step para qualquer evento e rejeitava a configuração antes de
simular ("'t_end' (0.5) deve ser maior que 't_step' (1.0)"). O teste falhava na
referência v7, nunca chegando a avaliar a v8 — falha da suíte, não da implementação.
Correção: acrescentar "--t-step", 0.2 aos argumentos comuns de C1 (mesmo valor usado
por curto() em R4/R5). O evento é "nenhum", logo t_step não altera a dinâmica; o
critério de C1 (termos v8 nulos e séries iguais às da v7) permanece o mesmo.

CONTRATO MATEMÁTICO CONGELADO
=============================
A. Bases e corrente limitada
----------------------------
  Z_b       = V_LL²/S_n
  I_hat,n   = sqrt(2) S_n/(sqrt(3) V_LL)       [pico]
  i_f,pu    = |I_f|/I_hat,n
A corrente limitada é I_f (indutor do conversor), nunca I_g nem o RMS atrasado.

B. Lei pura da impedância virtual
---------------------------------
A implementação deve expor a função pública:

  calcular_impedancia_virtual(if_complex_a, vc_complex_v, e_bruta_complex_v,
      rf_ohm, i_base_pico_a, z_base_ohm, imax_pu, i_on_pu,
      rv_max_pu, xv_rv) -> dict

O dicionário contém ao menos: if_env_pu, rv_pu, xv_pu, ativo,
teto_rv_atingido e rv_requerido_pu. Se imax_pu == 0, todos os termos de
atuação (R_v, X_v, estado ativo e anti-windup) são zero; if_env_pu permanece
como grandeza diagnóstica. i_on_pu=None seleciona 0,98*imax_pu. Para imax_pu>0:

  xi = clip((i_f,pu-i_on)/(imax-i_on), 0, 1)
  s  = xi²(3-2xi)
  u  = I_f/|I_f|                       (avaliado somente acima de i_on)
  r_bar = max(0, Re{conj(u)(e_bruta-V_c)}/(imax*I_hat,n)
                    - R_f + 0,02 Z_b)
  R_v,requerido = s*r_bar
  R_v = min(R_v,requerido, rv_max_pu*Z_b)
  X_v = xv_rv*R_v

A FEM aplicada é e_aplicada=e_bruta-(R_v+jX_v)I_f. Na fronteira, se o
teto não saturar, g-R_v|I_f|²<=0, onde
  g=Re{conj(I_f)(e_bruta-V_c)}-R_f|I_f|².
A função suave é exatamente zero até i_on, contínua e monotônica para direção
e tensão radial fixas. O padrão é imax_pu=0, i_on_pu=auto, rv_max_pu=2,0,
xv_rv=0 e k_aw=20 s^-1.

C. Anti-windup
--------------
  delta_sat = wrap(angle(e_aplicada)-angle(e_bruta)) em [-pi, pi)
  delta_aw_rad_s = k_aw*delta_sat
O termo é somado à derivada do ângulo interno. É exatamente zero quando não há
limitação ou k_aw=0. O estado angular não pode ser recortado nem enrolado.

D. Falta e segmentação
----------------------
evento=falta_3f aplica vg_falta*vg entre t_step (inclusive) e t_clear
(exclusive) e restaura vg em t_clear. t_step e t_clear são limites de segmentos;
o vetor de tempo não tem duplicatas e os estados são contínuos na eliminação.

E. Sincronismo
--------------
  delta_rel = delta_v-delta_g; a decisão usa delta_rel_unwrapped.
  pole_slips = soma das mudanças absolutas de
      floor((delta_rel_unwrapped+pi)/(2*pi)).
Na janela final de 0,5 s:
  perdido       se pole_slips>=1, ou max|f_vsg-f_g|>0,20 Hz,
                ou a excursão pico-a-pico de delta_rel excede 0,50 rad;
  mantido       se não perdido, max|f_vsg-f_g|<=0,05 Hz e excursão<=0,10 rad;
  indeterminado nos demais casos.

F. Saídas
---------
O NPZ acrescenta if_env_pu, rv_pu, xv_pu, limitador_ativo,
delta_aw_rad_s, delta_rel_rad, delta_rel_unwrapped_rad e vg_aplicada_pu.
O JSON acrescenta o bloco limitador definido em tasks_v8.md. Quando o limitador
atua em evento de frequência, inercial.limitado=true e o relatório não pode
tratar 2H*RoCoF/f0 como validado sem ressalva.

USO
===
  python testes_v8.py --autoteste
  python testes_v8.py
  python testes_v8.py --so R1,R4,C1,C2
  python testes_v8.py --script v8.py --v7 v7.py --v6 v6.py --v5 v5.py

A suíte completa executa R1-R5 e C1-C14, gera relatorio_testes_v8.json/.md e
só retorna código zero quando todos os testes selecionados passam e o hash está
íntegro. O autoteste não importa nem executa a implementação v8.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
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
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ======================================================================================
# Constantes e tolerâncias congeladas
# ======================================================================================
AUTO_I_ON_FRAC = 0.98
MARGEM_RV_PU = 0.02
RV_MAX_PADRAO_PU = 2.0
K_AW_PADRAO = 20.0
TOL = {
    "R1_rel": 1e-3,
    "R1_f_abs": 1e-3,
    "serie_pu": 1e-6,
    "serie_hz": 1e-6,
    "serie_rad": 1e-6,
    "modal_rel": 1e-5,
    "barreira_abs": 1e-10,
    "formula_rel": 2e-12,
    "teto_corrente_frac": 0.005,
    "saida_rel": 5e-4,
    "continuidade_pu": 2e-6,
    "convergencia_frac": 0.005,
    "sync_f_mantido_hz": 0.05,
    "sync_f_perdido_hz": 0.20,
    "sync_delta_mantido_rad": 0.10,
    "sync_delta_perdido_rad": 0.50,
    "sync_janela_s": 0.50,
    "fronteira_falta_s": 0.15,
    "desempenho_razao": 2.0,
}

BASE = {
    "sn": 100000.0, "vll": 380.0, "f0": 60.0,
    "rf": 0.002947423967, "xf": 0.1178969587, "bc": 0.05,
    "p1": 50000.0, "q1": 14967.025, "p2": 30000.0, "q2": 8980.215002,
    "H": 4.5, "mp": 0.05, "nq": 0.05, "fc": 10.0,
    "xr_rede": 10.0, "vg": 1.0,
}
DW_REF, TW_REF = 103.0, 1.0
HASH_TESTES_V7 = "d019398d74c4072f5233f4a96c94cb4926cd6648e8ffe07290626866c1797b4d"

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
REF_V5 = {
    "antes": {"Vll_rms_v": 380.00, "I_rms_a": 77.46, "f_hz": 60.0000},
    "final": {"Vll_rms_v": 373.54, "I_rms_a": 122.98, "f_hz": 59.1770},
    "rocof_max_hz_s": 1.608,
}
SERIES_COMUNS = {
    "Pf_pu": "serie_pu", "P_pu": "serie_pu", "P_rede_pu": "serie_pu",
    "Q_rede_pu": "serie_pu", "f_hz": "serie_hz", "delta_v_rad": "serie_rad",
}
SERIES_V8 = (
    "if_env_pu", "rv_pu", "xv_pu", "limitador_ativo", "delta_aw_rad_s",
    "delta_rel_rad", "delta_rel_unwrapped_rad", "vg_aplicada_pu",
)
CAMPOS_LIMITADOR = (
    "habilitado", "imax_pu", "i_on_pu", "I_pico_env_pu", "I_pico_rms_pu",
    "tempo_ativo_s", "rv_max_usado_pu", "xv_max_usado_pu", "teto_rv_atingido",
    "violacao_max_pu", "tempo_em_violacao_s", "delta_aw_max_rad_s", "pole_slips",
    "sincronismo",
)

# Casos congelados. O caso longo de falta fica afastado da fronteira do oráculo.
FALTA_CURTA = ["--modo", "rede", "--scr", 5, "--evento", "falta_3f",
               "--vg-falta", 0, "--t-step", 0.5, "--t-clear", 0.60,
               "--t-end", 3.0, "--dt-out", 1e-4]
RAMPA_LIMITADA = ["--modo", "rede", "--scr", 5, "--dw", DW_REF, "--tw", TW_REF,
                   "--evento", "freq_rampa", "--df-g", -0.5, "--rocof-g", 0.5,
                   "--t-step", 0.5, "--t-end", 3.0, "--dt-out", 2e-4]
FASE_RECUPERAVEL = ["--modo", "rede", "--scr", 5, "--evento", "fase", "--d-fase", 60,
                     "--t-step", 0.5, "--t-end", 4.0, "--dt-out", 1e-4]
FASE_AB = ["--modo", "rede", "--scr", 5, "--evento", "fase", "--d-fase", 80,
           "--t-step", 0.5, "--t-end", 5.0, "--dt-out", 1e-4]
FALTA_LONGA = ["--modo", "rede", "--scr", 5, "--evento", "falta_3f",
               "--vg-falta", 0, "--t-step", 0.5, "--t-clear", 1.40,
               "--t-end", 6.0, "--dt-out", 2e-4]
LIMITADOR = ["--imax-pu", 1.2, "--i-on-pu", 1.176, "--rv-max-pu", 3.0,
             "--xv-rv", 0.0, "--k-aw", K_AW_PADRAO]


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-15)


def checar(cond: Any, msg: str) -> None:
    if not bool(cond):
        raise AssertionError(msg)


def wrap(ang: Any) -> Any:
    return (np.asarray(ang) + np.pi) % (2 * np.pi) - np.pi


def sha256(caminho: str) -> str:
    h = hashlib.sha256()
    with open(caminho, "rb") as fh:
        for bloco in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(bloco)
    return h.hexdigest()


# ======================================================================================
# Oráculo local independente da barreira
# ======================================================================================
def barreira_oraculo(if_complex_a: complex, vc_complex_v: complex, e_bruta_complex_v: complex,
                     rf_ohm: float, i_base_pico_a: float, z_base_ohm: float,
                     imax_pu: float, i_on_pu: float | None, rv_max_pu: float,
                     xv_rv: float) -> dict[str, Any]:
    """Implementação literal do contrato B, sem importar o simulador."""
    I = complex(if_complex_a)
    if imax_pu == 0:
        return {"if_env_pu": abs(I) / i_base_pico_a, "rv_pu": 0.0, "xv_pu": 0.0,
                "ativo": False, "teto_rv_atingido": False, "rv_requerido_pu": 0.0}
    if imax_pu < 0 or i_base_pico_a <= 0 or z_base_ohm <= 0 or rv_max_pu <= 0 or xv_rv < 0:
        raise ValueError("parâmetros inválidos do limitador")
    ion = AUTO_I_ON_FRAC * imax_pu if i_on_pu is None else float(i_on_pu)
    if not 0 <= ion < imax_pu:
        raise ValueError("i_on_pu deve satisfazer 0 <= i_on_pu < imax_pu")
    ipu = abs(I) / i_base_pico_a
    if ipu <= ion or abs(I) <= np.finfo(float).tiny:
        return {"if_env_pu": ipu, "rv_pu": 0.0, "xv_pu": 0.0, "ativo": False,
                "teto_rv_atingido": False, "rv_requerido_pu": 0.0}
    xi = float(np.clip((ipu - ion) / (imax_pu - ion), 0.0, 1.0))
    s = xi * xi * (3.0 - 2.0 * xi)
    u = I / abs(I)
    radial_v = float(np.real(np.conj(u) * (complex(e_bruta_complex_v) - complex(vc_complex_v))))
    r_bar = max(0.0, radial_v / (imax_pu * i_base_pico_a) - rf_ohm + MARGEM_RV_PU * z_base_ohm)
    rv_req_pu = s * r_bar / z_base_ohm
    rv_pu = min(rv_req_pu, rv_max_pu)
    return {"if_env_pu": ipu, "rv_pu": rv_pu, "xv_pu": xv_rv * rv_pu,
            "ativo": rv_pu > 0.0, "teto_rv_atingido": rv_req_pu > rv_max_pu,
            "rv_requerido_pu": rv_req_pu}


def residual_radial(I: complex, Vc: complex, E: complex, rf: float, rv_ohm: float) -> float:
    return float(np.real(np.conj(I) * (E - Vc)) - (rf + rv_ohm) * abs(I) ** 2)


def antiwindup_oraculo(e_bruta: complex, e_aplicada: complex, k_aw: float, ativo: bool) -> float:
    if not ativo or k_aw == 0 or abs(e_bruta) == 0 or abs(e_aplicada) == 0:
        return 0.0
    return float(k_aw * wrap(np.angle(e_aplicada) - np.angle(e_bruta)))


def contar_pole_slips(delta_unwrapped: np.ndarray) -> int:
    bandas = np.floor((np.asarray(delta_unwrapped) + np.pi) / (2 * np.pi)).astype(np.int64)
    return int(np.sum(np.abs(np.diff(bandas)))) if len(bandas) > 1 else 0


def classificar_sincronismo(t: np.ndarray, delta: np.ndarray, f_diff: np.ndarray) -> tuple[str, int]:
    t, delta, f_diff = map(np.asarray, (t, delta, f_diff))
    slips = contar_pole_slips(delta)
    m = t >= t[-1] - TOL["sync_janela_s"]
    fmax = float(np.max(np.abs(f_diff[m])))
    dpp = float(np.ptp(delta[m]))
    if slips >= 1 or fmax > TOL["sync_f_perdido_hz"] or dpp > TOL["sync_delta_perdido_rad"]:
        return "perdido", slips
    if fmax <= TOL["sync_f_mantido_hz"] and dpp <= TOL["sync_delta_mantido_rad"]:
        return "mantido", slips
    return "indeterminado", slips


# ======================================================================================
# Oráculo reduzido de estabilidade transitória
# ======================================================================================
@dataclass
class ResultadoOraculo:
    t: np.ndarray
    delta: np.ndarray
    f_diff: np.ndarray
    i_pu: np.ndarray
    ativo: np.ndarray
    delta_aw: np.ndarray
    sincronismo: str
    pole_slips: int


class OraculoTransit:
    """Fonte de tensão atrás de Z equivalente; não reproduz o LCL da implementação."""

    def __init__(self, scr: float = 5.0, imax: float = 1.2, k_aw: float = K_AW_PADRAO):
        xr = BASE["xr_rede"]
        zg = (1.0 / scr) / math.sqrt(1.0 + 1.0 / xr ** 2)
        self.R = BASE["rf"] + zg / xr
        self.X = BASE["xf"] + zg
        self.E, self.Pm = 1.014920847, 0.5
        self.H, self.Dp = BASE["H"], 1.0 / BASE["mp"]
        self.w0, self.imax, self.k_aw = 2 * np.pi * BASE["f0"], imax, k_aw
        self.delta0 = brentq(lambda d: self.eletrica(d, 1.0, 0.0)[0] - self.Pm, 0.0, np.pi / 2)

    def eletrica(self, delta_v: float, vg: float, delta_g: float) -> tuple[float, float, float, bool]:
        e = self.E * np.exp(1j * delta_v)
        v = vg * np.exp(1j * delta_g)
        z = self.R + 1j * self.X
        i0 = (e - v) / z
        ativo = abs(i0) > self.imax
        i = self.imax * i0 / abs(i0) if ativo else i0
        e_aplicada = v + z * i
        p = float(np.real(v * np.conj(i)))
        aw = antiwindup_oraculo(e, e_aplicada, self.k_aw, ativo)
        return p, aw, abs(i), ativo

    def simular(self, tipo: str, severidade: float, t_end: float = 8.0,
                t_step: float = 0.5, dt: float = 0.002) -> ResultadoOraculo:
        if tipo == "fase":
            t_clear = None
            d_fase = math.radians(severidade)
        elif tipo == "falta_3f":
            t_clear = t_step + severidade
            d_fase = 0.0
        else:
            raise ValueError("tipo do oráculo deve ser fase ou falta_3f")

        def evento(t: float) -> tuple[float, float]:
            if tipo == "fase":
                return 1.0, 0.0 if t < t_step else d_fase
            return (0.0 if t_step <= t < t_clear else 1.0), 0.0

        def f(t: float, x: np.ndarray) -> np.ndarray:
            vg, dg = evento(t)
            pe, aw, _, _ = self.eletrica(x[0], vg, dg)
            return np.array([self.w0 * x[1] + aw, (self.Pm - pe - self.Dp * x[1]) / (2 * self.H)])

        quebras = [0.0, t_step]
        if t_clear is not None:
            quebras.append(t_clear)
        quebras.append(t_end)
        quebras = sorted(set(q for q in quebras if 0 <= q <= t_end))
        x, tt, xx = np.array([self.delta0, 0.0]), [], []
        for a, b in zip(quebras[:-1], quebras[1:]):
            s = solve_ivp(f, (a, b), x, rtol=3e-8, atol=1e-10, max_step=dt)
            if not s.success:
                raise RuntimeError(s.message)
            tt.append(s.t[:-1]); xx.append(s.y[:, :-1]); x = s.y[:, -1]
        tt.append(np.array([t_end])); xx.append(x[:, None])
        t = np.concatenate(tt); X = np.concatenate(xx, axis=1)
        dg = np.where((tipo == "fase") & (t >= t_step), d_fase, 0.0)
        delta = X[0] - dg
        f_diff = BASE["f0"] * X[1]
        vals = [self.eletrica(X[0, k], *evento(float(t[k]))) for k in range(len(t))]
        i_pu = np.array([v[2] for v in vals]); ativo = np.array([v[3] for v in vals], dtype=bool)
        daw = np.array([v[1] for v in vals])
        sincronismo, slips = classificar_sincronismo(t, delta, f_diff)
        return ResultadoOraculo(t, delta, f_diff, i_pu, ativo, daw, sincronismo, slips)

    def duracao_critica_falta(self, lo: float = 0.45, hi: float = 0.90) -> float:
        checar(self.simular("falta_3f", lo).sincronismo == "mantido", "limite inferior não é estável")
        checar(self.simular("falta_3f", hi).sincronismo == "perdido", "limite superior não é instável")
        for _ in range(11):
            mid = (lo + hi) / 2
            if self.simular("falta_3f", mid).sincronismo == "mantido":
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2


# ======================================================================================
# Execução externa
# ======================================================================================
class Ambiente:
    def __init__(self, script: str, gerador: str, v7: str, v6: str, v5: str,
                 testes_v7: str, workdir: str):
        self.script, self.gerador = os.path.abspath(script), os.path.abspath(gerador)
        self.v7, self.v6, self.v5 = map(os.path.abspath, (v7, v6, v5))
        self.testes_v7 = os.path.abspath(testes_v7)
        self.dir, self.cont, self.stdout = workdir, 0, ""
        self.base_json = os.path.join(workdir, "base.json")
        with open(self.base_json, "w", encoding="utf-8") as fh:
            json.dump(BASE, fh)
        self._modulo = None

    def _pref(self, nome: str) -> str:
        self.cont += 1
        return os.path.join(self.dir, f"{self.cont:03d}_{nome}")

    def sim(self, nome: str, args: Iterable[Any], config: str | None = None,
            graficos: bool = False, script: str | None = None, esperar_erro: bool = False,
            ler: bool = True) -> tuple[Any, Any, float] | subprocess.CompletedProcess[str]:
        pref = self._pref(nome)
        cmd = [sys.executable, script or self.script, "--config", config or self.base_json,
               *[str(a) for a in args], "--prefixo", pref]
        if not graficos:
            cmd.append("--sem-graficos")
        t0 = time.perf_counter()
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1800,
                            encoding="utf-8", errors="replace")
        dt = time.perf_counter() - t0
        self.stdout = (pr.stdout or "") + (pr.stderr or "")
        if esperar_erro:
            return pr
        if pr.returncode != 0:
            raise AssertionError(f"simulador falhou (código {pr.returncode}):\n{self.stdout[-1800:]}")
        if not ler:
            return None, None, dt
        with open(pref + "_resultados.json", encoding="utf-8") as fh:
            js = json.load(fh)
        with np.load(pref + "_series.npz") as arq:
            ser = {k: np.asarray(arq[k]) for k in arq.files}
        return js, ser, dt

    def gerar(self, nome: str, args: Iterable[Any], simular: bool = False) -> tuple[subprocess.CompletedProcess[str], str]:
        saida = self._pref(nome) + ".yaml"
        cmd = [sys.executable, self.gerador, *[str(a) for a in args], "-o", saida]
        if simular:
            cmd += ["--simular", "--prefixo", self._pref(nome + "_sim")]
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1200,
                            encoding="utf-8", errors="replace")
        return pr, saida

    def modulo_v8(self):
        if self._modulo is None:
            spec = importlib.util.spec_from_file_location("vsg_v8_sob_teste", self.script)
            checar(spec is not None and spec.loader is not None, "não foi possível criar spec da v8")
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            self._modulo = mod
        return self._modulo


def rede(scr: float, dw: float = 0.0, tw: float = 1.0) -> list[Any]:
    return ["--modo", "rede", "--scr", scr, "--dw", dw, "--tw", tw]


def curto(extra: Iterable[Any] = ()) -> list[Any]:
    return ["--evento", "nenhum", "--t-step", 0.2, "--t-end", 0.3, "--dt-out", 1e-4, *extra]


def comparar_series(s_ref: dict[str, np.ndarray], s_novo: dict[str, np.ndarray],
                    nomes: dict[str, str] = SERIES_COMUNS) -> float:
    checar(len(s_ref["t"]) == len(s_novo["t"]) and
           np.allclose(s_ref["t"], s_novo["t"], atol=1e-12, rtol=0),
           "vetores de tempo diferentes")
    pior = 0.0
    for k, tol_nome in nomes.items():
        checar(k in s_ref and k in s_novo, f"série comum ausente: {k}")
        d = float(np.nanmax(np.abs(s_ref[k] - s_novo[k])))
        pior = max(pior, d / TOL[tol_nome])
        checar(d <= TOL[tol_nome], f"{k}: diferença {d:.3e} > {TOL[tol_nome]:.1e}")
    return pior


def autovalores(js: dict[str, Any]) -> np.ndarray:
    return np.array([complex(r, i) for r, i in js["modal"]["autovalores"] if abs(complex(r, i)) >= 1e-6])


def bloco_limitador(js: dict[str, Any]) -> dict[str, Any]:
    lim = js.get("limitador")
    checar(isinstance(lim, dict), "bloco JSON 'limitador' ausente")
    for k in CAMPOS_LIMITADOR:
        checar(k in lim, f"bloco 'limitador' sem o campo '{k}'")
    checar(lim["sincronismo"] in ("mantido", "perdido", "indeterminado"),
           f"sincronismo inválido: {lim['sincronismo']!r}")
    return lim


def checar_series_v8(s: dict[str, np.ndarray]) -> None:
    n = len(s["t"])
    for k in SERIES_V8:
        checar(k in s, f"série v8 ausente: {k}")
        checar(len(s[k]) == n, f"série {k} tem {len(s[k])} amostras, esperado {n}")
        checar(np.all(np.isfinite(s[k])), f"série {k} contém NaN/Inf")


def verificar_teto(js: dict[str, Any], s: dict[str, np.ndarray], exigir_ativo: bool = True) -> str:
    lim = bloco_limitador(js); checar_series_v8(s)
    imax = float(lim["imax_pu"])
    pico = float(np.max(s["if_env_pu"]))
    checar(abs(float(lim["I_pico_env_pu"]) - pico) <= TOL["saida_rel"] * max(1.0, pico),
           f"I_pico_env_pu JSON {lim['I_pico_env_pu']} != NPZ {pico}")
    if exigir_ativo:
        checar(np.any(s["limitador_ativo"] > 0.5), "limitador não atuou no caso que deveria saturar")
    if not bool(lim["teto_rv_atingido"]):
        checar(pico <= imax * (1 + TOL["teto_corrente_frac"]),
               f"corrente {pico:.5f} pu > teto tolerado {imax * 1.005:.5f} pu")
    return f"I_pico={pico:.4f} pu; ativo={lim['tempo_ativo_s']:.4f} s; sync={lim['sincronismo']}"


# ======================================================================================
# Regressões R1-R5
# ======================================================================================
def R1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Ilha histórica da v5: v8 sem limite reproduz os valores publicados."""
    y = os.path.join(amb.dir, "caso_v5.yaml")
    with open(y, "w", encoding="utf-8") as fh:
        fh.write(YAML_V5)
    js, s, _ = amb.sim("R1", ["--imax-pu", 0], config=y)
    checar(js.get("versao") == "v8", f"versao={js.get('versao')!r}, esperado 'v8'")
    for bloco, chave in (("antes", "antes_evento"), ("final", "final")):
        for k, refv in REF_V5[bloco].items():
            v = float(js[chave][k])
            ok = abs(v - refv) <= TOL["R1_f_abs"] if k == "f_hz" else rel(v, refv) <= TOL["R1_rel"]
            checar(ok, f"{bloco}.{k}={v:.6g}, referência v5={refv}")
    checar(rel(js["rocof_max_hz_s"], REF_V5["rocof_max_hz_s"]) <= TOL["R1_rel"], "RoCoF difere da v5")
    lim = bloco_limitador(js); checar(not lim["habilitado"], "imax=0 deveria desabilitar o limitador")
    checar_series_v8(s)
    for k in ("rv_pu", "xv_pu", "limitador_ativo", "delta_aw_rad_s"):
        checar(np.all(s[k] == 0), f"{k} não é exatamente zero com imax_pu=0")
    return f"Vfinal={js['final']['Vll_rms_v']:.3f} V; ffinal={js['final']['f_hz']:.4f} Hz"


R2_CASOS = [
    ("pref_scr2", [*rede(2), "--evento", "pref", "--t-step", .5, "--t-end", 3, "--dt-out", 1e-4]),
    ("pref_scr5_dw", [*rede(5, 103), "--evento", "pref", "--t-step", .5, "--t-end", 3, "--dt-out", 1e-4]),
    ("carga_scr5", [*rede(5, 100), "--evento", "carga", "--t-step", .5, "--t-end", 3, "--dt-out", 2e-4]),
    ("fase_scr20", [*rede(20, 100), "--evento", "fase", "--d-fase", 5, "--t-step", .5, "--t-end", 2, "--dt-out", 5e-5]),
    ("carga_ilha", ["--modo", "ilhado", "--evento", "carga", "--dw", 100, "--tw", .1,
                     "--t-step", .5, "--t-end", 3, "--dt-out", 1e-4]),
]


def R2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Eventos da v6: séries comuns da v8 sem limite são iguais às da v7."""
    checar(os.path.isfile(amb.v7), f"v7 não encontrada: {amb.v7}")
    out = []
    for nome, args in R2_CASOS:
        _, s7, _ = amb.sim("R2_" + nome + "_v7", args, script=amb.v7)
        _, s8, _ = amb.sim("R2_" + nome + "_v8", [*args, "--imax-pu", 0])
        out.append(f"{nome} {comparar_series(s7, s8):.2f}xtol")
    return "; ".join(out)


def R3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """B1-B10 da suíte v7 passam na v8 com imax_pu=0 (adaptação apenas do rótulo de versão)."""
    checar(os.path.isfile(amb.testes_v7), f"testes_v7.py não encontrado: {amb.testes_v7}")
    h = sha256(amb.testes_v7)
    checar(h == HASH_TESTES_V7, f"testes_v7.py alterado: {h}")
    txt = open(amb.testes_v7, encoding="utf-8").read()
    txt = txt.replace('js.get("versao") == "v7"', 'js.get("versao") == "v8"')
    txt = txt.replace("campo 'versao' diferente de 'v7'", "campo 'versao' diferente de 'v8'")
    adapt = os.path.join(amb.dir, "testes_v7_contrato_v8.py")
    with open(adapt, "w", encoding="utf-8") as fh:
        fh.write(txt)
    with open(os.path.join(amb.dir, "testes_v7.sha256"), "w", encoding="ascii") as fh:
        fh.write(sha256(adapt) + "  testes_v7_contrato_v8.py\n")
    selecionados = "B1,B2,B3,B4,B5,B6,B7,B8,B9,B10"
    cmd = [sys.executable, adapt, "--so", selecionados, "--script", amb.script,
           "--gerador", amb.gerador, "--v6", amb.v6, "--v5", amb.v5]
    pr = subprocess.run(cmd, cwd=amb.dir, capture_output=True, text=True, timeout=2400,
                        encoding="utf-8", errors="replace")
    checar(pr.returncode == 0, "B1-B10 não passaram contra a v8:\n" + ((pr.stdout or "") + (pr.stderr or ""))[-2500:])
    return "10/10 contratos funcionais da v7 aprovados; fonte v7 íntegra"


def R4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Análise modal sem limitação: autovalores da v8 coincidem com os da v7."""
    out = []
    for scr in (2, 5, 20):
        for dw in (0, 100, 103):
            j7, _, _ = amb.sim(f"R4_{scr}_{dw}_v7", [*rede(scr, dw), *curto()], script=amb.v7)
            j8, _, _ = amb.sim(f"R4_{scr}_{dw}_v8", [*rede(scr, dw), *curto(["--imax-pu", 0])])
            l7, l8 = autovalores(j7), autovalores(j8)
            checar(len(l7) == len(l8), f"SCR{scr}/Dw{dw}: número de autovalores diferente")
            pior = 0.0
            for lam in l7:
                d = float(np.min(np.abs(l8 - lam)) / max(1.0, abs(lam)))
                pior = max(pior, d); checar(d <= TOL["modal_rel"], f"autovalor {lam} ausente na v8")
            out.append(f"{scr}/{dw}:{pior:.1e}")
    return ", ".join(out)


def R5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Contrato de CLI, JSON, NPZ, YAML histórico e gerador permanece compatível."""
    help_sim = subprocess.run([sys.executable, amb.script, "--help"], capture_output=True, text=True,
                              timeout=120, encoding="utf-8", errors="replace")
    checar(help_sim.returncode == 0, "--help do simulador falhou")
    for op in ("--evento", "--df-g", "--rocof-g", "--imax-pu", "--i-on-pu", "--rv-max-pu",
               "--xv-rv", "--k-aw", "--vg-falta", "--t-clear", "--gerar-config"):
        checar(op in help_sim.stdout, f"opção ausente no --help: {op}")
    j7, s7, _ = amb.sim("R5_v7", [*rede(5), *curto()], script=amb.v7)
    j8, s8, _ = amb.sim("R5_v8", [*rede(5), *curto(["--imax-pu", 0])])
    for k in j7:
        if k != "versao":
            checar(k in j8, f"campo JSON antigo removido: {k}")
    for k in s7:
        checar(k in s8, f"série NPZ antiga removida: {k}")
    pr, y = amb.gerar("R5", ["--modo", "rede", "--scr", 5, "--evento", "freq_rampa",
                                      "--df-rede", -.1, "--rocof-rede", .1])
    checar(pr.returncode == 0 and os.path.isfile(y), "gerador não produziu YAML histórico")
    amb.sim("R5_yaml", ["--imax-pu", 0], config=y)
    return f"{len(j7)} campos JSON e {len(s7)} séries antigas preservadas; YAML executável"


# ======================================================================================
# Aceitação funcional C1-C14
# ======================================================================================
def C1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Limitador habilitado mas inativo: termos v8 são zero e resultado coincide com a v7."""
    args = [*rede(5), "--evento", "nenhum", "--t-step", .2, "--t-end", .5, "--dt-out", 1e-4]
    _, s7, _ = amb.sim("C1_v7", args, script=amb.v7)
    js, s8, _ = amb.sim("C1_v8", [*args, "--imax-pu", 5, "--i-on-pu", 4.9,
                                   "--rv-max-pu", RV_MAX_PADRAO_PU, "--k-aw", K_AW_PADRAO])
    pior = comparar_series(s7, s8)
    lim = bloco_limitador(js); checar(lim["habilitado"], "limitador deveria estar habilitado")
    for k in ("rv_pu", "xv_pu", "limitador_ativo", "delta_aw_rad_s"):
        checar(np.all(s8[k] == 0), f"{k} não é zero no caso não limitado")
    return f"regressão {pior:.2f}xtol; todos os termos de atuação zero"


def C2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Função pura: fórmula, quadrantes e condição radial na fronteira de I_max."""
    mod = amb.modulo_v8()
    fn = getattr(mod, "calcular_impedancia_virtual", None)
    checar(callable(fn), "função pública calcular_impedancia_virtual ausente")
    ib, zb, rf, imax, ion = 100.0, 1.444, 0.004, 1.2, 1.176
    n = 0
    for ang in np.linspace(-np.pi, np.pi, 9)[:-1]:
        I = imax * ib * np.exp(1j * ang)
        for sinal in (-1, 1):
            U = sinal * 250.0 * np.exp(1j * ang) + 40j * np.exp(1j * ang)
            kw = dict(if_complex_a=I, vc_complex_v=30-20j, e_bruta_complex_v=30-20j+U,
                      rf_ohm=rf, i_base_pico_a=ib, z_base_ohm=zb, imax_pu=imax,
                      i_on_pu=ion, rv_max_pu=20.0, xv_rv=0.3)
            got, ref = fn(**kw), barreira_oraculo(**kw)
            for k in ("if_env_pu", "rv_pu", "xv_pu", "rv_requerido_pu"):
                checar(rel(float(got[k]), float(ref[k])) <= TOL["formula_rel"], f"{k} difere do contrato")
            res = residual_radial(I, kw["vc_complex_v"], kw["e_bruta_complex_v"], rf,
                                  float(got["rv_pu"]) * zb)
            checar(res <= TOL["barreira_abs"] * max(1.0, abs(I) ** 2),
                   f"derivada radial positiva na fronteira: {res:.4e}")
            n += 1
    pequeno = fn(0j, 0j, 1+0j, rf, ib, zb, imax, ion, 20.0, 0.0)
    checar(all(np.isfinite(float(pequeno[k])) for k in ("if_env_pu", "rv_pu", "xv_pu")),
           "NaN/Inf em corrente nula")
    return f"{n} estados na fronteira; residual radial não positivo"


def C3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Varredura: zero abaixo de I_on, continuidade e R_v monotônico com a severidade."""
    fn = getattr(amb.modulo_v8(), "calcular_impedancia_virtual")
    ib, zb, rf, imax, ion = 100.0, 1.444, 0.004, 1.2, 1.176
    mags = np.linspace(0, 1.25, 2501)
    rv = []
    for m in mags:
        I = m * ib * np.exp(0.37j)
        U = 250.0 * np.exp(0.37j)
        r = fn(I, 10-5j, 10-5j+U, rf, ib, zb, imax, ion, 20.0, 0.4)
        rv.append(float(r["rv_pu"]))
        checar(float(r["rv_pu"]) >= 0 and float(r["xv_pu"]) >= 0, "impedância negativa")
    rv = np.asarray(rv)
    checar(np.all(rv[mags <= ion] == 0), "R_v não é zero até I_on")
    checar(np.min(np.diff(rv)) >= -1e-12, "R_v diminui com severidade crescente")
    k = int(np.argmin(abs(mags-ion)))
    checar(rv[k] <= TOL["continuidade_pu"], f"salto em I_on: {rv[k]:.3e} pu")
    return f"{len(mags)} pontos; Rv_max_varredura={rv.max():.5f} pu"


def C4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """R_v,max insuficiente: saturação e violação são quantificadas e avisadas."""
    js, s, _ = amb.sim("C4", [*FALTA_CURTA, "--imax-pu", 1.05, "--i-on-pu", 1.029,
                              "--rv-max-pu", .01, "--xv-rv", 0, "--k-aw", K_AW_PADRAO])
    lim = bloco_limitador(js)
    checar(bool(lim["teto_rv_atingido"]), "teto_rv_atingido deveria ser true")
    checar(float(lim["violacao_max_pu"]) > 0 and float(lim["tempo_em_violacao_s"]) > 0,
           "violação não foi quantificada")
    avisos = [l for l in amb.stdout.splitlines() if l.strip().startswith("(!)")]
    checar(any("corrente" in l.lower() and "r_v,max" in l.lower() for l in avisos),
           "faltou aviso (!) contendo corrente e R_v,max")
    return f"violacao={lim['violacao_max_pu']:.4f} pu por {lim['tempo_em_violacao_s']:.4f} s"


def C5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Falta trifásica curta: corrente limitada, estados contínuos e recuperação sem pole slip."""
    js, s, _ = amb.sim("C5", [*FALTA_CURTA, *LIMITADOR])
    det = verificar_teto(js, s)
    lim = bloco_limitador(js)
    checar(lim["sincronismo"] == "mantido" and int(lim["pole_slips"]) == 0,
           f"falta curta não recuperou: {lim['sincronismo']}, slips={lim['pole_slips']}")
    t = s["t"]; checar(len(np.unique(t)) == len(t) and np.all(np.diff(t) > 0), "tempo duplicado/não crescente")
    for q in (.5, .6):
        checar(np.any(np.isclose(t, q, atol=1e-12)), f"ponto de quebra {q} s ausente")
        k = int(np.searchsorted(t, q))
        if 0 < k < len(t):
            checar(abs(float(s["delta_v_rad"][k] - s["delta_v_rad"][k-1])) < .05 and
                   abs(float(s["f_hz"][k] - s["f_hz"][k-1])) < .05,
                   f"estado eletromecânico descontínuo em {q} s")
    m = t >= t[-1] - .2
    checar(np.max(np.abs(s["f_hz"][m] - s["f_rede_hz"][m])) <= TOL["sync_f_mantido_hz"],
           "frequência não recuperou")
    checar(np.max(np.abs(s["vg_aplicada_pu"][(t >= .5) & (t < .6)])) <= 1e-12,
           "vg_falta=0 não foi aplicada exatamente")
    cache["C5"] = (js, s)
    return det


def C6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Falta curta em SCR 2, 5 e 20 com dimensionamento produzido pelo gerador v8."""
    out = []
    for scr in (2, 5, 20):
        args = ["--modo", "rede", "--scr", scr, "--evento", "falta_3f", "--vg-falta", 0,
                "--t-step", .5, "--t-clear", .6, "--t-end", 3, "--dt-out", 1e-4,
                "--imax-pu", 1.2, "--i-on-pu", 1.176, "--k-aw", K_AW_PADRAO]
        pr, y = amb.gerar(f"C6_scr{scr}", args)
        checar(pr.returncode == 0 and os.path.isfile(y), f"gerador falhou em SCR {scr}:\n{(pr.stderr or pr.stdout)[-1000:]}")
        checar("I_max" in pr.stdout and ("R_v,max" in pr.stdout or "R_v" in pr.stdout),
               f"gerador SCR {scr} não informou I_max e R_v,max")
        js, s, _ = amb.sim(f"C6_scr{scr}_sim", [], config=y)
        det = verificar_teto(js, s)
        checar(not bloco_limitador(js)["teto_rv_atingido"], f"gerador subdimensionou R_v,max em SCR {scr}")
        out.append(f"SCR{scr}: {det}")
    return "; ".join(out)


def C7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Rampa inercial da v7 com I_max=1,2: corrente limitada e métrica marcada como limitada."""
    js, s, _ = amb.sim("C7", [*RAMPA_LIMITADA, *LIMITADOR])
    det = verificar_teto(js, s)
    ine = js.get("inercial") or {}
    checar(ine.get("limitado") is True, "inercial.limitado deveria ser true")
    checar("limit" in amb.stdout.lower() and ("inercial" in amb.stdout.lower() or "2h" in amb.stdout.lower()),
           "relatório não ressalvou a métrica inercial limitada")
    cache["C7"] = (js, s)
    return det + f"; H_eff={ine.get('H_eff_s')} s; limitado=true"


def C8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Estresse angular recuperável: teto respeitado, zero pole slip e convergência final."""
    js, s, _ = amb.sim("C8", [*FASE_RECUPERAVEL, *LIMITADOR])
    det = verificar_teto(js, s)
    lim = bloco_limitador(js)
    checar(lim["sincronismo"] == "mantido" and int(lim["pole_slips"]) == 0,
           f"caso recuperável classificado {lim['sincronismo']}")
    cls, slips = classificar_sincronismo(s["t"], s["delta_rel_unwrapped_rad"],
                                        s["f_hz"] - s["f_rede_hz"])
    checar(cls == lim["sincronismo"] and slips == int(lim["pole_slips"]),
           "classificação JSON não coincide com as séries")
    cache["C8"] = (js, s)
    return det


def erro_final(s: dict[str, np.ndarray]) -> tuple[float, float]:
    m = s["t"] >= s["t"][-1] - TOL["sync_janela_s"]
    return float(np.ptp(s["delta_rel_unwrapped_rad"][m])), float(np.max(np.abs(s["f_hz"][m]-s["f_rede_hz"][m])))


def tempo_recuperacao(s: dict[str, np.ndarray], ts: float = .5) -> float:
    ok = (np.abs(s["f_hz"] - s["f_rede_hz"]) <= TOL["sync_f_mantido_hz"])
    for k in np.flatnonzero((s["t"] >= ts) & ok):
        if np.all(ok[k:]):
            return float(s["t"][k] - ts)
    return float("inf")


def C9(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Comparação A/B: anti-windup nominal não piora erro final nem recuperação."""
    j0, s0, _ = amb.sim("C9_k0", [*FASE_AB, *LIMITADOR[:-2], "--k-aw", 0])
    j1, s1, _ = amb.sim("C9_kaw", [*FASE_AB, *LIMITADOR])
    verificar_teto(j0, s0); verificar_teto(j1, s1)
    d0, f0 = erro_final(s0); d1, f1 = erro_final(s1)
    tr0, tr1 = tempo_recuperacao(s0), tempo_recuperacao(s1)
    checar(d1 <= d0 + .01 and f1 <= f0 + .005 and tr1 <= tr0 + .05,
           f"anti-windup piorou: delta {d0:.3g}->{d1:.3g}, f {f0:.3g}->{f1:.3g}, tr {tr0:.3g}->{tr1:.3g}")
    checar(np.all(s0["delta_aw_rad_s"] == 0), "delta_aw não é zero com k_aw=0")
    checar(np.any(np.abs(s1["delta_aw_rad_s"]) > 0), "anti-windup nominal não atuou")
    return f"delta_pp {d0:.4f}->{d1:.4f} rad; fmax {f0:.4f}->{f1:.4f} Hz; trec {tr0:.3f}->{tr1:.3f} s"


def C10(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Falta longa: perda de sincronismo permanece detectável com corrente limitada."""
    js, s, _ = amb.sim("C10", [*FALTA_LONGA, *LIMITADOR])
    det = verificar_teto(js, s)
    lim = bloco_limitador(js)
    m = s["t"] >= s["t"][-1] - TOL["sync_janela_s"]
    persistente = (np.max(np.abs(s["f_hz"][m] - s["f_rede_hz"][m])) > TOL["sync_f_perdido_hz"] or
                   np.ptp(s["delta_rel_unwrapped_rad"][m]) > TOL["sync_delta_perdido_rad"])
    checar(lim["sincronismo"] == "perdido" and (int(lim["pole_slips"]) >= 1 or persistente),
           f"caso instável mascarado: {lim['sincronismo']}, slips={lim['pole_slips']}")
    checar(not bool(lim["teto_rv_atingido"]), "C10 deve testar perda física, não R_v,max insuficiente")
    cache["C10"] = (js, s)
    return det + f"; slips={lim['pole_slips']}; divergencia_persistente={persistente}"


def C11(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fronteira crítica de falta: transição da v8 ocorre na faixa congelada do oráculo."""
    cct = OraculoTransit().duracao_critica_falta()
    baixo, alto = cct - TOL["fronteira_falta_s"], cct + TOL["fronteira_falta_s"]
    classes = []
    for nome, dur in (("baixo", baixo), ("centro", cct), ("alto", alto)):
        args = ["--modo", "rede", "--scr", 5, "--evento", "falta_3f", "--vg-falta", 0,
                "--t-step", .5, "--t-clear", .5 + dur, "--t-end", 8, "--dt-out", 2e-4, *LIMITADOR]
        js, s, _ = amb.sim("C11_" + nome, args)
        verificar_teto(js, s)
        classes.append(bloco_limitador(js)["sincronismo"])
    checar(classes[0] == "mantido", f"abaixo da fronteira: {classes[0]}")
    checar(classes[2] == "perdido", f"acima da fronteira: {classes[2]}")
    checar(classes[1] in ("mantido", "perdido", "indeterminado"), "classe central inválida")
    return f"CCT_oraculo={cct:.3f} s; classes -0,15/centro/+0,15={classes}"


def metrica_conv(js: dict[str, Any], s: dict[str, np.ndarray]) -> tuple[float, float, float, str]:
    lim = bloco_limitador(js)
    return (float(np.max(s["if_env_pu"])), float(lim["rv_max_usado_pu"]),
            float(lim["tempo_ativo_s"]), str(lim["sincronismo"]))


def C12(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Convergência: C5, C7, C8 e C10 com max_step e rtol dez vezes menores."""
    casos = {"C5": FALTA_CURTA, "C7": RAMPA_LIMITADA, "C8": FASE_RECUPERAVEL, "C10": FALTA_LONGA}
    out = []
    for nome, args in casos.items():
        if nome in cache:
            j1, s1 = cache[nome]
        else:
            j1, s1, _ = amb.sim("C12_" + nome + "_base", [*args, *LIMITADOR])
        j2, s2, _ = amb.sim("C12_" + nome + "_fino", [*args, *LIMITADOR, "--max-step", 2e-5, "--rtol", 1e-9])
        a, b = metrica_conv(j1, s1), metrica_conv(j2, s2)
        for k, rot in enumerate(("I_pico", "Rv_max", "tempo_ativo")):
            checar(rel(a[k], b[k]) <= TOL["convergencia_frac"], f"{nome} {rot}: {a[k]} vs {b[k]}")
        checar(a[3] == b[3], f"{nome}: sincronismo mudou {a[3]}->{b[3]}")
        out.append(f"{nome}:{max(rel(a[k], b[k]) for k in range(3))*100:.3f}%")
    return "; ".join(out)


def C13(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Entradas inválidas da v8 falham com Erro e sem traceback."""
    casos = {
        "imax<0": ["--imax-pu", -1],
        "ion=imax": ["--imax-pu", 1.2, "--i-on-pu", 1.2, "--rv-max-pu", 2],
        "ion>imax": ["--imax-pu", 1.2, "--i-on-pu", 1.21, "--rv-max-pu", 2],
        "rvmax=0": ["--imax-pu", 1.2, "--i-on-pu", 1.18, "--rv-max-pu", 0],
        "xv_rv<0": ["--imax-pu", 1.2, "--i-on-pu", 1.18, "--rv-max-pu", 2, "--xv-rv", -.1],
        "k_aw<0": ["--imax-pu", 1.2, "--i-on-pu", 1.18, "--rv-max-pu", 2, "--k-aw", -1],
        "clear<=step": [*rede(5), "--evento", "falta_3f", "--t-step", .5, "--t-clear", .5, "--t-end", 1],
        "clear>=end": [*rede(5), "--evento", "falta_3f", "--t-step", .5, "--t-clear", 1, "--t-end", 1],
        "falta_ilha": ["--modo", "ilhado", "--evento", "falta_3f", "--t-step", .5, "--t-clear", .6],
    }
    for nome, args in casos.items():
        pr = amb.sim("C13_" + nome, args, esperar_erro=True)
        txt = (pr.stdout or "") + (pr.stderr or "")
        checar(pr.returncode != 0, f"{nome}: deveria falhar")
        checar("Traceback" not in txt, f"{nome}: traceback exposto")
        checar(re.search(r"(^|\n)Erro", txt) is not None, f"{nome}: mensagem não inicia por Erro")
    return f"{len(casos)} combinações inválidas rejeitadas corretamente"


def C14(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Desempenho: caso v8 de 4 s com falta e gráficos <=2x referência v7 equivalente."""
    comum = [*rede(5), "--t-step", .5, "--t-end", 4, "--dt-out", 1e-4]
    _, _, t7 = amb.sim("C14_v7", [*comum, "--evento", "fase", "--d-fase", 5],
                       script=amb.v7, graficos=True, ler=False)
    _, _, t8 = amb.sim("C14_v8", [*comum, "--evento", "falta_3f", "--vg-falta", 0,
                                    "--t-clear", .6, *LIMITADOR], graficos=True, ler=False)
    razao = t8 / max(t7, 1e-9)
    checar(razao <= TOL["desempenho_razao"], f"v8 {t8:.1f}s / v7 {t7:.1f}s = {razao:.2f} > 2")
    return f"v7={t7:.1f} s; v8={t8:.1f} s; razao={razao:.2f}"


TESTES = [("R1", R1), ("R2", R2), ("R3", R3), ("R4", R4), ("R5", R5),
          ("C1", C1), ("C2", C2), ("C3", C3), ("C4", C4), ("C5", C5),
          ("C6", C6), ("C7", C7), ("C8", C8), ("C9", C9), ("C10", C10),
          ("C11", C11), ("C12", C12), ("C13", C13), ("C14", C14)]


# ======================================================================================
# Autoteste — sem implementação v8
# ======================================================================================
def autoteste() -> bool:
    ok = True

    def item(cond: Any, msg: str) -> None:
        nonlocal ok
        cond = bool(cond); ok &= cond
        print(f"  [{'ok' if cond else 'FALHA'}] {msg}")

    # Bases.
    zb = BASE["vll"] ** 2 / BASE["sn"]
    ib = math.sqrt(2) * BASE["sn"] / (math.sqrt(3) * BASE["vll"])
    item(abs(zb - 1.444) < 1e-12 and abs(ib - 214.867521) < 1e-6,
         f"bases: Zb={zb:.6f} ohm, Ibase_pico={ib:.6f} A")

    # Desabilitação exata, auto de I_on e ausência de NaN em I=0.
    z = barreira_oraculo(0j, 0j, 1+0j, .004, ib, zb, 0, None, 2, 0)
    item(z["rv_pu"] == z["xv_pu"] == 0 and not z["ativo"], "imax=0 zera exatamente a v8")
    z = barreira_oraculo(0j, 0j, 1+0j, .004, ib, zb, 1.2, None, 2, 0)
    item(np.isfinite(z["if_env_pu"]) and z["rv_pu"] == 0, "corrente nula é finita e não ativa")

    # Fronteira em quatro quadrantes e dois sentidos radiais.
    max_res = -np.inf
    for ang in np.linspace(-np.pi, np.pi, 9)[:-1]:
        I = 1.2 * ib * np.exp(1j * ang)
        for sinal in (-1, 1):
            U = sinal * 300 * np.exp(1j * ang) + 25j * np.exp(1j * ang)
            r = barreira_oraculo(I, 20-10j, 20-10j+U, .004, ib, zb, 1.2, 1.176, 20, .3)
            res = residual_radial(I, 20-10j, 20-10j+U, .004, r["rv_pu"] * zb)
            max_res = max(max_res, res)
    item(max_res <= TOL["barreira_abs"] * ib ** 2, f"barreira: maior residual radial={max_res:.3e}")

    # Suavidade e monotonicidade.
    mags = np.linspace(1.15, 1.21, 1201)
    rv = [barreira_oraculo(m*ib, 0j, 300+0j, .004, ib, zb, 1.2, 1.176, 20, 0)["rv_pu"] for m in mags]
    item(max(np.diff(rv)) > 0 and min(np.diff(rv)) >= -1e-13 and
         max(np.asarray(rv)[mags <= 1.176]) == 0, "ativação suave, nula e monotônica")

    # Anti-windup: sinal puxa o comando para a tensão aplicada e é zero fora da limitação.
    aw = antiwindup_oraculo(np.exp(.4j), np.exp(.2j), 20, True)
    item(aw < 0 and abs(aw + 4) < 1e-12 and antiwindup_oraculo(1+0j, 1j, 20, False) == 0,
         f"anti-windup: sinal={aw:.3f} rad/s e zero fora da limitação")

    # Classificador sintético.
    t = np.linspace(0, 2, 2001)
    d_ok = .2 + .01*np.exp(-4*t)*np.sin(10*t); f_ok = .02*np.exp(-4*t)
    d_bad = .2 + 2*np.pi*1.2*t; f_bad = np.ones_like(t)*.3
    c1, p1 = classificar_sincronismo(t, d_ok, f_ok)
    c2, p2 = classificar_sincronismo(t, d_bad, f_bad)
    item(c1 == "mantido" and p1 == 0 and c2 == "perdido" and p2 >= 1,
         f"classificador: estável={c1}/{p1}, instável={c2}/{p2}")

    # Oráculo reduzido: equilíbrio, caso estável, A/B e caso instável.
    o = OraculoTransit()
    p0, _, i0, a0 = o.eletrica(o.delta0, 1.0, 0.0)
    item(abs(p0-o.Pm) < 1e-10 and not a0 and i0 < o.imax,
         f"equilíbrio do oráculo: delta0={math.degrees(o.delta0):.3f} graus, I={i0:.3f} pu")
    curto_ = o.simular("falta_3f", .10, 4)
    longo = o.simular("falta_3f", .90, 8)
    item(curto_.sincronismo == "mantido" and curto_.pole_slips == 0,
         f"falta curta do oráculo: {curto_.sincronismo}, slips={curto_.pole_slips}")
    item(longo.sincronismo == "perdido" and longo.pole_slips >= 1,
         f"falta longa do oráculo: {longo.sincronismo}, slips={longo.pole_slips}")
    sem_aw = OraculoTransit(k_aw=0).simular("fase", 80, 6)
    com_aw = OraculoTransit(k_aw=K_AW_PADRAO).simular("fase", 80, 6)
    e0, e1 = np.ptp(sem_aw.delta[sem_aw.t >= 5.5]), np.ptp(com_aw.delta[com_aw.t >= 5.5])
    item(e1 <= e0 and np.any(np.abs(com_aw.delta_aw) > 0),
         f"A/B oráculo: excursão final {e0:.5f}->{e1:.5f} rad")
    cct = o.duracao_critica_falta()
    item(.60 < cct < .72, f"fronteira crítica congelada do oráculo: {cct:.4f} s")
    return ok


def verificar_integridade() -> tuple[str, str]:
    eu = os.path.abspath(__file__)
    atual = sha256(eu)
    arq = os.path.join(os.path.dirname(eu), "testes_v8_1.sha256")
    if not os.path.isfile(arq):
        return atual, "SEM REFERÊNCIA (testes_v8_1.sha256 ausente)"
    esperado = open(arq, encoding="utf-8").read().split()[0].strip().lower()
    return atual, "ÍNTEGRO" if esperado == atual else "ALTERADO — resultados sem validade"


def escrever_relatorios(resultados: list[dict[str, Any]], estado: str, hash_atual: str,
                         script: str, total_s: float) -> None:
    n_ok = sum(r["status"] == "APROVADO" for r in resultados)
    rel_json = {"sha256_testes": hash_atual, "integridade": estado, "script": script,
                "sha256_script": sha256(script), "aprovados": n_ok, "total": len(resultados),
                "contrato": {"i_on_auto_frac": AUTO_I_ON_FRAC, "margem_rv_pu": MARGEM_RV_PU,
                             "rv_max_padrao_pu": RV_MAX_PADRAO_PU, "k_aw_padrao_s_1": K_AW_PADRAO,
                             "tolerancia_teto_frac": TOL["teto_corrente_frac"]},
                "resultados": resultados}
    with open("relatorio_testes_v8.json", "w", encoding="utf-8") as fh:
        json.dump(rel_json, fh, ensure_ascii=False, indent=2)
    with open("relatorio_testes_v8.md", "w", encoding="utf-8") as fh:
        fh.write(f"# Relatório de aceitação — v8\n\n- testes_v8_1.py SHA-256: `{hash_atual}` ({estado})\n"
                 f"- implementação: `{script}` SHA-256 `{rel_json['sha256_script']}`\n"
                 f"- resultado: **{n_ok}/{len(resultados)} aprovados** em {total_s:.0f} s\n\n"
                 "| Teste | Status | Tempo (s) | Descrição | Detalhe |\n|---|---|---|---|---|\n")
        for r in resultados:
            det = str(r["detalhe"]).replace("\n", " ").replace("|", "/")[:500]
            fh.write(f"| {r['teste']} | {r['status']} | {r['tempo_s']} | {r['descricao']} | {det} |\n")


def main() -> None:
    for st in (sys.stdout, sys.stderr):
        try:
            st.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    aqui = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Suíte de aceitação congelada da v8.")
    ap.add_argument("--script", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v8.py"))
    ap.add_argument("--gerador", default=os.path.join(aqui, "gerar_caso_vsg.py"))
    ap.add_argument("--v7", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v7.py"))
    ap.add_argument("--v6", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v6.py"))
    ap.add_argument("--v5", default=os.path.join(aqui, "vsg_2a_ordem_degrau_carga_v5.py"))
    ap.add_argument("--testes-v7", default=os.path.join(aqui, "testes_v7.py"))
    ap.add_argument("--so", default=None, help="lista separada por vírgulas, por exemplo R1,C1,C2")
    ap.add_argument("--autoteste", action="store_true", help="valida contrato e oráculos sem usar a v8")
    ap.add_argument("--manter", action="store_true", help="mantém a pasta temporária")
    a = ap.parse_args()

    hash_atual, estado = verificar_integridade()
    print(f"testes_v8_1.py  SHA-256 {hash_atual}  [{estado}]")
    if a.autoteste:
        print("Autoteste do contrato e dos oráculos:")
        ok = autoteste()
        print("AUTOTESTE", "APROVADO" if ok else "REPROVADO")
        raise SystemExit(0 if ok else 1)

    sel = [x.strip().upper() for x in a.so.split(",") if x.strip()] if a.so else [n for n, _ in TESTES]
    conhecidos = {n for n, _ in TESTES}
    desconhecidos = set(sel) - conhecidos
    if desconhecidos:
        raise SystemExit("Erro: testes desconhecidos: " + ", ".join(sorted(desconhecidos)))
    for caminho, nome in ((a.script, "implementação v8"), (a.v7, "v7"), (a.v6, "v6"),
                          (a.v5, "v5"), (a.gerador, "gerador"), (a.testes_v7, "testes_v7.py")):
        if not os.path.isfile(caminho):
            raise SystemExit(f"Erro: {nome} não encontrado: {caminho}")

    work = tempfile.mkdtemp(prefix="testes_v8_")
    amb = Ambiente(a.script, a.gerador, a.v7, a.v6, a.v5, a.testes_v7, work)
    cache: dict[str, Any] = {}
    resultados: list[dict[str, Any]] = []
    inicio = time.perf_counter()
    try:
        for nome, fn in TESTES:
            if nome not in sel:
                continue
            t0 = time.perf_counter()
            try:
                detalhe, status = fn(amb, cache), "APROVADO"
            except AssertionError as exc:
                detalhe, status = str(exc), "REPROVADO"
            except Exception as exc:  # noqa: BLE001
                detalhe = f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=3)}"
                status = "ERRO"
            dt = time.perf_counter() - t0
            r = {"teste": nome, "status": status, "tempo_s": round(dt, 1),
                 "descricao": (fn.__doc__ or "").strip(), "detalhe": detalhe}
            resultados.append(r)
            print(f"[{status:^9}] {nome:<4} ({dt:6.1f} s) {r['descricao']}\n            {detalhe}")
        total = time.perf_counter() - inicio
        escrever_relatorios(resultados, estado, hash_atual, a.script, total)
        n_ok = sum(r["status"] == "APROVADO" for r in resultados)
        print(f"\nResultado: {n_ok}/{len(resultados)} aprovados em {total:.0f} s | integridade: {estado}")
        raise SystemExit(0 if n_ok == len(resultados) and estado == "ÍNTEGRO" else 1)
    finally:
        if a.manter:
            print(f"Pasta de trabalho: {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
testes_v9.py — suíte de aceitação CONGELADA da v9 do simulador VSG.

A v9 acrescenta à v8 um disjuntor de paralelismo, função 25 (sync-check),
antecipação do tempo mecânico de fechamento, pré-sincronizador ativo de fase e
tensão e retirada bumpless. Esta suíte é o contrato executável da T1 e deve ser
congelada antes da implementação de vsg_2a_ordem_degrau_carga_v9.py.

CONTRATO CONGELADO
==================
1. Compatibilidade
   - Recursos v9 são inertes por padrão; casos legados começam com o disjuntor
     fechado e devem reproduzir a v8.
   - evento=fechamento exige modo=rede e disjuntor_inicial=aberto.
   - A fonte v8 de referência e testes_v8_1.py são verificadas por SHA-256.

2. Convenções do sync-check
   delta_v = |V_pcc|/Vbase_pico - |V_rede|/Vbase_pico                 [pu]
   delta_f = f_vsg - f_rede_estimada                                  [Hz]
   delta_theta = wrap(angle(V_pcc)-angle(V_rede)) em [-pi, pi)         [rad]
   delta_theta_pred = wrap(delta_theta + 2*pi*delta_f*t_fechamento)    [rad]

   A medição é válida somente quando ambos os lados têm magnitude >=
   vmin_medicao_pu. As fronteiras são inclusivas. A autorização exige, ao
   mesmo tempo, qualidade válida, |delta_v|<=0,05 pu, |delta_f|<=0,10 Hz,
   |delta_theta|<=5 graus E |delta_theta_pred|<=5 graus. A janela atual é
   deliberadamente exigida junto com a prevista.

3. Supervisor e disjuntor
   - Amostragem fixa dt_rele=1 ms; hold contínuo de 0,10 s.
   - Qualquer amostra inválida ou fora da janela zera o timer.
   - Comando emitido é latched; contato ocorre exatamente t_fechamento=0,06 s
     depois, mesmo se a janela se perder após o comando.
   - timeout efetivo = min(t_step+t_sync_timeout_s, t_end). Sem comando até o
     limite, o estado final é timeout e o disjuntor permanece aberto.
   - No modo forçado o comando ocorre em t_step, sem aprovação do sync-check.
   - I_g(0-)=I_g(0+)=0. Após contato:
       dI_g/dt=(V_c-V_g-R_g I_g)/L_g - j*omega0*I_g.
     Nenhum outro estado contínuo é projetado no contato.

4. Pré-sincronizador
   e_theta=-delta_theta; delta_f_bias=sat(Kp_theta*e_theta+Ki_theta*xi_theta)
   e_v=-delta_v;         delta_e_bias=sat(Kp_v*e_v+Ki_v*xi_v)
   O integrador usa anti-windup condicional: integra quando a saída não está
   saturada ou quando o erro conduz a saída de volta à faixa. A atuação ocorre
   somente em estrategia_sync=ativo e com o disjuntor aberto.
   Após o contato, cada bias mantém o valor de contato e é retirado por rampa
   linear contínua até zero em t_release_sync_s=0,20 s.
   O PLL é dedicado à manobra e é inerte em casos legados. Banda padrão: 5 Hz.

5. Estados publicados
   0 aberto; 1 aquisicao; 2 pre_sincronizando; 3 janela_valida;
   4 comando_emitido; 5 fechando; 6 fechado; 7 fechado_fora_da_janela;
   8 timeout; 9 bloqueado.

6. Tolerâncias principais
   - teto de I_f: Imax*(1+0,005), herdado da v8;
   - contato: max(dt_rele,dt_out);
   - equivalência legada: 1e-6 pu/Hz/rad;
   - convergência contínua: 0,5%; supervisor: uma amostra do relay;
   - desempenho: <=2,5 vezes o caso v8 equivalente.

API PURA EXIGIDA DA IMPLEMENTAÇÃO V9
====================================
Para isolar erros, a implementação deve expor:
  avaliar_sync_check(dv_pu, df_hz, dtheta_rad, medicao_valida,
      dv_max_pu, df_max_hz, dtheta_max_rad, t_fechamento_s,
      antecipar=True, exigir_janela_atual=True) -> dict
  passo_pi_sync(erro, integral, kp, ki, limite, dt) -> dict
  rampa_retirada_sync(valor_contato, t, t_contato, t_release) -> float
  derivada_corrente_rede(vc_complex_v, vg_complex_v, ig_complex_a,
      rg_ohm, lg_h, omega0_rad_s) -> complex

Os dicionários das duas primeiras funções contêm, respectivamente, ao menos
{ok, dtheta_pred_rad} e {saida, integral, saturado}.

USO
===
  python testes_v9.py --autoteste
  python testes_v9.py
  python testes_v9.py --so R1,D1,S1,P1
  python testes_v9.py --script vsg_2a_ordem_degrau_carga_v9.py

O autoteste valida o contrato e os oráculos sem importar/executar a v9. A suíte
completa contém 41 testes (R1-R5, D1-D6, S1-S8, P1-P9, F1-F6 e E1-E7), gera
relatorio_testes_v9.json/.md e só retorna zero quando todos os testes
selecionados passam e o hash publicado está íntegro.
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

VERSAO_TESTES = "v9.0.0-contrato-t1"
HASH_TESTES_V8_1 = "f09ad750167d77e95ebe9215147cf0fb2322917941cbb4abb9e917aabb5d7f5d"
HASH_SCRIPT_V8 = "11cbf6f9861ed240a079d90ce07a428fbafed49280de85751d900a421dfc53e3"

CONTRATO = {
    "dv_sync_max_pu": 0.05,
    "df_sync_max_hz": 0.10,
    "dtheta_sync_max_graus": 5.0,
    "t_sync_hold_s": 0.10,
    "t_fechamento_s": 0.06,
    "t_sync_timeout_s": 10.0,
    "dt_rele_s": 0.001,
    "vmin_medicao_pu": 0.20,
    "df_sync_lim_hz": 0.50,
    "de_sync_lim_pu": 0.10,
    "t_release_sync_s": 0.20,
    "pll_sync_bw_hz": 5.0,
    "exigir_janela_atual": True,
    "timeout_limitado_por_t_end": True,
    "retirada": "rampa_linear",
    "antiwindup_pi": "integracao_condicional",
}
TOL = {
    "serie_pu": 1e-6,
    "serie_hz": 1e-6,
    "serie_rad": 1e-6,
    "modal_rel": 1e-5,
    "formula_abs": 2e-12,
    "teto_corrente_frac": 0.005,
    "saida_rel": 5e-4,
    "contato_estado": 0.05,
    "contato_corrente_ideal_pu": 0.05,
    "contato_tempo_extra_s": 1e-12,
    "convergencia_frac": 0.005,
    "supervisor_amostras": 1.0,
    "pll_f_hz": 0.02,
    "bumpless_bias": 2e-3,
    "bias_final": 1e-6,
    "desempenho_razao": 2.5,
    "sync_f_mantido_hz": 0.05,
    "sync_f_perdido_hz": 0.20,
    "sync_delta_mantido_rad": 0.10,
    "sync_delta_perdido_rad": 0.50,
    "sync_janela_s": 0.50,
}
ESTADOS_SYNC = {
    0: "aberto", 1: "aquisicao", 2: "pre_sincronizando", 3: "janela_valida",
    4: "comando_emitido", 5: "fechando", 6: "fechado",
    7: "fechado_fora_da_janela", 8: "timeout", 9: "bloqueado",
}
SERIES_V9 = (
    "disjuntor_fechado", "comando_fechamento", "sync_check_ok",
    "sync_estado_codigo", "sync_tempo_janela_s", "delta_sync_rad",
    "delta_sync_pred_rad", "df_sync_hz", "dv_sync_pu",
    "medicao_sync_valida", "pre_sync_ativo", "df_bias_sync_hz",
    "de_bias_sync_pu", "theta_pll_rede_rad", "f_pll_rede_hz", "ig_env_pu",
)
CAMPOS_SYNC = (
    "solicitada", "estrategia", "estado_final", "disjuntor_inicial",
    "disjuntor_final", "sync_check_habilitado", "antecipacao_habilitada",
    "tempo_inicio_s", "tempo_comando_s", "tempo_contato_s",
    "tempo_sincronizacao_s", "timeout", "motivo_bloqueio", "janelas",
    "no_comando", "no_contato", "transitorio", "mapa_estados",
)
SERIES_COMUNS = {
    "Pf_pu": "serie_pu", "P_pu": "serie_pu", "P_rede_pu": "serie_pu",
    "Q_rede_pu": "serie_pu", "f_hz": "serie_hz", "delta_v_rad": "serie_rad",
    "if_env_pu": "serie_pu", "rv_pu": "serie_pu", "xv_pu": "serie_pu",
    "delta_aw_rad_s": "serie_rad", "vg_aplicada_pu": "serie_pu",
}
CLI_V9 = (
    "--delta-g0-graus", "--df-g0-hz", "--disjuntor-inicial", "--estrategia-sync",
    "--dv-sync-max-pu", "--df-sync-max-hz", "--dtheta-sync-max-graus",
    "--t-sync-hold-s", "--t-fechamento-s", "--t-sync-timeout-s", "--dt-rele-s",
    "--antecipar-fechamento", "--vmin-medicao-pu", "--kp-theta-hz-rad",
    "--ki-theta-hz-rad-s", "--df-sync-lim-hz", "--kp-v", "--ki-v-s",
    "--de-sync-lim-pu", "--t-release-sync-s", "--pll-sync-bw-hz",
)


def rel(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-15)


def checar(cond: Any, msg: str) -> None:
    if not bool(cond):
        raise AssertionError(msg)


def wrap(ang: Any) -> Any:
    a = (np.asarray(ang) + np.pi) % (2 * np.pi) - np.pi
    return float(a) if np.ndim(a) == 0 else a


def sha256(caminho: str) -> str:
    h = hashlib.sha256()
    with open(caminho, "rb") as fh:
        for bloco in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(bloco)
    return h.hexdigest()


def predizer_delta_contato(delta_rad: float, df_hz: float, t_fechamento_s: float,
                           antecipar: bool = True) -> float:
    if not antecipar:
        return float(wrap(delta_rad))
    return float(wrap(delta_rad + 2 * np.pi * df_hz * t_fechamento_s))


def avaliar_sync_oraculo(dv_pu: float, df_hz: float, dtheta_rad: float,
                         medicao_valida: bool, dv_max_pu: float = CONTRATO["dv_sync_max_pu"],
                         df_max_hz: float = CONTRATO["df_sync_max_hz"],
                         dtheta_max_rad: float = math.radians(CONTRATO["dtheta_sync_max_graus"]),
                         t_fechamento_s: float = CONTRATO["t_fechamento_s"],
                         antecipar: bool = True,
                         exigir_janela_atual: bool = True) -> dict[str, Any]:
    pred = predizer_delta_contato(dtheta_rad, df_hz, t_fechamento_s, antecipar)
    atual = abs(wrap(dtheta_rad)) <= dtheta_max_rad + TOL["formula_abs"]
    previsto = abs(pred) <= dtheta_max_rad + TOL["formula_abs"]
    ok = (bool(medicao_valida) and abs(dv_pu) <= dv_max_pu + TOL["formula_abs"] and
          abs(df_hz) <= df_max_hz + TOL["formula_abs"] and previsto and
          (atual or not exigir_janela_atual))
    return {"ok": bool(ok), "dtheta_pred_rad": pred, "janela_atual": bool(atual),
            "janela_prevista": bool(previsto)}


def passo_pi_oraculo(erro: float, integral: float, kp: float, ki: float,
                      limite: float, dt: float) -> dict[str, Any]:
    if min(kp, ki, limite, dt) < 0 or limite == 0 or dt == 0:
        raise ValueError("parâmetros inválidos do PI")
    bruto = kp * erro + ki * integral
    saida = float(np.clip(bruto, -limite, limite))
    saturado = abs(bruto - saida) > TOL["formula_abs"]
    integrar = (not saturado or (saida >= limite and erro < 0) or
                (saida <= -limite and erro > 0))
    integ_novo = integral + dt * erro if integrar else integral
    bruto_novo = kp * erro + ki * integ_novo
    saida_nova = float(np.clip(bruto_novo, -limite, limite))
    return {"saida": saida_nova, "integral": float(integ_novo),
            "saturado": bool(abs(bruto_novo - saida_nova) > TOL["formula_abs"])}


def rampa_retirada_oraculo(valor_contato: float, t: float, t_contato: float,
                            t_release: float) -> float:
    if t_release <= 0:
        raise ValueError("t_release deve ser positivo")
    if t <= t_contato:
        return float(valor_contato)
    ganho = float(np.clip(1.0 - (t - t_contato) / t_release, 0.0, 1.0))
    return float(valor_contato * ganho)


def derivada_ig_oraculo(vc: complex, vg: complex, ig: complex, rg: float,
                        lg: float, omega0: float) -> complex:
    if lg <= 0:
        raise ValueError("L_g deve ser positivo")
    return (complex(vc) - complex(vg) - rg * complex(ig)) / lg - 1j * omega0 * complex(ig)


def contar_pole_slips(delta_unwrapped: np.ndarray) -> int:
    bandas = np.floor((np.asarray(delta_unwrapped) + np.pi) / (2 * np.pi)).astype(np.int64)
    return int(np.sum(np.abs(np.diff(bandas)))) if len(bandas) > 1 else 0


def classificar_pos_contato(t: np.ndarray, delta: np.ndarray, f_diff: np.ndarray,
                            t_contato: float) -> tuple[str, int]:
    t, delta, f_diff = map(np.asarray, (t, delta, f_diff))
    pos = t >= t_contato
    if not np.any(pos):
        return "indeterminado", 0
    d, f, tp = delta[pos], f_diff[pos], t[pos]
    slips = contar_pole_slips(d)
    m = tp >= tp[-1] - TOL["sync_janela_s"]
    fmax, dpp = float(np.max(np.abs(f[m]))), float(np.ptp(d[m]))
    if slips >= 1 or fmax > TOL["sync_f_perdido_hz"] or dpp > TOL["sync_delta_perdido_rad"]:
        return "perdido", slips
    if fmax <= TOL["sync_f_mantido_hz"] and dpp <= TOL["sync_delta_mantido_rad"]:
        return "mantido", slips
    return "indeterminado", slips


@dataclass
class EstadoRele:
    estado: str = "aberto"
    tempo_janela_s: float = 0.0
    comando: bool = False
    tempo_comando_s: float | None = None
    tempo_contato_s: float | None = None
    fechado: bool = False
    fora_da_janela: bool = False
    timeout: bool = False


def passo_rele_oraculo(st: EstadoRele, t: float, dv: float, df: float, dtheta: float,
                        valida: bool, estrategia: str = "passivo", t_step: float = 0.5,
                        t_end: float = 5.0, dt: float = CONTRATO["dt_rele_s"],
                        hold: float = CONTRATO["t_sync_hold_s"],
                        fechamento: float = CONTRATO["t_fechamento_s"],
                        timeout: float = CONTRATO["t_sync_timeout_s"],
                        antecipar: bool = True) -> EstadoRele:
    if st.fechado or st.timeout:
        return st
    if t < t_step:
        st.estado = "aberto"
        return st
    if st.comando:
        assert st.tempo_contato_s is not None
        if t + TOL["formula_abs"] >= st.tempo_contato_s:
            agora = avaliar_sync_oraculo(dv, df, dtheta, valida, antecipar=antecipar)
            st.fechado, st.fora_da_janela = True, not agora["ok"]
            st.estado = "fechado_fora_da_janela" if st.fora_da_janela else "fechado"
        else:
            st.estado = "fechando"
        return st
    prazo = min(t_step + timeout, t_end)
    if t + TOL["formula_abs"] >= prazo:
        st.estado, st.timeout = "timeout", True
        return st
    if estrategia == "forcado":
        st.comando, st.tempo_comando_s = True, float(t_step)
        st.tempo_contato_s = float(t_step + fechamento)
        st.estado = "comando_emitido"
        return st
    r = avaliar_sync_oraculo(dv, df, dtheta, valida, t_fechamento_s=fechamento,
                             antecipar=antecipar)
    if not valida:
        st.tempo_janela_s, st.estado = 0.0, "bloqueado"
    elif r["ok"]:
        st.tempo_janela_s += dt
        st.estado = "janela_valida"
        if st.tempo_janela_s + TOL["formula_abs"] >= hold:
            st.comando, st.tempo_comando_s = True, float(t)
            st.tempo_contato_s = float(t + fechamento)
            st.estado = "comando_emitido"
    else:
        st.tempo_janela_s = 0.0
        st.estado = "pre_sincronizando" if estrategia == "ativo" else "aquisicao"
    return st


@dataclass
class ResultadoReduzido:
    fechado: bool
    timeout: bool
    tempo_comando: float | None
    tempo_contato: float | None
    delta_contato: float | None
    dv_contato: float | None
    df_contato: float | None
    pico_bias_f: float
    pico_bias_v: float


def oraculo_reduzido(estrategia: str, delta0: float, df0: float, dv0: float,
                      t_end: float = 8.0, limite_f: float = 0.5,
                      limite_v: float = 0.10) -> ResultadoReduzido:
    dt, t_step = CONTRATO["dt_rele_s"], 0.5
    st, delta, dv = EstadoRele(), float(delta0), float(dv0)
    it, iv, bf, bv = 0.0, 0.0, 0.0, 0.0
    pbf, pbv = 0.0, 0.0
    contato = (None, None, None)
    for t in np.arange(0.0, t_end + dt/2, dt):
        if estrategia == "ativo" and t >= t_step and not st.comando:
            qf = passo_pi_oraculo(-delta, it, 0.60, 0.20, limite_f, dt)
            qv = passo_pi_oraculo(-dv, iv, 0.80, 0.30, limite_v, dt)
            bf, it, bv, iv = qf["saida"], qf["integral"], qv["saida"], qv["integral"]
        df = df0 + bf
        st = passo_rele_oraculo(st, float(t), dv, df, delta, True, estrategia,
                                t_step=t_step, t_end=t_end, timeout=t_end-t_step)
        pbf, pbv = max(pbf, abs(bf)), max(pbv, abs(bv))
        if st.fechado:
            contato = (delta, dv, df)
            break
        delta = float(wrap(delta + 2*np.pi*df*dt))
        dv += dt * (bv - dv) / 0.20
    return ResultadoReduzido(st.fechado, st.timeout, st.tempo_comando_s,
                             st.tempo_contato_s, contato[0], contato[1], contato[2], pbf, pbv)


class Ambiente:
    def __init__(self, script: str, v8: str, testes_v8: str, gerador: str,
                 v7: str, v6: str, v5: str, testes_v7: str, workdir: str):
        self.script, self.v8, self.testes_v8 = map(os.path.abspath, (script, v8, testes_v8))
        self.gerador = os.path.abspath(gerador)
        self.v7, self.v6, self.v5, self.testes_v7 = map(os.path.abspath, (v7, v6, v5, testes_v7))
        self.dir, self.cont, self.stdout, self._modulo = workdir, 0, "", None
        self.base_json = os.path.join(workdir, "base.json")
        with open(self.base_json, "w", encoding="utf-8") as fh:
            json.dump({"sistema": {"sn": 100000, "vll": 380, "f0": 60}}, fh)

    def _pref(self, nome: str) -> str:
        self.cont += 1
        return os.path.join(self.dir, f"{self.cont:03d}_{nome}")

    def sim(self, nome: str, args: Iterable[Any], script: str | None = None,
            config: str | None = None, graficos: bool = False, ler: bool = True,
            esperar_erro: bool = False) -> tuple[Any, Any, float] | subprocess.CompletedProcess[str]:
        pref = self._pref(nome)
        cmd = [sys.executable, script or self.script, "--config", config or self.base_json,
               *[str(x) for x in args], "--prefixo", pref]
        if not graficos:
            cmd.append("--sem-graficos")
        t0 = time.perf_counter()
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=2400,
                            encoding="utf-8", errors="replace")
        dt = time.perf_counter() - t0
        self.stdout = (pr.stdout or "") + (pr.stderr or "")
        if esperar_erro:
            return pr
        if pr.returncode:
            raise AssertionError(f"simulador falhou ({pr.returncode}):\n{self.stdout[-2400:]}")
        if not ler:
            return None, None, dt
        with open(pref + "_resultados.json", encoding="utf-8") as fh:
            js = json.load(fh)
        with np.load(pref + "_series.npz") as arq:
            ser = {k: np.asarray(arq[k]) for k in arq.files}
        return js, ser, dt

    def gerar(self, nome: str, args: Iterable[Any]) -> tuple[subprocess.CompletedProcess[str], str]:
        saida = self._pref(nome) + ".yaml"
        cmd = [sys.executable, self.gerador, *[str(x) for x in args], "-o", saida]
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1200,
                            encoding="utf-8", errors="replace")
        return pr, saida

    def modulo_v9(self):
        if self._modulo is None:
            spec = importlib.util.spec_from_file_location("vsg_v9_sob_teste", self.script)
            checar(spec is not None and spec.loader is not None, "não foi possível importar a v9")
            mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
            self._modulo = mod
        return self._modulo


def rede(scr: float = 5.0, dw: float = 100.0) -> list[Any]:
    return ["--modo", "rede", "--scr", scr, "--dw", dw, "--tw", 1.0]


def curto(extra: Iterable[Any] = ()) -> list[Any]:
    return ["--evento", "nenhum", "--t-step", .2, "--t-end", .4, "--dt-out", 1e-4, *extra]


def fechamento(estrategia: str = "ativo", delta_g: float = 20.0, df_g: float = .12,
               t_end: float = 8.0, scr: float = 5.0, extra: Iterable[Any] = ()) -> list[Any]:
    args = [*rede(scr), "--evento", "fechamento", "--disjuntor-inicial", "aberto",
            "--estrategia-sync", estrategia, "--delta-g0-graus", delta_g,
            "--df-g0-hz", df_g, "--t-step", .5, "--t-end", t_end,
            "--dt-out", 2e-4, "--dv-sync-max-pu", CONTRATO["dv_sync_max_pu"],
            "--df-sync-max-hz", CONTRATO["df_sync_max_hz"],
            "--dtheta-sync-max-graus", CONTRATO["dtheta_sync_max_graus"],
            "--t-sync-hold-s", CONTRATO["t_sync_hold_s"],
            "--t-fechamento-s", CONTRATO["t_fechamento_s"],
            "--t-sync-timeout-s", min(6.5, t_end-.5), "--dt-rele-s", CONTRATO["dt_rele_s"],
            "--vmin-medicao-pu", CONTRATO["vmin_medicao_pu"], *extra]
    if estrategia == "ativo":
        args += ["--kp-theta-hz-rad", .60, "--ki-theta-hz-rad-s", .20,
                 "--df-sync-lim-hz", .50, "--kp-v", .80, "--ki-v-s", .30,
                 "--de-sync-lim-pu", .10, "--t-release-sync-s", .20,
                 "--pll-sync-bw-hz", 5.0]
    return args


LIMITADOR = ["--imax-pu", 1.2, "--i-on-pu", 1.176, "--rv-max-pu", 3.0,
             "--xv-rv", 0.0, "--k-aw", 20.0]


def comparar_series(ref: dict[str, np.ndarray], novo: dict[str, np.ndarray],
                     nomes: dict[str, str] = SERIES_COMUNS) -> float:
    checar(len(ref["t"]) == len(novo["t"]) and np.allclose(ref["t"], novo["t"], atol=1e-12, rtol=0),
           "vetores de tempo diferentes")
    pior = 0.0
    for k, tk in nomes.items():
        checar(k in ref and k in novo, f"série comum ausente: {k}")
        d = float(np.nanmax(np.abs(ref[k] - novo[k])))
        pior = max(pior, d/TOL[tk]); checar(d <= TOL[tk], f"{k}: {d:.3e}>{TOL[tk]:.1e}")
    return pior


def bloco_sync(js: dict[str, Any]) -> dict[str, Any]:
    b = js.get("sincronizacao")
    checar(isinstance(b, dict), "bloco JSON 'sincronizacao' ausente")
    for k in CAMPOS_SYNC:
        checar(k in b, f"sincronizacao sem '{k}'")
    checar(b["estado_final"] in ESTADOS_SYNC.values(), f"estado final inválido: {b['estado_final']}")
    mapa = {int(k): v for k, v in b["mapa_estados"].items()}
    checar(mapa == ESTADOS_SYNC, f"mapa de estados divergente: {mapa}")
    return b


def checar_series_v9(s: dict[str, np.ndarray]) -> None:
    n = len(s["t"])
    for k in SERIES_V9:
        checar(k in s, f"série v9 ausente: {k}")
        checar(len(s[k]) == n, f"{k}: tamanho {len(s[k])}, esperado {n}")
        checar(np.all(np.isfinite(s[k])), f"{k} contém NaN/Inf")
    cod = np.unique(s["sync_estado_codigo"].astype(int))
    checar(all(int(x) in ESTADOS_SYNC for x in cod), f"código de estado desconhecido: {cod}")


def indice_tempo(t: np.ndarray, alvo: float) -> int:
    return int(np.argmin(np.abs(np.asarray(t)-alvo)))


def tempo_borda(s: dict[str, np.ndarray], chave: str) -> float | None:
    a = np.asarray(s[chave]) > .5
    k = np.flatnonzero(a)
    return float(s["t"][k[0]]) if len(k) else None


def exigir_funcao(mod: Any, nome: str):
    fn = getattr(mod, nome, None); checar(callable(fn), f"função pública {nome} ausente")
    return fn


def checar_teto_v8(js: dict[str, Any], s: dict[str, np.ndarray], exigir_ativo: bool = False) -> str:
    lim = js.get("limitador") or {}
    checar("I_pico_env_pu" in lim and "teto_rv_atingido" in lim, "bloco limitador v8 incompleto")
    pico = float(np.max(s["if_env_pu"])); imax = float(lim.get("imax_pu", 0))
    checar(abs(float(lim["I_pico_env_pu"])-pico) <= TOL["saida_rel"]*max(1.0, pico),
           "pico de corrente JSON/NPZ divergente")
    if exigir_ativo:
        checar(np.any(s["limitador_ativo"] > .5), "limitador deveria atuar")
    if imax > 0 and not bool(lim["teto_rv_atingido"]):
        checar(pico <= imax*(1+TOL["teto_corrente_frac"]), f"I_f={pico:.4f} pu excede Imax")
    return f"I_f,pico={pico:.4f} pu"


def autovalores(js: dict[str, Any]) -> np.ndarray:
    modal = js.get("modal_fechado", js.get("modal", {}))
    return np.array([complex(r, i) for r, i in modal["autovalores"] if abs(complex(r, i)) >= 1e-6])


# ======================================================================================
# Regressão R1-R5
# ======================================================================================
def R1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Baseline v8 integral: fonte íntegra e R1-R5/C1-C14 passam contra a v9."""
    checar(sha256(amb.testes_v8) == HASH_TESTES_V8_1, "testes_v8_1.py não corresponde à revisão congelada")
    checar(sha256(amb.v8) == HASH_SCRIPT_V8, "implementação v8 não corresponde à baseline congelada")
    txt = open(amb.testes_v8, encoding="utf-8").read()
    txt = txt.replace('js.get("versao") == "v8"', 'js.get("versao") == "v9"')
    txt = txt.replace("esperado 'v8'", "esperado 'v9'")
    txt = txt.replace("campo 'versao' diferente de 'v8'", "campo 'versao' diferente de 'v9'")
    txt = txt.replace("'js.get(\"versao\") == \"v8\"'", "'js.get(\"versao\") == \"v9\"'")
    adapt = os.path.join(amb.dir, "testes_v8_contrato_v9.py")
    with open(adapt, "w", encoding="utf-8") as fh:
        fh.write(txt)
    with open(os.path.join(amb.dir, "testes_v8_1.sha256"), "w", encoding="ascii") as fh:
        fh.write(sha256(adapt) + "  testes_v8_contrato_v9.py\n")
    cmd = [sys.executable, adapt, "--script", amb.script, "--gerador", amb.gerador,
           "--v7", amb.v7, "--v6", amb.v6, "--v5", amb.v5, "--testes-v7", amb.testes_v7]
    pr = subprocess.run(cmd, cwd=amb.dir, capture_output=True, text=True, timeout=7200,
                        encoding="utf-8", errors="replace")
    saida = (pr.stdout or "") + (pr.stderr or "")
    checar(pr.returncode == 0, "baseline v8 não passou integralmente contra a v9:\n" + saida[-4000:])
    checar(re.search(r"Resultado:\s*19/19 aprovados", saida) is not None,
           "execução adaptada não comprovou 19/19 testes")
    return "baseline v8 íntegra; 19/19 contratos aprovados contra a v9"


R2_CASOS = [
    ("ilha", ["--modo", "ilhado", "--evento", "carga", "--t-step", .5, "--t-end", 2, "--dt-out", 2e-4]),
    ("pref_scr2", [*rede(2), "--evento", "pref", "--t-step", .5, "--t-end", 2, "--dt-out", 2e-4]),
    ("fase_scr5", [*rede(5), "--evento", "fase", "--d-fase", 5, "--t-step", .5, "--t-end", 2, "--dt-out", 1e-4]),
    ("freq", [*rede(5, 103), "--evento", "freq_rampa", "--df-g", -.2, "--rocof-g", .2,
              "--t-step", .5, "--t-end", 2, "--dt-out", 2e-4]),
    ("falta", [*rede(20), "--evento", "falta_3f", "--vg-falta", 0, "--t-step", .5,
               "--t-clear", .6, "--t-end", 2, "--dt-out", 1e-4, *LIMITADOR]),
]


def R2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Equivalência numérica v8×v9 em modos, eventos, SCR e limitador."""
    out = []
    for nome, args in R2_CASOS:
        _, s8, _ = amb.sim("R2_"+nome+"_v8", args, script=amb.v8)
        j9, s9, _ = amb.sim("R2_"+nome+"_v9", args)
        checar(j9.get("versao") == "v9", f"{nome}: versão não é v9")
        out.append(f"{nome}:{comparar_series(s8, s9):.2f}xtol")
    return "; ".join(out)


def R3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """CLI e contratos JSON/NPZ da v8 permanecem e a interface v9 está completa."""
    hp = subprocess.run([sys.executable, amb.script, "--help"], capture_output=True, text=True,
                        timeout=120, encoding="utf-8", errors="replace")
    checar(hp.returncode == 0, "--help falhou")
    for op in CLI_V9:
        checar(op in hp.stdout, f"opção v9 ausente no --help: {op}")
    j8, s8, _ = amb.sim("R3_v8", [*rede(), *curto()], script=amb.v8)
    j9, s9, _ = amb.sim("R3_v9", [*rede(), *curto()])
    for k in j8:
        if k != "versao": checar(k in j9, f"campo JSON v8 removido: {k}")
    for k in s8: checar(k in s9, f"série NPZ v8 removida: {k}")
    checar_series_v9(s9); bloco_sync(j9)
    return f"{len(j8)} campos e {len(s8)} séries v8 preservados; {len(SERIES_V9)} séries v9 presentes"


def R4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Topologia legada fechada mantém os autovalores da v8."""
    out = []
    for scr in (2, 5, 20):
        j8, _, _ = amb.sim(f"R4_{scr}_v8", [*rede(scr), *curto()], script=amb.v8)
        j9, _, _ = amb.sim(f"R4_{scr}_v9", [*rede(scr), *curto()])
        a, b = autovalores(j8), autovalores(j9)
        checar(len(a) == len(b), f"SCR {scr}: número de autovalores mudou")
        pior = 0.0
        for x in a:
            d = float(np.min(np.abs(b-x))/max(1.0, abs(x)))
            pior = max(pior, d); checar(d <= TOL["modal_rel"], f"SCR {scr}: modo {x} ausente")
        out.append(f"SCR{scr}:{pior:.1e}")
    return ", ".join(out)


def R5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Portabilidade: caminho com espaços/Unicode, UTF-8 e NumPy corrente."""
    pasta = os.path.join(amb.dir, "caminho com espaços çã")
    os.makedirs(pasta); copia = os.path.join(pasta, "simulador v9.py")
    shutil.copy2(amb.script, copia)
    pref = os.path.join(pasta, "saída teste")
    cmd = [sys.executable, copia, "--config", amb.base_json, *map(str, [*rede(), *curto()]),
           "--prefixo", pref, "--sem-graficos"]
    pr = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                        encoding="utf-8", errors="strict")
    checar(pr.returncode == 0, "execução portátil falhou:\n" + ((pr.stdout or "")+(pr.stderr or ""))[-1500:])
    checar(os.path.isfile(pref+"_resultados.json") and os.path.isfile(pref+"_series.npz"),
           "arquivos não foram gravados em caminho Unicode")
    with open(pref+"_resultados.json", encoding="utf-8") as fh: json.load(fh)
    return f"Python {sys.version_info.major}.{sys.version_info.minor}; NumPy {np.__version__}; UTF-8 aprovado"


# ======================================================================================
# Disjuntor e topologia D1-D6
# ======================================================================================
def D1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Isolação aberta: I_g, P_rede e Q_rede exatamente zero antes do contato."""
    js, s, _ = amb.sim("D1", fechamento("forcado", 20, 0, 2, extra=LIMITADOR))
    b = bloco_sync(js); checar_series_v9(s)
    tc = float(b["tempo_contato_s"]); m = s["t"] < tc-TOL["formula_abs"]
    for k in ("ig_env_pu", "P_rede_pu", "Q_rede_pu", "disjuntor_fechado"):
        checar(np.all(s[k][m] == 0), f"{k} não é exatamente zero com disjuntor aberto")
    cache["D1"] = (js, s)
    return f"isolação exata em {int(np.sum(m))} amostras; contato={tc:.6f} s"


def D2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Equilíbrio inicial aberto coincide com o caso ilhado equivalente da v8."""
    j9, s9, _ = amb.sim("D2_v9", fechamento("forcado", 0, 0, 1.2))
    _, s8, _ = amb.sim("D2_v8", ["--modo", "ilhado", "--evento", "nenhum", "--t-step", .2,
                                       "--t-end", .4, "--dt-out", 2e-4], script=amb.v8)
    m = s9["t"] <= .4+1e-12
    checar(len(s8["t"]) == int(np.sum(m)), "grades pré-contato incompatíveis")
    pior = 0.0
    for k, tk in {"Pf_pu":"serie_pu", "P_pu":"serie_pu", "f_hz":"serie_hz",
                  "delta_v_rad":"serie_rad"}.items():
        d = float(np.max(np.abs(s8[k]-s9[k][m]))); pior=max(pior,d/TOL[tk])
        checar(d <= TOL[tk], f"equilíbrio ilhado divergente em {k}: {d}")
    return f"pré-fechamento igual à ilha v8; pior={pior:.2f}xtol"


def D3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fechamento ideal: estados contínuos, corrente pequena e limitador inativo."""
    js, s, _ = amb.sim("D3", fechamento("forcado", 0, 0, 2, extra=LIMITADOR))
    b = bloco_sync(js); tc=float(b["tempo_contato_s"]); k=indice_tempo(s["t"],tc)
    checar(bool(s["disjuntor_fechado"][k]), "disjuntor não fechou no contato")
    for q in ("delta_v_rad", "f_hz", "Vc_mag_pu"):
        if q in s and k>0:
            checar(abs(float(s[q][k]-s[q][k-1])) <= TOL["contato_estado"], f"salto artificial em {q}")
    checar(float(np.max(s["ig_env_pu"][k:k+20])) <= TOL["contato_corrente_ideal_pu"],
           "fechamento ideal gerou corrente excessiva")
    checar(not np.any(s["limitador_ativo"][k:k+20]>.5), "limitador atuou no fechamento ideal")
    cache["D3"]=(js,s)
    return f"contato ideal em {tc:.6f} s; Ig_pico_local={np.max(s['ig_env_pu'][k:k+20]):.4f} pu"


def D4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Derivada inicial do ramo após contato coincide com o oráculo Rg-Lg."""
    fn = exigir_funcao(amb.modulo_v9(), "derivada_corrente_rede")
    casos = [(230+0j, 230*np.exp(.2j), 0j), (210+20j, 225-10j, 3-2j)]
    pior=0.0
    for vc,vg,ig in casos:
        kw=dict(vc_complex_v=vc,vg_complex_v=vg,ig_complex_a=ig,rg_ohm=.02,lg_h=.001,
                omega0_rad_s=2*np.pi*60)
        got=complex(fn(**kw)); ref=derivada_ig_oraculo(vc,vg,ig,.02,.001,2*np.pi*60)
        d=abs(got-ref); pior=max(pior,d); checar(d<=1e-9*max(1,abs(ref)),f"derivada divergente: {got} != {ref}")
    return f"{len(casos)} estados; erro máximo={pior:.3e} A/s"


def D5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Contato não múltiplo de dt_out ocorre em breakpoint dedicado."""
    args=fechamento("forcado",0,0,2,extra=["--t-fechamento-s",.0637,"--dt-out",.002])
    js,s,_=amb.sim("D5",args); b=bloco_sync(js)
    esperado=.5637; tc=float(b["tempo_contato_s"])
    checar(abs(tc-esperado)<=TOL["formula_abs"],f"tempo JSON {tc} != {esperado}")
    checar(np.any(np.isclose(s["t"],esperado,atol=TOL["formula_abs"])),"breakpoint de contato ausente")
    return f"contato exato em {tc:.7f} s com dt_out=0,002 s"


def D6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Abertura/reclose e combinações de topologia não suportadas são rejeitadas."""
    casos={
        "fechamento_inicial_fechado":[*rede(),"--evento","fechamento","--disjuntor-inicial","fechado"],
        "aberto_em_carga":[*rede(),"--evento","carga","--disjuntor-inicial","aberto"],
        "fechamento_ilha":["--modo","ilhado","--evento","fechamento","--disjuntor-inicial","aberto"],
        "estrategia_invalida":[*rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","reclose"],
    }
    for nome,args in casos.items():
        pr=amb.sim("D6_"+nome,args,esperar_erro=True); txt=(pr.stdout or "")+(pr.stderr or "")
        checar(pr.returncode!=0 and "Traceback" not in txt and re.search(r"(^|\n)Erro",txt),f"{nome}: erro não limpo")
    return f"{len(casos)} configurações fora do escopo rejeitadas"


# ======================================================================================
# Sync-check S1-S8
# ======================================================================================
def S1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fronteiras das três janelas são inclusivas e pontos externos bloqueiam."""
    fn=exigir_funcao(amb.modulo_v9(),"avaliar_sync_check")
    dt=math.radians(CONTRATO["dtheta_sync_max_graus"])
    base=dict(medicao_valida=True,dv_max_pu=.05,df_max_hz=.10,dtheta_max_rad=dt,
              t_fechamento_s=.06,antecipar=True,exigir_janela_atual=True)
    casos=[(.05,0,0,True),(-.05,0,0,True),(0,.10,0,True),(0,-.10,0,True),
           (0,0,dt,True),(0,0,-dt,True),(.050001,0,0,False),(0,.100001,0,False),
           (0,0,dt+1e-6,False),(0,0,0,True)]
    for dv,df,dth,ok in casos:
        got=fn(dv_pu=dv,df_hz=df,dtheta_rad=dth,**base); ref=avaliar_sync_oraculo(dv,df,dth,**base)
        checar(bool(got["ok"])==ref["ok"]==ok,f"fronteira divergente: {(dv,df,dth)}")
        checar(abs(float(got["dtheta_pred_rad"])-ref["dtheta_pred_rad"])<=TOL["formula_abs"],"predição divergente")
    inval=fn(dv_pu=0,df_hz=0,dtheta_rad=0,**{**base,"medicao_valida":False})
    checar(not inval["ok"],"medição inválida autorizou")
    return f"{len(casos)+1} pontos da tabela verdade aprovados"


def S2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Timer exige permanência contínua e zera ao sair da janela."""
    dt=.001; st=EstadoRele(); seq=[]
    for n in range(900):
        t=n*dt
        dentro=(.5<=t<.56) or (.58<=t)
        st=passo_rele_oraculo(st,t,0,0,0 if dentro else .2,True,"passivo",t_end=2,dt=dt)
        seq.append(st.tempo_janela_s)
        if st.comando: break
    checar(abs(seq[int(.57/dt)])<=TOL["formula_abs"],"timer não zerou fora da janela")
    checar(st.comando and st.tempo_comando_s is not None and st.tempo_comando_s>=.679-1e-12,
           f"hold incompleto: comando={st.tempo_comando_s}")
    return f"timer zerado; comando após novo hold em {st.tempo_comando_s:.3f} s"


def S3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Tensão abaixo de vmin bloqueia medição e impede comando."""
    args=fechamento("passivo",0,0,2,extra=["--vg",.10,"--vmin-medicao-pu",.20,"--t-sync-timeout-s",1.0])
    js,s,_=amb.sim("S3",args); b=bloco_sync(js)
    checar(not np.any(s["comando_fechamento"]>.5),"houve comando com medição inválida")
    checar(np.any(s["medicao_sync_valida"]<.5),"qualidade inválida não foi publicada")
    checar(b["timeout"] and b["disjuntor_final"]=="aberto" and b["motivo_bloqueio"],"bloqueio/timeout não registrado")
    return f"estado={b['estado_final']}; motivo={b['motivo_bloqueio']}"


def S4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Antecipação angular coincide com delta+2π·df·t_fechamento."""
    fn=exigir_funcao(amb.modulo_v9(),"avaliar_sync_check"); pior=0.0
    for d in np.linspace(-np.pi,np.pi,13):
        for f in (-.1,-.03,.03,.1):
            g=fn(dv_pu=0,df_hz=f,dtheta_rad=d,medicao_valida=True,dv_max_pu=.05,
                 df_max_hz=.1,dtheta_max_rad=math.radians(5),t_fechamento_s=.063,
                 antecipar=True,exigir_janela_atual=True)
            ref=predizer_delta_contato(d,f,.063); e=abs(wrap(float(g["dtheta_pred_rad"])-ref)); pior=max(pior,e)
            checar(e<=TOL["formula_abs"],"predição angular incorreta")
    return f"52 combinações; erro máximo={pior:.3e} rad"


def S5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """A/B sem antecipação registra erro de contato não melhor que com antecipação."""
    base=fechamento("passivo",-4,.08,5)
    ja,_,_=amb.sim("S5_on",[*base,"--antecipar-fechamento","true"])
    jb,_,_=amb.sim("S5_off",[*base,"--antecipar-fechamento","false"])
    a,b=bloco_sync(ja),bloco_sync(jb)
    ea=abs(float(a["no_contato"]["dtheta_graus"])); eb=abs(float(b["no_contato"]["dtheta_graus"]))
    checar(ea<=eb+.1,f"antecipação piorou erro: {ea:.3f}>{eb:.3f} graus")
    return f"|erro contato| antecipado={ea:.3f}°, sem={eb:.3f}°"


def S6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Comando é latched até o contato e não é cancelado por perda posterior da janela."""
    js,s,_=amb.sim("S6",fechamento("passivo",-4,.08,5)); b=bloco_sync(js)
    tc,tm=float(b["tempo_contato_s"]),float(b["tempo_comando_s"])
    m=(s["t"]>=tm)&(s["t"]<tc)
    checar(np.all(s["comando_fechamento"][m]>.5),"comando deixou de ficar latched")
    checar(b["disjuntor_final"]=="fechado","contato não ocorreu após comando")
    checar(abs(tc-tm-CONTRATO["t_fechamento_s"])<=max(CONTRATO["dt_rele_s"],2e-4)+1e-12,"atraso mecânico incorreto")
    return f"comando={tm:.4f} s; contato={tc:.4f} s; latched em {int(np.sum(m))} amostras"


def S7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Condições impossíveis levam a timeout com disjuntor aberto."""
    args=fechamento("passivo",90,0,2,extra=["--t-sync-timeout-s",1.0])
    js,s,_=amb.sim("S7",args); b=bloco_sync(js)
    checar(b["timeout"] and b["estado_final"]=="timeout" and b["disjuntor_final"]=="aberto","timeout incorreto")
    checar(not np.any(s["comando_fechamento"]>.5) and not np.any(s["disjuntor_fechado"]>.5),"comando/contato indevido")
    return f"timeout em {s['t'][-1]:.3f} s; disjuntor permaneceu aberto"


def S8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Sincronização passiva fecha na primeira janela elegível após hold e predição."""
    js,s,_=amb.sim("S8",fechamento("passivo",-12,.03,5)); b=bloco_sync(js)
    checar(b["disjuntor_final"]=="fechado" and not b["timeout"],"caso passivo não fechou")
    tm=float(b["tempo_comando_s"]); k=indice_tempo(s["t"],tm)
    checar(float(s["sync_tempo_janela_s"][k])+CONTRATO["dt_rele_s"]>=CONTRATO["t_sync_hold_s"],"comando antes do hold")
    checar(bool(b["no_contato"]["dentro_das_janelas"]),"contato passivo fora das janelas")
    cache["S8"]=(js,s)
    return f"comando={tm:.3f} s; contato={float(b['tempo_contato_s']):.3f} s"


# ======================================================================================
# Pré-sincronizador P1-P9
# ======================================================================================
def _pi_impl(amb: Ambiente, erro: float, integral: float, kp: float, ki: float,
             limite: float, dt: float) -> dict[str, Any]:
    fn=exigir_funcao(amb.modulo_v9(),"passo_pi_sync")
    return fn(erro=erro,integral=integral,kp=kp,ki=ki,limite=limite,dt=dt)


def P1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Canal angular tem sinal correto para erros positivos e negativos."""
    for delta in (-.4,-.1,.1,.4):
        g=_pi_impl(amb,-delta,0,.6,.2,.5,.001)
        checar(float(g["saida"])*delta<0,f"bias {g['saida']} não reduz delta={delta}")
    return "quatro sinais angulares reduzem o erro"


def P2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """PI angular satura e não acumula windup contra o limite."""
    it=0.0
    for _ in range(5000):
        g=_pi_impl(amb,2.0,it,.6,.2,.5,.001); it=float(g["integral"])
    checar(abs(float(g["saida"]))<=.5+TOL["formula_abs"] and bool(g["saturado"]),"PI angular não saturou")
    checar(abs(it)<=.01,"integral angular acumulou windup")
    g2=_pi_impl(amb,-.2,it,.6,.2,.5,.001)
    checar(float(g2["integral"])<it,"integrador não recupera ao inverter o erro")
    return f"saída={g['saida']:.3f} Hz; integral={it:.4g}"


def P3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Canal de tensão produz bias que reduz ΔV nos dois sentidos."""
    for dv in (-.15,-.03,.03,.15):
        g=_pi_impl(amb,-dv,0,.8,.3,.1,.001)
        checar(float(g["saida"])*dv<0,f"bias E {g['saida']} não reduz dv={dv}")
    return "quatro sinais de tensão reduzem o mismatch"


def P4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """PI de tensão limita a saída e evita windup; caso inalcançável dá timeout."""
    it=0.0
    for _ in range(3000):
        g=_pi_impl(amb,1.0,it,.8,.3,.1,.001); it=float(g["integral"])
    checar(abs(float(g["saida"]))<=.1+TOL["formula_abs"] and abs(it)<=.01,"anti-windup de tensão falhou")
    js,_,_=amb.sim("P4_timeout",fechamento("ativo",0,0,2,extra=["--vg",.75,"--de-sync-lim-pu",.01,
                                                                 "--t-sync-timeout-s",1.0]))
    b=bloco_sync(js); checar(b["timeout"] and b["disjuntor_final"]=="aberto","caso inalcançável não deu timeout")
    return f"bias limitado a {g['saida']:.3f} pu; timeout honesto"


def P5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Caso ativo nominal converge e fecha dentro das três janelas."""
    js,s,_=amb.sim("P5",fechamento("ativo",20,.12,8,extra=LIMITADOR)); b=bloco_sync(js)
    checar(b["disjuntor_final"]=="fechado" and not b["timeout"],"pré-sync ativo não fechou")
    nc=b["no_contato"]
    checar(bool(nc["dentro_das_janelas"]),"contato ativo fora das janelas")
    checar(abs(float(nc["dv_pu"]))<=.05+1e-6 and abs(float(nc["df_hz"]))<=.10+1e-6 and
           abs(float(nc["dtheta_graus"]))<=5+1e-4,"mismatch no contato excede contrato")
    cache["P5"]=(js,s)
    return f"sync={b['tempo_sincronizacao_s']:.3f} s; contato={b['tempo_contato_s']:.3f} s"


def P6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Ativo supera o passivo no caso nominal sem piorar a corrente de contato."""
    ja,sa,_=amb.sim("P6_ativo",fechamento("ativo",20,.12,8,extra=LIMITADOR))
    jp,sp,_=amb.sim("P6_passivo",fechamento("passivo",20,.12,8,extra=LIMITADOR))
    a,p=bloco_sync(ja),bloco_sync(jp)
    checar(a["disjuntor_final"]=="fechado","ativo não fechou")
    melhora=(p["timeout"] or float(a["tempo_sincronizacao_s"])<=float(p["tempo_sincronizacao_s"])+CONTRATO["dt_rele_s"])
    checar(melhora,"ativo não reduziu tempo nem evitou timeout")
    ia=float(a["transitorio"]["ig_pico_pu"])
    if not p["timeout"]:
        checar(ia<=float(p["transitorio"]["ig_pico_pu"])+.02,"ativo piorou corrente de contato")
    return f"ativo={a['tempo_sincronizacao_s']} s; passivo={'timeout' if p['timeout'] else p['tempo_sincronizacao_s']}"


def P7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Biases são mantidos no contato e retirados linearmente, sem degrau."""
    js,s=cache.get("P5",(None,None))
    if js is None: js,s,_=amb.sim("P7",fechamento("ativo",20,.12,8,extra=LIMITADOR))
    b=bloco_sync(js); tc=float(b["tempo_contato_s"]); tr=CONTRATO["t_release_sync_s"]
    k=indice_tempo(s["t"],tc)
    for nome in ("df_bias_sync_hz","de_bias_sync_pu"):
        checar(k>0 and abs(float(s[nome][k]-s[nome][k-1]))<=TOL["bumpless_bias"],f"degrau em {nome} no contato")
        m=s["t"]>=tc+tr+max(CONTRATO["dt_rele_s"],2e-4)
        checar(np.max(np.abs(s[nome][m]))<=TOL["bias_final"],f"{nome} não chegou a zero")
    fn=exigir_funcao(amb.modulo_v9(),"rampa_retirada_sync")
    for q in (tc,tc+tr/2,tc+tr,tc+2*tr):
        got=float(fn(valor_contato=.4,t=q,t_contato=tc,t_release=tr))
        ref=rampa_retirada_oraculo(.4,q,tc,tr)
        checar(abs(got-ref)<=TOL["formula_abs"],"rampa pública diverge do contrato")
    return f"retirada linear concluída em {tc+tr:.3f} s"


def P8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """PLL é inerte em caso legado fechado, independentemente dos ganhos."""
    args=[*rede(),*curto()]
    _,a,_=amb.sim("P8_a",[*args,"--pll-sync-bw-hz",1])
    _,b,_=amb.sim("P8_b",[*args,"--pll-sync-bw-hz",20])
    pior=comparar_series(a,b)
    checar(np.all(a["pre_sync_ativo"]==0) and np.all(b["pre_sync_ativo"]==0),"pré-sync ativo em caso legado")
    return f"séries legadas iguais; pior={pior:.2f}xtol"


def P9(amb: Ambiente, cache: dict[str, Any]) -> str:
    """PLL dedicado rastreia pequeno desvio de frequência da rede."""
    js,s,_=amb.sim("P9",fechamento("ativo",5,.08,5)); b=bloco_sync(js)
    tc=float(b["tempo_contato_s"]); m=(s["t"]>=.5)&(s["t"]<tc)
    checar(np.any(m),"sem janela pré-contato para avaliar PLL")
    alvo=60+.08
    fim=np.flatnonzero(m)[-max(1,int(.1/(s["t"][1]-s["t"][0]))):]
    erro=float(np.max(np.abs(s["f_pll_rede_hz"][fim]-alvo)))
    checar(erro<=TOL["pll_f_hz"],f"PLL errou {erro:.4f} Hz")
    return f"erro máximo final do PLL={erro:.4f} Hz"


# ======================================================================================
# Fechamento, limite e estabilidade F1-F6
# ======================================================================================
def F1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fechamento ativo correto respeita Imax, não satura Rv e não perde sincronismo."""
    js,s=cache.get("P5",(None,None))
    if js is None: js,s,_=amb.sim("F1",fechamento("ativo",20,.12,8,extra=LIMITADOR))
    b=bloco_sync(js); det=checar_teto_v8(js,s)
    tr=b["transitorio"]
    checar(not bool(js["limitador"]["teto_rv_atingido"]),"R_v,max atingido no caso nominal")
    checar(int(tr["pole_slips"])==0 and tr["sincronismo_pos_fechamento"]=="mantido","sincronismo pós-fechamento falhou")
    return det+"; fechamento correto; zero pole slip"


def F2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fechamento forçado fora de fase atua o limitador sem falsa aprovação."""
    js,s,_=amb.sim("F2",fechamento("forcado",90,0,5,extra=LIMITADOR)); b=bloco_sync(js)
    checar(b["estado_final"] in ("fechado_fora_da_janela","fechado"),"modo forçado não fechou")
    checar(not bool(b["sync_check_habilitado"]),"forçado não deve habilitar sync-check")
    checar(not bool(b["no_contato"]["dentro_das_janelas"]),"fechamento severo marcado dentro das janelas")
    checar_teto_v8(js,s,exigir_ativo=True)
    txt=amb.stdout.lower(); checar("forç" in txt and ("janela" in txt or "sincron" in txt),"faltou alerta de fechamento forçado")
    cache["F2"]=(js,s)
    return f"contato a {b['no_contato']['dtheta_graus']:.2f}°; limitador atuou"


def F3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """R_v,max insuficiente mantém violação de corrente explícita."""
    lim=["--imax-pu",1.05,"--i-on-pu",1.029,"--rv-max-pu",.01,"--xv-rv",0,"--k-aw",20]
    js,s,_=amb.sim("F3",fechamento("forcado",90,0,4,extra=lim)); l=js.get("limitador",{})
    checar(bool(l.get("teto_rv_atingido")) and float(l.get("violacao_max_pu",0))>0,"saturação/violação não registrada")
    avisos=[x for x in amb.stdout.splitlines() if x.strip().startswith("(!)")]
    checar(any("corrente" in x.lower() and "r_v,max" in x.lower() for x in avisos),"aviso corrente/R_v,max ausente")
    return f"violação={l['violacao_max_pu']:.4f} pu"


def F4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Caso ativo opera em SCR 2, 5 e 20 sem oscilação crescente."""
    out=[]
    for scr in (2,5,20):
        js,s,_=amb.sim(f"F4_{scr}",fechamento("ativo",15,.08,8,scr,extra=LIMITADOR)); b=bloco_sync(js)
        checar(b["disjuntor_final"]=="fechado" and not b["timeout"],f"SCR {scr}: não fechou")
        checar_teto_v8(js,s)
        tc=float(b["tempo_contato_s"]); m=s["t"]>=max(tc,s["t"][-1]-.5)
        checar(np.ptp(s["delta_rel_unwrapped_rad"][m])<=TOL["sync_delta_perdido_rad"],f"SCR {scr}: oscilação crescente")
        out.append(f"SCR{scr}:{b['tempo_sincronizacao_s']:.2f}s")
    return ", ".join(out)


def F5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fechamento severo preserva a possibilidade física de perda de sincronismo."""
    js,s,_=amb.sim("F5",fechamento("forcado",150,.15,8,extra=LIMITADOR)); b=bloco_sync(js)
    tc=float(b["tempo_contato_s"])
    cls,slips=classificar_pos_contato(s["t"],s["delta_rel_unwrapped_rad"],s["f_hz"]-s["f_rede_hz"],tc)
    tr=b["transitorio"]
    checar(cls=="perdido" and (slips>=1 or int(tr["pole_slips"])>=1),f"perda física mascarada: {cls}/{slips}")
    checar(tr["sincronismo_pos_fechamento"]==cls,"classificação JSON/NPZ divergente")
    return f"classificação={cls}; pole slips={slips}"


def F6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Após retirada do pré-sync, frequência/ângulo convergem e biases zeram."""
    js,s=cache.get("P5",(None,None))
    if js is None: js,s,_=amb.sim("F6",fechamento("ativo",20,.12,8,extra=LIMITADOR))
    b=bloco_sync(js); tc=float(b["tempo_contato_s"]); m=s["t"]>=s["t"][-1]-.5
    checar(np.max(np.abs(s["df_bias_sync_hz"][m]))<=TOL["bias_final"] and
           np.max(np.abs(s["de_bias_sync_pu"][m]))<=TOL["bias_final"],"bias residual")
    cls,slips=classificar_pos_contato(s["t"],s["delta_rel_unwrapped_rad"],s["f_hz"]-s["f_rede_hz"],tc)
    checar(cls=="mantido" and slips==0,"regime pós-fechamento não convergiu")
    return "biases zero; sincronismo mantido"


# ======================================================================================
# Engenharia e produto E1-E7
# ======================================================================================
def _metricas_sync(js: dict[str, Any], s: dict[str, np.ndarray]) -> tuple[float,float,float,str]:
    b=bloco_sync(js)
    return (float(b["tempo_comando_s"]),float(b["tempo_contato_s"]),
            float(np.max(s["if_env_pu"])),str(b["estado_final"]))


def E1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Refino de max_step/rtol mantém picos, tempos e classe."""
    base=fechamento("ativo",20,.12,8,extra=LIMITADOR)
    j1,s1,_=amb.sim("E1_base",base)
    j2,s2,_=amb.sim("E1_fino",[*base,"--max-step",2e-5,"--rtol",1e-9])
    a,b=_metricas_sync(j1,s1),_metricas_sync(j2,s2)
    tol_t=max(CONTRATO["dt_rele_s"],2e-4)
    checar(abs(a[0]-b[0])<=tol_t and abs(a[1]-b[1])<=tol_t,"tempos mudaram com refino")
    checar(rel(a[2],b[2])<=TOL["convergencia_frac"],"pico de corrente não convergiu")
    checar(a[3]==b[3],"classe mudou com refino")
    return f"Δtcmd={abs(a[0]-b[0]):.4g}s; ΔIpico={100*rel(a[2],b[2]):.3f}%"


def E2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Refino de dt_rele preserva a decisão e tempos dentro de uma amostra."""
    out=[]
    for dt in (.001,.0005,.0001):
        js,_,_=amb.sim(f"E2_{dt}",fechamento("passivo",-12,.03,5,extra=["--dt-rele-s",dt]))
        b=bloco_sync(js); out.append((dt,float(b["tempo_comando_s"]),b["estado_final"]))
    ref=out[-1]
    for dt,t,st in out:
        checar(st==ref[2] and abs(t-ref[1])<=dt+ref[0]+1e-12,f"dt_rele {dt}: decisão/tempo divergente")
    return "; ".join(f"{d:g}->{t:.4f}s" for d,t,_ in out)


def E3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Entradas inválidas da v9 falham com 'Erro' e sem traceback."""
    raiz=[*rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--t-step",.5,"--t-end",2]
    casos={
      "dv0":[*raiz,"--dv-sync-max-pu",0],"df0":[*raiz,"--df-sync-max-hz",0],
      "dtheta0":[*raiz,"--dtheta-sync-max-graus",0],"hold0":[*raiz,"--t-sync-hold-s",0],
      "close0":[*raiz,"--t-fechamento-s",0],"timeout0":[*raiz,"--t-sync-timeout-s",0],
      "dt>hold":[*raiz,"--dt-rele-s",.2,"--t-sync-hold-s",.1],
      "vmin0":[*raiz,"--vmin-medicao-pu",0],"kpneg":[*raiz,"--kp-theta-hz-rad",-1],
      "kineg":[*raiz,"--ki-theta-hz-rad-s",-1],"limf0":[*raiz,"--df-sync-lim-hz",0],
      "kpvneg":[*raiz,"--kp-v",-1],"kivneg":[*raiz,"--ki-v-s",-1],
      "lime0":[*raiz,"--de-sync-lim-pu",0],"release0":[*raiz,"--t-release-sync-s",0],
      "pll0":[*raiz,"--pll-sync-bw-hz",0],
      "ativo_sem_ganho":[*raiz,"--estrategia-sync","ativo","--kp-theta-hz-rad",0,
                          "--ki-theta-hz-rad-s",0,"--kp-v",0,"--ki-v-s",0],
    }
    for nome,args in casos.items():
        pr=amb.sim("E3_"+nome,args,esperar_erro=True); txt=(pr.stdout or "")+(pr.stderr or "")
        checar(pr.returncode!=0 and "Traceback" not in txt and re.search(r"(^|\n)Erro",txt),f"{nome}: falha não limpa")
    return f"{len(casos)} entradas inválidas rejeitadas"


def E4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """JSON e NPZ concordam em tempos, erros, picos e estado final."""
    js,s=cache.get("P5",(None,None))
    if js is None: js,s,_=amb.sim("E4",fechamento("ativo",20,.12,8,extra=LIMITADOR))
    b=bloco_sync(js); checar_series_v9(s)
    tm,tc=float(b["tempo_comando_s"]),float(b["tempo_contato_s"])
    km,kc=indice_tempo(s["t"],tm),indice_tempo(s["t"],tc)
    nc=b["no_comando"]; nt=b["no_contato"]
    pares=[(nc["dv_pu"],s["dv_sync_pu"][km]),(nc["df_hz"],s["df_sync_hz"][km]),
           (nc["dtheta_graus"],np.degrees(s["delta_sync_rad"][km])),
           (nc["dtheta_pred_graus"],np.degrees(s["delta_sync_pred_rad"][km])),
           (nt["dv_pu"],s["dv_sync_pu"][kc]),(nt["df_hz"],s["df_sync_hz"][kc]),
           (nt["dtheta_graus"],np.degrees(s["delta_sync_rad"][kc]))]
    for a,z in pares: checar(abs(float(a)-float(z))<=TOL["saida_rel"]*max(1,abs(float(z))),"JSON/NPZ divergente")
    checar(abs(float(b["transitorio"]["if_pico_pu"])-float(np.max(s["if_env_pu"])))<=TOL["saida_rel"],"pico If divergente")
    cod=int(round(float(s["sync_estado_codigo"][-1]))); checar(ESTADOS_SYNC[cod]==b["estado_final"],"estado final divergente")
    return "tempos, sete erros, pico e estado final consistentes"


def E5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Gerador cria e executa os seis exemplos v9 sem edição."""
    casos=[
      ("ideal",["--modo","rede","--evento","fechamento","--estrategia-sync","forcado","--delta-g0-graus",0]),
      ("passivo",["--modo","rede","--evento","fechamento","--estrategia-sync","passivo","--delta-g0-graus",-12,"--df-g0-hz",.03]),
      ("ativo",["--modo","rede","--evento","fechamento","--estrategia-sync","ativo","--delta-g0-graus",20,"--df-g0-hz",.12]),
      ("timeout",["--modo","rede","--evento","fechamento","--estrategia-sync","passivo","--delta-g0-graus",90]),
      ("forcado",["--modo","rede","--evento","fechamento","--estrategia-sync","forcado","--delta-g0-graus",90]),
      ("scr2",["--modo","rede","--scr",2,"--evento","fechamento","--estrategia-sync","ativo","--delta-g0-graus",15,"--df-g0-hz",.08]),
    ]
    nomes=[]
    for nome,args in casos:
        pr,y=amb.gerar("E5_"+nome,args)
        checar(pr.returncode==0 and os.path.isfile(y),f"gerador falhou em {nome}:\n{(pr.stderr or pr.stdout)[-1000:]}")
        txt=open(y,encoding="utf-8").read(); checar("disjuntor_inicial" in txt and "estrategia_sync" in txt,f"YAML {nome} incompleto")
        amb.sim("E5_run_"+nome,[],config=y)
        nomes.append(nome)
    return "seis YAMLs gerados/executados: "+", ".join(nomes)


def E6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Desempenho do caso ativo é no máximo 2,5× o caso v8 equivalente."""
    comum=[*rede(),"--t-step",.5,"--t-end",4,"--dt-out",2e-4]
    _,_,t8=amb.sim("E6_v8",[*comum,"--evento","fase","--d-fase",5,*LIMITADOR],script=amb.v8,graficos=True,ler=False)
    _,_,t9=amb.sim("E6_v9",fechamento("ativo",15,.08,4,extra=LIMITADOR),graficos=True,ler=False)
    r=t9/max(t8,1e-9); checar(r<=TOL["desempenho_razao"],f"v9/v8={r:.2f}>2,5")
    return f"v8={t8:.1f}s; v9={t9:.1f}s; razão={r:.2f}"


def E7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Repetições no mesmo ambiente têm estado e tempos determinísticos."""
    vals=[]
    for n in range(2):
        js,_,_=amb.sim(f"E7_{n}",fechamento("ativo",20,.12,8)); b=bloco_sync(js)
        vals.append((b["estado_final"],float(b["tempo_comando_s"]),float(b["tempo_contato_s"])))
    checar(vals[0][0]==vals[1][0] and abs(vals[0][1]-vals[1][1])<=CONTRATO["dt_rele_s"] and
           abs(vals[0][2]-vals[1][2])<=CONTRATO["dt_rele_s"],f"resultados não determinísticos: {vals}")
    return f"estado={vals[0][0]}; tcmd={vals[0][1]:.4f}s; contato={vals[0][2]:.4f}s"


TESTES = [
    ("R1",R1),("R2",R2),("R3",R3),("R4",R4),("R5",R5),
    ("D1",D1),("D2",D2),("D3",D3),("D4",D4),("D5",D5),("D6",D6),
    ("S1",S1),("S2",S2),("S3",S3),("S4",S4),("S5",S5),("S6",S6),("S7",S7),("S8",S8),
    ("P1",P1),("P2",P2),("P3",P3),("P4",P4),("P5",P5),("P6",P6),("P7",P7),("P8",P8),("P9",P9),
    ("F1",F1),("F2",F2),("F3",F3),("F4",F4),("F5",F5),("F6",F6),
    ("E1",E1),("E2",E2),("E3",E3),("E4",E4),("E5",E5),("E6",E6),("E7",E7),
]


# ======================================================================================
# Autoteste independente, integridade, relatórios e CLI
# ======================================================================================
def autoteste() -> bool:
    ok=True
    def item(cond: Any,msg: str)->None:
        nonlocal ok
        cond=bool(cond); ok &= cond
        print(f"  [{'ok' if cond else 'FALHA'}] {msg}")

    item(len(TESTES)==41 and len({n for n,_ in TESTES})==41,"matriz congelada: 41 testes distintos")

    vals=np.array([-3*np.pi,-np.pi,-.1,0,.1,np.pi,3*np.pi])
    w=wrap(vals)
    item(np.all(w>=-np.pi) and np.all(w<np.pi) and w[1]==-np.pi and w[5]==-np.pi,
         "wrap em [-pi,pi), incluindo ±pi")

    th=math.radians(5)
    fronteiras=[(.05,0,0),(-.05,0,0),(0,.10,0),(0,-.10,0),(0,0,th),(0,0,-th)]
    item(all(avaliar_sync_oraculo(a,b,c,True)["ok"] for a,b,c in fronteiras),
         "fronteiras inclusivas das três janelas")
    externos=[(.050001,0,0),(0,.100001,0),(0,0,th+1e-6)]
    item(all(not avaliar_sync_oraculo(a,b,c,True)["ok"] for a,b,c in externos) and
         not avaliar_sync_oraculo(0,0,0,False)["ok"],"pontos externos e medição inválida bloqueiam")

    d=.02; f=.08; tf=.063
    pred=predizer_delta_contato(d,f,tf)
    item(abs(wrap(pred-(d+2*np.pi*f*tf)))<1e-14,"predição angular analítica")
    cruza=avaliar_sync_oraculo(0,.10,math.radians(-6),True,t_fechamento_s=.06)
    item(not cruza["ok"] and cruza["janela_prevista"] and not cruza["janela_atual"],
         "janela atual também é obrigatória")

    st=EstadoRele(); hist=[]
    for n in range(900):
        t=n*.001; dentro=(.5<=t<.56) or t>=.58
        st=passo_rele_oraculo(st,t,0,0,0 if dentro else .2,True,"passivo",t_end=2)
        hist.append(st.tempo_janela_s)
        if st.comando: break
    item(hist[570]==0 and st.comando and st.tempo_comando_s is not None and st.tempo_comando_s>=.679,
         f"timer contínuo/reset; comando={st.tempo_comando_s}")
    tc=st.tempo_contato_s
    for n in range(1,100):
        t=float(st.tempo_comando_s+n*.001)
        st=passo_rele_oraculo(st,t,.2,.2,.5,True,"passivo",t_end=2)
        if st.fechado: break
    item(st.fechado and st.fora_da_janela and tc is not None and abs(tc-st.tempo_contato_s)<1e-15,
         "comando latched e contato fora da janela registrado")

    g1=passo_pi_oraculo(2,0,.6,.2,.5,.001)
    it=g1["integral"]
    for _ in range(1000): it=passo_pi_oraculo(2,it,.6,.2,.5,.001)["integral"]
    g2=passo_pi_oraculo(-.2,it,.6,.2,.5,.001)
    item(g1["saturado"] and abs(g1["saida"]-.5)<1e-15 and abs(it)<1e-15 and g2["integral"]<it,
         "PI: saturação, anti-windup condicional e recuperação")

    r=[rampa_retirada_oraculo(.4,t,1,.2) for t in (.9,1,1.1,1.2,1.3)]
    item(np.allclose(r,[.4,.4,.2,0,0],atol=1e-15),"retirada bumpless por rampa linear")

    der=derivada_ig_oraculo(230+0j,230*np.exp(.2j),0j,.02,.001,2*np.pi*60)
    ref=(230+0j-230*np.exp(.2j))/.001
    item(abs(der-ref)<1e-10,"derivada do ramo no contato e continuidade de I_g")

    passivo=oraculo_reduzido("passivo",math.radians(-12),.03,0,t_end=5)
    ativo=oraculo_reduzido("ativo",math.radians(20),.12,.08,t_end=8)
    timeout=oraculo_reduzido("passivo",math.radians(90),0,0,t_end=2)
    item(passivo.fechado and not passivo.timeout and passivo.delta_contato is not None and
         abs(passivo.delta_contato)<=th,"oráculo passivo fecha dentro da janela")
    item(ativo.fechado and not ativo.timeout and ativo.delta_contato is not None and
         abs(ativo.delta_contato)<=th and abs(ativo.df_contato)<=.10 and abs(ativo.dv_contato)<=.05,
         "oráculo ativo converge nas três janelas")
    item(timeout.timeout and not timeout.fechado,"oráculo de condições impossíveis termina em timeout")
    item(ativo.pico_bias_f<=.5+1e-12 and ativo.pico_bias_v<=.1+1e-12,
         "oráculo ativo respeita limites dos biases")

    t=np.linspace(0,3,3001); d_ok=.1+.01*np.exp(-3*t)*np.sin(8*t); f_ok=.02*np.exp(-3*t)
    d_bad=.1+2*np.pi*.4*t; f_bad=np.full_like(t,.4)
    c1,p1=classificar_pos_contato(t,d_ok,f_ok,.5); c2,p2=classificar_pos_contato(t,d_bad,f_bad,.5)
    item(c1=="mantido" and p1==0 and c2=="perdido" and p2>=1,
         f"classificador pós-contato: {c1}/{p1}, {c2}/{p2}")
    return bool(ok)


def verificar_integridade() -> tuple[str,str]:
    eu=os.path.abspath(__file__); atual=sha256(eu)
    arq=os.path.join(os.path.dirname(eu),"testes_v9.sha256")
    if not os.path.isfile(arq): return atual,"SEM REFERÊNCIA (testes_v9.sha256 ausente)"
    esperado=open(arq,encoding="utf-8").read().split()[0].strip().lower()
    return atual,"ÍNTEGRO" if esperado==atual else "ALTERADO — resultados sem validade"


def escrever_relatorios(resultados: list[dict[str,Any]], estado: str, hash_atual: str,
                         script: str, total_s: float) -> None:
    n_ok=sum(r["status"]=="APROVADO" for r in resultados)
    relatorio={"versao_suite":VERSAO_TESTES,"sha256_testes":hash_atual,"integridade":estado,
               "script":script,"sha256_script":sha256(script),"aprovados":n_ok,
               "total":len(resultados),"contrato":CONTRATO,"tolerancias":TOL,
               "resultados":resultados}
    with open("relatorio_testes_v9.json","w",encoding="utf-8") as fh:
        json.dump(relatorio,fh,ensure_ascii=False,indent=2)
    with open("relatorio_testes_v9.md","w",encoding="utf-8") as fh:
        fh.write(f"# Relatório de aceitação — v9\n\n- suíte: `{VERSAO_TESTES}`\n"
                 f"- testes_v9.py SHA-256: `{hash_atual}` ({estado})\n"
                 f"- implementação: `{script}` SHA-256 `{relatorio['sha256_script']}`\n"
                 f"- resultado: **{n_ok}/{len(resultados)} aprovados** em {total_s:.0f} s\n\n"
                 "| Teste | Status | Tempo (s) | Descrição | Detalhe |\n|---|---|---|---|---|\n")
        for r in resultados:
            det=str(r["detalhe"]).replace("\n"," ").replace("|","/")[:600]
            fh.write(f"| {r['teste']} | {r['status']} | {r['tempo_s']} | {r['descricao']} | {det} |\n")


def main() -> None:
    for st in (sys.stdout,sys.stderr):
        try: st.reconfigure(errors="replace")
        except (AttributeError,ValueError): pass
    aqui=os.path.dirname(os.path.abspath(__file__))
    ap=argparse.ArgumentParser(description="Suíte de aceitação congelada da v9.")
    ap.add_argument("--script",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v9.py"))
    ap.add_argument("--v8",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v8.py"))
    ap.add_argument("--testes-v8",default=os.path.join(aqui,"testes_v8_1.py"))
    ap.add_argument("--gerador",default=os.path.join(aqui,"gerar_caso_vsg.py"))
    ap.add_argument("--v7",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v7.py"))
    ap.add_argument("--v6",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v6.py"))
    ap.add_argument("--v5",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v5.py"))
    ap.add_argument("--testes-v7",default=os.path.join(aqui,"testes_v7.py"))
    ap.add_argument("--so",default=None,help="lista separada por vírgulas, ex.: R1,D1,S1,P1")
    ap.add_argument("--autoteste",action="store_true",help="valida contrato/oráculos sem usar a v9")
    ap.add_argument("--manter",action="store_true",help="mantém a pasta temporária")
    a=ap.parse_args()

    hash_atual,estado=verificar_integridade()
    print(f"testes_v9.py  SHA-256 {hash_atual}  [{estado}]")
    if a.autoteste:
        print("Autoteste do contrato e dos oráculos:")
        aprovado=autoteste(); print("AUTOTESTE", "APROVADO" if aprovado else "REPROVADO")
        raise SystemExit(0 if aprovado else 1)

    sel=[x.strip().upper() for x in a.so.split(",") if x.strip()] if a.so else [n for n,_ in TESTES]
    conhecidos={n for n,_ in TESTES}; desconhecidos=set(sel)-conhecidos
    if desconhecidos: raise SystemExit("Erro: testes desconhecidos: "+", ".join(sorted(desconhecidos)))
    obrigatorios=[(a.script,"implementação v9")]
    if any(x in sel for x in ("R1","R2","R3","R4","R5","D2","P8","E6")):
        obrigatorios += [(a.v8,"v8")]
    if "R1" in sel:
        obrigatorios += [(a.testes_v8,"testes_v8_1.py"),(a.gerador,"gerador"),(a.v7,"v7"),
                         (a.v6,"v6"),(a.v5,"v5"),(a.testes_v7,"testes_v7.py")]
    if "E5" in sel: obrigatorios += [(a.gerador,"gerador")]
    for caminho,nome in obrigatorios:
        if not os.path.isfile(caminho): raise SystemExit(f"Erro: {nome} não encontrado: {caminho}")

    work=tempfile.mkdtemp(prefix="testes_v9_")
    amb=Ambiente(a.script,a.v8,a.testes_v8,a.gerador,a.v7,a.v6,a.v5,a.testes_v7,work)
    cache: dict[str,Any]={}; resultados=[]; inicio=time.perf_counter()
    try:
        for nome,fn in TESTES:
            if nome not in sel: continue
            t0=time.perf_counter()
            try: detalhe,status=fn(amb,cache),"APROVADO"
            except AssertionError as exc: detalhe,status=str(exc),"REPROVADO"
            except Exception as exc: detalhe,status=f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=3)}","ERRO"
            dt=time.perf_counter()-t0
            r={"teste":nome,"status":status,"tempo_s":round(dt,1),
               "descricao":(fn.__doc__ or "").strip(),"detalhe":detalhe}
            resultados.append(r)
            print(f"[{status:^9}] {nome:<3} ({dt:6.1f} s) {r['descricao']}\n            {detalhe}")
        total=time.perf_counter()-inicio
        escrever_relatorios(resultados,estado,hash_atual,a.script,total)
        n_ok=sum(r["status"]=="APROVADO" for r in resultados)
        print(f"\nResultado: {n_ok}/{len(resultados)} aprovados em {total:.0f} s | integridade: {estado}")
        raise SystemExit(0 if n_ok==len(resultados) and estado=="ÍNTEGRO" else 1)
    finally:
        if a.manter: print(f"Pasta de trabalho: {work}")
        else: shutil.rmtree(work,ignore_errors=True)


if __name__ == "__main__":
    main()

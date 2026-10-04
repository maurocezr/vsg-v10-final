#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
testes_v10.py — contrato executável CONGELADO da v10 do simulador VSG.

A v10 acrescenta, de forma opt-in, bateria Thévenin OCV(SOC)+R0, capacitor de
barramento CC, estados VDC/SOC, acoplamento energético CA/CC e limite SVPWM
dependente de VDC. O modo modelo_cc=ideal preserva exatamente a v9.

CONTRATO CONGELADO — v10.0.0-contrato-t1
=========================================
1. Convenções e equações
   I_bat > 0: descarga. P_conv > 0: fluxo CC -> CA.
   V_oc = interp_linear(SOC, ocv_soc_pu, ocv_v)
   I_bat = (V_oc - V_dc) / R0
   P_bat_quim = V_oc I_bat
   P_bat_term = V_dc I_bat
   P_perda_R0 = R0 I_bat^2
   dSOC/dt = -I_bat/(3600 capacidade_Ah)
   P_conv = 1,5 Re(E_aplicada conj(I_f))
   I_dc_conv = P_conv/V_dc
   C_dc dV_dc/dt = I_bat - I_dc_conv
   P_cap = V_dc (I_bat - I_dc_conv)
   E_cap = 0,5 C_dc V_dc^2

2. Modulação
   E_max_pu = m_max_pu V_dc/(sqrt(2) V_LL_base)
   e_pos_rv = e_req - Z_v I_f
   O limitador é radial, preserva fase e só marca saturação quando
   |e_pos_rv| > E_max. Igualdade não é saturação. m_utilizado_pu mede a
   solicitação e pode ser maior que 1.

3. Inicialização e validade
   D = V_oc0^2 - 4 R0 P0; a raiz de alta tensão é obrigatória.
   VDC manual desequilibrada é rejeitada, salvo opt-in explícito.
   OCV não extrapola. Saída do domínio OCV e VDC <= piso são terminais.
   piso auto = 1% de VDC inicial e nunca é clamp do denominador.

4. Estados
   0..14 CA; 15 XI_SYNC_P; 16 XI_SYNC_V; 17 VDC; 18 SOC.
   ideal/comum=[0..14], ideal/sync=[0..16],
   thevenin/comum=[0..14,17,18], thevenin/sync=[0..18].

5. API pura exigida da implementação v10
   interpolar_ocv(soc_pu, ocv_soc_pu, ocv_v) -> escalar/array
   bateria_thevenin(soc_pu, vdc_v, r0_ohm, ocv_soc_pu, ocv_v) -> dict com
       {ocv_v, i_bat_a, p_bat_quim_w, p_bat_term_w, p_perda_r0_w}
   potencia_ponte(e_aplicada_complex_v, if_complex_a) -> W
   limitar_modulacao(e_pos_rv_pu, vdc_v, m_max_pu, vll_base_v) -> dict com
       {e_aplicada_pu, e_max_mod_pu, k_mod, m_utilizado_pu,
        modulacao_saturada}
   equilibrio_cc(p_conv_w, soc_pu, r0_ohm, ocv_soc_pu, ocv_v) -> dict com
       {ocv_v, discriminante_v2, vdc_v, i_bat_a, p_max_descarga_w}
   derivadas_cc(vdc_v, soc_pu, p_conv_w, c_dc_f, capacidade_ah, r0_ohm,
       ocv_soc_pu, ocv_v) -> dict com ao menos
       {dvdc_dt_v_s, dsoc_dt_pu_s, i_bat_a, i_dc_conv_a, p_cap_w,
        residuo_potencia_cc_w}
   indices_integracao(modelo_cc, sincronizacao_ativa) -> sequência de índices

6. Interface de produto
   CLI usa hífens; YAML usa sublinhados. JSON sempre contém lado_cc. As séries
   CC existem sempre no NPZ; em ideal, contínuas são NaN e flags são zero.
   Tensões internas e corrente de filtro necessárias à recomposição de P_conv
   são complex128.

7. Integridade e aceite
   A matriz tem exatamente 66 IDs congelados. Oráculos não importam/executam a
   v10. O hash fica em testes_v10.sha256. Execução seletiva auxilia diagnóstico,
   mas não fecha release; aceite final exige 66/66 e regressões internas.

USO
===
  python testes_v10.py --autoteste
  python testes_v10.py
  python testes_v10.py --so M1,M4,M8,I2
  python testes_v10.py --script vsg_2a_ordem_degrau_carga_v10.py
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
from typing import Any, Iterable, Sequence

import numpy as np

VERSAO_TESTES = "v10.0.0-contrato-t1"
DATA_CONGELAMENTO = "2026-10-02"
HASH_SCRIPT_V9 = "077b3e17509cb813bc230326bb7dc4a523050182dd4ad1777f6419921bfad713"
HASH_TESTES_V9 = "0f0d6265eff12b27a8e44d533643fe8cb5ea1085649b55c3ac713b67030a9b8d"
HASH_SCRIPT_V8 = "11cbf6f9861ed240a079d90ce07a428fbafed49280de85751d900a421dfc53e3"
HASH_TESTES_V8 = "f09ad750167d77e95ebe9215147cf0fb2322917941cbb4abb9e917aabb5d7f5d"

CONTRATO = {
    "modelo_cc_padrao": "ideal",
    "modelo_cc_ativo": "thevenin",
    "corrente_positiva": "descarga",
    "pconv_positiva": "cc_para_ca",
    "ponte": "ideal_bidirecional",
    "ocv": "linear_sem_extrapolacao",
    "vdc_piso_auto_frac": 0.01,
    "m_max_pu_padrao": 1.0,
    "soc_min_oper_pu_padrao": 0.10,
    "soc_max_oper_pu_padrao": 0.90,
    "i_desc_max_a_padrao": 0.0,
    "i_carga_max_a_padrao": 0.0,
    "igualdade_svpwm_satura": False,
    "vdc_manual_desequilibrado_padrao": "rejeitar",
    "estado_terminal_vdc": "colapso_vdc",
    "estado_terminal_ocv": "fora_dominio_ocv",
    "indices": {
        "ideal_comum": list(range(15)),
        "ideal_sync": list(range(17)),
        "thevenin_comum": list(range(15)) + [17, 18],
        "thevenin_sync": list(range(19)),
    },
}

TOL = {
    "formula_rel": 1e-10,
    "formula_abs": 1e-12,
    "equilibrio_corrente_abs_a": 1e-9,
    "equilibrio_corrente_rel": 1e-8,
    "residuo_potencia_rel": 1e-8,
    "energia_abs_j": 0.5,
    "energia_rel": 0.005,
    "continuidade_vdc_v": 1e-5,
    "continuidade_soc_pu": 1e-9,
    "svpwm_rel": 1e-12,
    "vdc_atol_v": 1e-6,
    "soc_atol_pu": 1e-10,
    "serie_pu": 1e-6,
    "serie_hz": 1e-6,
    "serie_rad": 1e-6,
    "modal_rel": 1e-5,
    "convergencia_frac": 0.005,
    "desempenho_ideal_razao": 1.05,
    "desempenho_thevenin_razao": 2.0,
}

OCV_SOC = np.array([0.0, 0.25, 0.50, 0.75, 1.0], dtype=float)
OCV_V = np.array([660.0, 680.0, 700.0, 720.0, 740.0], dtype=float)

CLI_V10 = (
    "--modelo-cc", "--c-dc-f", "--vdc-inicial-v",
    "--permitir-desequilibrio-inicial", "--vdc-min-oper-v",
    "--vdc-max-oper-v", "--m-max-pu", "--vdc-piso-numerico-v",
    "--capacidade-ah", "--soc-inicial", "--r0-ohm", "--ocv-soc-pu",
    "--ocv-v", "--soc-min-oper-pu", "--soc-max-oper-pu",
    "--i-desc-max-a", "--i-carga-max-a",
)

SERIES_CC_CONTINUAS = (
    "vdc_v", "vdc_pu", "soc_pu", "ocv_v", "i_bat_a",
    "i_dc_conv_a", "p_conv_w", "p_bat_term_w", "p_bat_quim_w",
    "p_perda_r0_w", "p_cap_w", "e_cap_j", "residuo_potencia_cc_w",
    "e_req_pu", "e_pos_rv_pu", "e_aplicada_pu", "e_max_mod_pu",
    "m_utilizado_pu", "if_complex_a",
)
SERIES_CC_FLAGS = (
    "modulacao_saturada", "vdc_subtensao", "vdc_sobretensao",
    "soc_fora_faixa", "i_bat_fora_limite",
)
SERIES_COMPLEXAS = ("e_req_pu", "e_pos_rv_pu", "e_aplicada_pu", "if_complex_a")
CAMPOS_LADO_CC = (
    "modelo", "ativo", "convencao_corrente", "estado_final",
    "soc_inicial_pu", "soc_final_pu", "soc_min_pu", "soc_max_pu",
    "vdc_inicial_v", "vdc_final_v", "vdc_min_v", "vdc_max_v",
    "i_bat_min_a", "i_bat_max_a", "p_bat_quim_pico_w",
    "p_bat_term_pico_w", "p_perda_r0_pico_w", "p_conv_pico_w",
    "energia_bateria_quimica_j", "energia_bateria_terminal_j",
    "perda_r0_j", "energia_ponte_j", "energia_cap_delta_j",
    "residuo_energia_abs_j", "residuo_energia_rel", "modulacao",
    "limites", "equilibrio_inicial", "residuo_corrente_inicial_a",
    "dvdc_dt_inicial_v_s",
)

GRUPOS = {
    "R": [f"R{i}" for i in range(1, 7)],
    "M": [f"M{i}" for i in range(1, 11)],
    "I": [f"I{i}" for i in range(1, 9)],
    "D": [f"D{i}" for i in range(1, 11)],
    "B": [f"B{i}" for i in range(1, 7)],
    "C": [f"C{i}" for i in range(1, 7)],
    "X": [f"X{i}" for i in range(1, 11)],
    "E": [f"E{i}" for i in range(1, 11)],
}


def rel(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-15)


def escala(*valores: Any) -> float:
    return max(1.0, *(abs(float(v)) for v in valores))


def perto(a: float, b: float, reltol: float = TOL["formula_rel"],
          abstol: float = TOL["formula_abs"]) -> bool:
    return abs(float(a) - float(b)) <= max(abstol, reltol * escala(a, b))


def checar(cond: Any, msg: str) -> None:
    if not bool(cond):
        raise AssertionError(msg)


def sha256(caminho: str) -> str:
    h = hashlib.sha256()
    with open(caminho, "rb") as fh:
        for bloco in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(bloco)
    return h.hexdigest()


def integral_trapezio(y: Sequence[float], t: Sequence[float]) -> float:
    yy, tt = np.asarray(y, dtype=float), np.asarray(t, dtype=float)
    checar(yy.shape == tt.shape and yy.ndim == 1 and len(yy) >= 2,
           "vetores inválidos para integração")
    return float(np.trapezoid(yy, tt) if hasattr(np, "trapezoid") else np.trapz(yy, tt))


def ocv_oraculo(soc: Any, xs: Sequence[float] = OCV_SOC,
                ys: Sequence[float] = OCV_V) -> Any:
    x, y, q = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float), np.asarray(soc, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 2:
        raise ValueError("curva OCV exige listas unidimensionais de mesmo comprimento >=2")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(np.diff(x) > 0):
        raise ValueError("curva OCV inválida")
    if np.any(x < 0) or np.any(x > 1) or np.any(y <= 0):
        raise ValueError("domínio/valores da curva OCV inválidos")
    if np.any(q < x[0]) or np.any(q > x[-1]):
        raise ValueError("SOC fora do domínio OCV; extrapolação proibida")
    out = np.interp(q, x, y)
    return float(out) if out.ndim == 0 else out


def bateria_oraculo(soc: float, vdc: float, r0: float,
                    xs: Sequence[float] = OCV_SOC,
                    ys: Sequence[float] = OCV_V) -> dict[str, float]:
    if not np.isfinite(vdc) or vdc <= 0 or not np.isfinite(r0) or r0 <= 0:
        raise ValueError("VDC e R0 devem ser positivos")
    voc = float(ocv_oraculo(soc, xs, ys))
    ibat = (voc - float(vdc)) / float(r0)
    return {
        "ocv_v": voc,
        "i_bat_a": ibat,
        "p_bat_quim_w": voc * ibat,
        "p_bat_term_w": float(vdc) * ibat,
        "p_perda_r0_w": float(r0) * ibat * ibat,
    }


def potencia_ponte_oraculo(e: complex, corrente: complex) -> float:
    return 1.5 * float(np.real(complex(e) * np.conj(complex(corrente))))


def modulacao_oraculo(e_pos_rv: complex, vdc: float, mmax: float,
                      vll_base: float) -> dict[str, Any]:
    if min(float(vdc), float(mmax), float(vll_base)) <= 0:
        raise ValueError("VDC, m_max e VLL devem ser positivos")
    e = complex(e_pos_rv)
    emax = float(mmax) * float(vdc) / (math.sqrt(2.0) * float(vll_base))
    mag = abs(e)
    uso = mag / emax
    saturada = bool(mag > emax * (1.0 + TOL["svpwm_rel"]))
    k = emax / mag if saturada and mag > 0 else 1.0
    return {"e_aplicada_pu": k * e, "e_max_mod_pu": emax, "k_mod": k,
            "m_utilizado_pu": uso, "modulacao_saturada": saturada}


def equilibrio_oraculo(pconv: float, soc: float, r0: float,
                       xs: Sequence[float] = OCV_SOC,
                       ys: Sequence[float] = OCV_V) -> dict[str, float]:
    if r0 <= 0:
        raise ValueError("R0 deve ser positivo")
    voc = float(ocv_oraculo(soc, xs, ys))
    disc = voc * voc - 4.0 * float(r0) * float(pconv)
    if disc < 0:
        pmax = voc * voc / (4.0 * float(r0))
        raise ValueError(f"potência impossível; potência máxima de descarga={pmax:.12g} W")
    vdc = 0.5 * (voc + math.sqrt(max(0.0, disc)))
    ibat = float(pconv) / vdc
    return {"ocv_v": voc, "discriminante_v2": disc, "vdc_v": vdc,
            "i_bat_a": ibat, "p_max_descarga_w": voc * voc / (4.0 * float(r0))}


def derivadas_cc_oraculo(vdc: float, soc: float, pconv: float, cdc: float,
                         qah: float, r0: float, xs: Sequence[float] = OCV_SOC,
                         ys: Sequence[float] = OCV_V) -> dict[str, float]:
    if min(float(vdc), float(cdc), float(qah), float(r0)) <= 0:
        raise ValueError("parâmetros CC devem ser positivos")
    b = bateria_oraculo(soc, vdc, r0, xs, ys)
    iconv = float(pconv) / float(vdc)
    dv = (b["i_bat_a"] - iconv) / float(cdc)
    dsoc = -b["i_bat_a"] / (3600.0 * float(qah))
    pcap = float(vdc) * float(cdc) * dv
    resid = b["p_bat_quim_w"] - b["p_perda_r0_w"] - float(pconv) - pcap
    return {**b, "i_dc_conv_a": iconv, "dvdc_dt_v_s": dv,
            "dsoc_dt_pu_s": dsoc, "p_cap_w": pcap,
            "residuo_potencia_cc_w": resid}


def limite_energia(residuo_j: float, *energias_j: float) -> float:
    return max(TOL["energia_abs_j"], TOL["energia_rel"] * escala(*energias_j))


def exigir_funcao(mod: Any, nome: str):
    fn = getattr(mod, nome, None)
    checar(callable(fn), f"função pública obrigatória ausente: {nome}")
    return fn


def exigir_dict(d: Any, campos: Iterable[str], contexto: str) -> dict[str, Any]:
    checar(isinstance(d, dict), f"{contexto} deve retornar dict")
    for k in campos:
        checar(k in d, f"{contexto} sem campo '{k}'")
    return d


def texto_erro_limpo(pr: subprocess.CompletedProcess[str]) -> str:
    txt = ((pr.stdout or "") + (pr.stderr or "")).strip()
    checar(pr.returncode != 0, "entrada inválida deveria falhar")
    checar(re.search(r"(^|\n)Erro", txt) is not None, "erro não começa por 'Erro'")
    checar("Traceback" not in txt, "entrada inválida expôs traceback")
    return txt[-500:]


@dataclass
class Ambiente:
    script: str
    v9: str
    testes_v9: str
    v8: str
    testes_v8: str
    gerador: str
    v7: str
    v6: str
    v5: str
    testes_v7: str
    workdir: str
    _modulo: Any = None
    contador: int = 0

    def __post_init__(self) -> None:
        for nome in ("script", "v9", "testes_v9", "v8", "testes_v8", "gerador",
                     "v7", "v6", "v5", "testes_v7"):
            setattr(self, nome, os.path.abspath(getattr(self, nome)))
        self.base_json = os.path.join(self.workdir, "base.json")
        with open(self.base_json, "w", encoding="utf-8") as fh:
            json.dump({"sistema": {"sn": 100000, "vll": 380, "f0": 60}}, fh)

    def prefixo(self, nome: str) -> str:
        self.contador += 1
        return os.path.join(self.workdir, f"{self.contador:03d}_{nome}")

    def modulo(self):
        if self._modulo is None:
            spec = importlib.util.spec_from_file_location("vsg_v10_sob_teste", self.script)
            checar(spec is not None and spec.loader is not None, "não foi possível importar v10")
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            self._modulo = mod
        return self._modulo

    def sim(self, nome: str, args: Iterable[Any] = (), *, script: str | None = None,
            config: str | None = None, esperar_erro: bool = False,
            graficos: bool = False, timeout: int = 2400) -> Any:
        pref = self.prefixo(nome)
        cmd = [sys.executable, script or self.script, "--config", config or self.base_json,
               *[str(x) for x in args], "--prefixo", pref]
        if not graficos:
            cmd.append("--sem-graficos")
        t0 = time.perf_counter()
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                            encoding="utf-8", errors="replace")
        dur = time.perf_counter() - t0
        if esperar_erro:
            return pr
        saida = (pr.stdout or "") + (pr.stderr or "")
        checar(pr.returncode == 0, f"simulador falhou ({pr.returncode}):\n{saida[-3000:]}")
        with open(pref + "_resultados.json", encoding="utf-8") as fh:
            js = json.load(fh)
        with np.load(pref + "_series.npz") as arq:
            ser = {k: np.asarray(arq[k]) for k in arq.files}
        return js, ser, dur, pref


def args_rede(scr: float = 5.0, dw: float = 100.0) -> list[Any]:
    return ["--modo", "rede", "--scr", scr, "--dw", dw, "--tw", 1.0]


def args_cc(*extra: Any, cdc: float = 0.10, r0: float = 0.05,
            capacidade: float = 100.0, soc: float = 0.60,
            mmax: float = 1.0) -> list[Any]:
    return ["--modelo-cc", "thevenin", "--c-dc-f", cdc,
            "--vdc-inicial-v", "auto", "--vdc-min-oper-v", 500,
            "--vdc-max-oper-v", 800, "--m-max-pu", mmax,
            "--vdc-piso-numerico-v", "auto", "--capacidade-ah", capacidade,
            "--soc-inicial", soc, "--r0-ohm", r0,
            "--ocv-soc-pu", "0,0.25,0.5,0.75,1",
            "--ocv-v", "660,680,700,720,740", "--soc-min-oper-pu", 0.10,
            "--soc-max-oper-pu", 0.90, "--i-desc-max-a", 0,
            "--i-carga-max-a", 0, *extra]


def args_curto(*extra: Any) -> list[Any]:
    return ["--evento", "nenhum", "--t-step", .2, "--t-end", .4,
            "--dt-out", 1e-4, *extra]


def caso(amb: Ambiente, cache: dict[str, Any], nome: str, args: Iterable[Any]) -> Any:
    if nome not in cache:
        cache[nome] = amb.sim(nome, args)
    return cache[nome]


def caso_cc_nominal(amb: Ambiente, cache: dict[str, Any]) -> Any:
    """Caso curto, nominal e sem evento, compartilhado sem depender da ordem da suíte."""
    return caso(amb, cache, "cc_nominal", args_curto(*args_cc()))


def caso_degrau_carga(amb: Ambiente, cache: dict[str, Any]) -> Any:
    """Degrau de carga nominal compartilhado sem depender da execução prévia de D1."""
    args = ["--modo", "ilhado", "--evento", "carga", "--t-step", .5,
            "--t-end", 1.2, "--dt-out", 1e-4, *args_cc()]
    return caso(amb, cache, "D1_carga", args)


def bloco_cc(js: dict[str, Any], ativo: bool | None = None) -> dict[str, Any]:
    b = js.get("lado_cc")
    checar(isinstance(b, dict), "bloco JSON lado_cc ausente")
    if ativo is not None:
        checar(bool(b.get("ativo")) is ativo, f"lado_cc.ativo deveria ser {ativo}")
    return b


def checar_series_cc(s: dict[str, np.ndarray], *, ideal: bool = False) -> None:
    checar("t" in s and np.ndim(s["t"]) == 1, "série t ausente/inválida")
    n = len(s["t"])
    for k in SERIES_CC_CONTINUAS + SERIES_CC_FLAGS:
        checar(k in s, f"série v10 ausente: {k}")
        checar(np.shape(s[k]) == (n,), f"{k}: shape {np.shape(s[k])}, esperado {(n,)}")
    for k in SERIES_COMPLEXAS:
        checar(np.issubdtype(s[k].dtype, np.complexfloating), f"{k} deve ser complex")
    if ideal:
        for k in SERIES_CC_CONTINUAS:
            checar(np.all(np.isnan(s[k])), f"{k} deve ser NaN no modo ideal")
        for k in SERIES_CC_FLAGS:
            checar(np.all(np.asarray(s[k]) == 0), f"{k} deve ser zero no modo ideal")
    else:
        for k in SERIES_CC_CONTINUAS:
            checar(np.all(np.isfinite(s[k])), f"{k} contém NaN/Inf no modo ativo")


def comparar_series_v9(ref: dict[str, np.ndarray], novo: dict[str, np.ndarray]) -> float:
    nomes = {"Pf_pu": "serie_pu", "P_pu": "serie_pu", "P_rede_pu": "serie_pu",
             "Q_rede_pu": "serie_pu", "f_hz": "serie_hz",
             "delta_v_rad": "serie_rad", "if_env_pu": "serie_pu",
             "rv_pu": "serie_pu", "xv_pu": "serie_pu",
             "delta_aw_rad_s": "serie_rad", "vg_aplicada_pu": "serie_pu"}
    checar(len(ref["t"]) == len(novo["t"]) and np.allclose(ref["t"], novo["t"], atol=1e-12, rtol=0),
           "vetores de tempo v9/v10 diferentes")
    pior = 0.0
    for k, tk in nomes.items():
        checar(k in ref and k in novo, f"série comum ausente: {k}")
        d = float(np.nanmax(np.abs(ref[k] - novo[k])))
        pior = max(pior, d / TOL[tk])
        checar(d <= TOL[tk], f"{k}: {d:.3e}>{TOL[tk]:.1e}")
    return pior


def saltos_em(t: np.ndarray, x: np.ndarray, alvo: float) -> float:
    i = int(np.searchsorted(t, alvo, side="left"))
    if i <= 0 or i >= len(t):
        raise AssertionError(f"evento {alvo} fora da malha")
    return float(abs(x[i] - x[i - 1]))


def reconstruir_energia(s: dict[str, np.ndarray]) -> dict[str, float]:
    t = s["t"]
    eq = integral_trapezio(s["p_bat_quim_w"], t)
    et = integral_trapezio(s["p_bat_term_w"], t)
    er = integral_trapezio(s["p_perda_r0_w"], t)
    ep = integral_trapezio(s["p_conv_w"], t)
    ec = float(s["e_cap_j"][-1] - s["e_cap_j"][0])
    resid = eq - er - ep - ec
    return {"energia_bateria_quimica_j": eq, "energia_bateria_terminal_j": et,
            "perda_r0_j": er, "energia_ponte_j": ep,
            "energia_cap_delta_j": ec, "residuo_energia_abs_j": abs(resid)}


# Regressão R1-R6

def R1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Executa a suíte v9 congelada completa contra v10 em modelo_cc=ideal."""
    checar(sha256(amb.testes_v9) == HASH_TESTES_V9, "testes_v9.py diverge da baseline")
    cmd = [sys.executable, amb.testes_v9, "--script", amb.script, "--v8", amb.v8,
           "--testes-v8", amb.testes_v8, "--gerador", amb.gerador, "--v7", amb.v7,
           "--v6", amb.v6, "--v5", amb.v5, "--testes-v7", amb.testes_v7]
    pr = subprocess.run(cmd, cwd=amb.workdir, capture_output=True, text=True, timeout=14400,
                        encoding="utf-8", errors="replace")
    txt = (pr.stdout or "") + (pr.stderr or "")
    checar(pr.returncode == 0 and re.search(r"Resultado:\s*41/41 aprovados", txt),
           "v9 integral não aprovou contra v10/ideal:\n" + txt[-4000:])
    return "41/41 testes v9 aprovados contra v10/ideal"


def R2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Confirma explicitamente os 19 contratos herdados da v8."""
    checar(sha256(amb.testes_v9) == HASH_TESTES_V9, "testes_v9.py diverge da baseline")
    cmd = [sys.executable, amb.testes_v9, "--so", "R1", "--script", amb.script,
           "--v8", amb.v8, "--testes-v8", amb.testes_v8, "--gerador", amb.gerador,
           "--v7", amb.v7, "--v6", amb.v6, "--v5", amb.v5, "--testes-v7", amb.testes_v7]
    pr = subprocess.run(cmd, cwd=amb.workdir, capture_output=True, text=True, timeout=10800,
                        encoding="utf-8", errors="replace")
    txt = (pr.stdout or "") + (pr.stderr or "")
    checar(pr.returncode == 0 and "19/19" in txt, "19 contratos v8 não comprovados:\n" + txt[-4000:])
    return "19/19 contratos herdados aprovados"


def R3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Compara numericamente v9 e v10/ideal em modos, eventos, SCR e limitador."""
    checar(sha256(amb.v9) == HASH_SCRIPT_V9, "simulador v9 diverge da baseline")
    casos = [
        ("ilha", ["--modo", "ilhado", "--evento", "carga", "--t-step", .3, "--t-end", 1.2, "--dt-out", 2e-4]),
        ("pref_scr2", [*args_rede(2), "--evento", "pref", "--t-step", .3, "--t-end", 1.2, "--dt-out", 2e-4]),
        ("fase_scr5", [*args_rede(5), "--evento", "fase", "--d-fase", 5, "--t-step", .3, "--t-end", 1.2, "--dt-out", 1e-4]),
        ("freq_scr20", [*args_rede(20), "--evento", "freq_rampa", "--df-g", -.2, "--rocof-g", .2, "--t-step", .3, "--t-end", 1.2, "--dt-out", 2e-4]),
        ("falta", [*args_rede(20), "--evento", "falta_3f", "--vg-falta", 0, "--t-step", .3, "--t-clear", .4, "--t-end", 1.2, "--dt-out", 1e-4, "--imax-pu", 1.2, "--i-on-pu", 1.176, "--rv-max-pu", 3]),
    ]
    out = []
    for nome, a in casos:
        _, s9, _, _ = amb.sim("R3_" + nome + "_v9", a, script=amb.v9)
        j10, s10, _, _ = amb.sim("R3_" + nome + "_v10", [*a, "--modelo-cc", "ideal"])
        checar(j10.get("versao") == "v10", "campo versão não é v10")
        out.append(f"{nome}:{comparar_series_v9(s9, s10):.2f}xtol")
    return "; ".join(out)


def R4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Valida mapas congelados 15/17/17/19 e índices VDC=17, SOC=18."""
    mod = amb.modulo()
    fn = exigir_funcao(mod, "indices_integracao")
    obt = {
        "ideal_comum": list(fn("ideal", False)), "ideal_sync": list(fn("ideal", True)),
        "thevenin_comum": list(fn("thevenin", False)),
        "thevenin_sync": list(fn("thevenin", True)),
    }
    checar(obt == CONTRATO["indices"], f"mapas de estados divergentes: {obt}")
    checar(getattr(mod, "VDC", 17) == 17 and getattr(mod, "SOC", 18) == 18,
           "índices VDC/SOC divergentes")
    return "mapas 15/17/17/19 preservam índices históricos"


def R5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Valida interface histórica aditiva e representação inerte do lado CC."""
    hp = subprocess.run([sys.executable, amb.script, "--help"], capture_output=True,
                        text=True, timeout=120, encoding="utf-8", errors="replace")
    checar(hp.returncode == 0, "--help falhou")
    for op in CLI_V10:
        checar(op in hp.stdout, f"opção v10 ausente: {op}")
    js, s, _, _ = caso(amb, cache, "ideal_curto", args_curto("--modelo-cc", "ideal"))
    b = bloco_cc(js, False)
    checar(b.get("modelo") == "ideal" and b.get("estado_final") == "ideal",
           "semântica JSON ideal divergente")
    checar_series_cc(s, ideal=True)
    return "CLI aditiva; JSON/NPZ ideal inertes conforme contrato"


def R6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Preserva autovalores históricos no modo ideal."""
    a = [*args_rede(5), "--evento", "nenhum", "--t-step", .2, "--t-end", .4, "--dt-out", 1e-4]
    j9, _, _, _ = amb.sim("R6_v9", a, script=amb.v9)
    j10, _, _, _ = amb.sim("R6_v10", [*a, "--modelo-cc", "ideal"])
    def eig(js: dict[str, Any]) -> np.ndarray:
        md = js.get("modal_fechado", js.get("modal", {}))
        return np.sort_complex(np.array([complex(r, i) for r, i in md["autovalores"]]))
    e9, e10 = eig(j9), eig(j10)
    checar(e9.shape == e10.shape, "quantidade de autovalores mudou no modo ideal")
    erro = float(np.max(np.abs(e9 - e10) / np.maximum(1.0, np.abs(e9))))
    checar(erro <= TOL["modal_rel"], f"autovalores divergiram: {erro:.3e}")
    return f"erro modal relativo máximo={erro:.3e}"


# Modelo físico e funções puras M1-M10

def M1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """OCV coincide exatamente com todos os nós da tabela."""
    fn = exigir_funcao(amb.modulo(), "interpolar_ocv")
    got = np.asarray(fn(OCV_SOC, OCV_SOC, OCV_V), dtype=float)
    checar(np.array_equal(got, OCV_V), f"OCV nos nós divergente: {got}")
    return "todos os nós OCV exatos"


def M2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """OCV usa interpolação linear em pontos médios e determinísticos."""
    fn = exigir_funcao(amb.modulo(), "interpolar_ocv")
    q = np.array([.125, .375, .625, .875, .03125, .9375])
    ref = ocv_oraculo(q)
    got = np.asarray(fn(q, OCV_SOC, OCV_V), dtype=float)
    checar(np.allclose(got, ref, rtol=TOL["formula_rel"], atol=TOL["formula_abs"]),
           f"interpolação OCV divergente: {got} != {ref}")
    return "interpolação linear aprovada"


def M3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """OCV rejeita extrapolação e tabelas inválidas."""
    fn = exigir_funcao(amb.modulo(), "interpolar_ocv")
    casos = [(-.01, OCV_SOC, OCV_V), (1.01, OCV_SOC, OCV_V),
             (.5, [0, .5, .5, 1], [660, 690, 700, 740]),
             (.5, [0, 1], [660]), (.5, [0, 1], [660, -1])]
    for c in casos:
        try: fn(*c)
        except (ValueError, SystemExit): pass
        else: raise AssertionError(f"curva/domínio inválido aceito: {c}")
    return f"{len(casos)} entradas OCV inválidas rejeitadas"


def M4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Thévenin em descarga fecha corrente e identidades de potência."""
    fn = exigir_funcao(amb.modulo(), "bateria_thevenin")
    ref = bateria_oraculo(.5, 690.0, .05)
    got = exigir_dict(fn(.5, 690.0, .05, OCV_SOC, OCV_V), ref, "bateria_thevenin")
    for k, v in ref.items(): checar(perto(got[k], v), f"{k}: {got[k]} != {v}")
    checar(got["i_bat_a"] > 0 and perto(got["p_bat_quim_w"],
           got["p_bat_term_w"] + got["p_perda_r0_w"]), "identidade de descarga falhou")
    return "descarga e identidades algébricas aprovadas"


def M5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Thévenin em carga preserva sinais e identidade de potência."""
    fn = exigir_funcao(amb.modulo(), "bateria_thevenin")
    ref = bateria_oraculo(.5, 710.0, .05)
    got = exigir_dict(fn(.5, 710.0, .05, OCV_SOC, OCV_V), ref, "bateria_thevenin")
    checar(got["i_bat_a"] < 0 and got["p_bat_term_w"] < 0 and got["p_perda_r0_w"] > 0,
           "sinais de regeneração incorretos")
    for k, v in ref.items(): checar(perto(got[k], v), f"{k}: {got[k]} != {v}")
    return "carga/regeneração e sinais aprovados"


def M6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """SOC segue reta analítica para correntes positivas e negativas."""
    fn = exigir_funcao(amb.modulo(), "derivadas_cc")
    for vdc in (690.0, 710.0):
        got = fn(vdc, .5, 0.0, .2, 100.0, .05, OCV_SOC, OCV_V)
        ib = bateria_oraculo(.5, vdc, .05)["i_bat_a"]
        ref = -ib / (3600 * 100.0)
        checar(perto(got["dsoc_dt_pu_s"], ref), "derivada SOC divergente")
        soc10 = .5 + 10 * got["dsoc_dt_pu_s"]
        checar((soc10 < .5) == (ib > 0), "sentido do SOC incorreto")
    return "integração coulômbica bidirecional aprovada"


def M7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Capacitor reproduz inclinação analítica com corrente líquida constante."""
    fn = exigir_funcao(amb.modulo(), "derivadas_cc")
    got = fn(700.0, .5, 35000.0, .2, 100.0, .05, OCV_SOC, OCV_V)
    ref = derivadas_cc_oraculo(700.0, .5, 35000.0, .2, 100.0, .05)
    checar(perto(got["dvdc_dt_v_s"], ref["dvdc_dt_v_s"]), "dVDC/dt divergente")
    checar(perto(.5 * .2 * 700**2, 49000.0), "oráculo de energia do capacitor inválido")
    return f"dVDC/dt={ref['dvdc_dt_v_s']:.9g} V/s"


def M8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Modulação radial cobre zero, limite, fora e quatro quadrantes."""
    fn = exigir_funcao(amb.modulo(), "limitar_modulacao")
    emax = 1.0
    vdc = math.sqrt(2) * 380 * emax
    vetores = [0j, 1+0j, -1+0j, 1j, -1j, 2*np.exp(1j*.4), 2*np.exp(1j*2.3)]
    for e in vetores:
        g = exigir_dict(fn(e, vdc, 1.0, 380.0),
            ("e_aplicada_pu", "e_max_mod_pu", "k_mod", "m_utilizado_pu", "modulacao_saturada"),
            "limitar_modulacao")
        r = modulacao_oraculo(e, vdc, 1.0, 380.0)
        checar(abs(complex(g["e_aplicada_pu"]) - r["e_aplicada_pu"]) <= 1e-12, "vetor aplicado divergente")
        checar(bool(g["modulacao_saturada"]) == r["modulacao_saturada"], "flag SVPWM divergente")
        if e and abs(g["e_aplicada_pu"]) > 0:
            checar(abs(np.angle(g["e_aplicada_pu"] / e)) <= 1e-12, "fase não preservada")
    return "zero, fronteira e quadrantes aprovados; igualdade não satura"


def M9(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Potência da ponte usa 1,5 Re(E aplicada vezes conjugado de If)."""
    fn = exigir_funcao(amb.modulo(), "potencia_ponte")
    pares = [(230+0j, 100+0j), (230j, 100+0j), (100+80j, 40-20j), (230+0j, -50+0j)]
    for e, i in pares:
        ref = potencia_ponte_oraculo(e, i); got = float(fn(e, i))
        checar(perto(got, ref), f"Pconv divergente: {got} != {ref}")
    return "fasores resistivo, reativo e regenerativo aprovados"


def M10(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Balanço instantâneo CC fecha para combinações sintéticas."""
    fn = exigir_funcao(amb.modulo(), "derivadas_cc")
    for v, soc, p in [(690, .5, 20e3), (710, .5, -15e3), (700, .5, 0)]:
        got = fn(v, soc, p, .2, 100, .05, OCV_SOC, OCV_V)
        mag = escala(got["p_bat_quim_w"], got["p_perda_r0_w"], p, got["p_cap_w"])
        checar(abs(got["residuo_potencia_cc_w"]) <= TOL["formula_rel"] * mag,
               f"resíduo instantâneo excessivo: {got['residuo_potencia_cc_w']}")
    return "balanço instantâneo fecha nos três sentidos"


# Inicialização I1-I8

def I1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Equilíbrio automático nominal fecha corrente e derivada inicial."""
    js, s, _, _ = caso_cc_nominal(amb, cache)
    b = bloco_cc(js, True)
    checar(bool(b["equilibrio_inicial"]), "equilíbrio inicial não reconhecido")
    lim = max(TOL["equilibrio_corrente_abs_a"], TOL["equilibrio_corrente_rel"] * escala(s["i_bat_a"][0], s["i_dc_conv_a"][0]))
    checar(abs(float(s["i_bat_a"][0] - s["i_dc_conv_a"][0])) <= lim, "KCL inicial não fecha")
    checar(abs(float(b["dvdc_dt_inicial_v_s"])) <= lim / .1, "dVDC/dt inicial excessivo")
    return f"resíduo inicial={b['residuo_corrente_inicial_a']:.3e} A"


def I2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Equilíbrio escolhe a raiz positiva de alta tensão."""
    fn = exigir_funcao(amb.modulo(), "equilibrio_cc")
    p = 50e3; ref = equilibrio_oraculo(p, .6, .05)
    got = exigir_dict(fn(p, .6, .05, OCV_SOC, OCV_V), ref, "equilibrio_cc")
    for k, v in ref.items(): checar(perto(got[k], v), f"{k} divergente")
    raiz_baixa = .5 * (ref["ocv_v"] - math.sqrt(ref["discriminante_v2"]))
    checar(got["vdc_v"] > raiz_baixa, "raiz baixa escolhida")
    return f"raiz alta={got['vdc_v']:.6f} V; baixa={raiz_baixa:.6f} V"


def I3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Regeneração inicial produz corrente negativa e equilíbrio válido."""
    fn = exigir_funcao(amb.modulo(), "equilibrio_cc")
    got = fn(-20e3, .6, .05, OCV_SOC, OCV_V)
    checar(got["i_bat_a"] < 0 and got["vdc_v"] > got["ocv_v"], "regeneração inicial inválida")
    checar(perto(got["vdc_v"] * got["i_bat_a"], -20e3), "potência regenerativa não fecha")
    return f"VDC={got['vdc_v']:.3f} V; Ibat={got['i_bat_a']:.3f} A"


def I4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Potência acima do máximo Thévenin é rejeitada com causa explícita."""
    fn = exigir_funcao(amb.modulo(), "equilibrio_cc")
    pmax = float(ocv_oraculo(.5))**2/(4*.05)
    try: fn(pmax*(1+1e-9), .5, .05, OCV_SOC, OCV_V)
    except (ValueError, SystemExit) as e:
        txt = str(e).lower(); checar("pot" in txt and ("máx" in txt or "max" in txt or "discrimin" in txt), "causa pouco clara")
    else: raise AssertionError("potência impossível aceita")
    return f"Pmax={pmax:.6g} W rejeita excedente"


def I5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Headroom inicial insuficiente é erro/aviso explícito, nunca oculto."""
    pr = amb.sim("I5_headroom", args_curto(*args_cc(mmax=.2)), esperar_erro=True)
    txt = (pr.stdout or "") + (pr.stderr or "")
    if pr.returncode:
        checar("Traceback" not in txt and re.search(r"Erro", txt), "falha de headroom não foi limpa")
        return "headroom insuficiente rejeitado"
    pref = os.path.join(amb.workdir, f"{amb.contador:03d}_I5_headroom")
    with open(pref + "_resultados.json", encoding="utf-8") as fh: js = json.load(fh)
    b = bloco_cc(js, True)
    checar(bool(b["modulacao"]["saturou"]), "headroom insuficiente ficou oculto")
    return "headroom insuficiente explicitamente classificado"


def I6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """VDC manual igual ao oráculo reproduz o estado automático."""
    fn = exigir_funcao(amb.modulo(), "equilibrio_cc")
    v = float(fn(50e3, .6, .05, OCV_SOC, OCV_V)["vdc_v"])
    auto = caso(amb, cache, "I6_auto", args_curto(*args_cc()))[1]
    man = amb.sim("I6_manual", args_curto(*args_cc("--vdc-inicial-v", f"{v:.15g}")))[1]
    checar(abs(float(auto["vdc_v"][0]-man["vdc_v"][0])) <= TOL["continuidade_vdc_v"], "VDC manual não reproduz auto")
    return f"VDC manual equilibrada={v:.9f} V"


def I7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """VDC manual desequilibrada é rejeitada por padrão e aceita apenas por opt-in."""
    pr = amb.sim("I7_rejeita", args_curto(*args_cc("--vdc-inicial-v", 650)), esperar_erro=True)
    texto_erro_limpo(pr)
    js, _, _, _ = amb.sim("I7_optin", args_curto(*args_cc("--vdc-inicial-v", 650,
        "--permitir-desequilibrio-inicial", "true")))
    b = bloco_cc(js, True)
    checar(b["equilibrio_inicial"] is False and abs(float(b["dvdc_dt_inicial_v_s"])) > 0,
           "opt-in não publicou desequilíbrio")
    return "rejeição padrão e opt-in explícito aprovados"


def I8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """SOC, OCV e limites inconsistentes falham limpidamente."""
    casos = [
        ["--soc-inicial", 1.1], ["--r0-ohm", 0], ["--c-dc-f", 0],
        ["--ocv-soc-pu", "0,.5,.5,1"], ["--ocv-v", "660,700"],
        ["--soc-min-oper-pu", .7, "--soc-max-oper-pu", .6],
    ]
    for n, extra in enumerate(casos):
        pr = amb.sim(f"I8_{n}", args_curto(*args_cc(*extra)), esperar_erro=True)
        texto_erro_limpo(pr)
    return f"{len(casos)} combinações inválidas rejeitadas"


# Dinâmica CC D1-D10

def D1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """VDC e SOC permanecem contínuos em degraus de carga e Pref."""
    casos = [
        caso_degrau_carga(amb, cache),
        caso(amb, cache, "D1_pref", [*args_rede(), "--evento", "pref", "--t-step", .5, "--t-end", 1.2, "--dt-out", 1e-4, *args_cc()]),
    ]
    for _, s, _, _ in casos:
        checar(saltos_em(s["t"], s["vdc_v"], .5) <= TOL["continuidade_vdc_v"], "salto em VDC")
        checar(saltos_em(s["t"], s["soc_pu"], .5) <= TOL["continuidade_soc_pu"], "salto em SOC")
    return "continuidade em carga e Pref aprovada"


def D2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Derivada imediatamente após evento coincide com a KCL independente."""
    _, s, _, _ = caso_degrau_carga(amb, cache)
    i = int(np.searchsorted(s["t"], .5, side="right"))
    ref = (s["i_bat_a"][i] - s["i_dc_conv_a"][i])/.1
    got = (s["vdc_v"][i+1]-s["vdc_v"][i])/(s["t"][i+1]-s["t"][i])
    checar(abs(got-ref) <= max(.5, .02*abs(ref)), f"derivada evento {got} != {ref}")
    return f"dVDC/dt num={got:.3g}; KCL={ref:.3g} V/s"


def D3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Afundamento inicial reduz monotonicamente para C, 2C e 4C."""
    quedas=[]
    for c in (.05,.10,.20):
        _,s,_,_=caso(amb,cache,f"D3_{c}",["--modo","ilhado","--evento","carga","--t-step",.5,"--t-end",.7,"--dt-out",1e-4,*args_cc(cdc=c)])
        i0=int(np.searchsorted(s["t"],.5)); i1=int(np.searchsorted(s["t"],.55))
        quedas.append(float(s["vdc_v"][i0]-np.min(s["vdc_v"][i0:i1+1])))
    checar(quedas[0] > quedas[1] > quedas[2] >= 0, f"escala Cdc incoerente: {quedas}")
    return f"quedas C/2C/4C={quedas} V"


def D4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Maior R0 aumenta perda e queda de VDC em descarga."""
    vals=[]
    for r in (.02,.05,.10):
        js,s,_,_=caso(amb,cache,f"D4_{r}",["--modo","ilhado","--evento","carga","--t-step",.5,"--t-end",1,"--dt-out",2e-4,*args_cc(r0=r)])
        vals.append((float(np.min(s["vdc_v"])), integral_trapezio(s["p_perda_r0_w"],s["t"])))
    checar(vals[0][0] > vals[1][0] > vals[2][0], f"queda não cresce com R0: {vals}")
    checar(vals[0][1] < vals[1][1] < vals[2][1], f"perda não cresce com R0: {vals}")
    return f"Vmin/perdas={vals}"


def D5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Integral de Pcap concorda com 0,5 C delta(VDC²)."""
    _,s,_,_=caso_degrau_carga(amb,cache)
    integ=integral_trapezio(s["p_cap_w"],s["t"])
    delta=float(s["e_cap_j"][-1]-s["e_cap_j"][0])
    checar(abs(integ-delta)<=max(TOL["energia_abs_j"],TOL["energia_rel"]*escala(integ,delta)), f"energia capacitor: {integ} != {delta}")
    return f"integral={integ:.6g} J; delta={delta:.6g} J"


def D6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Balanço energético integrado fecha dentro do limite congelado."""
    _,s,_,_=caso_degrau_carga(amb,cache)
    e=reconstruir_energia(s); lim=limite_energia(e["residuo_energia_abs_j"],*e.values())
    checar(e["residuo_energia_abs_j"]<=lim, f"resíduo energético={e['residuo_energia_abs_j']} J > {lim} J")
    return f"resíduo={e['residuo_energia_abs_j']:.6g} J; limite={lim:.6g} J"


def D7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Degrau de carga ilhado afunda VDC e termina classificado."""
    js,s,_,_=caso_degrau_carga(amb,cache); b=bloco_cc(js,True)
    i=int(np.searchsorted(s["t"],.5)); checar(np.min(s["vdc_v"][i:])<s["vdc_v"][i],"VDC não afundou")
    checar(b["estado_final"] in ("normal","colapso_vdc","fora_dominio_ocv"),"estado final desconhecido")
    return f"VDC={s['vdc_v'][i]:.3f}->{np.min(s['vdc_v'][i:]):.3f} V; {b['estado_final']}"


def D8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Degraus positivo e negativo de Pref produzem sinais CC coerentes."""
    saida=[]
    for dp in (.10,-.10):
        _,s,_,_=caso(amb,cache,f"D8_{dp}",[*args_rede(),"--evento","pref","--d-pref",dp,"--t-step",.5,"--t-end",1.2,"--dt-out",2e-4,*args_cc()])
        i=int(np.searchsorted(s["t"],.55)); j=int(np.searchsorted(s["t"],.8)); med=float(np.mean(s["p_conv_w"][i:j]))
        saida.append(med)
    checar(saida[0]>saida[1],f"Pref não refletiu em Pconv: {saida}")
    return f"Pconv médio alto/baixo={saida} W"


def D9(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Colapso deliberado termina limpidamente, sem NaN ou traceback."""
    js,s,_,_=amb.sim("D9",["--modo","ilhado","--evento","carga","--p2",250000,"--t-step",.2,"--t-end",2,"--dt-out",2e-4,*args_cc(cdc=.001,r0=.2),"--vdc-piso-numerico-v",500])
    b=bloco_cc(js,True); checar(b["estado_final"] in ("colapso_vdc","fora_dominio_ocv"),f"caso inviável não terminou: {b['estado_final']}")
    checar(all(np.all(np.isfinite(s[k])) for k in SERIES_CC_CONTINUAS),"NaN/Inf antes do término")
    checar(b.get("tempo_terminal_s") is not None and b.get("causa_terminal"),"causa/tempo terminal ausente")
    return f"{b['estado_final']} em {b['tempo_terminal_s']} s"


def D10(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Horizonte longo mostra evolução coerente de SOC e OCV."""
    _,s,_,_=amb.sim("D10",["--modo","ilhado","--evento","nenhum","--t-end",20,"--dt-out",.01,*args_cc(capacidade=2)])
    checar(s["i_bat_a"].mean()>0 and s["soc_pu"][-1]<s["soc_pu"][0],"SOC não caiu em descarga")
    checar(s["ocv_v"][-1]<=s["ocv_v"][0],"OCV não acompanhou SOC")
    return f"delta SOC={s['soc_pu'][-1]-s['soc_pu'][0]:.6g} pu"


# Bateria B1-B6

def B1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Potência CA positiva sustentada descarrega a bateria."""
    _,s,_,_=caso_cc_nominal(amb,cache)
    checar(np.mean(s["i_bat_a"])>0 and s["soc_pu"][-1]<s["soc_pu"][0],"descarga incoerente")
    return "Ibat>0 e SOC decrescente"


def B2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Potência de ponte negativa sustentada carrega a bateria."""
    js,s,_,_=amb.sim("B2",[*args_rede(),"--p1",1000,"--q1",0,"--pref",-.15,"--evento","nenhum","--t-end",1,"--dt-out",2e-4,*args_cc()])
    checar(np.mean(s["p_conv_w"][-100:])<0 and np.mean(s["i_bat_a"][-100:])<0,"regeneração não estabelecida")
    checar(s["soc_pu"][-1]>s["soc_pu"][0],"SOC não aumentou")
    return "Pconv<0, Ibat<0 e SOC crescente"


def B3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Cruzamento descarga-carga não cria descontinuidade espúria."""
    _,s,_,_=amb.sim("B3",[*args_rede(),"--p1",1000,"--q1",0,"--pref",.10,"--evento","pref","--d-pref",-.25,"--t-step",.5,"--t-end",1.2,"--dt-out",1e-4,*args_cc()])
    checar(np.any(s["p_conv_w"]>0) and np.any(s["p_conv_w"]<0),"Pconv não cruzou zero")
    checar(saltos_em(s["t"],s["vdc_v"],.5)<=TOL["continuidade_vdc_v"],"salto em VDC no cruzamento")
    return "Pconv cruzou zero com VDC contínua"


def B4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Limite de corrente é diagnóstico e não recorta Ibat."""
    js,s,_,_=amb.sim("B4",args_curto(*args_cc("--i-desc-max-a",10)))
    b=bloco_cc(js,True); checar(np.max(s["i_bat_a"])>10,"corrente não excedeu limite de teste")
    checar(np.any(s["i_bat_fora_limite"]>0),"flag de corrente ausente")
    checar(b["limites"]["corrente_descarga"],"JSON não registrou limite")
    return f"Ibat pico={np.max(s['i_bat_a']):.3f} A sem clamp"


def B5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Limite operacional de SOC gera flag sem BMS implícito."""
    js,s,_,_=amb.sim("B5",["--modo","ilhado","--evento","nenhum","--t-end",3,"--dt-out",1e-3,*args_cc(capacidade=.2,soc=.101)])
    b=bloco_cc(js,True); checar(np.any(s["soc_fora_faixa"]>0),"flag SOC ausente")
    checar(np.min(s["soc_pu"])<.10 and b["limites"]["soc"],"SOC foi recortado ou não reportado")
    return f"SOC mínimo={np.min(s['soc_pu']):.6f} sem clamp"


def B6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Perda R0 no JSON coincide com integral independente do NPZ."""
    js,s,_,_=caso_degrau_carga(amb,cache); b=bloco_cc(js,True)
    ref=integral_trapezio(s["p_perda_r0_w"],s["t"])
    checar(abs(float(b["perda_r0_j"])-ref)<=max(.5,.005*escala(ref)),"perda R0 JSON/NPZ diverge")
    return f"perda R0={ref:.6g} J"


# Acoplamento e modulação C1-C6

def C1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Com headroom amplo, modulação não altera a tensão pós-Rv."""
    _,s,_,_=amb.sim("C1",args_curto(*args_cc(mmax=2.0)))
    checar(not np.any(s["modulacao_saturada"]),"saturação inesperada")
    checar(np.max(np.abs(s["e_aplicada_pu"]-s["e_pos_rv_pu"]))<=TOL["serie_pu"],"clamp alterou trajetória")
    return "trajetória sem saturação preservada"


def C2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fronteira SVPWM é inclusiva no comando e exclusiva na flag."""
    fn=exigir_funcao(amb.modulo(),"limitar_modulacao"); v=math.sqrt(2)*380
    igual=fn(1+0j,v,1,380); fora=fn((1+2e-12)+0j,v,1,380)
    checar(not igual["modulacao_saturada"] and fora["modulacao_saturada"],"fronteira SVPWM divergente")
    checar(abs(igual["e_aplicada_pu"]-1)<=1e-12 and abs(fora["e_aplicada_pu"]-1)<=3e-12,"limite aplicado incorreto")
    return "igualdade não satura; ponto externo satura"


def C3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Emax acompanha VDC e a saturação é registrada."""
    js,s,_,_=amb.sim("C3",["--modo","ilhado","--evento","carga","--p2",100000,"--t-step",.3,"--t-end",1,"--dt-out",1e-4,*args_cc(cdc=.005,mmax=.75)])
    ref=.75*s["vdc_v"]/(math.sqrt(2)*380)
    checar(np.allclose(s["e_max_mod_pu"],ref,rtol=1e-10,atol=1e-12),"Emax não acompanha VDC")
    checar(np.any(s["modulacao_saturada"]),"saturação esperada ausente")
    checar(bloco_cc(js,True)["modulacao"]["saturou"],"JSON não registrou saturação")
    return f"utilização máxima={np.max(s['m_utilizado_pu']):.4f} pu"


def C4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Ordem comando→Rv→modulação é observável nas séries complexas."""
    _,s,_,_=amb.sim("C4",[*args_rede(20),"--evento","falta_3f","--vg-falta",0,"--t-step",.3,"--t-clear",.4,"--t-end",.8,"--dt-out",1e-4,"--imax-pu",1.1,"--i-on-pu",1.0,"--rv-max-pu",3,*args_cc(mmax=.7)])
    checar(np.any(np.abs(s["e_req_pu"]-s["e_pos_rv_pu"])>1e-8),"efeito de Rv não observado")
    sat=s["modulacao_saturada"]>.5; checar(np.any(sat),"saturação não observada")
    checar(np.all(np.abs(s["e_aplicada_pu"][sat])<=s["e_max_mod_pu"][sat]*(1+1e-10)),"modulação aplicada antes/incorretamente")
    return "tensões requerida, pós-Rv e aplicada obedecem à ordem"


def C5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Pré-sync ativo não acumula integral contra limite físico."""
    js,s,_,_=amb.sim("C5",[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","ativo","--delta-g0-graus",20,"--df-g0-hz",.12,"--t-step",.5,"--t-end",4,"--dt-out",2e-4,*args_cc(mmax=.5)])
    checar(np.any(s["modulacao_saturada"]),"caso não saturou")
    checar(np.all(np.isfinite(s["de_bias_sync_pu"])),"bias contém NaN/Inf")
    checar(np.max(np.abs(s["de_bias_sync_pu"]))<=.100001,"PI ultrapassou limite")
    return f"bias máximo={np.max(np.abs(s['de_bias_sync_pu'])):.6g} pu"


def C6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Remoção da sobrecarga encerra saturação sem salto em VDC."""
    _,s,_,_=amb.sim("C6",[*args_rede(),"--evento","pref","--d-pref",-.3,"--pref",.3,"--t-step",.6,"--t-end",1.5,"--dt-out",1e-4,*args_cc(cdc=.01,mmax=.65)])
    sat=s["modulacao_saturada"]>.5
    checar(np.any(sat) and not bool(sat[-1]),"saturação não cessou")
    checar(saltos_em(s["t"],s["vdc_v"],.6)<=TOL["continuidade_vdc_v"],"salto não físico em VDC")
    return "saturação cessou com estado contínuo"


# Interações X1-X10

def X1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Limitador de corrente v9 mantém teto ou declara Rvmax."""
    js,s,_,_=amb.sim("X1",["--modo","ilhado","--evento","carga","--t-step",.3,"--t-end",1,"--dt-out",1e-4,"--imax-pu",1.1,"--i-on-pu",1.0,"--rv-max-pu",3,*args_cc()])
    lim=js["limitador"]; pico=float(np.max(s["if_env_pu"])); imax=float(lim["imax_pu"])
    checar(bool(lim["teto_rv_atingido"]) or pico<=imax*1.005,"contrato Imax violado sem Rvmax")
    return f"If pico={pico:.4f} pu"


def X2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Falta trifásica mantém energia coerente e estados CC contínuos para três Cdc."""
    vals=[]
    for c in (.02,.05,.10):
        _,s,_,_=amb.sim(f"X2_{c}",[*args_rede(20),"--evento","falta_3f","--vg-falta",0,"--t-step",.3,"--t-clear",.4,"--t-end",1,"--dt-out",1e-4,*args_cc(cdc=c)])
        for te in (.3,.4):
            checar(saltos_em(s["t"],s["vdc_v"],te)<=TOL["continuidade_vdc_v"],"salto VDC na falta")
            checar(saltos_em(s["t"],s["soc_pu"],te)<=TOL["continuidade_soc_pu"],"salto SOC na falta")
        e=reconstruir_energia(s)
        vals.append((e["residuo_energia_abs_j"], limite_energia(e["residuo_energia_abs_j"], *e.values())))
    checar(all(v<=lim for v,lim in vals),"energia da falta não fecha")
    return f"resíduos/limites={vals} J"


def X3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Degrau e rampa de frequência aparecem na energia do lado CC."""
    energ=[]
    for ev in ("freq_degrau","freq_rampa"):
        _,s,_,_=amb.sim("X3_"+ev,[*args_rede(),"--evento",ev,"--df-g",-.2,"--rocof-g",.5,"--t-step",.3,"--t-end",1.2,"--dt-out",2e-4,*args_cc()])
        energ.append(abs(integral_trapezio(s["p_conv_w"]-s["p_conv_w"][0],s["t"])))
    checar(all(e>0 for e in energ),"evento de frequência não refletiu energia")
    return f"energias incrementais={energ} J"


def X4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Maior inércia virtual exige energia/potência BESS não menor no mesmo evento."""
    vals=[]
    for h in (2,5,10):
        _,s,_,_=amb.sim(f"X4_{h}",[*args_rede(),"--H",h,"--evento","freq_rampa","--df-g",-.2,"--rocof-g",.5,"--t-step",.3,"--t-end",1.2,"--dt-out",2e-4,*args_cc()])
        vals.append(float(np.max(np.abs(s["p_conv_w"]-s["p_conv_w"][0]))))
    checar(vals[0]<=vals[1]*(1+.01) and vals[1]<=vals[2]*(1+.01),f"relação H-energia incoerente: {vals}")
    return f"picos incrementais H=2/5/10: {vals} W"


def X5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Casos válidos em SCR 2/5/20 não apresentam instabilidade crescente."""
    out=[]
    for scr in (2,5,20):
        js,s,_,_=amb.sim(f"X5_{scr}",[*args_rede(scr),"--evento","pref","--d-pref",.05,"--t-step",.3,"--t-end",1.5,"--dt-out",2e-4,*args_cc()])
        checar(np.all(np.isfinite(s["vdc_v"])) and bloco_cc(js,True)["estado_final"]=="normal",f"SCR {scr} inválido")
        out.append(float(np.ptp(s["vdc_v"][-max(10,len(s["vdc_v"])//10):])))
    return f"amplitudes finais SCR 2/5/20={out} V"


def X6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fechamento ativo mantém VDC/SOC contínuos e respeita as janelas."""
    js,s,_,_=amb.sim("X6",[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","ativo","--delta-g0-graus",20,"--df-g0-hz",.12,"--t-step",.5,"--t-end",8,"--dt-out",2e-4,*args_cc()])
    sync=js["sincronizacao"]; tc=float(sync["tempo_contato_s"])
    checar(saltos_em(s["t"],s["vdc_v"],tc)<=TOL["continuidade_vdc_v"],"salto VDC no contato")
    checar(saltos_em(s["t"],s["soc_pu"],tc)<=TOL["continuidade_soc_pu"],"salto SOC no contato")
    checar(sync["estado_final"] in ("fechado","fechado_fora_da_janela"),"não fechou")
    return f"contato={tc:.6f} s; estado={sync['estado_final']}"


def X7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Fechamento forçado fora de fase registra pico sem falsa aprovação."""
    js,s,_,_=amb.sim("X7",[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","forcado","--delta-g0-graus",90,"--t-step",.5,"--t-end",1.5,"--dt-out",1e-4,*args_cc()])
    sync=js["sincronizacao"]
    checar(not bool(sync["sync_check_habilitado"]) or not np.any(s["sync_check_ok"]>.5),"falsa aprovação sync-check")
    checar(np.max(np.abs(s["i_bat_a"]))>0 and np.max(np.abs(s["p_conv_w"]))>0,"pico energético não registrado")
    return f"Ibat pico={np.max(np.abs(s['i_bat_a'])):.3f} A"


def X8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Timeout impossível mantém disjuntor aberto enquanto o BESS alimenta a ilha."""
    js,s,_,_=amb.sim("X8",[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","passivo","--delta-g0-graus",90,"--df-g0-hz",0,"--t-step",.3,"--t-sync-timeout-s",.5,"--t-end",1,"--dt-out",2e-4,*args_cc()])
    sync=js["sincronizacao"]
    checar(sync["estado_final"]=="timeout" and not np.any(s["disjuntor_fechado"]>.5),"timeout/disjuntor incorreto")
    checar(np.mean(s["i_bat_a"])>0,"BESS não alimentou ilha")
    return "timeout com disjuntor aberto e descarga ativa"


def X9(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Headroom insuficiente não produz janela falsa no sync-check."""
    js,s,_,_=amb.sim("X9",[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","ativo","--delta-g0-graus",20,"--df-g0-hz",.12,"--t-step",.3,"--t-sync-timeout-s",1,"--t-end",1.5,"--dt-out",2e-4,*args_cc(mmax=.45)])
    sync=js["sincronizacao"]
    checar(np.any(s["modulacao_saturada"]),"caso não saturou")
    checar(sync["estado_final"] in ("timeout","bloqueado"),"headroom insuficiente recebeu aprovação falsa")
    return f"estado={sync['estado_final']}"


def X10(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Após fechamento, biases zeram e o balanço CC permanece válido."""
    js,s,_,_=amb.sim("X10",[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","ativo","--delta-g0-graus",20,"--df-g0-hz",.12,"--t-step",.5,"--t-end",8,"--dt-out",2e-4,*args_cc()])
    sync=js["sincronizacao"]; checar(sync["tempo_contato_s"] is not None,"sem contato")
    checar(abs(float(s["df_bias_sync_hz"][-1]))<=1e-6 and abs(float(s["de_bias_sync_pu"][-1]))<=1e-6,"bias final não zerou")
    e=reconstruir_energia(s); checar(e["residuo_energia_abs_j"]<=limite_energia(e["residuo_energia_abs_j"],*e.values()),"balanço CC pós-fechamento falhou")
    return f"bias final zerado; resíduo={e['residuo_energia_abs_j']:.6g} J"


# Engenharia e produto E1-E10

def E1(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Refino de max_step/rtol/atol por estado preserva classe, picos e energia."""
    base=["--modo","ilhado","--evento","carga","--t-step",.3,"--t-end",1,"--dt-out",2e-4,*args_cc()]
    j1,s1,_,_=amb.sim("E1a",[*base,"--max-step",1e-3,"--rtol",1e-8])
    j2,s2,_,_=amb.sim("E1b",[*base,"--max-step",5e-4,"--rtol",2.5e-9])
    checar(bloco_cc(j1,True)["estado_final"]==bloco_cc(j2,True)["estado_final"],"classe mudou")
    for k in ("vdc_v","i_bat_a","p_conv_w"):
        a,b=float(np.max(np.abs(s1[k]))),float(np.max(np.abs(s2[k])))
        checar(rel(a,b)<=TOL["convergencia_frac"],f"{k} não convergiu: {a}/{b}")
    return "classe e picos convergentes"


def E2(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Refino de dt_out preserva métricas integrais e picos."""
    base=["--modo","ilhado","--evento","carga","--t-step",.3,"--t-end",1,*args_cc()]
    _,s1,_,_=amb.sim("E2a",[*base,"--dt-out",2e-4]); _,s2,_,_=amb.sim("E2b",[*base,"--dt-out",1e-4])
    for k in ("p_conv_w","p_perda_r0_w"):
        checar(rel(integral_trapezio(s1[k],s1["t"]),integral_trapezio(s2[k],s2["t"]))<=TOL["convergencia_frac"],f"integral {k} não convergiu")
    return "integrais convergentes com dt_out"


def E3(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Refino de dt_rele preserva decisão e tempo dentro de uma amostra."""
    base=[*args_rede(),"--evento","fechamento","--disjuntor-inicial","aberto","--estrategia-sync","ativo","--delta-g0-graus",20,"--df-g0-hz",.12,"--t-step",.5,"--t-end",8,"--dt-out",2e-4,*args_cc()]
    j1,_,_,_=amb.sim("E3a",[*base,"--dt-rele-s",.001]); j2,_,_,_=amb.sim("E3b",[*base,"--dt-rele-s",.0005])
    a,b=j1["sincronizacao"],j2["sincronizacao"]
    checar(a["estado_final"]==b["estado_final"],"decisão do relé mudou")
    checar(abs(float(a["tempo_contato_s"])-float(b["tempo_contato_s"]))<=.001000000001,"tempo variou mais de uma amostra")
    return "decisão e tempo do relé convergentes"


def E4(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Repetições no mesmo ambiente são determinísticas."""
    a=args_curto(*args_cc()); j1,s1,_,_=amb.sim("E4a",a); j2,s2,_,_=amb.sim("E4b",a)
    checar(j1["lado_cc"]==j2["lado_cc"],"JSON CC não determinístico")
    for k in SERIES_CC_CONTINUAS+SERIES_CC_FLAGS:
        checar(np.array_equal(s1[k],s2[k],equal_nan=True),f"série não determinística: {k}")
    return "JSON e séries CC determinísticos"


def E5(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Matriz de enums, listas, unidades e combinações inválidas falha sem traceback."""
    casos=[["--modelo-cc","x"],["--m-max-pu",0],["--capacidade-ah",0],["--ocv-soc-pu","a,b"],["--vdc-piso-numerico-v",900]]
    for n,x in enumerate(casos): texto_erro_limpo(amb.sim(f"E5_{n}",args_curto(*args_cc(*x)),esperar_erro=True))
    return f"{len(casos)} entradas inválidas rejeitadas limpidamente"


def E6(amb: Ambiente, cache: dict[str, Any]) -> str:
    """JSON recompõe as métricas do NPZ nas tolerâncias congeladas."""
    js,s,_,_=caso_degrau_carga(amb,cache); b=bloco_cc(js,True); e=reconstruir_energia(s)
    for k,v in e.items():
        checar(k in b,f"JSON sem {k}")
        checar(abs(float(b[k])-v)<=max(.5,.005*escala(v)),f"{k} JSON/NPZ divergente")
    parametros=js.get("parametros",{})
    vll=float(parametros.get("vll",parametros.get("vll_base_v",380.0)))
    base=float(parametros.get("vbase_pico_v",math.sqrt(2.0)*vll/math.sqrt(3.0)))
    p=1.5*np.real((s["e_aplicada_pu"]*base)*np.conj(s["if_complex_a"]))
    checar(np.max(np.abs(p-s["p_conv_w"]))<=TOL["residuo_potencia_rel"]*escala(*s["p_conv_w"]),"Pconv não recomposto")
    return "energias e Pconv recompostos do NPZ"


def E7(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Gerador canônico cria e executa os dez exemplos v10 sem edição."""
    checar(os.path.isfile(amb.gerador),f"gerador não encontrado: {amb.gerador}")
    exemplos=("caso_cc_ideal_regressao_v10.yaml","caso_bess_degrau_carga_v10.yaml","caso_bess_cdc_baixo_v10.yaml","caso_bess_cdc_alto_v10.yaml","caso_bess_regeneracao_v10.yaml","caso_bess_saturacao_modulacao_v10.yaml","caso_bess_falta_3f_v10.yaml","caso_bess_freq_rampa_v10.yaml","caso_bess_fechamento_ativo_v10.yaml","caso_bess_inviavel_v10.yaml")
    out=os.path.join(amb.workdir,"exemplos"); os.makedirs(out)
    pr=subprocess.run([sys.executable,amb.gerador,"--gerar-exemplos",out],capture_output=True,text=True,timeout=2400,encoding="utf-8",errors="replace")
    checar(pr.returncode==0,"gerador falhou:\n"+((pr.stdout or "")+(pr.stderr or ""))[-3000:])
    for nome in exemplos:
        p=os.path.join(out,nome); checar(os.path.isfile(p),f"exemplo ausente: {nome}")
        sim=amb.sim("E7_"+nome[:-5],config=p,esperar_erro=nome.endswith("inviavel_v10.yaml"))
        if nome.endswith("inviavel_v10.yaml"): texto_erro_limpo(sim)
    return "dez YAMLs gerados e executados sem edição"


def E8(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Mediana de três execuções respeita metas de desempenho congeladas."""
    a=args_curto(); t9=[]; ti=[]; tt=[]
    for n in range(3):
        t9.append(amb.sim(f"E8_v9_{n}",a,script=amb.v9)[2])
        ti.append(amb.sim(f"E8_i_{n}",[*a,"--modelo-cc","ideal"])[2])
        tt.append(amb.sim(f"E8_t_{n}",args_curto(*args_cc()))[2])
    m9,mi,mt=map(float,map(np.median,(t9,ti,tt)))
    checar(mi<=TOL["desempenho_ideal_razao"]*m9,f"v10 ideal={mi/m9:.3f}x v9")
    checar(mt<=TOL["desempenho_thevenin_razao"]*mi,f"thevenin={mt/mi:.3f}x ideal")
    return f"ideal/v9={mi/m9:.3f}; thevenin/ideal={mt/mi:.3f}"


def E9(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Executa em caminho com espaços/Unicode e preserva UTF-8/NumPy suportado."""
    pasta=os.path.join(amb.workdir,"espaço Δ v10"); os.makedirs(pasta)
    copia=os.path.join(pasta,"simulador v10.py"); shutil.copy2(amb.script,copia)
    js,s,_,_=amb.sim("E9_unicode",args_curto(*args_cc()),script=copia)
    checar(js.get("versao")=="v10" and np.__version__,"portabilidade falhou")
    checar_series_cc(s,ideal=False)
    return f"UTF-8/Unicode aprovado com NumPy {np.__version__}"


def E10(amb: Ambiente, cache: dict[str, Any]) -> str:
    """Alteração da suíte torna a integridade explicitamente inválida."""
    origem=os.path.abspath(__file__); copia=os.path.join(amb.workdir,"testes_v10_alterado.py")
    shutil.copy2(origem,copia); open(copia,"a",encoding="utf-8").write("\n# alteração controlada E10\n")
    ref=open(os.path.join(os.path.dirname(origem),"testes_v10.sha256"),encoding="utf-8").read().split()[0]
    checar(sha256(copia)!=ref,"alteração não mudou hash")
    return "cópia alterada perde correspondência SHA-256"


TESTES = [
    ("R1",R1),("R2",R2),("R3",R3),("R4",R4),("R5",R5),("R6",R6),
    ("M1",M1),("M2",M2),("M3",M3),("M4",M4),("M5",M5),("M6",M6),("M7",M7),("M8",M8),("M9",M9),("M10",M10),
    ("I1",I1),("I2",I2),("I3",I3),("I4",I4),("I5",I5),("I6",I6),("I7",I7),("I8",I8),
    ("D1",D1),("D2",D2),("D3",D3),("D4",D4),("D5",D5),("D6",D6),("D7",D7),("D8",D8),("D9",D9),("D10",D10),
    ("B1",B1),("B2",B2),("B3",B3),("B4",B4),("B5",B5),("B6",B6),
    ("C1",C1),("C2",C2),("C3",C3),("C4",C4),("C5",C5),("C6",C6),
    ("X1",X1),("X2",X2),("X3",X3),("X4",X4),("X5",X5),("X6",X6),("X7",X7),("X8",X8),("X9",X9),("X10",X10),
    ("E1",E1),("E2",E2),("E3",E3),("E4",E4),("E5",E5),("E6",E6),("E7",E7),("E8",E8),("E9",E9),("E10",E10),
]


def autoteste() -> bool:
    ok=True
    def item(cond: Any, msg: str) -> None:
        nonlocal ok
        try: valor=bool(cond)
        except Exception: valor=False
        ok &= valor; print(f"  [{'ok' if valor else 'FALHA'}] {msg}")

    ids=[n for n,_ in TESTES]; esperado=sum(GRUPOS.values(),[])
    item(len(ids)==66 and len(set(ids))==66 and ids==esperado,"matriz congelada: 66 IDs distintos e ordenados")
    item(all(callable(fn) and (fn.__doc__ or "").strip() for _,fn in TESTES),"todos os testes têm função e contrato descritivo")
    item(CONTRATO["indices"]=={"ideal_comum":list(range(15)),"ideal_sync":list(range(17)),"thevenin_comum":list(range(15))+[17,18],"thevenin_sync":list(range(19))},"mapas de estados 15/17/17/19")

    item(np.array_equal(ocv_oraculo(OCV_SOC),OCV_V),"OCV exata nos nós")
    q=np.array([.125,.375,.625,.875]); item(np.allclose(ocv_oraculo(q),[670,690,710,730],rtol=0,atol=1e-12),"OCV linear nos pontos médios")
    rejeitou=0
    for x,y,q0 in [(OCV_SOC,OCV_V,-.1),([0,.5,.5,1],OCV_V,.5),([0,1],[700],.5)]:
        try: ocv_oraculo(q0,x,y)
        except ValueError: rejeitou+=1
    item(rejeitou==3,"OCV rejeita domínio/tabelas inválidos sem extrapolar")

    bd=bateria_oraculo(.5,690,.05); bc=bateria_oraculo(.5,710,.05)
    item(bd["i_bat_a"]>0 and bc["i_bat_a"]<0,"convenção bidirecional de Ibat")
    item(perto(bd["p_bat_quim_w"],bd["p_bat_term_w"]+bd["p_perda_r0_w"]) and perto(bc["p_bat_quim_w"],bc["p_bat_term_w"]+bc["p_perda_r0_w"]),"identidade de potência Thévenin")

    d=derivadas_cc_oraculo(690,.5,20e3,.2,100,.05)
    item(perto(d["dsoc_dt_pu_s"],-d["i_bat_a"]/(3600*100)),"oráculo coulômbico")
    item(perto(d["dvdc_dt_v_s"],(d["i_bat_a"]-d["i_dc_conv_a"])/.2),"oráculo KCL do capacitor")
    item(abs(d["residuo_potencia_cc_w"])<=1e-10*escala(d["p_bat_quim_w"],d["p_conv_w"] if "p_conv_w" in d else 20e3),"oráculo de balanço instantâneo")

    pares=[(230+0j,100+0j,34500),(230j,100+0j,0),(230+0j,-50+0j,-17250)]
    item(all(perto(potencia_ponte_oraculo(e,i),p) for e,i,p in pares),"oráculo de potência da ponte e sinais")

    v=math.sqrt(2)*380
    z=modulacao_oraculo(0j,v,1,380); eq=modulacao_oraculo(1+0j,v,1,380); fo=modulacao_oraculo(2*np.exp(1j*.7),v,1,380)
    item(z["k_mod"]==1 and not z["modulacao_saturada"],"SVPWM trata vetor nulo")
    item(not eq["modulacao_saturada"] and perto(eq["e_max_mod_pu"],1),"igualdade SVPWM não satura")
    item(fo["modulacao_saturada"] and perto(abs(fo["e_aplicada_pu"]),1) and abs(np.angle(fo["e_aplicada_pu"]/(2*np.exp(1j*.7))))<1e-12,"SVPWM radial preserva fase")

    for p in (50e3,0,-20e3):
        e=equilibrio_oraculo(p,.6,.05)
        item(perto(e["vdc_v"]*(e["ocv_v"]-e["vdc_v"])/.05,p),f"equilíbrio Thévenin P={p:g} W")
    pmax=float(ocv_oraculo(.5))**2/(4*.05)
    try: equilibrio_oraculo(pmax*(1+1e-12),.5,.05); impossivel=False
    except ValueError: impossivel=True
    item(impossivel,"discriminante negativo rejeita potência impossível")

    t=np.linspace(0,2,2001); p=np.full_like(t,1234.5)
    item(perto(integral_trapezio(p,t),2469.0,reltol=1e-12),"integração trapezoidal independente")
    item(perto(TOL["continuidade_vdc_v"],10*TOL["vdc_atol_v"],reltol=1e-15) and perto(TOL["continuidade_soc_pu"],10*TOL["soc_atol_pu"],reltol=1e-15),"continuidade = 10 vezes atol por estado")
    item(TOL["desempenho_ideal_razao"]==1.05 and TOL["desempenho_thevenin_razao"]==2.0,"metas de desempenho 1,05x/2,0x")
    item(HASH_SCRIPT_V9=="077b3e17509cb813bc230326bb7dc4a523050182dd4ad1777f6419921bfad713" and HASH_TESTES_V9=="0f0d6265eff12b27a8e44d533643fe8cb5ea1085649b55c3ac713b67030a9b8d","hashes canônicos da baseline v9")
    item(len(SERIES_CC_CONTINUAS)==19 and len(SERIES_CC_FLAGS)==5 and set(SERIES_COMPLEXAS)<=set(SERIES_CC_CONTINUAS),"contrato NPZ completo e tipos complexos explícitos")
    item(set(CLI_V10)=={"--"+x.replace("_","-") for x in ("modelo_cc","c_dc_f","vdc_inicial_v","permitir_desequilibrio_inicial","vdc_min_oper_v","vdc_max_oper_v","m_max_pu","vdc_piso_numerico_v","capacidade_ah","soc_inicial","r0_ohm","ocv_soc_pu","ocv_v","soc_min_oper_pu","soc_max_oper_pu","i_desc_max_a","i_carga_max_a")},"nomes CLI congelados")
    return bool(ok)


def verificar_integridade() -> tuple[str, str]:
    atual=sha256(os.path.abspath(__file__))
    arq=os.path.join(os.path.dirname(os.path.abspath(__file__)),"testes_v10.sha256")
    if not os.path.isfile(arq): return atual,"SEM REFERÊNCIA (testes_v10.sha256 ausente)"
    esperado=open(arq,encoding="utf-8").read().split()[0].strip().lower()
    return atual,("ÍNTEGRO" if esperado==atual else "ALTERADO — resultados sem validade")


def escrever_relatorios(resultados: list[dict[str,Any]], estado: str, hash_atual: str,
                         script: str, total_s: float, integral: bool) -> None:
    n_ok=sum(r["status"]=="APROVADO" for r in resultados)
    relatorio={"versao_suite":VERSAO_TESTES,"data_congelamento":DATA_CONGELAMENTO,
        "sha256_testes":hash_atual,"integridade":estado,"execucao_integral":integral,
        "script":script,"sha256_script":sha256(script),"aprovados":n_ok,
        "total":len(resultados),"matriz_total":66,"contrato":CONTRATO,
        "tolerancias":TOL,"resultados":resultados}
    with open("relatorio_testes_v10.json","w",encoding="utf-8") as fh:
        json.dump(relatorio,fh,ensure_ascii=False,indent=2)
    with open("relatorio_testes_v10.md","w",encoding="utf-8") as fh:
        fh.write(f"# Relatório de aceitação — v10\n\n- suíte: `{VERSAO_TESTES}`\n"
                 f"- testes_v10.py SHA-256: `{hash_atual}` ({estado})\n"
                 f"- implementação: `{script}` SHA-256 `{relatorio['sha256_script']}`\n"
                 f"- execução integral: **{'sim' if integral else 'não'}**\n"
                 f"- resultado: **{n_ok}/{len(resultados)} aprovados** em {total_s:.1f} s\n\n"
                 "| Teste | Status | Tempo (s) | Contrato | Detalhe |\n|---|---|---:|---|---|\n")
        for r in resultados:
            det=str(r["detalhe"]).replace("\n"," ").replace("|","/")[:800]
            desc=str(r["descricao"]).replace("|","/")
            fh.write(f"| {r['teste']} | {r['status']} | {r['tempo_s']} | {desc} | {det} |\n")


def main() -> None:
    for st in (sys.stdout,sys.stderr):
        try: st.reconfigure(errors="replace")
        except (AttributeError,ValueError): pass
    aqui=os.path.dirname(os.path.abspath(__file__))
    ap=argparse.ArgumentParser(description="Contrato executável congelado da v10 (66 testes).")
    ap.add_argument("--script",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v10.py"))
    ap.add_argument("--v9",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v9.py"))
    ap.add_argument("--testes-v9",default=os.path.join(aqui,"testes_v9.py"))
    ap.add_argument("--v8",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v8.py"))
    ap.add_argument("--testes-v8",default=os.path.join(aqui,"testes_v8_1.py"))
    ap.add_argument("--gerador",default=os.path.join(aqui,"gerar_caso_vsg_v10.py"))
    ap.add_argument("--v7",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v7.py"))
    ap.add_argument("--v6",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v6.py"))
    ap.add_argument("--v5",default=os.path.join(aqui,"vsg_2a_ordem_degrau_carga_v5.py"))
    ap.add_argument("--testes-v7",default=os.path.join(aqui,"testes_v7.py"))
    ap.add_argument("--so",default=None,help="IDs/grupos separados por vírgula; ex.: M,I2,C")
    ap.add_argument("--autoteste",action="store_true",help="valida contrato/oráculos sem importar/executar v10")
    ap.add_argument("--manter",action="store_true",help="mantém pasta temporária")
    a=ap.parse_args()

    hash_atual,estado=verificar_integridade()
    print(f"testes_v10.py  SHA-256 {hash_atual}  [{estado}]")
    if a.autoteste:
        print("Autoteste independente do contrato e dos oráculos:")
        aprovado=autoteste()
        print("AUTOTESTE", "APROVADO" if aprovado and estado=="ÍNTEGRO" else "REPROVADO")
        raise SystemExit(0 if aprovado and estado=="ÍNTEGRO" else 1)

    mapa=dict(TESTES); selecionados=[]
    tokens=[x.strip().upper() for x in a.so.split(",") if x.strip()] if a.so else list(mapa)
    for tok in tokens:
        if tok in GRUPOS: selecionados.extend(GRUPOS[tok])
        elif tok in mapa: selecionados.append(tok)
        else: raise SystemExit(f"Erro: teste/grupo desconhecido: {tok}")
    selecionados=list(dict.fromkeys(selecionados)); integral=selecionados==list(mapa)

    obrigatorios=[(a.script,"implementação v10")]
    if any(x in selecionados for x in ("R1","R2","R3","R6","E8")):
        obrigatorios += [(a.v9,"simulador v9"),(a.testes_v9,"testes_v9.py")]
    if any(x in selecionados for x in ("R1","R2")):
        obrigatorios += [(a.v8,"simulador v8"),(a.testes_v8,"testes_v8_1.py"),
            (a.v7,"simulador v7"),(a.v6,"simulador v6"),(a.v5,"simulador v5"),
            (a.testes_v7,"testes_v7.py"),(a.gerador,"gerador v10")]
    if "E7" in selecionados: obrigatorios.append((a.gerador,"gerador v10"))
    for caminho,nome in obrigatorios:
        if not os.path.isfile(caminho): raise SystemExit(f"Erro: {nome} não encontrado: {caminho}")
    if os.path.isfile(a.v9): checar(sha256(a.v9)==HASH_SCRIPT_V9,"simulador v9 não é a baseline canônica")
    if os.path.isfile(a.testes_v9): checar(sha256(a.testes_v9)==HASH_TESTES_V9,"testes v9 não são a baseline canônica")

    work=tempfile.mkdtemp(prefix="testes_v10_")
    amb=Ambiente(a.script,a.v9,a.testes_v9,a.v8,a.testes_v8,a.gerador,a.v7,a.v6,a.v5,a.testes_v7,work)
    cache:dict[str,Any]={}; resultados=[]; inicio=time.perf_counter()
    try:
        for nome,fn in TESTES:
            if nome not in selecionados: continue
            t0=time.perf_counter()
            try: detalhe,status=fn(amb,cache),"APROVADO"
            except AssertionError as exc: detalhe,status=str(exc),"REPROVADO"
            except Exception as exc: detalhe,status=f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=3)}","ERRO"
            dt=time.perf_counter()-t0
            r={"teste":nome,"status":status,"tempo_s":round(dt,3),
               "descricao":(fn.__doc__ or "").strip(),"detalhe":detalhe}
            resultados.append(r)
            print(f"[{status:^9}] {nome:<3} ({dt:7.2f} s) {r['descricao']}\n            {detalhe}")
        total=time.perf_counter()-inicio
        escrever_relatorios(resultados,estado,hash_atual,a.script,total,integral)
        n_ok=sum(r["status"]=="APROVADO" for r in resultados)
        print(f"\nResultado: {n_ok}/{len(resultados)} aprovados em {total:.1f} s | integridade: {estado}")
        if not integral: print("AVISO: execução seletiva; não constitui aceite de release.")
        raise SystemExit(0 if n_ok==len(resultados) and estado=="ÍNTEGRO" else 1)
    finally:
        if a.manter: print(f"Pasta de trabalho: {work}")
        else: shutil.rmtree(work,ignore_errors=True)


if __name__ == "__main__":
    main()

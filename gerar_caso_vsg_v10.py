#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gerador de casos para vsg_2a_ordem_degrau_carga_v10.py (versão v10).

A partir de dados de projeto em unidades de engenharia, gera um YAML pronto para o
simulador, aplicando:

  1) Cargas (inicial e degrau) em W, kW, MW, GW ou pu + fator de potência (indutivo).
  2) Filtro R_f, L_f, C_f pela regra de projeto (ripple, reativo de C_f, X/R).
  3) Constante de inércia pelo critério de RoCoF, com verificação de H_max.
  4) Ponto de operação inicial (P_ref, Q_ref, E_ref).
  5) [v6] Modo "rede": SCR, X/R da rede, ressonância LCL verificada nos extremos de SCR,
     e ganho de washout D_w dimensionado pelo MODELO LINEARIZADO COMPLETO do simulador
     para que o modo eletromecânico tenha ζ ≥ ζ_alvo no SCR mais forte da faixa.

  6) [v7] Eventos de frequência da rede (freq_degrau, freq_rampa): grava df_g e rocof_g,
     relata a inércia efetiva H_eff = H + D_w·T_w/2 e o pico de potência previsto
        P_pico = P1 + (|R|/f0)·[2H + D_w·T_w·(1 − e^(−T_r/T_w))] + D_p·|Δf|/f0,
     com aviso "(!)" quando P_pico > --pmax-pcs.
  7) [v8] Limitação de corrente, anti-windup e falta_3f: grava os parâmetros novos,
     dimensiona R_v,max de forma conservadora e alerta para recorte da resposta inercial.
  8) [v9] Disjuntor, sync-check e pré-sincronizador: grava as janelas do relé 25,
      tempo de qualificação, tempo mecânico e ganhos PI de ângulo/frequência e tensão.
v6.0.1: saída de console segura em qualquer codificação (Windows cp1252 incluso);
        o YAML continua gravado em UTF-8.
v6.0.2: strings do YAML entre aspas simples (caminhos Windows com '\\' eram lidos
        como sequências de escape).

No modo "rede", P_ref e Q_ref são gravados como null: o simulador parte do equilíbrio
com troca nula de potência com a rede (a carga local é suprida pelo VSG).

Unidades: potência W/kW/MW/GW/pu; aparente VA/kVA/MVA/GVA; tensão V/kV/pu; freq. Hz/kHz.

Uso:
  python gerar_caso_vsg.py                                  # ilhado, padrões
  python gerar_caso_vsg.py --modo rede --scr 5 --scr-faixa 2,20 --zeta-alvo 0.3
  python gerar_caso_vsg.py --modo rede --evento freq_rampa --df-rede -0.5 --rocof-rede 0.5
  python gerar_caso_vsg_v9.py --modo rede --evento sincronizacao --d-fase 20 --df-rede 0.15
  python vsg_2a_ordem_degrau_carga_v9.py --config meu_caso.yaml
"""

import argparse
import importlib.util
import math
import os
import re
import sys

import numpy as np

SIMULADOR = "vsg_2a_ordem_degrau_carga_v10.py"
EVENTOS = ("nenhum", "carga", "pref", "fase", "freq_degrau", "freq_rampa", "falta_3f", "fechamento", "sincronizacao")
EVENTOS_FREQ = ("freq_degrau", "freq_rampa")


def console_seguro():
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


# ---------------------------------------------------------------------------------------
# Unidades
# ---------------------------------------------------------------------------------------
_MULT = {"": 1.0, "k": 1e3, "m": 1e6, "g": 1e9}
_NUM = r"^\s*([-+]?\d*[.,]?\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z/]*)\s*$"


def erro(msg):
    sys.exit(f"Erro: {msg}")


def _split(txt, nome):
    m = re.match(_NUM, str(txt))
    if not m:
        erro(f"não foi possível interpretar {nome} = '{txt}'")
    return float(m.group(1).replace(",", ".")), m.group(2).lower()


def parse_power(txt, sn, nome):
    v, u = _split(txt, nome)
    if u == "pu":
        return v * sn
    if u == "":
        return v
    if u.endswith("w") and u[:-1] in _MULT:
        return v * _MULT[u[:-1]]
    erro(f"unidade '{u}' inválida para {nome} (use W, kW, MW, GW ou pu)")


def parse_apparent(txt, nome):
    v, u = _split(txt, nome)
    if u in ("", "va"):
        return v
    if u.endswith("va") and u[:-2] in _MULT:
        return v * _MULT[u[:-2]]
    erro(f"unidade '{u}' inválida para {nome} (use VA, kVA, MVA ou GVA)")


def parse_voltage(txt, nome, base=None):
    v, u = _split(txt, nome)
    if u in ("", "v"):
        return v
    if u == "kv":
        return v * 1e3
    if u == "pu" and base is not None:
        return v * base
    erro(f"unidade '{u}' inválida para {nome} (use V, kV{' ou pu' if base else ''})")


def parse_freq(txt, nome):
    v, u = _split(txt, nome)
    if u in ("", "hz"):
        return v
    if u == "khz":
        return v * 1e3
    erro(f"unidade '{u}' inválida para {nome} (use Hz ou kHz)")


def fmt_p(w):
    for lim, pre in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(w) >= lim:
            return f"{w / lim:.4g} {pre}"
    return f"{w:.4g} "


def carregar_simulador():
    caminho = os.path.join(os.path.dirname(os.path.abspath(__file__)), SIMULADOR)
    if not os.path.isfile(caminho):
        erro(f"{SIMULADOR} não encontrado na pasta do gerador (necessário para a análise modal)")
    spec = importlib.util.spec_from_file_location("vsg_sim", caminho)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------------------
# Regime fasorial ilhado (convenção do simulador)
# ---------------------------------------------------------------------------------------
def rl_from_pq(P, Q, Vll, w):
    S2 = P ** 2 + Q ** 2
    return Vll ** 2 * P / S2, (Vll ** 2 * Q / S2) / w


def circuit(E_pu, dw, d, cargas):
    w = d["w0"] * (1 + dw)
    Zf = d["Rf"] + 1j * w * d["Lf"]
    Y = 1j * w * d["Cf"] + sum(1.0 / (R + 1j * w * L) for R, L in cargas)
    e = -1j * math.sqrt(2) * E_pu * d["Vph"]
    Vc = e / (1 + Zf * Y)
    If = (e - Vc) / Zf
    S = 1.5 * Vc * np.conj(If) / d["Sb"]
    return Vc, If, S


def steady_state(d, cargas, Eref, Pref, Qref, Dp, nq):
    E, dw = Eref, 0.0
    for _ in range(1000):
        Vc, If, S = circuit(E, dw, d, cargas)
        E_new = Eref - nq * (S.imag - Qref)
        dw_new = (Pref - S.real) / Dp
        if abs(E_new - E) < 1e-13 and abs(dw_new - dw) < 1e-13:
            break
        E, dw = E_new, dw_new
    return E, dw, Vc, If, S


# ---------------------------------------------------------------------------------------
# Projeto
# ---------------------------------------------------------------------------------------
def projetar(a):
    avisos, r = [], {}
    if a.modo not in ("ilhado", "rede"):
        erro(f"--modo deve ser 'ilhado' ou 'rede' (valor: {a.modo!r})")
    Sn = parse_apparent(a.sn, "--sn")
    Vll = parse_voltage(a.vll, "--vll")
    f0 = float(a.f0)
    w0 = 2 * math.pi * f0
    Vph = Vll / math.sqrt(3)
    Zb = Vll ** 2 / Sn
    Lb, Cb = Zb / w0, 1 / (w0 * Zb)
    Ipk = math.sqrt(2) * Sn / (math.sqrt(3) * Vll)
    r.update(Sn=Sn, Vll=Vll, f0=f0, w0=w0, Vph=Vph, Zb=Zb, Lb=Lb, Cb=Cb, Ipk=Ipk, Sb=Sn)
    if a.imax_pu < 0:
        erro("--imax-pu deve ser >= 0")
    ion = 0.98 * a.imax_pu if a.i_on_pu is None else a.i_on_pu
    rvmax = (3.0 if a.rv_max_pu is None and a.imax_pu > 0 else
             (2.0 if a.rv_max_pu is None else a.rv_max_pu))
    if a.imax_pu > 0 and not 0 <= ion < a.imax_pu:
        erro("--i-on-pu deve satisfazer 0 <= i_on < imax")
    if a.imax_pu > 0 and rvmax <= 0:
        erro("--rv-max-pu deve ser > 0 com o limitador ativo")
    if a.xv_rv < 0 or a.k_aw < 0:
        erro("--xv-rv e --k-aw devem ser >= 0")
    r.update(imax_pu=a.imax_pu, i_on_pu=ion, rv_max_pu=rvmax,
             xv_rv=a.xv_rv, k_aw=a.k_aw)
    if a.modelo_cc not in ("ideal", "thevenin"):
        erro("--modelo-cc deve ser ideal ou thevenin")
    if a.m_max_pu <= 0:
        erro("--m-max-pu deve ser > 0")
    if a.modelo_cc == "thevenin":
        obrig = {"c_dc_f": a.c_dc_f, "capacidade_ah": a.capacidade_ah,
                 "soc_inicial": a.soc_inicial, "r0_ohm": a.r0_ohm,
                 "ocv_soc_pu": a.ocv_soc_pu, "ocv_v": a.ocv_v}
        faltam = [k for k, v_ in obrig.items() if v_ is None]
        if faltam:
            erro("parâmetros obrigatórios do modelo Thévenin ausentes: " + ", ".join(faltam))
    if a.pre_sync not in (0.0, 1.0):
        erro("--pre-sync deve ser 0 ou 1")
    for nome in ("kp_sync_p", "ki_sync_p", "p_sync_max_pu", "kp_sync_v", "ki_sync_v", "e_sync_max_pu"):
        if getattr(a, nome) < 0:
            erro(f"--{nome.replace('_', '-')} deve ser >= 0")
    for nome in ("sync_dv_max_pu", "sync_df_max_hz", "sync_delta_max_deg", "sync_hold_s"):
        if getattr(a, nome) <= 0:
            erro(f"--{nome.replace('_', '-')} deve ser > 0")
    if a.breaker_delay_s < 0:
        erro("--breaker-delay-s deve ser >= 0")
    if a.evento == "fechamento":
        if a.modo != "rede": erro("o evento 'fechamento' exige --modo rede")
        if a.disjuntor_inicial != "aberto": erro("fechamento exige --disjuntor-inicial aberto")
        if a.estrategia_sync not in ("forcado", "passivo", "ativo"):
            erro("--estrategia-sync deve ser forcado, passivo ou ativo")

    def carga(txt, fp, nome):
        P = parse_power(txt, Sn, nome)
        if P <= 0:
            erro(f"{nome} deve ser > 0")
        if not 0 < fp <= 1:
            erro(f"fator de potência de {nome} deve estar em (0, 1]")
        return P, P * math.tan(math.acos(fp))

    fp_deg = a.fp if a.fp_degrau is None else a.fp_degrau
    P1, Q1 = carga(a.carga, a.fp, "--carga")
    P2, Q2 = carga(a.degrau, fp_deg, "--degrau")
    r.update(P1=P1, Q1=Q1, P2=P2, Q2=Q2, fp1=a.fp, fp2=fp_deg)
    if math.hypot(P1 + P2, Q1 + Q2) > Sn:
        avisos.append(f"Carga total após o degrau excede S_n = {fmt_p(Sn)}VA — sem limite de corrente no modelo.")

    # --- Filtro --------------------------------------------------------------------------
    fsw = parse_freq(a.fsw, "--fsw")
    Vdc = (1.3 * math.sqrt(2) * Vll) if str(a.vdc).lower() == "auto" else parse_voltage(a.vdc, "--vdc")
    dI = a.ripple * Ipk
    Lf = Vdc / (6 * fsw * dI)
    Cf = a.qc_max * Cb
    Xf = w0 * Lf
    Rf = Xf / a.xr
    fr = 1 / (2 * math.pi * math.sqrt(Lf * Cf))
    att = 1 / abs(1 - (2 * math.pi * fsw) ** 2 * Lf * Cf)
    Rd = 1 / (3 * 2 * math.pi * fr * Cf)
    r.update(fsw=fsw, Vdc=Vdc, vdc_auto=str(a.vdc).lower() == "auto", dI=dI, Lf=Lf, Cf=Cf, Rf=Rf,
             Xf=Xf, fr=fr, att=att, Rd=Rd, rf_pu=Rf / Zb, xf_pu=Xf / Zb, bc_pu=Cf / Cb)
    if not 10 * f0 < fr < fsw / 2:
        avisos.append(f"Ressonância do filtro LC f_r = {fr:.0f} Hz fora da faixa "
                      f"({10 * f0:.0f} Hz < f_r < {fsw / 2:.0f} Hz). Revise f_sw, ripple ou q_C.")
    if Xf / Zb > 0.15:
        avisos.append(f"L_f = {Xf / Zb:.3f} pu > 0,15 pu: queda de tensão elevada.")

    # --- Inércia -------------------------------------------------------------------------
    dPmax = P2 if a.dpmax is None else parse_power(a.dpmax, Sn, "--dpmax")
    dPmax_pu = dPmax / Sn
    H_min = dPmax_pu * f0 / (2 * a.rocof_max)
    H = math.ceil(H_min * a.margem_h * 10 - 1e-9) / 10
    if a.h_faixa:
        h_lo, h_hi = (float(x) for x in a.h_faixa.split(","))
        if H < h_lo:
            avisos.append(f"H calculado ({H:.1f} s) abaixo da faixa do PCS; usando {h_lo:.1f} s.")
            H = h_lo
        elif H > h_hi:
            avisos.append(f"H necessário ({H:.1f} s) acima da faixa ajustável do PCS ({h_hi:.1f} s): "
                          f"o RoCoF máximo NÃO será atendido com H = {h_hi:.1f} s.")
            H = h_hi
    rocof_rede = a.rocof_rede if a.rocof_rede is not None else a.rocof_max
    P1_pu = P1 / Sn
    H_max = (a.pmax_pcs - P1_pu) * f0 / (2 * rocof_rede)
    if H > H_max:
        avisos.append(f"H = {H:.1f} s > H_max = {H_max:.2f} s: resposta inercial limitada pela corrente.")
    r.update(dPmax=dPmax, dPmax_pu=dPmax_pu, H_min=H_min, H=H, H_max=H_max,
             rocof_rede=rocof_rede, rocof_prev=dPmax_pu * f0 / (2 * H))

    # --- Ponto de operação ilhado ---------------------------------------------------------
    Dp = 1 / a.mp
    RL1, LL1 = rl_from_pq(P1, Q1, Vll, w0)
    RL2, LL2 = rl_from_pq(P2, Q2, Vll, w0)
    r.update(RL1=RL1, LL1=LL1, RL2=RL2, LL2=LL2, Dp=Dp)
    if str(a.v_inversor).lower() == "auto":
        Vc1, _, _ = circuit(1.0, 0.0, r, [(RL1, LL1)])
        Eref = 1.0 / (abs(Vc1) / (math.sqrt(2) * Vph))
        modo_e = "auto (V_PCC nominal na carga inicial)"
    else:
        Eref = parse_voltage(a.v_inversor, "--v-inversor", base=Vph) / Vph
        modo_e = f"definido pelo usuário ({a.v_inversor})"
    Vc1, If1, S1 = circuit(Eref, 0.0, r, [(RL1, LL1)])
    r.update(Eref=Eref, modo_e=modo_e, Pref=S1.real, Qref=S1.imag,
             Vpcc1=abs(Vc1) / math.sqrt(2) * math.sqrt(3), I1=abs(If1) / math.sqrt(2))
    E2, dw2, Vc2, If2, S2 = steady_state(r, [(RL1, LL1), (RL2, LL2)], Eref, S1.real, S1.imag, Dp, a.nq)
    r.update(E2=E2, f2=f0 * (1 + dw2), Vpcc2=abs(Vc2) / math.sqrt(2) * math.sqrt(3),
             I2=abs(If2) / math.sqrt(2), P2_out=S2.real)
    Sr = 0.9 - 1j * math.sqrt(1 - 0.81)
    E_pior = max(E2, Eref, abs(1 + (Rf + 1j * Xf) / Zb * Sr))
    Vdc_min = math.sqrt(2) * Vll * E_pior * (1 + a.margem_vdc)
    r.update(E_pior=E_pior, Vdc_min=Vdc_min)
    if Vdc < Vdc_min:
        avisos.append(f"V_dc = {Vdc:.0f} V < mínimo de {Vdc_min:.0f} V (SVPWM, margem "
                      f"{100 * a.margem_vdc:.0f}%): haverá sobremodulação.")
    if a.t_end <= a.t_step:
        erro("--t-end deve ser maior que --t-step")
    if a.evento == "falta_3f":
        if a.t_clear is None or not a.t_step < a.t_clear < a.t_end:
            erro("--t-clear deve satisfazer t_step < t_clear < t_end no evento falta_3f")
        if a.vg_falta < 0:
            erro("--vg-falta deve ser >= 0")

    # --- Modo rede (v6) --------------------------------------------------------------------
    r["modo"] = a.modo
    r["Dw"] = 0.0
    if a.modo == "rede":
        if a.scr <= 0 or a.xr_rede <= 0 or a.tw <= 0:
            erro("--scr, --xr-rede e --tw devem ser > 0")
        if a.scr_faixa:
            try:
                lo, hi = sorted(float(x) for x in a.scr_faixa.split(","))
            except ValueError:
                erro("--scr-faixa deve ter o formato 'min,max'")
        else:
            lo = hi = a.scr
        scrs = sorted({lo, hi, a.scr})
        r["scrs"] = scrs
        r["lcl"] = []
        for s in scrs:
            Xg = (Zb / s) / math.sqrt(1 + 1 / a.xr_rede ** 2)
            Lg = Xg / w0
            frl = math.sqrt((Lf + Lg) / (Lf * Lg * Cf)) / (2 * math.pi)
            r["lcl"].append((s, Lg, frl))
            if not 10 * f0 < frl < fsw / 2:
                avisos.append(f"Ressonância LCL com SCR {s:g}: f_r = {frl:.0f} Hz fora da faixa "
                              f"({10 * f0:.0f} Hz < f_r < {fsw / 2:.0f} Hz). Revise f_sw, L_f ou C_f; "
                              f"considere amortecimento (R_d ≈ {Rd * 1e3:.0f} mΩ).")
        sim = carregar_simulador()
        base = dict(sim.DEFAULTS)
        base.update(sn=Sn, vll=Vll, f0=f0, rf=Rf / Zb, xf=Xf / Zb, bc=Cf / Cb, p1=P1, q1=Q1,
                    p2=P2, q2=Q2, H=H, mp=a.mp, nq=a.nq, fc=a.fc, modo="rede",
                    xr_rede=a.xr_rede, vg=1.0, tw=a.tw, pref=None, qref=None)

        def zeta(dw, s):
            c = dict(base, dw=dw, scr=s)
            lam, _, _ = sim.analise_modal(c)
            dom = sim.modo_dominante(lam)
            if dom is None:
                return -1.0 if np.any(lam.real > 1e-9) else 1.0
            return dom[1]

        alvo = a.zeta_alvo + 0.005
        s_forte = max(scrs)
        if a.dw is not None:
            Dw = float(a.dw)
        elif zeta(0.0, s_forte) >= alvo:
            Dw = 0.0
        else:
            ant, Dw = 0.0, None
            for cand in np.arange(10.0, 2001.0, 10.0):
                if zeta(cand, s_forte) >= alvo:
                    Dw = cand
                    break
                ant = cand
            if Dw is None:
                erro(f"não foi possível atingir ζ = {a.zeta_alvo} com D_w <= 2000 no SCR {s_forte:g}")
            lo_, hi_ = ant, Dw
            for _ in range(25):
                mid = 0.5 * (lo_ + hi_)
                if zeta(mid, s_forte) >= alvo:
                    hi_ = mid
                else:
                    lo_ = mid
            Dw = math.ceil(hi_ * 10) / 10
        r["Dw"] = Dw
        r["zetas"] = [(s, zeta(0.0, s), zeta(Dw, s)) for s in scrs]
        for s, z0, z1 in r["zetas"]:
            if z1 < a.zeta_alvo:
                avisos.append(f"Com D_w = {Dw:g}, ζ = {z1:.3f} no SCR {s:g} (< alvo {a.zeta_alvo}).")
    # --- Eventos de frequência da rede (v7) ---------------------------------------------
    r["df_g"] = float(a.df_rede) if a.df_rede is not None else 0.0
    r["rocof_g"] = float(rocof_rede)
    if a.evento in EVENTOS_FREQ:
        if r["df_g"] == 0:
            erro(f"o evento '{a.evento}' exige --df-rede diferente de zero")
        if r["rocof_g"] <= 0:
            erro("--rocof-rede (ou --rocof-max) deve ser > 0")
        tw = a.tw
        Heff = H + r["Dw"] * tw / 2
        dP_droop = Dp * abs(r["df_g"]) / f0
        if a.evento == "freq_rampa":
            Tr = abs(r["df_g"]) / r["rocof_g"]
            dP_in = (r["rocof_g"] / f0) * (2 * H + r["Dw"] * tw * (1 - math.exp(-Tr / tw)))
        else:
            Tr, dP_in = 0.0, 0.0
        pico = P1 / Sn + dP_in + dP_droop
        r.update(Heff=Heff, Tr=Tr, dP_in=dP_in, dP_droop=dP_droop, P_pico=pico)
        if pico > a.pmax_pcs:
            avisos.append(f"P_pico previsto {pico:.3f} pu > p_max,pcs = {a.pmax_pcs} pu com H_eff = "
                          f"{Heff:.1f} s (H = {H:g} s + D_w·T_w/2 = {r['Dw'] * tw / 2:.1f} s): "
                          f"reduza RoCoF/Δf, H ou D_w·T_w, ou aumente a capacidade do PCS.")
        if r["Dw"] * tw / 2 > H:
            avisos.append(f"O washout domina a inércia efetiva (D_w·T_w/2 = {r['Dw'] * tw / 2:.1f} s "
                          f"> H = {H:g} s): numa rampa, a inércia efetiva é {Heff:.1f} s.")
    if a.evento == "sincronizacao":
        if abs(r["df_g"]) > 2.0:
            avisos.append(f"Desvio inicial de frequência de {r['df_g']:+.3f} Hz é elevado para uma manobra de paralelo.")
        if abs(a.d_fase) > 90:
            avisos.append(f"Ângulo inicial de {a.d_fase:+.1f}° pode exigir mais tempo de pré-sincronização.")
    if a.imax_pu > 0:
        # Pré-verificação conservadora: fonte de aproximadamente 1 pu atrás da impedância
        # física e virtual. A aceitação final permanece a simulação EMT.
        zg_pu = 1.0 / a.scr if a.modo == "rede" else 0.0
        zfis = math.hypot(r["rf_pu"] + zg_pu / max(a.xr_rede, 1e-12), r["xf_pu"] + zg_pu)
        if a.evento == "falta_3f":
            dv_prev = abs(1.0 - a.vg_falta)
            i_prev = dv_prev / max(zfis, 1e-12)
        elif a.evento == "fase":
            dv_prev = 2 * math.sin(abs(math.radians(a.d_fase)) / 2)
            i_prev = dv_prev / max(zfis, 1e-12)
        elif a.evento in EVENTOS_FREQ and "P_pico" in r:
            dv_prev, i_prev = 0.0, r["P_pico"]
        else:
            dv_prev, i_prev = 0.0, P1 / Sn
        rv_est = max(0.0, dv_prev / max(a.imax_pu, 1e-12) - zfis)
        r.update(i_prev_pu=i_prev, rv_est_pu=rv_est)
        if rv_est > rvmax:
            avisos.append(f"R_v,max = {rvmax:.3f} pu pode ser insuficiente; pré-estimativa conservadora "
                          f"requer {rv_est:.3f} pu.")
        if a.evento in EVENTOS_FREQ and i_prev > a.imax_pu:
            avisos.append(f"Resposta inercial prevista ({i_prev:.3f} pu) será recortada por I_max = "
                          f"{a.imax_pu:.3f} pu.")
    return r, avisos


# ---------------------------------------------------------------------------------------
# Saída
# ---------------------------------------------------------------------------------------
def resumo(r, a, avisos):
    L = ["=== Bases ===",
         f"S_n = {fmt_p(r['Sn'])}VA | V_LL = {r['Vll']:.4g} V | f0 = {r['f0']:g} Hz | "
         f"Z_b = {r['Zb']:.4g} Ω | Î_n = {r['Ipk']:.4g} A | modo: {r['modo']}",
         "=== Cargas ===",
         f"Inicial: {fmt_p(r['P1'])}W + j{fmt_p(r['Q1'])}var (FP {r['fp1']:.3f}) = {r['P1'] / r['Sn']:.4f} pu",
         f"Degrau : {fmt_p(r['P2'])}W + j{fmt_p(r['Q2'])}var (FP {r['fp2']:.3f}) = {r['P2'] / r['Sn']:.4f} pu",
         "=== Filtro (regra de projeto) ===",
         f"V_dc = {r['Vdc']:.0f} V{' (auto)' if r['vdc_auto'] else ''} | f_sw = {r['fsw']:g} Hz | "
         f"ripple = {100 * a.ripple:.0f}%",
         f"L_f = {r['Lf'] * 1e3:.4g} mH ({r['xf_pu']:.4f} pu) | C_f = {r['Cf'] * 1e6:.4g} µF "
         f"({r['bc_pu']:.4f} pu) | R_f = {r['Rf'] * 1e3:.4g} mΩ ({r['rf_pu']:.5f} pu)",
         f"f_r LC = {r['fr']:.0f} Hz | atenuação em f_sw = {20 * math.log10(r['att']):.1f} dB | "
         f"R_d sugerido = {r['Rd'] * 1e3:.4g} mΩ | V_dc mínimo = {r['Vdc_min']:.0f} V",
         "=== Inércia ===",
         f"ΔP_max = {r['dPmax_pu']:.4f} pu | RoCoF_max = {a.rocof_max} Hz/s -> H_min = {r['H_min']:.3f} s | "
         f"H adotado = {r['H']:.1f} s | H_max = {r['H_max']:.2f} s"]
    if r["modo"] == "ilhado":
        L += ["=== Ponto de operação inicial (ilhado) ===",
              f"E_ref = {r['Eref']:.5f} pu -> V_fase,inv = {r['Eref'] * r['Vph']:.2f} V ({r['modo_e']})",
              f"V_LL,PCC = {r['Vpcc1']:.2f} V | P_ref = {r['Pref']:.5f} pu | Q_ref = {r['Qref']:.5f} pu",
              "=== Regime previsto após o degrau ===",
              f"f = {r['f2']:.4f} Hz | V_LL,PCC = {r['Vpcc2']:.2f} V | I = {r['I2']:.4g} A"]
    else:
        L += ["=== Rede (v6) ===",
              f"X/R da rede = {a.xr_rede:g} | P_ref, Q_ref, E_ref: automáticos (troca nula com a rede)"]
        for s, Lg, frl in r["lcl"]:
            L.append(f"SCR {s:g}: L_g = {Lg * 1e3:.4g} mH | f_r LCL = {frl:.0f} Hz")
        L.append(f"=== Amortecimento (modelo linearizado) === D_w = {r['Dw']:g} pu, T_w = {a.tw:g} s, "
                 f"ζ_alvo = {a.zeta_alvo}")
        for s, z0, z1 in r["zetas"]:
            L.append(f"SCR {s:g}: ζ sem washout = {z0:+.3f}{' (INSTÁVEL)' if z0 < 0 else ''} | "
                     f"ζ com washout = {z1:.3f}")
        if r["imax_pu"] > 0:
            L += ["=== Limitador de corrente (v8) ===",
                  f"I_max = {r['imax_pu']:.3f} pu = {r['imax_pu']*r['Ipk']:.1f} A pico "
                  f"({r['imax_pu']*r['Ipk']/math.sqrt(2):.1f} A rms) | I_on = {r['i_on_pu']:.3f} pu",
                  f"R_v,max = {r['rv_max_pu']:.3f} pu = {r['rv_max_pu']*r['Zb']:.3f} Ω | "
                  f"X_v,max = {r['xv_rv']*r['rv_max_pu']:.3f} pu | k_aw = {r['k_aw']:g} 1/s",
                  f"Pré-verificação: I ≈ {r['i_prev_pu']:.3f} pu | R_v adicional ≈ {r['rv_est_pu']:.3f} pu"]
        if a.evento in EVENTOS_FREQ:
            perfil = (f"rampa de {r['df_g']:+g} Hz a {r['rocof_g']:g} Hz/s (T_r = {r['Tr']:g} s)"
                      if a.evento == "freq_rampa" else f"degrau de {r['df_g']:+g} Hz")
            L += ["=== Frequência da rede (v7) ===",
                  f"Evento: {perfil}",
                  f"H_eff = {r['Heff']:.2f} s",
                  f"ΔP inercial previsto = {r['dP_in']:.4f} pu | ΔP de droop = {r['dP_droop']:.4f} pu",
                  f"P_pico previsto = {r['P_pico']:.4f} pu (p_max,pcs = {a.pmax_pcs} pu)"]
        if a.evento == "sincronizacao":
            L += ["=== Sincronização e disjuntor (v9) ===",
                  f"Condição inicial: δ_g = {a.d_fase:+g}° | Δf_g = {r['df_g']:+g} Hz | "
                  f"pré-sincronizador: {'ativo' if a.pre_sync else 'desabilitado'}",
                  f"Janelas do relé 25: |ΔV| ≤ {a.sync_dv_max_pu:g} pu | "
                  f"|Δf| ≤ {a.sync_df_max_hz:g} Hz | |δ previsto| ≤ {a.sync_delta_max_deg:g}°",
                  f"Qualificação = {a.sync_hold_s:g} s | tempo mecânico = {a.breaker_delay_s:g} s"]
    if avisos:
        L.append("=== AVISOS ===")
        L += [f"(!) {x}" for x in avisos]
    return L


def yq(s):
    """String YAML entre aspas simples: não interpreta barras invertidas (caminhos Windows)."""
    return "'" + str(s).replace("'", "''") + "'"


def yaml_caso(r, a, linhas):
    v = lambda x: "null" if x is None else f"{x:.10g}"
    rede = r["modo"] == "rede"
    cab = ["# Caso gerado por gerar_caso_vsg_v10.py (v10)", "# Comando: " + " ".join(sys.argv), "#"] + \
          [f"# {x}" for x in linhas] + [""]
    corpo = f"""sistema:
  sn:       {v(r['Sn'])}
  vll:      {v(r['Vll'])}
  f0:       {v(r['f0'])}
  modo:     {yq(r['modo'])}

rede:
  scr:      {v(a.scr)}
  xr_rede:  {v(a.xr_rede)}
  vg:       1

filtro:
  rf:       {v(r['rf_pu'])}      # R_f = {r['Rf'] * 1e3:.4g} mΩ
  xf:       {v(r['xf_pu'])}      # L_f = {r['Lf'] * 1e3:.4g} mH
  bc:       {v(r['bc_pu'])}      # C_f = {r['Cf'] * 1e6:.4g} µF

carga:
  p1:       {v(r['P1'])}
  q1:       {v(r['Q1'])}
  p2:       {v(r['P2'])}
  q2:       {v(r['Q2'])}

controle:
  pref:     {v(None if rede else r['Pref'])}
  qref:     {v(None if rede else r['Qref'])}
  eref:     {v(r['Eref'])}
  H:        {v(r['H'])}
  mp:       {v(a.mp)}
  nq:       {v(a.nq)}
  fc:       {v(a.fc)}
  dw:       {v(r['Dw'])}
  tw:       {v(a.tw)}
  imax_pu:  {v(r['imax_pu'])}
  i_on_pu:  {v(r['i_on_pu'])}
  rv_max_pu: {v(r['rv_max_pu'])}
  xv_rv:    {v(r['xv_rv'])}
  k_aw:     {v(r['k_aw'])}

sincronismo:
  delta_g0_graus: {v(a.delta_g0_graus)}
  df_g0_hz: {v(a.df_g0_hz)}
  disjuntor_inicial: {yq(a.disjuntor_inicial)}
  estrategia_sync: {yq(a.estrategia_sync)}
  dv_sync_max_pu: {v(a.dv_sync_max_pu)}
  df_sync_max_hz: {v(a.df_sync_max_hz)}
  dtheta_sync_max_graus: {v(a.dtheta_sync_max_graus)}
  t_sync_hold_s: {v(a.t_sync_hold_s)}
  t_fechamento_s: {v(a.t_fechamento_s)}
  t_sync_timeout_s: {v(a.t_sync_timeout_s)}
  dt_rele_s: {v(a.dt_rele_s)}
  antecipar_fechamento: {yq(a.antecipar_fechamento)}
  vmin_medicao_pu: {v(a.vmin_medicao_pu)}
  kp_theta_hz_rad: {v(a.kp_theta_hz_rad)}
  ki_theta_hz_rad_s: {v(a.ki_theta_hz_rad_s)}
  df_sync_lim_hz: {v(a.df_sync_lim_hz)}
  kp_v: {v(a.kp_v)}
  ki_v_s: {v(a.ki_v_s)}
  de_sync_lim_pu: {v(a.de_sync_lim_pu)}
  t_release_sync_s: {v(a.t_release_sync_s)}
  pll_sync_bw_hz: {v(a.pll_sync_bw_hz)}
  pre_sync: {v(a.pre_sync)}
  kp_sync_p: {v(a.kp_sync_p)}
  ki_sync_p: {v(a.ki_sync_p)}
  p_sync_max_pu: {v(a.p_sync_max_pu)}
  kp_sync_v: {v(a.kp_sync_v)}
  ki_sync_v: {v(a.ki_sync_v)}
  e_sync_max_pu: {v(a.e_sync_max_pu)}
  sync_dv_max_pu: {v(a.sync_dv_max_pu)}
  sync_df_max_hz: {v(a.sync_df_max_hz)}
  sync_delta_max_deg: {v(a.sync_delta_max_deg)}
  sync_hold_s: {v(a.sync_hold_s)}
  breaker_delay_s: {v(a.breaker_delay_s)}

lado_cc:
  modelo_cc: {yq(a.modelo_cc)}
  c_dc_f: {v(a.c_dc_f)}
  vdc_inicial_v: {yq('auto') if str(a.vdc_inicial_v).lower() == 'auto' else v(float(a.vdc_inicial_v))}
  permitir_desequilibrio_inicial: {'true' if a.permitir_desequilibrio_inicial else 'false'}
  vdc_min_oper_v: {v(a.vdc_min_oper_v)}
  vdc_max_oper_v: {v(a.vdc_max_oper_v)}
  m_max_pu: {v(a.m_max_pu)}
  vdc_piso_numerico_v: {yq('auto') if str(a.vdc_piso_numerico_v).lower() == 'auto' else v(float(a.vdc_piso_numerico_v))}

bateria:
  capacidade_ah: {v(a.capacidade_ah)}
  soc_inicial: {v(a.soc_inicial)}
  r0_ohm: {v(a.r0_ohm)}
  ocv_soc_pu: {('null' if a.ocv_soc_pu is None else '[' + ', '.join(a.ocv_soc_pu.split(',')) + ']')}
  ocv_v: {('null' if a.ocv_v is None else '[' + ', '.join(a.ocv_v.split(',')) + ']')}
  soc_min_oper_pu: {v(a.soc_min_oper_pu)}
  soc_max_oper_pu: {v(a.soc_max_oper_pu)}
  i_desc_max_a: {v(a.i_desc_max_a)}
  i_carga_max_a: {v(a.i_carga_max_a)}

evento:
  evento:   {yq(a.evento)}
  d_fase:   {v(a.d_fase)}
  df_g:     {v(r['df_g'])}
  rocof_g:  {v(r['rocof_g'])}
  vg_falta: {v(a.vg_falta)}
  t_clear:  {v(a.t_clear)}

simulacao:
  t_step:   {v(a.t_step)}
  t_end:    {v(a.t_end)}
  dt_out:   {v(a.dt_out)}

saida:
  prefixo:  {yq(a.prefixo)}
"""
    return "\n".join(cab) + corpo


def gerar_exemplos(pasta):
    """Gera a matriz reproduzível de dez exemplos definida pelo contrato v10."""
    sim = carregar_simulador()
    os.makedirs(pasta, exist_ok=True)
    cc = {
        "modelo_cc": "thevenin", "c_dc_f": 0.10, "vdc_inicial_v": "auto",
        "permitir_desequilibrio_inicial": False, "vdc_min_oper_v": 500.0,
        "vdc_max_oper_v": 800.0, "m_max_pu": 1.0,
        "vdc_piso_numerico_v": "auto", "capacidade_ah": 100.0,
        "soc_inicial": 0.60, "r0_ohm": 0.05,
        "ocv_soc_pu": [0.0, 0.25, 0.50, 0.75, 1.0],
        "ocv_v": [660.0, 680.0, 700.0, 720.0, 740.0],
        "soc_min_oper_pu": 0.10, "soc_max_oper_pu": 0.90,
        "i_desc_max_a": 0.0, "i_carga_max_a": 0.0,
    }
    curto = {"evento": "nenhum", "t_step": 0.2, "t_end": 0.4, "dt_out": 1e-4}
    casos = {
        "caso_cc_ideal_regressao_v10.yaml": {**curto, "modelo_cc": "ideal"},
        "caso_bess_degrau_carga_v10.yaml": {**cc, "modo": "ilhado", "evento": "carga",
                                              "t_step": 0.5, "t_end": 1.2, "dt_out": 1e-4},
        "caso_bess_cdc_baixo_v10.yaml": {**cc, **curto, "c_dc_f": 0.02},
        "caso_bess_cdc_alto_v10.yaml": {**cc, **curto, "c_dc_f": 0.20},
        "caso_bess_regeneracao_v10.yaml": {**cc, "modo": "rede", "p1": 1000.0,
                                             "q1": 0.0, "pref": -0.15, "evento": "nenhum",
                                             "t_step": 0.2, "t_end": 1.0, "dt_out": 2e-4},
        "caso_bess_saturacao_modulacao_v10.yaml": {**cc, "modo": "ilhado",
                                                     "evento": "carga", "p2": 100000.0,
                                                     "m_max_pu": 0.75, "t_step": 0.3,
                                                     "t_end": 1.0, "dt_out": 1e-4},
        "caso_bess_falta_3f_v10.yaml": {**cc, "modo": "rede", "scr": 20.0,
                                         "dw": 100.0, "evento": "falta_3f",
                                         "vg_falta": 0.0, "t_step": 0.3,
                                         "t_clear": 0.4, "t_end": 1.0, "dt_out": 1e-4},
        "caso_bess_freq_rampa_v10.yaml": {**cc, "modo": "rede", "dw": 100.0,
                                           "evento": "freq_rampa", "df_g": -0.2,
                                           "rocof_g": 0.5, "t_step": 0.3,
                                           "t_end": 1.2, "dt_out": 2e-4},
        "caso_bess_fechamento_ativo_v10.yaml": {**cc, "modo": "rede",
                                                  "evento": "fechamento",
                                                  "disjuntor_inicial": "aberto",
                                                  "estrategia_sync": "ativo",
                                                  "delta_g0_graus": 20.0,
                                                  "df_g0_hz": 0.12, "t_step": 0.5,
                                                  "t_end": 8.0, "dt_out": 2e-4},
        "caso_bess_inviavel_v10.yaml": {**cc, **curto, "vdc_inicial_v": 650.0},
    }
    for nome, alteracoes in casos.items():
        cfg = dict(sim.DEFAULTS)
        cfg.update(alteracoes)
        caminho = os.path.join(pasta, nome)
        with open(caminho, "w", encoding="utf-8") as fh:
            fh.write(sim.config_to_yaml(cfg, f"Exemplo v10 — {nome}"))
        print(f"Exemplo gravado em: {caminho}")


def main():
    console_seguro()
    ap = argparse.ArgumentParser(description="Gera YAMLs de caso para o simulador VSG v10.")
    ap.add_argument("--gerar-exemplos", metavar="PASTA",
                    help="Gera os dez casos de aceitação v10 e encerra")
    g = ap.add_argument_group("sistema")
    g.add_argument("--sn", default="100 kVA")
    g.add_argument("--vll", default="380 V")
    g.add_argument("--f0", type=float, default=60.0)
    g.add_argument("--modo", default="ilhado", help="ilhado | rede (padrão: ilhado)")
    g = ap.add_argument_group("rede (v6)")
    g.add_argument("--scr", type=float, default=5.0, help="SCR do caso simulado (padrão: 5)")
    g.add_argument("--scr-faixa", default=None, help="Faixa 'min,max' para verificação e D_w")
    g.add_argument("--xr-rede", type=float, default=10.0)
    g.add_argument("--tw", type=float, default=1.0, help="Constante de tempo do washout [s]")
    g.add_argument("--zeta-alvo", type=float, default=0.3, help="ζ mínimo do modo eletromecânico")
    g.add_argument("--dw", type=float, default=None, help="Fixa D_w manualmente (desliga o cálculo)")
    g.add_argument("--evento", default="carga",
                    help="nenhum | carga | pref | fase | freq_degrau | freq_rampa | falta_3f | sincronizacao")
    g.add_argument("--df-rede", type=float, default=None,
                   help="[v7] Desvio final da frequência da rede [Hz], com sinal (eventos freq_*)")
    g.add_argument("--d-fase", type=float, default=5.0, help="Salto angular da rede [graus]")
    g.add_argument("--vg-falta", type=float, default=0.0, help="Tensão relativa da rede durante falta_3f")
    g.add_argument("--t-clear", type=float, default=None, help="Instante de eliminação da falta_3f [s]")
    g = ap.add_argument_group("cargas")
    g.add_argument("--carga", default="50 kW")
    g.add_argument("--fp", type=float, default=0.958)
    g.add_argument("--degrau", default="30 kW")
    g.add_argument("--fp-degrau", type=float, default=None)
    g = ap.add_argument_group("filtro")
    g.add_argument("--vdc", default="auto")
    g.add_argument("--fsw", default="8 kHz")
    g.add_argument("--ripple", type=float, default=0.15)
    g.add_argument("--qc-max", type=float, default=0.05)
    g.add_argument("--xr", type=float, default=40.0)
    g.add_argument("--margem-vdc", type=float, default=0.10)
    g = ap.add_argument_group("inércia")
    g.add_argument("--rocof-max", type=float, default=2.0)
    g.add_argument("--dpmax", default=None)
    g.add_argument("--margem-h", type=float, default=1.0)
    g.add_argument("--h-faixa", default=None)
    g.add_argument("--pmax-pcs", type=float, default=1.1)
    g.add_argument("--rocof-rede", type=float, default=None,
                   help="RoCoF da rede [Hz/s]: H_max e, na v7, taxa da rampa (padrão: --rocof-max)")
    g = ap.add_argument_group("controle")
    g.add_argument("--v-inversor", default="auto")
    g.add_argument("--mp", type=float, default=0.05)
    g.add_argument("--nq", type=float, default=0.05)
    g.add_argument("--fc", type=float, default=10.0)
    g.add_argument("--imax-pu", type=float, default=0.0)
    g.add_argument("--i-on-pu", type=float, default=None)
    g.add_argument("--rv-max-pu", type=float, default=None,
                   help="Teto de R_v [pu]; auto usa 3 pu com limitação ativa")
    g.add_argument("--xv-rv", type=float, default=0.0)
    g.add_argument("--k-aw", type=float, default=20.0)
    g = ap.add_argument_group("lado CC e bateria (v10)")
    g.add_argument("--modelo-cc", choices=("ideal", "thevenin"), default="ideal")
    g.add_argument("--c-dc-f", type=float, default=None)
    g.add_argument("--vdc-inicial-v", default="auto")
    g.add_argument("--permitir-desequilibrio-inicial", action="store_true")
    g.add_argument("--vdc-min-oper-v", type=float, default=None)
    g.add_argument("--vdc-max-oper-v", type=float, default=None)
    g.add_argument("--m-max-pu", type=float, default=1.0)
    g.add_argument("--vdc-piso-numerico-v", default="auto")
    g.add_argument("--capacidade-ah", type=float, default=None)
    g.add_argument("--soc-inicial", type=float, default=None)
    g.add_argument("--r0-ohm", type=float, default=None)
    g.add_argument("--ocv-soc-pu", default=None)
    g.add_argument("--ocv-v", default=None)
    g.add_argument("--soc-min-oper-pu", type=float, default=0.10)
    g.add_argument("--soc-max-oper-pu", type=float, default=0.90)
    g.add_argument("--i-desc-max-a", type=float, default=0.0)
    g.add_argument("--i-carga-max-a", type=float, default=0.0)
    g = ap.add_argument_group("sincronismo e disjuntor (v9)")
    g.add_argument("--delta-g0-graus", type=float, default=0.0)
    g.add_argument("--df-g0-hz", type=float, default=0.0)
    g.add_argument("--disjuntor-inicial", choices=("aberto","fechado"), default=None)
    g.add_argument("--estrategia-sync", choices=("forcado","passivo","ativo"), default="passivo")
    g.add_argument("--dv-sync-max-pu", type=float, default=0.05)
    g.add_argument("--df-sync-max-hz", type=float, default=0.10)
    g.add_argument("--dtheta-sync-max-graus", type=float, default=5.0)
    g.add_argument("--t-sync-hold-s", type=float, default=0.10)
    g.add_argument("--t-fechamento-s", type=float, default=0.06)
    g.add_argument("--t-sync-timeout-s", type=float, default=10.0)
    g.add_argument("--dt-rele-s", type=float, default=0.001)
    g.add_argument("--antecipar-fechamento", choices=("true","false"), default="true")
    g.add_argument("--vmin-medicao-pu", type=float, default=0.20)
    g.add_argument("--kp-theta-hz-rad", type=float, default=0.60)
    g.add_argument("--ki-theta-hz-rad-s", type=float, default=0.20)
    g.add_argument("--df-sync-lim-hz", type=float, default=0.50)
    g.add_argument("--kp-v", type=float, default=0.80)
    g.add_argument("--ki-v-s", type=float, default=0.30)
    g.add_argument("--de-sync-lim-pu", type=float, default=0.10)
    g.add_argument("--t-release-sync-s", type=float, default=0.20)
    g.add_argument("--pll-sync-bw-hz", type=float, default=5.0)
    g.add_argument("--pre-sync", type=float, default=1.0, help="1 habilita; 0 desabilita")
    g.add_argument("--kp-sync-p", type=float, default=4.0)
    g.add_argument("--ki-sync-p", type=float, default=2.0)
    g.add_argument("--p-sync-max-pu", type=float, default=0.50)
    g.add_argument("--kp-sync-v", type=float, default=0.50)
    g.add_argument("--ki-sync-v", type=float, default=1.0)
    g.add_argument("--e-sync-max-pu", type=float, default=0.20)
    g.add_argument("--sync-dv-max-pu", type=float, default=0.05)
    g.add_argument("--sync-df-max-hz", type=float, default=0.10)
    g.add_argument("--sync-delta-max-deg", type=float, default=10.0)
    g.add_argument("--sync-hold-s", type=float, default=0.10)
    g.add_argument("--breaker-delay-s", type=float, default=0.05)
    g = ap.add_argument_group("simulação e saída")
    g.add_argument("--t-step", type=float, default=1.0)
    g.add_argument("--t-end", type=float, default=4.0)
    g.add_argument("--dt-out", type=float, default=5e-5)
    g.add_argument("--prefixo", default="vsg")
    g.add_argument("-o", "--saida", default="meu_caso.yaml")
    try:
        a = ap.parse_args()
    except SystemExit as e:
        if e.code not in (0, None):
            erro("argumentos inválidos (veja --help)")
        raise
    if a.gerar_exemplos:
        gerar_exemplos(a.gerar_exemplos)
        return
    if a.evento not in EVENTOS:
        erro(f"--evento inválido: {a.evento!r}")
    if a.evento in ("fase", "falta_3f", "sincronizacao", "fechamento", *EVENTOS_FREQ) and a.modo != "rede":
        erro(f"o evento '{a.evento}' exige --modo rede")
    if a.disjuntor_inicial is None:
        a.disjuntor_inicial = "aberto" if a.evento == "fechamento" else "fechado"

    r, avisos = projetar(a)
    linhas = resumo(r, a, avisos)
    print("\n".join(linhas))
    pasta = os.path.dirname(a.saida)
    if pasta:
        os.makedirs(pasta, exist_ok=True)
    with open(a.saida, "w", encoding="utf-8") as fh:
        fh.write(yaml_caso(r, a, linhas))
    print(f"\nArquivo de caso gravado em: {a.saida}")
    print(f"Para simular: python {SIMULADOR} --config {a.saida}")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001
        sys.exit(f"Erro inesperado: {type(e).__name__}: {e}")

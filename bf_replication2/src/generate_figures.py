#!/usr/bin/env python3
"""
src/generate_figures.py — Phase 6: Paper-quality figures from calibration output

Reads CSVs from data/results/ and generates PNG figures matching the paper:
  Figure 2: Benchmark (complementarity) — RGDP decomposition bars
  Figure 3: CD regime — RGDP decomposition bars
  Figure 4: HtM sweep — RGDP as function of HtM share
  Appendix A1: Histogram of sectoral price changes
  Appendix A2: Scatter plot: model prices vs PPI data
  Appendix A3: Scatter plot: model prices vs PCE shock
  Appendix A4: Sector unemployment vs BLS shocks

Usage:
  /opt/hermes/.venv/bin/python src/generate_figures.py

Output: figures/*.png

Paper targets (Baqaee & Farhi 2020):
  Benchmark: RGDP −8.1%, supply-only −5.7%, demand-only −5.1%
  CD: RGDP −8.2% (approx), supply-only −5.9%, demand-only −6.0%
"""

import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

# ─── paths ────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / 'data' / 'results'
FIGS = REPO / 'figures'
FIGS.mkdir(exist_ok=True)

# Use Hermes venv python explicitly
print(f"→ Data dir: {DATA}")
print(f"→ Output:   {FIGS}")

# ─── load data ────────────────────────────────────────────────
summary_loop1 = pd.read_csv(DATA / 'summary_loop1.csv')
summary_htm   = pd.read_csv(DATA / 'summary_loop2_htm.csv')
baseline      = pd.read_csv(DATA / 'baseline_fit.csv')

SHOCK_LABELS = ['Baseline', 'Supply\nonly', 'Demand\nonly',
                'Agg.\ndemand', 'Supply+\nsectoral']

print("✓ Data loaded")
print(f"  Loop1: {summary_loop1.shape[0]} shock types")
print(f"  HtM:   {summary_htm.shape[0]} HtM shares")
print(f"  Baseline: {baseline.shape[0]} sectors")

# ─── colour palette (accessible, print-friendly) ──────────────
BLUE  = '#2166AC'
RED   = '#B2182B'
GREEN = '#4DAF4A'
GREY  = '#757575'
ORANGE = '#E69F00'
PALETTE = ['#2166AC', '#D6604D', '#4DAF4A', '#F4A582', '#92C5DE']

# ─── helper: styled bar axis ──────────────────────────────────
def style_bar_ax(ax, title, ylabel='Real GDP change (%)', ylim=None):
    ax.axhline(0, color='black', linewidth=0.6)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xlabel('')
    ax.tick_params(axis='both', labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ylim:
        ax.set_ylim(ylim)

# ═══════════════════════════════════════════════════════════════
# Figure 2: Benchmark (complementarity) — 5 shock types
#           Bars: RGDP Δ, plus markers for unemployment & inflation
# ═══════════════════════════════════════════════════════════════
def fig2_benchmark():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    
    rgdp = summary_loop1['RGDP_benchmark'].values  # negative = decline
    unemp = summary_loop1['Unemp_benchmark'].values
    
    x = np.arange(5)
    width = 0.55
    
    # RGDP bar chart (negative = decline, faithful to B&F)
    bars = ax1.bar(x, rgdp, width, color=PALETTE, edgecolor='white', linewidth=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(SHOCK_LABELS, fontsize=9, rotation=20, ha='right')
    style_bar_ax(ax1, 'Figure 2: Real GDP change (benchmark)', ylim=(-11, 1))
    
    # Value labels on bars (sign-aware: below tip for declines)
    for bar, val in zip(bars, rgdp):
        h = bar.get_height()
        off = 0.15 if h >= 0 else -0.5
        va = 'bottom' if h >= 0 else 'top'
        ax1.text(bar.get_x() + bar.get_width()/2, h + off,
                f'{val:.1f}%', ha='center', va=va, fontsize=8.5, fontweight='bold')
    
    # Unemployment chart
    ax2.bar(x, unemp, width, color=PALETTE, edgecolor='white', linewidth=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(SHOCK_LABELS, fontsize=9, rotation=20, ha='right')
    style_bar_ax(ax2, 'Unemployment change (benchmark)',
                 ylabel='Unemployment (pp)', ylim=(-2, 12))
    for i, v in enumerate(unemp):
        ax2.text(i, v + 0.2, f'{v:.1f}', ha='center', va='bottom', fontsize=8.5, fontweight='bold')
    
    plt.tight_layout()
    path = FIGS / 'fig2_benchmark.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Figure 2 saved: {path}")

# ═══════════════════════════════════════════════════════════════
# Figure 3: CD regime
# ═══════════════════════════════════════════════════════════════
def fig3_cd():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    
    rgdp = summary_loop1['RGDP_cd'].values  # negative = decline
    unemp = summary_loop1['Unemp_cd'].values
    
    x = np.arange(5)
    width = 0.55
    
    bars = ax1.bar(x, rgdp, width, color=PALETTE, edgecolor='white', linewidth=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(SHOCK_LABELS, fontsize=9, rotation=20, ha='right')
    style_bar_ax(ax1, 'Figure 3: Real GDP change (Cobb–Douglas)', ylim=(-11, 1))
    for bar, val in zip(bars, rgdp):
        h = bar.get_height()
        off = 0.15 if h >= 0 else -0.5
        va = 'bottom' if h >= 0 else 'top'
        ax1.text(bar.get_x() + bar.get_width()/2, h + off,
                f'{val:.1f}%', ha='center', va=va, fontsize=8.5, fontweight='bold')
    
    # Handle negative zero
    unemp_display = np.where(np.abs(unemp) < 0.01, 0.0, unemp)
    ax2.bar(x, unemp_display, width, color=PALETTE, edgecolor='white', linewidth=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(SHOCK_LABELS, fontsize=9, rotation=20, ha='right')
    style_bar_ax(ax2, 'Unemployment change (Cobb–Douglas)',
                 ylabel='Unemployment (pp)', ylim=(-2, 14))
    for i, v in enumerate(unemp_display):
        ax2.text(i, v + 0.2, f'{v:.1f}', ha='center', va='bottom', fontsize=8.5, fontweight='bold')
    
    plt.tight_layout()
    path = FIGS / 'fig3_cd.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Figure 3 saved: {path}")

# ═══════════════════════════════════════════════════════════════
# Figure 4: HtM sweep — line chart
# ═══════════════════════════════════════════════════════════════
def fig4_htm():
    # Check if we have real data (not all zeros)
    htm_rgdp = summary_htm['RGDP_benchmark'].values
    htm_cd   = summary_htm['RGDP_cd'].values
    
    # Even if summary is zeros, try reading from cell files
    # to reconstruct from individual timeseries
    if np.all(np.abs(htm_rgdp) < 1e-10) and np.all(np.abs(htm_cd) < 1e-10):
        print("⚠  summary_loop2_htm.csv contains zeros — reconstructing from cell files...")
        htm_data = reconstruct_htm()
        if htm_data is not None:
            htm_rgdp, htm_cd, htm_shares = htm_data
        else:
            print("⚠  Cannot reconstruct — showing placeholder")
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.text(0.5, 0.5, 'HtM sweep data not yet available\n(re-run loop=2 on Mac)',
                    transform=ax.transAxes, ha='center', va='center', fontsize=13)
            ax.set_title('Figure 4: HtM sweep — data pending')
            path = FIGS / 'fig4_htm_sweep.png'
            fig.savefig(path, dpi=200, bbox_inches='tight')
            plt.close(fig)
            print(f"⚠ Placeholder saved: {path}")
            return
    
    htm_shares = summary_htm['htm_share'].values
    
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(htm_shares, htm_rgdp * -100, 'o-', color=BLUE, linewidth=2.2,
            markersize=7, label='Benchmark (complementarity)')
    ax.plot(htm_shares, htm_cd * -100, 's--', color=RED, linewidth=2.2,
            markersize=7, label='Cobb–Douglas')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xlabel('Hand-to-mouth share', fontsize=11)
    ax.set_ylabel('Real GDP change (%)', fontsize=11)
    ax.set_title('Figure 4: HtM share and real GDP response', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    
    # Annotate values
    for s, v_b, v_c in zip(htm_shares, htm_rgdp * -100, htm_cd * -100):
        ax.annotate(f'{v_b:.1f}', (s, v_b), textcoords='offset points',
                   xytext=(0, -14), ha='center', fontsize=7.5, color=BLUE)
    
    plt.tight_layout()
    path = FIGS / 'fig4_htm_sweep.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Figure 4 saved: {path}")

def reconstruct_htm():
    """Reconstruct HtM data from individual cell timeseries when summary is missing."""
    htm_shares = np.arange(0, 1.01, 0.2)
    htm_rgdp_b = np.full(6, np.nan)
    htm_rgdp_c = np.full(6, np.nan)
    
    for s_idx in range(1, 7):
        # Try benchmark regime (el=1)
        cell_dir = DATA / 'loop2' / 'el1' / 'st1' / f's{s_idx}'
        ts_path = cell_dir / 'timeseries.csv'
        if ts_path.exists():
            ts = pd.read_csv(ts_path)
            if len(ts) > 0:
                # Chain-weighted real GDP
                gdp = ts['GDP'].values
                nom = ts['nominal_GDP'].values
                htm_rgdp_b[s_idx - 1] = (1 - gdp[-1]) * 100
        
        # Try CD regime (el=2)
        cell_dir = DATA / 'loop2' / 'el2' / 'st1' / f's{s_idx}'
        ts_path = cell_dir / 'timeseries.csv'
        if ts_path.exists():
            ts = pd.read_csv(ts_path)
            if len(ts) > 0:
                gdp = ts['GDP'].values
                htm_rgdp_c[s_idx - 1] = (1 - gdp[-1]) * 100
    
    if np.all(np.isnan(htm_rgdp_b)) and np.all(np.isnan(htm_rgdp_c)):
        return None
    
    return htm_rgdp_b / 100 * -1, htm_rgdp_c / 100 * -1, htm_shares

# ═══════════════════════════════════════════════════════════════
# Appendix A1: Histogram of sectoral price changes
# ═══════════════════════════════════════════════════════════════
def fig_a1_price_histogram():
    fig, ax = plt.subplots(figsize=(7, 4.5))
    
    prices = baseline['price_model'].values * 100  # to percent
    
    # Mark sticky vs flexible
    sticky = baseline['sticky'].values
    
    ax.hist(prices, bins=20, color=BLUE, edgecolor='white', alpha=0.8,
            label=f'All sectors (n={len(prices)})')
    
    # Add sticky-sector markers
    sticky_prices = prices[sticky]
    if len(sticky_prices) > 0:
        ax.scatter(sticky_prices, np.full_like(sticky_prices, 1),
                  color=RED, s=30, zorder=5, marker='v',
                  label=f'Sticky-price sectors (n={len(sticky_prices)})')
    
    ax.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    ax.set_xlabel('Model-implied sectoral price change (%)', fontsize=11)
    ax.set_ylabel('Number of sectors', fontsize=11)
    ax.set_title('Appendix A1: Distribution of sectoral price changes', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    path = FIGS / 'fig_a1_price_histogram.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Appendix A1 saved: {path}")

# ═══════════════════════════════════════════════════════════════
# Appendix A2: Model prices vs PPI data
# ═══════════════════════════════════════════════════════════════
def fig_a2_model_vs_ppi():
    fig, ax = plt.subplots(figsize=(7, 6))
    
    model_p = baseline['price_model'].values * 100
    ppi     = baseline['PPI_data'].values * 100
    names   = baseline['name'].values
    sticky  = baseline['sticky'].values
    
    # Colour by sticky/flexible
    colors = np.where(sticky, RED, BLUE)
    
    ax.scatter(ppi, model_p, c=colors, s=25, alpha=0.7, edgecolors='none')
    
    # 45-degree line
    lims = [min(ppi.min(), model_p.min()) - 5, max(ppi.max(), model_p.max()) + 5]
    ax.plot(lims, lims, 'k--', linewidth=0.8, alpha=0.5, label='45° line')
    
    # OLS fit
    from numpy.polynomial.polynomial import polyfit
    coeffs = np.polyfit(ppi, model_p, 1)
    x_fit = np.linspace(lims[0], lims[1], 100)
    ax.plot(x_fit, np.polyval(coeffs, x_fit), 'grey', linewidth=1.0, alpha=0.6,
            label=f'OLS (slope={coeffs[0]:.2f})')
    
    # Label outliers
    resid = np.abs(model_p - np.polyval(coeffs, ppi))
    outliers = np.argsort(resid)[-5:]
    for idx in outliers:
        ax.annotate(names[idx][:20], (ppi[idx], model_p[idx]),
                   textcoords='offset points', xytext=(5, 5), fontsize=6.5, alpha=0.7)
    
    ax.set_xlabel('PPI inflation data (%)', fontsize=11)
    ax.set_ylabel('Model-implied price change (%)', fontsize=11)
    ax.set_title('Appendix A2: Model vs data sectoral price changes', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(0, color='grey', linewidth=0.4)
    ax.axvline(0, color='grey', linewidth=0.4)
    ax.tick_params(labelsize=10)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    plt.tight_layout()
    path = FIGS / 'fig_a2_model_vs_ppi.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Appendix A2 saved: {path}")

# ═══════════════════════════════════════════════════════════════
# Appendix A3: Sector unemployment vs BLS shocks
# ═══════════════════════════════════════════════════════════════
def fig_a3_unemp_vs_bls():
    fig, ax = plt.subplots(figsize=(7, 6))
    
    model_unemp = baseline['unemp_model'].values * 100
    bls_shock   = baseline['BLS_shock'].values * 100
    names       = baseline['name'].values
    sticky      = baseline['sticky'].values
    
    colors = np.where(sticky, RED, BLUE)
    
    ax.scatter(bls_shock, model_unemp, c=colors, s=25, alpha=0.7, edgecolors='none')
    
    lims = [min(bls_shock.min(), model_unemp.min()) - 10,
            max(bls_shock.max(), model_unemp.max()) + 10]
    ax.plot(lims, lims, 'k--', linewidth=0.8, alpha=0.5, label='45° line')
    
    # OLS
    mask = np.isfinite(bls_shock) & np.isfinite(model_unemp)
    if mask.sum() > 2:
        coeffs = np.polyfit(bls_shock[mask], model_unemp[mask], 1)
        x_fit = np.linspace(lims[0], lims[1], 100)
        ax.plot(x_fit, np.polyval(coeffs, x_fit), 'grey', linewidth=1.0, alpha=0.6,
                label=f'OLS (slope={coeffs[0]:.2f})')
    
    # Label outliers
    resid = np.abs(model_unemp - np.polyval(coeffs, bls_shock))
    outliers = np.argsort(resid)[-5:]
    for idx in outliers:
        ax.annotate(names[idx][:20], (bls_shock[idx], model_unemp[idx]),
                   textcoords='offset points', xytext=(5, 5), fontsize=6.5, alpha=0.7)
    
    ax.set_xlabel('BLS shock (PP change in log employment)', fontsize=11)
    ax.set_ylabel('Model-implied sectoral employment change (%)', fontsize=11)
    ax.set_title('Appendix A3: Model employment vs BLS shocks', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(0, color='grey', linewidth=0.4)
    ax.axvline(0, color='grey', linewidth=0.4)
    ax.tick_params(labelsize=10)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    plt.tight_layout()
    path = FIGS / 'fig_a3_unemp_vs_bls.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Appendix A3 saved: {path}")

# ═══════════════════════════════════════════════════════════════
# Summary bar chart: Benchmark vs CD side by side (Figure 2+3 combined)
# ═══════════════════════════════════════════════════════════════
def fig_combined():
    """Side-by-side comparison of benchmark vs CD across shock types."""
    fig, ax = plt.subplots(figsize=(9, 5))
    
    rgdp_b = summary_loop1['RGDP_benchmark'].values  # negative = decline
    rgdp_c = summary_loop1['RGDP_cd'].values
    
    x = np.arange(5)
    width = 0.35
    
    ax.bar(x - width/2, rgdp_b, width, color=BLUE, edgecolor='white',
           linewidth=0.5, label='Benchmark (complementarity)', alpha=0.85)
    ax.bar(x + width/2, rgdp_c, width, color=RED, edgecolor='white',
           linewidth=0.5, label='Cobb–Douglas', alpha=0.85)
    
    ax.set_xticks(x)
    ax.set_xticklabels(['Baseline', 'Supply only', 'Demand only',
                        'Agg. demand', 'Supply +\nsectoral dem.'], fontsize=9)
    ax.set_ylabel('Real GDP change (%)', fontsize=11)
    ax.set_title('Figure 2+3: Real GDP decomposition — Benchmark vs CD', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.9)
    ax.axhline(0, color='black', linewidth=0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.set_ylim(-11, 1)
    
    # Value labels (sign-aware: below tip for declines)
    for i in range(5):
        hb, hc = rgdp_b[i], rgdp_c[i]
        ob = 0.15 if hb >= 0 else -0.5
        oc = 0.15 if hc >= 0 else -0.5
        vb = 'bottom' if hb >= 0 else 'top'
        vc = 'bottom' if hc >= 0 else 'top'
        ax.text(i - width/2, hb + ob, f'{hb:.1f}%',
               ha='center', va=vb, fontsize=7.5, fontweight='bold', color=BLUE)
        ax.text(i + width/2, hc + oc, f'{hc:.1f}%',
               ha='center', va=vc, fontsize=7.5, fontweight='bold', color=RED)
    
    plt.tight_layout()
    path = FIGS / 'fig_combined.png'
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Combined figure saved: {path}")

# ═══════════════════════════════════════════════════════════════
# Run all
# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("\n" + "=" * 50)
    print("GENERATING FIGURES FOR BAQAEE–FARHI (2020) REPLICATION")
    print("=" * 50 + "\n")
    
    fig2_benchmark()
    fig3_cd()
    fig4_htm()
    fig_a1_price_histogram()
    fig_a2_model_vs_ppi()
    fig_a3_unemp_vs_bls()
    fig_combined()
    
    print(f"\n{'=' * 50}")
    print(f"All figures saved to {FIGS.resolve()}/")
    print(f"Files: {sorted([f.name for f in FIGS.glob('*.png')])}")
    print(f"{'=' * 50}")
#!/usr/bin/env python3
"""
Visualise flux_test.dat - the single-timestep comparison of every radiative
transfer scheme against n-stream DISORT.

Produces
  flux_test_lw_error.png   longwave  error maps over (ssa, g)
  flux_test_sw_error.png   shortwave error maps over (ssa, g)
  flux_test_profiles.png   vertical profiles of heating-rate error

Run flux_test.py first to generate flux_test.dat.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm

ROOT = os.path.dirname(os.path.abspath(__file__))
DAT = os.path.join(ROOT, 'flux_test.dat')

# ---------------------------------------------------------------- palette
# Validated categorical slots, assigned in fixed order and held constant
# across every figure so a scheme keeps its colour everywhere.
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
         '#e87ba4', '#008300', '#4a3aa7', '#e34948', '#8c5a2b', '#5f8f99',
         '#7b61a8', '#00a5a8', '#b45f06', '#6a8f2a']
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#dcdbd7', '#fcfcfb'

# Sequential: one hue, light to dark (documented blue ramp, steps 100-700).
SEQ = LinearSegmentedColormap.from_list(
    'seq_blue', ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5',
                 '#256abf', '#184f95', '#0d366b'])
# Diverging: warm/cool poles with a neutral grey midpoint.
DIV = LinearSegmentedColormap.from_list(
    'div_br', ['#104281', '#2a78d6', '#9ec5f4', '#f0efec',
               '#f6b3b3', '#e34948', '#8c2020'])

LW_ORDER = ['lw_sc_scat', 'lw_sc_scat_ric',
            'lw_sc_iter_ali', 'lw_sc_iter_ng', 'lw_sc_iter_gs',
            'lw_sc_iter_krylov', 'lw_VIM_E', 'lw_Toon', 'lw_SH2', 'lw_SH4',
            'lw_disort_ts', 'lw_DFE']
SW_ORDER = ['sw_direct', 'sw_adding', 'sw_SDA', 'sw_Toon', 'sw_disort_ts']
COLOUR = {s: SLOTS[i] for i, s in enumerate(LW_ORDER)}
COLOUR.update({s: SLOTS[i] for i, s in enumerate(SW_ORDER)})

# Secondary encoding. Several schemes coincide almost exactly (the direct,
# Riccati, and converged iterative SC variants; sw_Toon and sw_disort_ts), so without it one curve simply hides
# another. It also keeps the lines separable for colour-vision-deficient
# readers and in print.
DASH = [(0, ()), (0, (5, 1.6)), (0, (1.4, 1.4)), (0, (7, 1.6, 1.4, 1.6)),
        (0, (3, 1.3)), (0, (1.4, 1.2, 4, 1.2)), (0, (9, 2))]
STYLE = {s: DASH[i % len(DASH)] for i, s in enumerate(LW_ORDER)}
STYLE.update({s: DASH[i % len(DASH)] for i, s in enumerate(SW_ORDER)})

LW_REF, SW_REF = 'lw_disort', 'sw_disort'

# Profile cases chosen per band: the shortwave has nothing to show at ssa = 0
# (no scattering to get wrong), so it gets its own more informative set.
LW_CASES = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.6), (0.9, 0.9)]
SW_CASES = [(0.25, 0.0), (0.5, 0.0), (0.75, 0.6), (0.9, 0.9)]


def load(path):
    """Read the long-format table into {(band,scheme,ssa,g): arrays}."""
    meta, rec = {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                for key in ('grav', 'cp_air'):
                    if key + ' = ' in line:
                        meta[key] = float(line.split(key + ' = ')[1].split()[0])
                continue
            p = line.split()
            if len(p) < 10:
                continue
            key = (p[0], p[1], float(p[2]), float(p[3]))
            rec.setdefault(key, []).append(
                [float(p[5]), float(p[6]), float(p[7]), float(p[8]), float(p[9])])
    out = {}
    for k, v in rec.items():
        a = np.array(v)
        out[k] = dict(pe=a[:, 0], tau=a[:, 1], up=a[:, 2],
                      down=a[:, 3], net=a[:, 4])
    return out, meta


def heating(pe, net, grav, cp):
    """Radiative heating rate per layer [K/day] from the net flux profile."""
    return (grav / cp) * np.diff(net) / np.diff(pe) * 86400.0


def layer_p(pe):
    return np.diff(pe) / np.log(pe[1:] / pe[:-1])


def symlog_ticks(lim, thr):
    """Ticks for a symmetric symlog colour scale."""
    t, d = [0.0], thr
    while d <= lim:
        t = [-d] + t + [d]
        d *= 10
    return t


def error_maps(data, meta, band, order, ref, fname, toa_label, title, note):
    ssas = sorted({k[2] for k in data if k[0] == band})
    gs = sorted({k[3] for k in data if k[0] == band})
    schemes = [s for s in order
               if sum(1 for k in data if k[0] == band and k[1] == s)
               == len(ssas) * len(gs)]

    toa, dq = {}, {}
    for s in schemes:
        A = np.full((len(ssas), len(gs)), np.nan)
        B = np.full((len(ssas), len(gs)), np.nan)
        for i, a in enumerate(ssas):
            for j, g in enumerate(gs):
                k, kr = (band, s, a, g), (band, ref, a, g)
                if k not in data or kr not in data:
                    continue
                ub, uk = data[kr]['up'][0], data[k]['up'][0]
                # A relative error in the reflected flux is undefined when the
                # benchmark reflects nothing at all (shortwave, ssa = 0).
                if abs(ub) > 1e-6:
                    A[i, j] = (uk - ub) / ub * 100.0
                qk = heating(data[k]['pe'], data[k]['net'], meta['grav'], meta['cp_air'])
                qr = heating(data[kr]['pe'], data[kr]['net'], meta['grav'], meta['cp_air'])
                B[i, j] = np.max(np.abs(qk - qr))
        toa[s], dq[s] = A, B

    alltoa = np.concatenate([toa[s].ravel() for s in schemes])
    lim = float(np.nanmax(np.abs(alltoa)))
    qmax = float(np.nanmax(np.concatenate([dq[s].ravel() for s in schemes])))
    # Errors span three orders of magnitude, so a linear scale would hide all
    # structure in the accurate schemes. Symlog keeps small errors legible.
    nrm_t = SymLogNorm(linthresh=1.0, vmin=-lim, vmax=lim, base=10)
    nrm_q = SymLogNorm(linthresh=1.0, vmin=0.0, vmax=qmax, base=10)

    n = len(schemes)
    fig, axes = plt.subplots(2, n, figsize=(2.15 * n + 1.8, 6.9))
    fig.patch.set_facecolor(SURF)
    if n == 1:
        axes = axes.reshape(2, 1)

    for row, (vals, cmap, norm, lab) in enumerate([
            (toa, DIV, nrm_t, toa_label),
            (dq, SEQ, nrm_q, 'max |heating-rate error|  [K/day]')]):
        for c, s in enumerate(schemes):
            ax = axes[row, c]
            ax.imshow(vals[s], cmap=cmap, norm=norm, aspect='auto',
                      origin='lower', interpolation='nearest')
            # Numeric cell labels are the direct-label relief for the three
            # low-contrast slots, and double as a readable table.
            for i in range(len(ssas)):
                for j in range(len(gs)):
                    v = vals[s][i, j]
                    if np.isnan(v):
                        ax.text(j, i, 'n/a', ha='center', va='center',
                                fontsize=6.2, color=INK2, style='italic')
                        continue
                    rgba = cmap(norm(v))
                    dark = (0.299 * rgba[0] + 0.587 * rgba[1]
                            + 0.114 * rgba[2]) < 0.55
                    txt = ('%.0f' if abs(v) >= 100 else '%.1f') % v
                    ax.text(j, i, txt, ha='center', va='center', fontsize=6.4,
                            color='#ffffff' if dark else INK)
            ax.set_xticks(range(len(gs)))
            ax.set_xticklabels(['%g' % g for g in gs])
            ax.set_yticks(range(len(ssas)))
            ax.set_yticklabels(['%g' % a for a in ssas] if c == 0 else [])
            ax.tick_params(colors=INK2, labelsize=8, length=0)
            for sp in ax.spines.values():
                sp.set_color(GRID)
            if row == 0:
                ax.set_title(s.replace(band + '_', ''), fontsize=9.5,
                             color=INK, pad=6)
            if c == 0:
                ax.set_ylabel('single scattering albedo', fontsize=8.5,
                              color=INK2)
            ax.set_xlabel('asymmetry  g', fontsize=8.5, color=INK2)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cax = fig.add_axes([0.915, 0.565 - 0.475 * row, 0.014, 0.34])
        ticks = symlog_ticks(lim, 1.0) if row == 0 else \
            [0.0] + [10.0 ** e for e in range(0, int(np.ceil(np.log10(qmax))) + 1)]
        cb = fig.colorbar(sm, cax=cax, ticks=ticks)
        cb.set_label(lab, fontsize=8.5, color=INK2)
        cb.ax.tick_params(colors=INK2, labelsize=7.5, length=3)
        cb.ax.set_yticklabels(['%g' % t for t in ticks])
        cb.outline.set_edgecolor(GRID)

    fig.suptitle(title, fontsize=13, color=INK, x=0.45, y=0.985)
    fig.text(0.45, 0.938, note, fontsize=8.5, color=INK2, ha='center')
    fig.subplots_adjust(left=0.055, right=0.895, top=0.885, bottom=0.075,
                        wspace=0.09, hspace=0.30)
    fig.savefig(os.path.join(ROOT, fname), dpi=160, facecolor=SURF)
    plt.close(fig)
    print('wrote', fname)


def profiles(data, meta, fname):
    rows = [('lw', LW_ORDER, LW_REF, LW_CASES),
            ('sw', SW_ORDER, SW_REF, SW_CASES)]
    ncol = max(len(r[3]) for r in rows)
    fig, axes = plt.subplots(2, ncol, figsize=(3.15 * ncol + 1.5, 8.8),
                             sharey=True)
    fig.patch.set_facecolor(SURF)

    for row, (band, order, ref, cases) in enumerate(rows):
        handles, labels = [], []
        for c, (ssa, g) in enumerate(cases):
            ax = axes[row, c]
            ax.set_facecolor('none')
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)
            for sp in ('left', 'bottom'):
                ax.spines[sp].set_color(GRID)
            ax.tick_params(colors=INK2, labelsize=8, length=3, width=0.8)
            ax.grid(True, color=GRID, lw=0.6, alpha=0.9)
            ax.set_axisbelow(True)

            kr = (band, ref, ssa, g)
            if kr not in data:
                ax.set_visible(False)
                continue
            pr = layer_p(data[kr]['pe']) / 1e5
            qr = heating(data[kr]['pe'], data[kr]['net'],
                         meta['grav'], meta['cp_air'])
            ax.axvline(0.0, color=INK2, lw=1.0, alpha=0.55)
            worst = 0.0
            for s in order:
                k = (band, s, ssa, g)
                if k not in data:
                    continue
                q = heating(data[k]['pe'], data[k]['net'],
                            meta['grav'], meta['cp_air'])
                ln, = ax.plot(q - qr, pr, lw=2.0, color=COLOUR[s],
                              linestyle=STYLE[s],
                              label=s.replace(band + '_', ''),
                              solid_capstyle='round')
                worst = max(worst, float(np.max(np.abs(q - qr))))
                if c == 0 or s not in labels:
                    if s.replace(band + '_', '') not in labels:
                        handles.append(ln)
                        labels.append(s.replace(band + '_', ''))
            ax.set_yscale('log')
            ax.set_xscale('symlog', linthresh=0.01)
            L = 1.6 * max(worst, 0.05)
            ax.set_xlim(-L, L)
            # One tick per two decades, otherwise the negative labels collide
            # around zero on a symlog axis.
            dec = [10.0 ** e for e in range(0, int(np.ceil(np.log10(L))) + 1, 2)]
            ax.set_xticks([-d for d in reversed(dec)] + [0.0] + dec)
            ax.invert_yaxis()
            ax.set_title('%s   ssa = %g,  g = %g' % (band.upper(), ssa, g),
                         fontsize=10, color=INK, pad=7)
            if c == 0:
                ax.set_ylabel('pressure  [bar]', fontsize=9, color=INK2)
            ax.set_xlabel('heating-rate error  [K/day]', fontsize=9, color=INK2)
        # One legend per row rather than one per panel.
        fig.legend(handles, labels, fontsize=8.5, frameon=False,
                   labelcolor=INK2, handlelength=1.8, ncol=1,
                   loc='center left', bbox_to_anchor=(0.885, 0.755 - 0.487 * row))

    fig.suptitle('Heating-rate error against n-stream DISORT',
                 fontsize=13, color=INK, x=0.45, y=0.985)
    fig.text(0.45, 0.952,
             'single timestep from an identical T-p profile, so every '
             'difference is radiative transfer error',
             fontsize=8.5, color=INK2, ha='center')
    fig.subplots_adjust(left=0.062, right=0.875, top=0.912, bottom=0.062,
                        wspace=0.16, hspace=0.26)
    fig.savefig(os.path.join(ROOT, fname), dpi=160, facecolor=SURF)
    plt.close(fig)
    print('wrote', fname)


def main():
    if not os.path.exists(DAT):
        sys.exit('%s not found - run flux_test.py first' % DAT)
    data, meta = load(DAT)
    meta.setdefault('grav', 10.0)
    meta.setdefault('cp_air', 1.3e4)
    error_maps(data, meta, 'lw', LW_ORDER, LW_REF, 'flux_test_lw_error.png',
               'OLR error  [%]',
               'Longwave scheme error against n-stream DISORT',
               'symlog colour scale, linear within +/-1 - single timestep, '
               'identical T-p profile for every scheme')
    error_maps(data, meta, 'sw', SW_ORDER, SW_REF, 'flux_test_sw_error.png',
               'reflected flux error  [%]',
               'Shortwave scheme error against n-stream DISORT',
               'symlog colour scale, linear within +/-1 - at ssa = 0 nothing '
               'is reflected, so a relative error is undefined')
    profiles(data, meta, 'flux_test_profiles.png')


if __name__ == '__main__':
    main()

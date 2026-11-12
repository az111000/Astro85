"""
Project 20, stretch goal: is there a large-scale "dipole" in cosmic rays above 8 EeV?

A dipole means slightly MORE cosmic rays arrive from one half of the sky than from the opposite half.
In 2017 the Auger collaboration announced exactly this, using ~32,000 events above 8 EeV: about 6.5%
more particles from one direction, which points away from the centre of our galaxy, so those
particles are probably coming from other galaxies (Science 357, 1266).

Here we use Auger's PUBLIC data release (roughly 10% of their events) to see how much of that signal
we can find ourselves.

Method ("Rayleigh analysis", the standard first test):
  Auger sees every right ascension (RA) equally over many years, so if the sky had no dipole,
  the RAs of events would be spread evenly around 0-360 degrees. We fit
        rate(RA) = average * (1 + r * cos(RA - phase))
  r = how strong the east-west wobble is, phase = which RA gets the most events.
  The chance of getting an r this big from a perfectly smooth sky is  P = exp(-N r^2 / 4).

Run from the repo root:  python3 project20/dipole_8EeV.py
Data: Pierre Auger Observatory public data release, SD1500 events, table auger3.ev_sd1500,
      served by GAVO's TAP service https://dc.g-vo.org/tap (see also https://opendata.auger.org).
"""
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parent.parent
DATA, PLOTS = ROOT / 'data', ROOT / 'plots'
DATA.mkdir(exist_ok=True); PLOTS.mkdir(exist_ok=True)
BLUE, ORANGE, GRAY, INK = '#2a78d6', '#eb6834', '#9a9993', '#52514e'

# Auger's published result for E >= 8 EeV (Science 2017), first harmonic in right ascension
AUGER_R, AUGER_PHASE = 0.047, 100.0      # amplitude, phase in degrees
AUGER_DIPOLE_GAL = (233.0, -13.0)        # 3-D dipole direction, Galactic (l, b) in degrees

# ---- 1. Download events above 4 EeV once (about 1 MB), then reuse the copy in data/ ----
cache = DATA / 'auger_sd1500_above4EeV.csv'
if not cache.exists():
    import pyvo
    tap = pyvo.dal.TAPService('https://dc.g-vo.org/tap')
    query = ('SELECT id, gpstime, sd_energy, sd_theta, sd_ra, sd_dec, sd_l, sd_b '
             'FROM auger3.ev_sd1500 WHERE sd_energy >= 4')
    tap.run_sync(query, maxrec=100000).to_table().to_pandas().to_csv(cache, index=False)
ev = pd.read_csv(cache)
print(f'{len(ev)} events >= 4 EeV, zenith angles up to {ev.sd_theta.max():.0f} deg')


def rayleigh(ra_deg):
    """First-harmonic (Rayleigh) analysis in right ascension. Returns amplitude r, phase (deg), chance P."""
    a = np.radians(ra_deg)
    n = len(a)
    A, B = 2 / n * np.cos(a).sum(), 2 / n * np.sin(a).sum()
    r = np.hypot(A, B)
    return r, np.degrees(np.arctan2(B, A)) % 360, np.exp(-n * r**2 / 4)


bins_cfg = {'8 EeV and above': ev.sd_energy >= 8, '4 to 8 EeV (comparison)': ev.sd_energy < 8}
results = {}
for label, sel in bins_cfg.items():
    ra = ev.sd_ra[sel].to_numpy()
    r, phase, P = rayleigh(ra)
    sigma = np.sqrt(2 / len(ra))           # typical size of r from random noise alone
    results[label] = (ra, r, phase, P, sigma)
    print(f'{label:>24}: N={len(ra):5d}  r={r:.3f} +/- {sigma:.3f}  phase={phase:5.0f} deg  chance P={P:.3f}')

# ---- 2. Figure 1: event rate vs right ascension, with the fitted cosine ----
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
edges = np.arange(0, 361, 30)
centers = (edges[:-1] + edges[1:]) / 2
x = np.linspace(0, 360, 361)
for ax, (label, (ra, r, phase, P, sigma)) in zip(axes, results.items()):
    counts, _ = np.histogram(ra, bins=edges)
    mean = counts.mean()
    ax.errorbar(centers, counts / mean, yerr=np.sqrt(counts) / mean, fmt='o', color=INK, ms=6,
                capsize=0, label=f'events per 30° of RA (N = {len(ra):,})')
    ax.plot(x, 1 + r * np.cos(np.radians(x - phase)), color=ORANGE, lw=2,
            label=f'fit: r = {r*100:.1f}% ± {sigma*100:.1f}%, phase {phase:.0f}°')
    if label.startswith('8'):
        ax.plot(x, 1 + AUGER_R * np.cos(np.radians(x - AUGER_PHASE)), color=BLUE, lw=1.5, ls='--',
                label='Auger 2017 (full data): r = 4.7%, phase 100°')
    ax.axhline(1, color=GRAY, lw=0.8)
    ax.set_title(f'{label}:  chance probability P = {P:.2f}', fontsize=11)
    ax.set_xlabel('Right ascension [degrees]')
    ax.set_xticks(range(0, 361, 60))
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5, loc='lower left')
axes[0].set_ylabel('Event rate relative to average')
axes[0].set_ylim(0.78, 1.2)
fig.text(0.99, -0.01, 'Data: Pierre Auger Observatory public data (auger3.ev_sd1500 via GAVO TAP, dc.g-vo.org)',
         ha='right', fontsize=7, color=GRAY)
fig.tight_layout()
out = PLOTS / 'cr_dipole_ra_harmonic.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print('saved', out.relative_to(ROOT))

# ---- 3. Figure 2: smoothed sky map of the excess above 8 EeV ----
# For each spot on the sky, count events within 45 degrees and compare with the average of
# shuffled skies (same trick as the Cen A prototype: keep declinations, randomise RA).
SMOOTH, N_SKIES = 45.0, 200
rng = np.random.default_rng(85)
hi = ev[ev.sd_energy >= 8]
ra_ev, dec_ev = np.radians(hi.sd_ra.to_numpy()), np.radians(hi.sd_dec.to_numpy())

# grid of sky positions in Galactic coordinates, converted to RA/Dec for the counting
lon = np.linspace(-180, 180, 145); lat = np.linspace(-90, 90, 73)
LON, LAT = np.meshgrid(lon, lat)
grid = SkyCoord(l=LON.ravel(), b=LAT.ravel(), unit='deg', frame='galactic').icrs
g_ra, g_dec = grid.ra.radian, grid.dec.radian
cos_r = np.cos(np.radians(SMOOTH))


def count_near(ev_ra, ev_dec):
    """For every grid point, the number of events within SMOOTH degrees."""
    c = (np.sin(g_dec)[:, None] * np.sin(ev_dec)[None, :] +
         np.cos(g_dec)[:, None] * np.cos(ev_dec)[None, :] * np.cos(g_ra[:, None] - ev_ra[None, :]))
    return (c > cos_r).sum(axis=1)


observed = count_near(ra_ev, dec_ev)
expected = np.mean([count_near(rng.uniform(0, 2 * np.pi, len(ra_ev)), dec_ev) for _ in range(N_SKIES)], axis=0)
with np.errstate(divide='ignore', invalid='ignore'):
    # hide parts of the sky Auger barely sees (less than 30% of its best coverage): too few events there to trust
    excess = np.where(expected > 0.3 * expected.max(), observed / expected - 1, np.nan)
excess = excess.reshape(LAT.shape)

fig = plt.figure(figsize=(10, 5.8))
ax = fig.add_subplot(projection='mollweide')
lim = 0.15   # fixed colour range of +/-15%
mesh = ax.pcolormesh(-np.radians(LON), np.radians(LAT), excess * 100, cmap='RdBu_r', vmin=-lim * 100,
                     vmax=lim * 100, shading='auto')
l0, b0 = AUGER_DIPOLE_GAL
ax.plot(-np.radians((l0 + 180) % 360 - 180), np.radians(b0), marker='X', ms=13, color='black',
        markeredgecolor='white', ls='none', label='Auger 2017 dipole direction (full data)')
ax.plot(0, 0, marker='+', ms=14, mew=2, color='black', ls='none', label='Galactic centre')
ax.set_xticks(np.radians(range(-150, 151, 30)), labels=[f'{-t % 360}°' for t in range(-150, 151, 30)], fontsize=8)
ax.tick_params(axis='y', labelsize=8)
ax.grid(True, color='#bbbbbb', lw=0.5)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False, fontsize=9)
ax.set_title(f'Cosmic rays ≥ 8 EeV: excess over a smooth sky ({SMOOTH:.0f}° smoothing, N = {len(hi):,})\n'
             'Galactic coordinates; white = sky Auger sees too little of', fontsize=11, pad=14)
fig.colorbar(mesh, orientation='horizontal', shrink=0.5, pad=0.1, extend='both', label='Excess relative to smooth sky [%]')
out = PLOTS / 'cr_dipole_skymap.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
i = np.nanargmax(excess)
print(f'hottest spot on map: l={LON.ravel()[i] % 360:.0f}, b={LAT.ravel()[i]:.0f}, excess {np.nanmax(excess)*100:.1f}%')
print('saved', out.relative_to(ROOT))

"""
Project 20, Thursday prototype: are Auger's highest-energy cosmic rays bunched up around Centaurus A?

The idea, in plain words:
  1. Count how many of Auger's top events landed within 25 degrees of Centaurus A.
  2. Make 1,000 "fake skies" where each event keeps its declination (how far north/south it is)
     but gets a random right ascension (east/west position). Auger sees the sky by letting the Earth
     spin it past, so over many years every right ascension gets the same amount of viewing time,
     while declination does not (Auger is in Argentina and mostly sees the south). Keeping the
     declinations means the fake skies have exactly Auger's view of the sky, just with no sources.
  3. The p-value is the fraction of fake skies with at least as many events near Cen A as the real sky.
     Small p-value = hard to get the real count by luck.

Run from the repo root:  python3 project20/cena_significance.py
Data: Pierre Auger Observatory, catalog of its highest-energy cosmic rays (2004-2020),
      https://opendata.auger.org/catalog/
"""
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parent.parent
DATA, PLOTS = ROOT / 'data', ROOT / 'plots'
DATA.mkdir(exist_ok=True); PLOTS.mkdir(exist_ok=True)

RADIUS = 25.0                     # degrees; chosen BEFORE looking (close to what Auger itself uses for Cen A)
N_SKIES = 1000                    # number of shuffled skies, as in the project brief
CEN_A = (201.37, -43.02)          # Centaurus A (RA, Dec) in degrees
rng = np.random.default_rng(85)   # fixed seed so you get the same answer every run

# ---- 1. Load the catalog (download once, then reuse the copy in data/) ----
cache = DATA / 'auger_top_events.csv'
if not cache.exists():
    pd.read_csv('https://opendata.auger.org/catalog/list/auger_catalogSD.csv').to_csv(cache, index=False)
df = pd.read_csv(cache)
ra, dec = df['ra'].to_numpy(), df['dec'].to_numpy()
print(f'{len(df)} events, energies {df.energy.min():.0f}-{df.energy.max():.0f} EeV')


def angle_to(ra, dec, ra0, dec0):
    """Angle on the sky (degrees) between points (ra, dec) and one point (ra0, dec0). Works on arrays."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cos_sep = np.sin(dec) * np.sin(dec0) + np.cos(dec) * np.cos(dec0) * np.cos(ra - ra0)
    return np.degrees(np.arccos(np.clip(cos_sep, -1, 1)))


# ---- 2. The real count ----
sep_real = angle_to(ra, dec, *CEN_A)
n_real = int((sep_real < RADIUS).sum())

# ---- 3. 1,000 shuffled skies: same declinations, random right ascensions ----
fake_ra = rng.uniform(0, 360, size=(N_SKIES, len(df)))          # one row per fake sky
sep_fake = angle_to(fake_ra, dec[None, :], *CEN_A)
n_fake = (sep_fake < RADIUS).sum(axis=1)

expected = n_fake.mean()
p_value = (np.sum(n_fake >= n_real) + 1) / (N_SKIES + 1)        # +1s: the real sky counts as one possible sky
print(f'Within {RADIUS:.0f} deg of Cen A: observed {n_real}, expected {expected:.1f} '
      f'(fake skies range {n_fake.min()}-{n_fake.max()})')
print(f'{np.sum(n_fake >= n_real)} of {N_SKIES} fake skies had >= {n_real} -> p = {p_value:.3f}')

# ---- 4. Figure: sky map (left) + "how often does luck do this?" histogram (right) ----
BLUE, ORANGE, GRAY = '#2a78d6', '#eb6834', '#9a9993'


def galactic(ra, dec):   # RA/Dec -> Galactic lon/lat in radians, longitude increasing to the left (astronomy style)
    g = SkyCoord(ra, dec, unit='deg').galactic
    return -g.l.wrap_at('180d').radian, g.b.radian


fig = plt.figure(figsize=(15, 5.6))
ax = fig.add_subplot(1, 2, 1, projection='mollweide')
inside = sep_real < RADIUS
ax.scatter(*galactic(ra[~inside], dec[~inside]), s=df.energy[~inside] / 3, color=GRAY, alpha=0.7,
           edgecolor='white', linewidth=0.5, label=f'other events ({(~inside).sum()})')
ax.scatter(*galactic(ra[inside], dec[inside]), s=df.energy[inside] / 3, color=ORANGE,
           edgecolor='white', linewidth=0.5, label=f'within {RADIUS:.0f}° of Cen A ({n_real})')
circle = SkyCoord(*CEN_A, unit='deg').directional_offset_by(np.linspace(0, 360, 200) * u.deg, RADIUS * u.deg)
ax.plot(*galactic(circle.ra.deg, circle.dec.deg), color=ORANGE, lw=1.5)
ax.plot(*galactic(*CEN_A), marker='*', ms=15, color=BLUE, markeredgecolor='white', ls='none', label='Centaurus A')
ax.set_xticks(np.radians(range(-150, 151, 30)), labels=[f'{-t % 360}°' for t in range(-150, 151, 30)], fontsize=8)
ax.tick_params(axis='y', labelsize=8)
ax.grid(True, color='#dddddd')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=9, frameon=False)
ax.set_title(f'Auger\'s {len(df)} highest-energy events (Galactic coordinates; dot size = energy)', pad=16, fontsize=11)

ax2 = fig.add_subplot(1, 2, 2)
bins = np.arange(n_fake.min(), max(n_fake.max(), n_real) + 2) - 0.5
ax2.hist(n_fake, bins=bins, color=BLUE, edgecolor='white', linewidth=2, label=f'{N_SKIES:,} shuffled skies')
ax2.axvline(n_real, color=ORANGE, lw=2.5, label=f'real sky: {n_real} events')
ax2.axvline(expected, color='#52514e', lw=1, ls='--', label=f'average shuffled sky: {expected:.1f}')
ax2.set_xlabel(f'Number of events within {RADIUS:.0f}° of Centaurus A')
ax2.set_ylabel('Number of shuffled skies')
ax2.set_title(f'{np.sum(n_fake >= n_real)} of {N_SKIES:,} shuffled skies have {n_real} or more events near Cen A  (p ≈ {p_value:.2f})',
              fontsize=11)
ax2.spines[['top', 'right']].set_visible(False)
ax2.legend(frameon=False, fontsize=9)
fig.text(0.99, 0.01, 'Data: Pierre Auger Observatory open data, opendata.auger.org', ha='right', fontsize=7, color=GRAY)
fig.tight_layout()
out = PLOTS / 'cr_cena_excess.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print('saved', out.relative_to(ROOT))

# ---- 5. Bonus figure: what if we had picked a different circle size? ----
# Warning: trying many radii and keeping the best one makes the result look stronger than it is
# ("look-elsewhere effect"). The 25 deg answer above is the honest one; this plot is for exploring.
radii = np.arange(5, 61, 1)
obs = np.array([(sep_real < r).sum() for r in radii])
fake_counts = np.array([(sep_fake < r).sum(axis=1) for r in radii])        # shape (n_radii, N_SKIES)
lo, mid, hi = np.percentile(fake_counts, [2.5, 50, 97.5], axis=1)
p_scan = ((fake_counts >= obs[:, None]).sum(axis=1) + 1) / (N_SKIES + 1)

fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
a1.fill_between(radii, lo, hi, color=BLUE, alpha=0.2, linewidth=0, label='95% of shuffled skies')
a1.plot(radii, mid, color=BLUE, lw=2, label='typical shuffled sky')
a1.plot(radii, obs, color=ORANGE, lw=2, label='real sky')
a1.set_ylabel('Events within radius of Cen A')
a1.legend(frameon=False, fontsize=9)
a1.set_title('Event count around Centaurus A vs. circle size', fontsize=11)
a2.semilogy(radii, p_scan, color=ORANGE, lw=2)
a2.set_yticks([0.01, 0.05, 0.1, 0.2, 0.5, 1], labels=['0.01', '0.05', '0.1', '0.2', '0.5', '1']); a2.minorticks_off()
a2.axhline(0.05, color=GRAY, ls=':', lw=1); a2.text(59, 0.055, 'p = 0.05', ha='right', fontsize=8, color=GRAY)
a2.axvline(RADIUS, color='#52514e', ls='--', lw=1); a2.text(RADIUS + 0.5, 0.5, 'pre-chosen 25°', fontsize=8)
a2.set_xlabel('Circle radius around Centaurus A [degrees]')
a2.set_ylabel('p-value (one radius)')
for a in (a1, a2):
    a.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
out = PLOTS / 'cr_cena_radius_scan.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
best = radii[np.argmin(p_scan)]
print(f'smallest p in scan: {p_scan.min():.4f} at {best} deg (not a fair significance, see comment)')
print('saved', out.relative_to(ROOT))

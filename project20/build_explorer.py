"""
Build the interactive results page for project 20.

It puts the real Auger data (from data/, downloaded by the other two project20 scripts) into
explorer_template.html and writes project20/results_explorer.html. The page then runs both tests
live in the browser.

Run from the repo root, after cena_significance.py and dipole_8EeV.py:
    python3 project20/build_explorer.py
"""
import json
from pathlib import Path
import pandas as pd
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parent.parent
HERE = ROOT / 'project20'

top = pd.read_csv(ROOT / 'data' / 'auger_top_events.csv')
g = SkyCoord(top.ra, top.dec, unit='deg').galactic
big = pd.read_csv(ROOT / 'data' / 'auger_sd1500_above4EeV.csv')

data = {
    'top': {'ra': top.ra.round(2).tolist(), 'dec': top.dec.round(2).tolist(),
            'l': [round(x, 2) for x in g.l.deg], 'b': [round(x, 2) for x in g.b.deg],
            'e': top.energy.astype(int).tolist(), 't': top.utctime.astype(int).tolist()},
    'big': {'ra': big.sd_ra.round(2).tolist(), 'dec': big.sd_dec.round(2).tolist(),
            'e': big.sd_energy.round(3).tolist()},
}
html = (HERE / 'explorer_template.html').read_text().replace('__DATA__', json.dumps(data, separators=(',', ':')))
out = HERE / 'results_explorer.html'
out.write_text(html)
print(f'wrote {out.relative_to(ROOT)} ({len(html) / 1024:.0f} KB): {len(top)} top events, {len(big)} SD1500 events')

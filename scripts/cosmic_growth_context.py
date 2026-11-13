"""Conditional cosmic clocks; these scenarios are not measured object histories."""

from pathlib import Path

import astropy
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
from matplotlib.figure import Figure

from cosmic_anomaly_lab.open_search import growth_efold_myr
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/open_search"
out.mkdir(parents=True, exist_ok=True)
formation_z = 20
formation_age = float(Planck18.age(formation_z).value)
rows = []
for z in [3, 5, 7, 10, 12]:
    age = float(Planck18.age(z).value)
    for seed in [100, 1e4, 1e5]:
        for efficiency in [0.057, 0.1, 0.2]:
            for duty in [0.3, 1.0]:
                for eddington in [0.3, 1.0]:
                    efold = growth_efold_myr(efficiency, eddington, duty)
                    available = (age - formation_age) * 1000
                    required = efold * np.log(1e8 / seed)
                    rows.append(
                        dict(
                            redshift=z,
                            age_Gyr=age,
                            formation_redshift=formation_z,
                            seed_mass_Msun=seed,
                            target_mass_Msun=1e8,
                            efficiency=efficiency,
                            duty_cycle=duty,
                            eddington_ratio=eddington,
                            efold_Myr=efold,
                            available_Myr=available,
                            required_Myr=float(required),
                            reaches_hypothetical_target=bool(required <= available),
                        )
                    )
pd.DataFrame(rows).to_csv(out / "conditional_growth_grid.csv", index=False)
summary = dict(
    executed_utc=utcnow(),
    astropy_version=astropy.__version__,
    cosmology=str(Planck18),
    formation_redshift=formation_z,
    formation_age_Gyr=formation_age,
    default_efold_Myr=growth_efold_myr(),
    formula="t_e = (c sigma_T / (4 pi G m_p)) * epsilon / ((1-epsilon) lambda duty); "
    "t_required = t_e ln(M_target/M_seed)",
    trials=len(rows),
    physical_classification="NOT_INFERRED",
    limitations=[
        "No measured mass inferred for any source",
        "Constant efficiency and averaged duty",
        "No finite gas supply, feedback, mergers or structure evolution",
        "No inference of future fate or new cosmological physics",
        "Classical radiatively efficient growth is not a quasi-star envelope solver",
    ],
    script_sha256=text_digest(Path(__file__)),
    analysis_sha256=text_digest(root / "src/cosmic_anomaly_lab/open_search.py"),
)
write_json(out / "cosmic_context.json", summary)
fig = Figure(figsize=(12, 5), layout="constrained")
ax, bx = fig.subplots(1, 2)
z = np.linspace(0.5, 20, 300)
ax.plot(z, Planck18.age(z).value, color="#216e79")
ax.scatter([3, 5, 7, 10, 12], Planck18.age([3, 5, 7, 10, 12]).value, color="#a65b35")
ax.set(xlabel="Redshift", ylabel="Cosmic age (Gyr)", title="Cosmic clock · Planck18 assumption")
ax.grid(alpha=0.2)
age = np.linspace(formation_age, 1.6, 400)
for seed, duty, color in [
    (100, 1.0, "#236d93"),
    (1e5, 1.0, "#a75a35"),
    (100, 0.3, "#81a9bd"),
    (1e5, 0.3, "#dbad87"),
]:
    log_mass = np.log10(seed) + (age - formation_age) * 1000 / growth_efold_myr(
        duty_cycle=duty
    ) / np.log(10)
    mask = log_mass <= 9.5
    bx.plot(
        age[mask],
        log_mass[mask],
        color=color,
        ls="-" if duty == 1 else "--",
        label=f"Seed {seed:g} solar masses; duty {duty:g}",
    )
bx.axhline(
    8, color="#444444", lw=0.8, ls=":", label="Illustrative target: 100 million solar masses"
)
bx.set(
    xlabel="Cosmic age (Gyr)",
    ylabel="log₁₀ mass / solar mass",
    title="SCENARIOS · ε=0.1, L/L_Edd=1, seeds at z=20",
    ylim=(1.8, 9.5),
)
bx.legend(fontsize=8, loc="lower right")
bx.grid(alpha=0.2)
fig.suptitle(
    "Time available for growth · no masses measured, no object fate predicted", fontsize=13
)
save_figure(fig, root / "figures/cosmic_growth_context")
print(summary)

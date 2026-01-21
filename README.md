# Sentinel

**Sentinel** is an automated daily situational-awareness system for
**viromics, metagenomics, and evolutionary biology**.

It scans the scientific literature every day, extracts newly published
signals relevant to global virome surveillance, and publishes a clean,
public-facing summary via GitHub Pages.

---

## 🌐 Live Dashboard
https://<your-username>.github.io/Sentinel/

---

## 🎯 What Sentinel Does
- Queries PubMed daily for:
  - viromes & viral metagenomics
  - wastewater & urban surveillance
  - influenza & H5N1 evolution
  - human & environmental viromes
  - AMR metagenomics
- Produces:
  - `latest.json` — machine-readable
  - `latest.md` — human-readable
  - `index.html` — public dashboard
- Updates automatically via GitHub Actions

---

## 🧬 Design Philosophy
- **Predictive evolutionary mindset**
- **Reproducible & auditable**
- **Static > fragile dashboards**
- **No backend, no servers, no keys**

This repo is intended as scientific infrastructure, not a demo.

---

## ⚙️ Automation
Runs daily using GitHub Actions.  
All outputs live in `/docs` and are immediately published via GitHub Pages.

---

## ⚠️ Disclaimer
For informational and research purposes only.  
Not medical, public-health, or policy advice.

---

## 👤 Author
Alexander G. Lucaci, PhD  
Computational Evolution • Viromics • Surveillance

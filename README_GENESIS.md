# Genesis probe plotting

This script clones the plotting approach from `nrho_plot.py` but samples a BSP for the Genesis probe.

Steps to use:

1. Install dependencies (preferably in a virtual environment):

```bash
pip install -r requirements.txt
```

2. Run the script:

```bash
python genesis_halo_plot.py
```

To write static HTML and a PNG preview (requires `kaleido`):

```bash
python genesis_halo_plot.py --static
```
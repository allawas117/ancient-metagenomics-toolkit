# Toolkit for Ancient Metagenomic Data (Currently in progress)

A modular Python toolkit for analyzing ancient metagenomic data — enabling streamlined workflows for damage analysis, contamination detection, database reference mapping, and visualization.

---

## Features

- 🔬 **Damage Analysis** — Identify and profile DNA damage patterns specific to ancient samples.  
- 🛡 **Authentication & Contamination Checks** — Validate sample authenticity and detect modern contamination.  
- 📚 **Reference Database Handling** — Effortlessly map to taxonomic references and manage databases.  
- 🛠 **Assembly** — Assemble ancient metagenomic sequences with specialized algorithms.  
- 📊 **Visualization** — Temporal exploration and visualization of ancient DNA data.  
- 🚀 **C++ Extensions** — High-performance computation for key tasks.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/AncientMetagenomicsToolkit.git
cd AncientMetagenomicsToolkit

# Install dependencies
pip install -r requirements.txt

# Optional: Install C++ extensions
cd amt/cpp_extensions
python setup.py install
```

Alternatively, use Docker:
```bash
docker build -t ancientmetagenomics .
docker run -it ancientmetagenomics
```

---

## Usage

Check out the examples and notebooks:
```bash
cd examples/notebooks
jupyter notebook damage_profiling.ipynb
```

Or use it programmatically:
```python
from amt.core import damage_analysis
damage_analysis.run_analysis("your_data_file")
```

---

## Directory Structure

- `amt/` — Main source code  
- `tests/` — Unit tests  
- `docs/` — Documentation and tutorials  
- `examples/` — Example datasets and Jupyter notebooks  

---

## Documentation

- [Installation Guide](docs/installation.md)  
- [Tutorials](docs/tutorials)  

---

## License

This project is licensed under the [MIT License](LICENSE).

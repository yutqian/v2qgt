# v2qgt 🚀  

This code **v2qgt** provides an efficient tool for calculating the **quantum geometric tensor (QGT)** directly from DFT calculations. The code is based on [VASPBERRY](https://github.com/Infant83/VASPBERRY) but includes important fixes and enhancements.  

## ✨ Features  
- Computes **Berry curvature** and **quantum metric tensor** directly from VASP's `WAVECAR`.  
- Implements the **Kubo formula** for accurate calculations.  
- Fixes bugs in **VASPBERRY** and extends functionality to support **both 2D and 3D materials**.  

## 🔧 Installation  
Simply compile using:  
```sh
make
```

## 📖 Citation
If you use **v2qgt** in your research, please cite our paper:

> **Sunje Kim**, **Yoonah Chung**, **Yuting Qian**, *et al.*,  
> **Direct measurement of the quantum metric tensor in solids**,  
> *Science*, **388**, 1050–1054 (2025).  
> [DOI: 10.1126/science.ado6049](https://www.science.org/doi/abs/10.1126/science.ado6049)

### BibTeX 
```bibtex
@article{doi:10.1126/science.ado6049,
  author = {Sunje Kim and Yoonah Chung and Yuting Qian and Soobin Park and Chris Jozwiak and Eli Rotenberg and Aaron Bostwick and Keun Su Kim and Bohm-Jung Yang},
  title = {Direct measurement of the quantum metric tensor in solids},
  journal = {Science},
  volume = {388},
  number = {6751},
  pages = {1050-1054},
  year = {2025},
  doi = {10.1126/science.ado6049},
  url = {https://www.science.org/doi/abs/10.1126/science.ado6049},
  eprint = {https://www.science.org/doi/pdf/10.1126/science.ado6049}
}

@software{yutqian_2025_14927897,
  author       = {yutqian},
  title        = {yutqian/v2qgt: v2qgt v1.0.0-alpha.1},
  month        = feb,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v1.0.0-alpha.1},
  doi          = {10.5281/zenodo.14927897},
  url          = {https://doi.org/10.5281/zenodo.14927897}
}


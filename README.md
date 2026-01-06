# Fast Recursive Projections (FRP)

This repository contains a C++ implementation of the Fast Recursive Projections (FRP) method for pricing American-style options under a variety of asset dynamics.

The code accompanies the paper:

Cosma, A., Galluccio, S., Pederzoli, P., & Scaillet, O.
Early Exercise Decision in American Options with Dividends, Stochastic Volatility, and Jumps
Journal of Financial and Quantitative Analysis, 55(1), 331–356 (2020)
DOI: 10.1017/S0022109018001229

The implementation supports American option pricing under:
- Black–Scholes dynamics
- Heston stochastic volatility
- Merton jump-diffusion

The repository provides both a shared library (libFRP.so) and a simple executable example illustrating how to use the pricing routines.

---

## Repository structure

Main files:
- FRP_*.{h,cpp}: model-specific FRP implementations (Black–Scholes, Heston, Merton)
- FRP_Pricing.{h,cpp}: pricing interface and orchestration logic
- main_pricing.cpp: example program used to run pricing experiments (edit the /INPUT PARAMETERS/ section in main())
- include/: headers for third-party numerical components
- lib/: directory where compiled artifacts are generated

---

## Dependencies

The code relies on the following external libraries:

- newmat (matrix algebra)
  http://www.robertnz.net/nm_intro.htm

- FFTW3 (Fast Fourier Transform)
  http://www.fftw.org

Some legacy references mention librecipes.a; in practice, the headers shipped in the include/ directory are sufficient for the provided build setup.

This is research code and is provided as is, without warranty.

---

## Build and run

The code is built using make.

Compile the shared library:

    make lib

Compile the example executable:

    make exec

Run the example:

    make run

The executable is built from main_pricing.cpp and linked against libFRP.so.

---

## Usage

To run your own experiments, edit the /INPUT PARAMETERS/ section at the beginning of main() in main_pricing.cpp.
That section is where you specify:
- the model (Black–Scholes, Heston, or Merton),
- option characteristics,
- numerical and grid parameters.

After modifying the parameters, rebuild and run using the commands above.

---

## Citation

If you use this code in academic work, please cite the associated paper:

BibTeX:

    @article{CosmaGalluccioPederzoliScaillet2020,
      title   = {Early Exercise Decision in American Options with Dividends, Stochastic Volatility, and Jumps},
      author  = {Cosma, Antonio and Galluccio, Stefano and Pederzoli, Paola and Scaillet, Olivier},
      journal = {Journal of Financial and Quantitative Analysis},
      year    = {2020},
      volume  = {55},
      number  = {1},
      pages   = {331--356},
      doi     = {10.1017/S0022109018001229}
    }

---

## License

This repository is released under the CC0-1.0 license.

---

## Acknowledgements

This code relies on open-source numerical libraries, in particular newmat and FFTW.

---

## Contact

Questions, comments, or issues can be raised via the GitHub issue tracker for this repository.

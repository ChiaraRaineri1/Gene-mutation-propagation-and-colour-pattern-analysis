
# Gene Mutation Propagation and Colour Pattern Analysis

This repository contains the course project for the subject **Mathematical and Numerical Methods in Engineering** at Politecnico di Milano, developed by Chiara Raineri during the academic year 2022–2023. The project is split into two core scientific and numerical modeling activities:

1. **Simulation of gene mutation propagation in a linear habitat using the Fisher Equation.**
2. **Pattern formation in biological systems using the Turing model for chemical morphogens.**

## 📘 Project Overview

This project aims to explore real-world biological phenomena through Partial Differential Equation (PDE) modeling and numerical simulations.

### 1️⃣ Gene Mutation Propagation

Models the spread of a mutant gene in a 1D linear habitat using the **Fisher Equation**:

- Domain: Ω = (−4, 4)
- Initial localized mutation: u(x,0) = Gaussian profile
- λ: mutation selection intensity
- u(x,t): gene frequency at position x and time t

#### Numerical Methods Used:
- **BEC** (Backward Euler Centered)
- **FEC** (Forward Euler Centered)

#### Key Parameters:
- Time step: Δt = 0.05 or 0.3
- Space step: h = 0.5
- λ tested from 0.1 to 10

#### Results:
- Simulation of gene diffusion and stability
- Impact of λ on mutation propagation
- Time evolution of the gene frequency u(x,t)
- Computation of p(t): fraction of affected individuals in domain Ω₀ = (−2, 2)

---

### 2️⃣ Pattern Formation – Turing Model

Investigates stable, spatially non-homogeneous chemical distributions responsible for coat colour patterns (e.g. in animals).

#### Model:
- Two-species reaction-diffusion system (u,v) over 2D domain
- β = 5, α and γ varied
- Initial values randomly distributed
- Simulates interactions of **morphogens** (chemical substances)

#### Output:
- Generated patterns based on different α and γ values
- Observations on how tuning parameters affects resulting shapes and structures

---

## 📁 Repository Structure

- `project.m` – Main MATLAB file for activity 1
- functions for activity 1: `FEnonlin.m`, `nonlinsolv.m`
- `turing2d.edp` – FreeFEM analysis for activity 2

## 🧪 Tools & Methods

- **MATLAB**: Partial Differential Equation Modeling (activity 1)
- **FREEFEM**: Finite Element Modeling (activity 2)

## 👩‍💻 Author

**Chiara Raineri**  

## 📜 License

This project is intended for educational and academic use.  
For licensing terms, see the [LICENSE](LICENSE) file.

---

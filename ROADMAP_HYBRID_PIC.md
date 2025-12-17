# Roadmap for Hybrid PIC Implementation in FLEKS

This roadmap outlines the steps to implement a Hybrid PIC (Particle-in-Cell) model in FLEKS. The Hybrid PIC model treats ions as kinetic particles and electrons as a massless, charge-neutralizing fluid.

## Phase 1: Basic Convection Term (Current Implementation)
**Goal:** Implement the fundamental Hybrid PIC loop with the convective term of Ohm's Law.

*   **Assumption:** Electrons are massless fluid; Quasi-neutrality ($n_e \approx n_i$).
*   **Ohm's Law:** $\mathbf{E} = -\mathbf{U}_i \times \mathbf{B}$
    *   $\mathbf{U}_i$: Ion bulk velocity (calculated from particle moments).
    *   $\mathbf{B}$: Magnetic field.
*   **Implementation Details:**
    *   Add `useHybridPIC` flag.
    *   Disable standard Maxwell solver (`solveEM = false`).
    *   Implement `update_E_hybrid` to calculate $\mathbf{E}$ at nodes.
    *   Use existing `update_B` (Faraday's Law) to evolve $\mathbf{B}$.

## Phase 2: Hall Term
**Goal:** Include the Hall effect, important for scales smaller than the ion inertial length.

*   **Ohm's Law:** $\mathbf{E} = -\mathbf{U}_i \times \mathbf{B} + \frac{\mathbf{J} \times \mathbf{B}}{n_e q_e}$
    *   $\mathbf{J} = \frac{1}{\mu_0} \nabla \times \mathbf{B}$ (Ampere's Law, neglecting displacement current).
*   **Implementation Steps:**
    *   Calculate current density $\mathbf{J}$ from $\nabla \times \mathbf{B}$ at cell centers/nodes.
    *   Add the Hall term to the Electric field calculation.
    *   *Note:* The Hall term introduces whistler waves, which puts a strict constraint on the time step ($\Delta t \sim (\Delta x)^2$). Sub-cycling or implicit methods may be required for stability.

## Phase 3: Electron Pressure Gradient
**Goal:** Account for electron pressure effects.

*   **Ohm's Law:** $\mathbf{E} = -\mathbf{U}_i \times \mathbf{B} + \dots - \frac{\nabla P_e}{n_e q_e}$
*   **Electron Model:**
    *   Option A: Isothermal ($P_e = n_e k T_e$, constant $T_e$).
    *   Option B: Adiabatic ($P_e \propto n_e^\gamma$).
*   **Implementation Steps:**
    *   Calculate $\nabla P_e$ (gradient of electron pressure).
    *   Add the pressure gradient term to $\mathbf{E}$.

## Phase 4: Resistivity
**Goal:** Model collisional or anomalous resistivity.

*   **Ohm's Law:** $\mathbf{E} = -\mathbf{U}_i \times \mathbf{B} + \dots + \eta \mathbf{J}$
*   **Implementation Steps:**
    *   Define resistivity $\eta$ (constant or spatially dependent).
    *   Add resistive term to $\mathbf{E}$.

## Phase 5: Electron Thermodynamics (Heating)
**Goal:** Evolve electron temperature dynamically.

*   **Equation:** Solve an energy equation for electrons (e.g., fluid entropy or internal energy equation).
    *   $\frac{\partial P_e}{\partial t} + \mathbf{U}_e \cdot \nabla P_e + \gamma P_e \nabla \cdot \mathbf{U}_e = \dots$
*   **Implementation Steps:**
    *   Implement a fluid solver/update step for electron scalar quantities ($T_e$ or $P_e$).

# Interactive Polarization 3D Simulator

Simulate light polarization through optical systems using Jones Calculus directly in your browser. No installation required. Features an interactive 3D visualization and real-time updates.

➡️ **[Live Demo & Introduction](https://physicify-simulab.github.io/Polarization/)** ⬅️
*(Click Play on the intro page to launch)*

---

## Key Features

*   **Interactive 3D Canvas:** Visualize E-field evolution using Three.js; navigate with OrbitControls (rotate/pan/zoom).
*   **Jones Calculus Engine:** Simulates polarization state changes with Jones vectors and matrices.
*   **Browser-Based:** Runs on desktops, tablets, and mobile devices.
*   **Multiple Visualizations:** Includes 3D view, intensity `I(z)` plot (Plotly), and data table.
*   **Supported Elements:** Linear Polarizers, HWP, QWP, General/Arbitrary Waveplates, Mirrors, Faraday Rotators.
*   **Real-time Updates:** Instantly see effects of adding, removing, or modifying elements and the initial state.
*   **Visualization Controls:** Toggle E-field, envelope, labels; play/pause animation; switch 3D theme (light/dark).
*   **Feedback Integrated:** Giscus comments on intro page; GitHub Issues linked.
*   **Open Source (MIT License).**

---

## Screenshot

![Simulator Screenshot](https://raw.githubusercontent.com/physicify-simulab/Polarization/refs/heads/main/Screenshot.png)

---

## Technical Overview

The simulation uses Jones Calculus, representing light with a complex Jones vector `[Ex, Ey]` and optical elements with 2x2 Jones matrices. The output vector is calculated by sequential matrix multiplication. **Jones matrices are often based on standard definitions, primarily sourced from Wikipedia's Jones Calculus page.** Uses phase convention `φ = kz - ωt`.

---

## Quick Start

1.  **Visit the Live Demo** and click the **Play** button.
2.  **Initial State:** Modify the first row in the "Optical Path" table for the input beam.
3.  **Add Elements:** Use the "Add Element" dropdown and button.
4.  **Modify:** Edit element properties (position, angle, retardation) directly in the table. Use "Remove" buttons as needed.
5.  **Observe:** See results instantly updated in the 3D view, intensity plot, and table. Use mouse/touch to explore the 3D view.

---

## Technical Stack

*   HTML5 / CSS3 / Vanilla JavaScript (ES6+)
*   Plotly.js (Plotting)
*   complex.min.js (for complex number arithmetic)
*   Three.js (3D Visualization)
*   Giscus (Comments)

---

## Contributing & Feedback

This simulator is in beta. Contributions and feedback are welcome!

➡️ **[Report Issues or Suggest Features Here](https://github.com/physicify-simulab/Polarization/issues)** ⬅️

---

## License

Distributed under the [MIT License](LICENSE).

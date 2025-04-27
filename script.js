// script.js (Complete - Includes Double-Pass Mirror Logic & Direct Waveplate Formulas)

// --- Global State ---
let opticalElements = [];
let initialJonesVector = [new Complex(1, 0), new Complex(0, 0)]; // Default: H-pol
let initialStateType = 'H'; // Keep track of the selected type/name
let nextElementId = 0;
const PI = Math.PI;
const sqrt2Inv = 1 / Math.sqrt(2);

// Animation State
let animationState = {
    running: false, time: 0, omega: 2 * PI * 0.5, k: 2 * PI / 1,
    requestId: null, lastTimestamp: 0
};

// SVG Namespace and Constants
const SVG_NS = "http://www.w3.org/2000/svg";
const SVG_MAX_RADIUS = 20; // Max radius for drawing in SVG cell
const RAD_TO_DEG = 180 / PI;
const DEG_TO_RAD = PI / 180;
const EPSILON = 1e-9; // Small number for float comparisons

// Reflection Coordinate Transformation Matrix
const REFLECTION_COORD_TRANSFORM = [[complex(1), complex(0)], [complex(0), complex(-1)]];


// --- DOM Elements ---
const addElementSelect = document.getElementById('add-element');
const addElementBtn = document.getElementById('add-element-btn');
const playPauseBtn = document.getElementById('play-pause-btn');
const intensityPlotDiv = document.getElementById('intensityPlot');
const showElementsCheckbox = document.getElementById('show-elements-checkbox');
// --- ADDED Export/Import Elements ---
const exportBtn = document.getElementById('export-btn');
const importBtn = document.getElementById('import-btn');
const importFileInput = document.getElementById('import-file-input');
// --- END Added Elements ---

// --- Explicitly Expose Key State and Functions to Window ---
window.opticalElements = opticalElements;
window.initialJonesVector = initialJonesVector;
window.animationState = animationState;
window.updateCanvasVisualization = null; // Placeholder, will be assigned by interactiveCanvas.js
window.initCanvas = null; // Placeholder
window.syncCanvasElements = null; // Placeholder
window.removeCanvasElement = null; // Placeholder
window.Complex = Complex; // Expose Complex library
window.complex = complex; // Expose helper function
window.multiplyMatrixVector = multiplyMatrixVector; // Expose needed function for canvas
window.REFLECTION_COORD_TRANSFORM = REFLECTION_COORD_TRANSFORM;


// --- Complex Number Helper Function ---
function complex(real, imag = 0) {
    if (typeof Complex === 'undefined') {
        console.error("Complex library not loaded!");
        return { re: NaN, im: NaN, mul: () => complex(NaN), add: () => complex(NaN), sub: () => complex(NaN), div: () => complex(NaN), exp: () => complex(NaN), abs: () => NaN, arg: () => NaN, conjugate: () => complex(NaN)};
    }
    return new Complex(real, imag);
}

// --- Matrix and Vector Operations ---
function multiplyMatrixVector(matrix, vector) {
    if (!Array.isArray(matrix) || matrix.length !== 2 || !Array.isArray(matrix[0]) || matrix[0].length !== 2 || !Array.isArray(matrix[1]) || matrix[1].length !== 2) { console.error("Invalid matrix format for multiplication:", matrix); return [complex(NaN), complex(NaN)]; }
    if (!Array.isArray(vector) || vector.length !== 2 || !(vector[0] instanceof Complex) || !(vector[1] instanceof Complex)) { console.error("Invalid vector format for multiplication:", vector); return [complex(NaN), complex(NaN)]; }
    const [[m11, m12], [m21, m22]] = matrix; const [v1, v2] = vector;
    // Ensure Complex instances
    const M11 = (m11 instanceof Complex) ? m11 : complex(NaN); const M12 = (m12 instanceof Complex) ? m12 : complex(NaN); const M21 = (m21 instanceof Complex) ? m21 : complex(NaN); const M22 = (m22 instanceof Complex) ? m22 : complex(NaN);
    const V1 = (v1 instanceof Complex) ? v1 : complex(NaN); const V2 = (v2 instanceof Complex) ? v2 : complex(NaN);
    if (isNaN(M11.re) || isNaN(M12.re) || isNaN(M21.re) || isNaN(M22.re) || isNaN(V1.re) || isNaN(V2.re)) {
        // console.warn("NaN detected in matrix/vector multiplication inputs"); // Reduce noise
        return [complex(NaN), complex(NaN)];
    }
    const r1 = M11.mul(V1).add(M12.mul(V2)); const r2 = M21.mul(V1).add(M22.mul(V2));
    return [r1, r2];
}

function multiplyMatrices(A, B) {
    if (!A || !B || !A[0] || !B[0] || !A[0][0] || !B[0][0]) { console.error("Invalid matrix for multiplication:", A, B); return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]]; }
    const [[a11, a12], [a21, a22]] = A; const [[b11, b12], [b21, b22]] = B;
    if ([a11,a12,a21,a22,b11,b12,b21,b22].some(c => !(c instanceof Complex))) { console.error("Non-complex element in matrix mult:", A, B); return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]]; }
    return [ [ a11.mul(b11).add(a12.mul(b21)), a11.mul(b12).add(a12.mul(b22)) ], [ a21.mul(b11).add(a22.mul(b21)), a21.mul(b12).add(a22.mul(b22)) ] ];
}

function vecMagSq(vector) {
    if (!Array.isArray(vector) || vector.length < 2 || !(vector[0] instanceof Complex) || !(vector[1] instanceof Complex)) { return NaN; }
    // Check for NaN components before calculation
    if (isNaN(vector[0].re) || isNaN(vector[0].im) || isNaN(vector[1].re) || isNaN(vector[1].im)) { return NaN; }
    const magSq = vector[0].conjugate().mul(vector[0]).add(vector[1].conjugate().mul(vector[1])).re;
    return isNaN(magSq) ? NaN : magSq;
}

function normalizeVector(vector) {
    if (!Array.isArray(vector) || vector.length < 2 || !(vector[0] instanceof Complex) || !(vector[1] instanceof Complex)) { return [complex(NaN), complex(NaN)]; }
     // Check for NaN components before calculation
    if (isNaN(vector[0].re) || isNaN(vector[0].im) || isNaN(vector[1].re) || isNaN(vector[1].im)) { return [complex(NaN), complex(NaN)]; }
    const magSq = vecMagSq(vector);
    if (isNaN(magSq) || magSq < EPSILON) return [complex(0), complex(0)]; // Avoid division by zero/small number
    const mag = Math.sqrt(magSq);
    return [vector[0].div(mag), vector[1].div(mag)];
}

// --- Jones Matrix Calculation (Uses Direct Rotated Waveplate Formulas) ---
function calculateJonesMatrix(element) {
    const physicalAngleRad = (element.angle || 0) * DEG_TO_RAD; // Angle theta
    const etaRad = (element.parameters?.eta || 0) * DEG_TO_RAD; // Retardation eta
    const deltaRad = (element.parameters?.delta || 0) * DEG_TO_RAD; // Ellipticity delta (Arbitrary Birefringent)
    const thetaRotRad = (element.parameters?.theta_rot || 0) * DEG_TO_RAD; // Faraday rotation theta_rot

    const c1 = complex(1); const c0 = complex(0); const ci = complex(0, 1);
    const c = Math.cos(physicalAngleRad); const s = Math.sin(physicalAngleRad);
    const c2 = c*c; const s2 = s*s; const sc = s*c;

    try {
        switch (element.type) {
            case 'Linear Polarizer': {
                return [ [complex(c2), complex(sc)], [complex(sc), complex(s2)] ];
            }
            case 'Mirror':
                 return [[complex(1), c0], [c0, complex(-1)]];

            case 'HWP': // Retardation eta = PI
            case 'QWP': // Retardation eta = PI/2
            case 'General Waveplate': {
                let specificEtaRad;
                if (element.type === 'HWP') specificEtaRad = PI;
                else if (element.type === 'QWP') specificEtaRad = PI / 2;
                else specificEtaRad = etaRad;

                const expEta = complex(0, specificEtaRad).exp(); // e^(i*eta)
                const oneMinusExpEta = c1.sub(expEta);         // (1 - e^(i*eta))

                const m11 = complex(c2).add(expEta.mul(s2));
                const m12 = oneMinusExpEta.mul(sc);
                const m22 = complex(s2).add(expEta.mul(c2));

                return [[m11, m12], [m12, m22]];
            }
            case 'Arbitrary Birefringent': {
                 const expEta = complex(0, etaRad).exp();
                 const expDelta = complex(0, deltaRad).exp();
                 const expNegDelta = complex(0, -deltaRad).exp();
                 const oneMinusExpEta = c1.sub(expEta);
                 const m11_arb = complex(c2).add(expEta.mul(s2));
                 const m12_arb = oneMinusExpEta.mul(expNegDelta).mul(sc);
                 const m21_arb = oneMinusExpEta.mul(expDelta).mul(sc);
                 const m22_arb = complex(s2).add(expEta.mul(c2));
                 return [[m11_arb, m12_arb], [m21_arb, m22_arb]];
             }
            case 'Faraday Rotator': {
                const cosTh = Math.cos(thetaRotRad); const sinTh = Math.sin(thetaRotRad);
                return [ [complex(cosTh), complex(sinTh)], [complex(-sinTh), complex(cosTh)] ];
            }
            default:
                console.warn("Unknown element type for Jones Matrix:", element.type);
                return [[c1, c0], [c0, c1]]; // Identity matrix
        }
    } catch (error) {
        console.error("Error calculating Jones matrix for:", element, error);
        return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]];
    }
}

// --- Calculate Backward Jones Matrix (for Reflected Path) ---
function calculateBackwardJonesMatrix(forwardMatrix, elementType) {
    if (!Array.isArray(forwardMatrix) || forwardMatrix.length !== 2 ||
        !Array.isArray(forwardMatrix[0]) || forwardMatrix[0].length !== 2 ||
        !Array.isArray(forwardMatrix[1]) || forwardMatrix[1].length !== 2 ||
        !(forwardMatrix[0][0] instanceof Complex)) {
        console.error("Invalid forward matrix format for backward calculation:", forwardMatrix);
        return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]];
    }
    if (forwardMatrix.flat().some(c => !(c instanceof Complex) || isNaN(c.re))) {
        console.warn(`NaN detected in forward matrix input for backward calculation (type: ${elementType})`, forwardMatrix);
        return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]];
    }

    const [[m11, m12], [m21, m22]] = forwardMatrix;
    let backwardMatrix;

    try {
        if (elementType === 'Faraday Rotator') {
            backwardMatrix = [ [m11, m21], [m12, m22] ];
        } else {
            const neg_m12 = m12.mul(-1);
            const neg_m21 = m21.mul(-1);
            backwardMatrix = [ [m11, neg_m12], [neg_m21, m22] ];
        }

        if (backwardMatrix.flat().some(c => !(c instanceof Complex) || isNaN(c.re))) {
            console.warn(`NaN generated during backward matrix calculation for type ${elementType}`, forwardMatrix, backwardMatrix);
            return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]];
        }
        return backwardMatrix;

    } catch (error) {
         console.error(`Error during backward matrix calculation for type ${elementType}:`, error, forwardMatrix);
         return [[complex(NaN), complex(NaN)], [complex(NaN), complex(NaN)]];
    }
}
// --- System Recalculation (Handles Forward and Backward Pass) ---
function recalculateSystem() {
    // 1. Sort Elements by Position
    opticalElements.sort((a, b) => a.position - b.position);
    window.opticalElements = opticalElements; // Update global reference

    // 2. Get Initial Jones Vector and Validate
    let currentVector = initialJonesVector;
    if (!Array.isArray(currentVector) || currentVector.length !== 2 ||
        !(currentVector[0] instanceof Complex) || !(currentVector[1] instanceof Complex) ||
        isNaN(currentVector[0].re) || isNaN(currentVector[1].re)) {
        console.error("Invalid initial Jones vector detected during recalculation, resetting to H.", currentVector);
        setInitialPolarization('H');
        currentVector = initialJonesVector;
    }

    // 3. Initialize Reflection State Variables
    let doublePassMirrorFound = false;
    let doublePassMirrorIndex = -1;
    let vectorAtMirrorInput = null;

    // --- 4. Forward Pass Calculation ---
    opticalElements.forEach((el, index) => {
        el.reflectedInputVector = undefined; el.reflectedOutputVector = undefined;
        el.reflectedIntensity = undefined; el.isReflectedStage = false;

        if (doublePassMirrorFound) {
            el.inputVector = [complex(0), complex(0)];
            el.jonesMatrix = calculateJonesMatrix(el);
            el.outputVector = [complex(0), complex(0)];
            el.intensity = 0;
        } else {
            el.inputVector = currentVector;
            el.jonesMatrix = calculateJonesMatrix(el);

            if (el.type === 'Mirror' && el.parameters?.isDoublePass) {
                doublePassMirrorFound = true;
                doublePassMirrorIndex = index;
                vectorAtMirrorInput = el.inputVector;
                el.outputVector = [complex(0), complex(0)];
                el.intensity = 0;
                currentVector = [complex(0), complex(0)];
            } else {
                el.outputVector = multiplyMatrixVector(el.jonesMatrix, currentVector);
                el.intensity = vecMagSq(el.outputVector);

                if (isNaN(el.outputVector[0].re) || isNaN(el.outputVector[1].re)) {
                    console.warn(`NaN detected in forward path output of element ${el.id} (${el.type})`);
                    currentVector = [complex(NaN), complex(NaN)];
                } else {
                    currentVector = el.outputVector;
                }
            }
        }
    });

    // --- 5. Backward Pass Calculation ---
    if (doublePassMirrorFound && vectorAtMirrorInput &&
        !isNaN(vectorAtMirrorInput[0]?.re) && !isNaN(vectorAtMirrorInput[1]?.re)) {

        let currentReflectedVector = multiplyMatrixVector(REFLECTION_COORD_TRANSFORM, vectorAtMirrorInput);

        if (isNaN(currentReflectedVector[0].re) || isNaN(currentReflectedVector[1].re)) {
            console.warn("NaN detected immediately after reflection transform. Skipping backward pass.");
            for (let i = doublePassMirrorIndex - 1; i >= 0; i--) {
                 opticalElements[i].isReflectedStage = true; opticalElements[i].reflectedInputVector = [complex(NaN), complex(NaN)]; opticalElements[i].reflectedOutputVector = [complex(NaN), complex(NaN)]; opticalElements[i].reflectedIntensity = NaN;
            }
        } else {
             for (let i = doublePassMirrorIndex - 1; i >= 0; i--) {
                 const el = opticalElements[i];
                 el.isReflectedStage = true;
                 el.reflectedInputVector = currentReflectedVector;
                 const forwardMatrix = el.jonesMatrix;
                 const reflectedJonesMatrix = calculateBackwardJonesMatrix(forwardMatrix, el.type);
                 el.reflectedJonesMatrix = reflectedJonesMatrix;
                 el.reflectedOutputVector = multiplyMatrixVector(reflectedJonesMatrix, currentReflectedVector);
                 el.reflectedIntensity = vecMagSq(el.reflectedOutputVector);

                 if (isNaN(el.reflectedOutputVector[0].re) || isNaN(el.reflectedOutputVector[1].re)) {
                     console.warn(`NaN detected in reflected path output of element ${el.id} (${el.type})`);
                     currentReflectedVector = [complex(NaN), complex(NaN)];
                     for (let j = i - 1; j >= 0; j--) {
                          opticalElements[j].isReflectedStage = true; opticalElements[j].reflectedInputVector = [complex(NaN), complex(NaN)]; opticalElements[j].reflectedOutputVector = [complex(NaN), complex(NaN)]; opticalElements[j].reflectedIntensity = NaN;
                     }
                     break;
                 } else {
                     currentReflectedVector = el.reflectedOutputVector;
                 }
             }
        }
    } else if (doublePassMirrorFound && ( !vectorAtMirrorInput || isNaN(vectorAtMirrorInput[0]?.re) || isNaN(vectorAtMirrorInput[1]?.re))) {
        console.warn("Double pass mirror active, but input vector to mirror is invalid. Skipping reflection calculation.");
         for (let i = doublePassMirrorIndex - 1; i >= 0; i--) {
              opticalElements[i].isReflectedStage = true; opticalElements[i].reflectedInputVector = [complex(NaN), complex(NaN)]; opticalElements[i].reflectedOutputVector = [complex(NaN), complex(NaN)]; opticalElements[i].reflectedIntensity = NaN;
         }
    }

    // --- 6. Update UI and Visualization ---
    updateTable();
    updateIntensityPlot();
    if (window.updateCanvasVisualization) {
        const visVector = (!isNaN(initialJonesVector[0].re) && !isNaN(initialJonesVector[1].re))
                          ? initialJonesVector
                          : [complex(1), complex(0)];
        window.updateCanvasVisualization(opticalElements, visVector, animationState);
    } else {
         console.warn("window.updateCanvasVisualization is not available or not yet assigned.");
    }
}


// --- Calculate Ellipse Parameters ---
function calculateEllipseParameters(jonesVector) {
    if (!Array.isArray(jonesVector) || jonesVector.length !== 2 || !(jonesVector[0] instanceof Complex) || !(jonesVector[1] instanceof Complex)) {
        return { psiDeg: NaN, chiDeg: NaN, majorAxis: 0, minorAxis: 0, isValid: false };
    }
    const Ex = jonesVector[0];
    const Ey = jonesVector[1];
    if (isNaN(Ex.re) || isNaN(Ex.im) || isNaN(Ey.re) || isNaN(Ey.im)) {
        return { psiDeg: NaN, chiDeg: NaN, majorAxis: 0, minorAxis: 0, isValid: false };
    }

    const E0x = Ex.abs();
    const E0y = Ey.abs();
    const delta = Ey.arg() - Ex.arg(); // Relative phase phi_y - phi_x
    const E0xSq = E0x * E0x;
    const E0ySq = E0y * E0y;
    const intensity = E0xSq + E0ySq;

    if (intensity < EPSILON) {
        return { psiDeg: 0, chiDeg: 0, majorAxis: 0, minorAxis: 0, isValid: true }; // Represent as a point
    }

    let psiRad;
    if (Math.abs(E0xSq - E0ySq) < EPSILON * intensity) {
         if (Math.abs(Math.cos(delta)) < EPSILON) { psiRad = PI / 4; }
         else if (Math.cos(delta) > 0) { psiRad = PI / 4; }
         else { psiRad = 3 * PI / 4; }
    } else {
         psiRad = 0.5 * Math.atan2(2 * E0x * E0y * Math.cos(delta), E0xSq - E0ySq);
    }
    if (psiRad < -EPSILON) { psiRad += PI; }
    else if (psiRad >= PI - EPSILON) { psiRad = 0; }

    let sin2ChiArg = (2 * E0x * E0y * Math.sin(delta)) / intensity;
    sin2ChiArg = Math.max(-1.0, Math.min(1.0, sin2ChiArg));
    const chiRad = 0.5 * Math.asin(sin2ChiArg);

    const majorAxis = SVG_MAX_RADIUS;
    const minorAxis = majorAxis * Math.abs(Math.tan(chiRad));

    return {
        psiDeg: psiRad * RAD_TO_DEG, chiDeg: chiRad * RAD_TO_DEG,
        majorAxis: majorAxis, minorAxis: minorAxis, isValid: true
    };
}


// --- UI Update Functions ---

// Helper to create SVG ellipse visualization
function createEllipseSVG(params) {
    const svg = document.createElementNS(SVG_NS, "svg");
    svg.setAttribute("viewBox", "-25 -25 50 50");

    if (!params || !params.isValid) {
        const text = document.createElementNS(SVG_NS, "text"); text.setAttribute("x", "0"); text.setAttribute("y", "5"); text.setAttribute("font-size", "10"); text.setAttribute("text-anchor", "middle"); text.setAttribute("dominant-baseline", "middle"); text.setAttribute("fill", "red"); text.textContent = "N/A"; svg.appendChild(text);
    } else if (params.majorAxis < EPSILON && params.minorAxis < EPSILON) {
        const circle = document.createElementNS(SVG_NS, "circle"); circle.setAttribute("cx", "0"); circle.setAttribute("cy", "0"); circle.setAttribute("r", "1.5"); circle.setAttribute("fill", "grey"); svg.appendChild(circle);
    } else {
        const ellipse = document.createElementNS(SVG_NS, "ellipse"); ellipse.setAttribute("cx", "0"); ellipse.setAttribute("cy", "0"); ellipse.setAttribute("rx", params.majorAxis.toFixed(2)); ellipse.setAttribute("ry", params.minorAxis.toFixed(2)); ellipse.setAttribute("transform", `rotate(${-params.psiDeg.toFixed(1)})`); ellipse.setAttribute("stroke", "blue"); ellipse.setAttribute("stroke-width", "1.5"); ellipse.setAttribute("fill", "none"); svg.appendChild(ellipse);
        const line = document.createElementNS(SVG_NS, "line"); line.setAttribute("x1", (-params.majorAxis).toFixed(2)); line.setAttribute("y1", "0"); line.setAttribute("x2", params.majorAxis.toFixed(2)); line.setAttribute("y2", "0"); line.setAttribute("stroke", "rgb(0, 191, 255)"); line.setAttribute("stroke-width", "1"); line.setAttribute("transform", `rotate(${-params.psiDeg.toFixed(1)})`); svg.appendChild(line);
    }
    return svg;
}

// Helper to create Psi/Chi text block
function createPsiChiText(params) {
    const textDiv = document.createElement('div');
    textDiv.classList.add('psi-chi-text');
    if (!params || !params.isValid || isNaN(params.psiDeg) || isNaN(params.chiDeg)) { textDiv.innerHTML = `ψ: NaN<br>χ: NaN`; }
    else { textDiv.innerHTML = `ψ: ${params.psiDeg.toFixed(1)}°<br>χ: ${params.chiDeg.toFixed(1)}°`; }
    textDiv.title = "Orientation Angle (ψ)\nEllipticity Angle (χ)";
    return textDiv;
}

// Function to calculate Jones Vector from Psi and Chi (Inverse Calculation)
// psiRad: Orientation angle (0 to PI)
// chiRad: Ellipticity angle (-PI/4 to PI/4)
// Returns a NORMALIZED Jones vector [Ex, Ey]
function calculateJonesVectorFromPsiChi(psiRad, chiRad) {
    // Clamp angles to avoid potential issues outside standard ranges, although math might handle it
    // psiRad = Math.max(0, Math.min(PI, psiRad)); // Allow full range potentially
    // chiRad = Math.max(-PI / 4, Math.min(PI / 4, chiRad));

    const cosPsi = Math.cos(psiRad);
    const sinPsi = Math.sin(psiRad);
    const cosChi = Math.cos(chiRad);
    const sinChi = Math.sin(chiRad);

    // Formulas based on common definition (e.g., Wikipedia Jones Calculus page)
    const Ex_re = cosPsi * cosChi;
    const Ex_im = -sinPsi * sinChi;
    const Ey_re = sinPsi * cosChi;
    const Ey_im = cosPsi * sinChi;

    const Ex = complex(Ex_re, Ex_im);
    const Ey = complex(Ey_re, Ey_im);

    // The formulas inherently produce a normalized vector if psi/chi are correct
    const vec = [Ex, Ey];
    // Optional: Explicitly normalize for safety, though shouldn't be strictly necessary
    // return normalizeVector(vec);
    return vec;
}

// Main Table Update Function (Handles Forward and Reflected Paths)
function updateTable() {
    const tableBody = document.querySelector('#elements-table tbody');
    if (!tableBody) { console.error("Cannot find table body #elements-table tbody"); return; }
    tableBody.innerHTML = ''; // Clear existing rows

    // --- Add Initial Beam Row ---
    const initialRow = tableBody.insertRow();
    initialRow.insertCell().textContent = 'Input Beam';
    initialRow.insertCell().textContent = 'z = 0';
    const initialPropsCell = initialRow.insertCell();
    initialPropsCell.classList.add('initial-properties-cell');
    const initialEllipseCell = initialRow.insertCell(); // Ellipse cell
    initialEllipseCell.classList.add('ellipse-cell', 'input-ellipse-cell'); // Added class for specific styling
    const initialVecCell = initialRow.insertCell();        // Output Vector
    const initialIntCell = initialRow.insertCell();        // Intensity
    initialRow.insertCell().textContent = '—';            // Actions

    // --- Initial Polarization Properties (Dropdown + Custom Jones) ---
    const initialLabel = document.createElement('label'); initialLabel.htmlFor = 'initial-polarization-inline'; initialLabel.textContent = 'Polarization:';
    const initialSelect = document.createElement('select'); initialSelect.id = 'initial-polarization-inline';
    initialSelect.innerHTML = `<option value="H">Horizontal |H⟩ [1, 0]</option><option value="V">Vertical |V⟩ [0, 1]</option><option value="D">Diagonal L+45 |D⟩ [1/√2, 1/√2]</option><option value="A">Anti-diagonal L-45 |A⟩ [1/√2, -1/√2]</option><option value="R">Circular CCW |R⟩ [1/√2, +i/√2]</option><option value="L">Circular CW |L⟩ [1/√2, -i/√2]</option><option value="Custom">Custom...</option>`;
    initialSelect.value = initialStateType;
    initialSelect.addEventListener('change', () => handleInitialPolarizationChange(true));

    const customDiv = document.createElement('div'); customDiv.id = 'custom-initial-polarization-inline'; customDiv.classList.add('custom-inputs'); customDiv.style.display = (initialStateType === 'Custom') ? 'block' : 'none';
    const exReal = initialJonesVector[0]?.re?.toFixed(3) ?? '0'; const exImag = initialJonesVector[0]?.im?.toFixed(3) ?? '0'; const eyReal = initialJonesVector[1]?.re?.toFixed(3) ?? '0'; const eyImag = initialJonesVector[1]?.im?.toFixed(3) ?? '0';
    customDiv.innerHTML = `<span>Ex:</span> <input type="text" id="initial-ex-real-inline" value="${exReal}" size="4"> + <input type="text" id="initial-ex-imag-inline" value="${exImag}" size="4"> i <br><span>Ey:</span> <input type="text" id="initial-ey-real-inline" value="${eyReal}" size="4"> + <input type="text" id="initial-ey-imag-inline" value="${eyImag}" size="4"> i <button id="set-custom-initial-inline">Set</button>`;

    initialPropsCell.appendChild(initialLabel); initialPropsCell.appendChild(initialSelect); initialPropsCell.appendChild(customDiv);
    const customBtn = initialPropsCell.querySelector('#set-custom-initial-inline'); if(customBtn) { customBtn.addEventListener('click', () => handleSetCustomInitial(true)); }

    if (initialStateType === 'Custom') {
        const customOption = initialSelect.querySelector('option[value="Custom"]');
        if (customOption) {
            const normVec = normalizeVector(initialJonesVector); // Use current (potentially unnormalized if just set from psi/chi)
            if (!isNaN(normVec[0]?.re)) {
                customOption.textContent = `Custom [${normVec[0].re.toFixed(2)}${normVec[0].im >= 0 ? '+' : ''}${normVec[0].im.toFixed(2)}i, ${normVec[1].re.toFixed(2)}${normVec[1].im >= 0 ? '+' : ''}${normVec[1].im.toFixed(2)}i]`;
            } else {
                customOption.textContent = `Custom [Invalid]`;
            }
        }
     }

    // --- Initial Ellipse Cell (SVG + EDITABLE Psi/Chi Inputs) ---
    const initialParams = calculateEllipseParameters(initialJonesVector);
    initialEllipseCell.appendChild(createEllipseSVG(initialParams)); // Keep the SVG

    // Create container for Psi/Chi inputs
    const psiChiInputsDiv = document.createElement('div');
    psiChiInputsDiv.classList.add('psi-chi-inputs'); // Add class for styling

    // Psi Input
    const psiLabel = document.createElement('label'); psiLabel.htmlFor = 'initial-psi-input'; psiLabel.textContent = 'ψ:';
    const psiInput = document.createElement('input');
    psiInput.type = 'number'; psiInput.id = 'initial-psi-input'; psiInput.step = 0.1;
    psiInput.value = (!initialParams || isNaN(initialParams.psiDeg)) ? '' : initialParams.psiDeg.toFixed(1);
    psiInput.title = 'Orientation Angle ψ (0° to 180°)';
    psiInput.addEventListener('keydown', handleInitialPsiChiEnter);
    psiInput.addEventListener('blur', (e) => handleInitialPsiChiChange(e, false)); // Handle blur

    // Chi Input
    const chiLabel = document.createElement('label'); chiLabel.htmlFor = 'initial-chi-input'; chiLabel.textContent = 'χ:';
    const chiInput = document.createElement('input');
    chiInput.type = 'number'; chiInput.id = 'initial-chi-input'; chiInput.step = 0.1;
    chiInput.value = (!initialParams || isNaN(initialParams.chiDeg)) ? '' : initialParams.chiDeg.toFixed(1);
    chiInput.title = 'Ellipticity Angle χ (-45° to +45°)';
    chiInput.addEventListener('keydown', handleInitialPsiChiEnter);
    chiInput.addEventListener('blur', (e) => handleInitialPsiChiChange(e, false)); // Handle blur

    // Append inputs and labels
    psiChiInputsDiv.appendChild(psiLabel);
    psiChiInputsDiv.appendChild(psiInput);
    psiChiInputsDiv.appendChild(document.createElement('br')); // Line break
    psiChiInputsDiv.appendChild(chiLabel);
    psiChiInputsDiv.appendChild(chiInput);

    initialEllipseCell.appendChild(psiChiInputsDiv); // Add inputs below SVG

    // --- Initial Vector and Intensity (Display Only) ---
    const iVec = initialJonesVector;
    if (!isNaN(iVec[0]?.re) && !isNaN(iVec[1]?.re)) {
        initialVecCell.textContent = `[${iVec[0].re.toFixed(3)}${iVec[0].im >= 0 ? '+' : ''}${iVec[0].im.toFixed(3)}i, ${iVec[1].re.toFixed(3)}${iVec[1].im >= 0 ? '+' : ''}${iVec[1].im.toFixed(3)}i]`;
        initialIntCell.textContent = vecMagSq(iVec).toFixed(4); // Show intensity of current vector
    } else {
        initialVecCell.textContent = '[NaN, NaN]';
        initialIntCell.textContent = 'NaN';
    }
    initialIntCell.style.textAlign = 'right';

    // --- Add Rows for Optical Elements (Forward Pass Rendering) ---
    // ... (rest of the element rows generation code remains the same) ...
    let doublePassMirrorFound = false;
    let doublePassMirrorIndex = -1;
    opticalElements.forEach((el, index) => {
        const isAfterMirror = doublePassMirrorFound;
        const isActiveMirror = !isAfterMirror && el.type === 'Mirror' && el.parameters?.isDoublePass;
        if (isActiveMirror) { doublePassMirrorFound = true; doublePassMirrorIndex = index; }

        const row = tableBody.insertRow(); row.dataset.id = el.id;
        if (isAfterMirror) { row.classList.add('zero-intensity-row'); }
        if (isActiveMirror) { row.classList.add('active-mirror-row'); }

        row.insertCell().textContent = el.type; // Type
        const posCell = row.insertCell(); const posInput = document.createElement('input'); posInput.type = 'number'; posInput.value = el.position.toFixed(2); posInput.step = 0.1; posInput.min = 0; posInput.dataset.param = 'position'; posCell.appendChild(posInput); posCell.title = `Absolute Position (z = ${el.position.toFixed(2)})`; // Position
        const propertiesCell = row.insertCell(); propertiesCell.classList.add('properties-cell'); populateElementProperties(propertiesCell, el); // Properties

        // Ellipse (Display only for elements)
        const ellipseCell = row.insertCell(); ellipseCell.classList.add('ellipse-cell'); let ellipseParams; let ellipseTitle = "Polarization ellipse after this element (forward path)";
        if (isActiveMirror && el.inputVector) {
            // ... (mirror ellipse logic) ...
            const inputVec = el.inputVector; let initialReflectedVectorForDisplay = [complex(NaN), complex(NaN)]; if (Array.isArray(inputVec) && inputVec.length === 2 && inputVec[0] instanceof Complex && inputVec[1] instanceof Complex && !isNaN(inputVec[0].re) && !isNaN(inputVec[1].re)) { const transformedVec = multiplyMatrixVector(REFLECTION_COORD_TRANSFORM, inputVec); if (!isNaN(transformedVec[0]?.re) && !isNaN(transformedVec[1]?.re)) { initialReflectedVectorForDisplay = transformedVec; } } ellipseParams = calculateEllipseParameters(initialReflectedVectorForDisplay); ellipseTitle = "Polarization ellipse of the *reflected* beam immediately after the mirror (in reflected coords)";
        } else {
            ellipseParams = calculateEllipseParameters(el.outputVector); if (el.type === 'Mirror' && !el.parameters?.isDoublePass) { ellipseTitle = "Polarization ellipse after non-double-pass mirror (forward path)"; }
        }
        ellipseCell.appendChild(createEllipseSVG(ellipseParams));
        // Display text for element rows
        const textDiv = document.createElement('div');
        textDiv.classList.add('psi-chi-text'); // Use existing class for display
        if (!ellipseParams || !ellipseParams.isValid || isNaN(ellipseParams.psiDeg) || isNaN(ellipseParams.chiDeg)) { textDiv.innerHTML = `ψ: NaN<br>χ: NaN`; }
        else { textDiv.innerHTML = `ψ: ${ellipseParams.psiDeg.toFixed(1)}°<br>χ: ${ellipseParams.chiDeg.toFixed(1)}°`; }
        textDiv.title = ellipseTitle; // Use the calculated title
        ellipseCell.appendChild(textDiv);


        // Output Vector (Forward Path Output)
        const outVecCell = row.insertCell(); let outExStr = 'NaN', outEyStr = 'NaN', intensityMag = NaN;
        if (Array.isArray(el.outputVector) && el.outputVector.length === 2 && el.outputVector[0] instanceof Complex && el.outputVector[1] instanceof Complex && !isNaN(el.outputVector[0].re) && !isNaN(el.outputVector[1].re)) { const [outEx, outEy] = el.outputVector; outExStr = `${outEx.re.toFixed(3)}${outEx.im >= 0 ? '+' : ''}${outEx.im.toFixed(3)}i`; outEyStr = `${outEy.re.toFixed(3)}${outEy.im >= 0 ? '+' : ''}${outEy.im.toFixed(3)}i`; intensityMag = Math.sqrt(el.intensity || 0); }
        outVecCell.textContent = `[${outExStr}, ${outEyStr}]`; outVecCell.title = `Forward Path Output Vector. Magnitude: ${isNaN(intensityMag) ? 'NaN' : intensityMag.toFixed(3)}`;

        // Intensity (Forward Path Output)
        const intensityCell = row.insertCell(); intensityCell.textContent = (el.intensity !== undefined && !isNaN(el.intensity)) ? el.intensity.toFixed(4) : 'NaN'; intensityCell.style.textAlign = 'right'; intensityCell.title = "Forward Path Intensity";

        // Actions
        const actionCell = row.insertCell(); const removeBtn = document.createElement('button'); removeBtn.textContent = 'X'; removeBtn.classList.add('remove-btn'); removeBtn.onclick = () => removeElement(el.id); actionCell.appendChild(removeBtn);

        row.querySelectorAll('input[type="number"]').forEach(input => { input.addEventListener('keydown', handleTableInputEnter); input.addEventListener('blur', handleTableInputBlurOrEnter); });
    });

    // --- Add Rows for Reflected Path (if applicable) ---
     if (doublePassMirrorFound) {
        const mirrorElement = opticalElements[doublePassMirrorIndex]; let mirrorZ = NaN; if (mirrorElement && typeof mirrorElement.position === 'number' && !isNaN(mirrorElement.position)) { mirrorZ = mirrorElement.position; } else { console.error("Could not find valid mirror position for relative calculation."); }

        for (let i = doublePassMirrorIndex - 1; i >= 0; i--) {
            const el = opticalElements[i];
            const hasReflectedData = el.reflectedOutputVector !== undefined;

            const row = tableBody.insertRow(); row.dataset.id = `${el.id}-reflected`; row.classList.add('reflected-row');
            row.insertCell().textContent = `${el.type} (Reflected)`; // Type

            // Position (Relative to Mirror)
            const posCell = row.insertCell(); let relativeZ = NaN; let posTitle = 'Cannot calculate relative position.';
            if (!isNaN(mirrorZ) && typeof el.position === 'number' && !isNaN(el.position)) { relativeZ = Math.abs(el.position - mirrorZ); posCell.textContent = `z' = ${relativeZ.toFixed(2)}`; posTitle = `Position relative to mirror (at z=${mirrorZ.toFixed(2)}). Absolute z=${el.position.toFixed(2)}`; }
            else { posCell.textContent = 'z\' = ?'; posTitle += ` Abs z=${el.position?.toFixed(2) ?? 'N/A'}`; } posCell.title = posTitle;

            // Properties (Disabled, show effective angle)
            const propertiesCell = row.insertCell(); propertiesCell.classList.add('properties-cell', 'reflected-properties'); populateElementProperties(propertiesCell, el); propertiesCell.querySelectorAll('input, button, select').forEach(input => input.disabled = true);
            const angleInput = propertiesCell.querySelector('input[data-param="angle"]');
            if (angleInput && ['Linear Polarizer', 'HWP', 'QWP', 'General Waveplate', 'Arbitrary Birefringent'].includes(el.type)) { const originalAngle = el.angle || 0; const effectiveAngle = -originalAngle; const infoSpan = document.createElement('span'); infoSpan.className = 'reflected-angle-info'; infoSpan.textContent = ` (Eff: ${effectiveAngle.toFixed(1)}°)`; angleInput.parentNode.appendChild(infoSpan); }

            // Ellipse (Reflected - Display Only)
            const ellipseCell = row.insertCell(); ellipseCell.classList.add('ellipse-cell');
            const reflectedParams = calculateEllipseParameters(el.reflectedOutputVector);
            ellipseCell.appendChild(createEllipseSVG(reflectedParams));
            // Display text for reflected rows
            const reflectedTextDiv = document.createElement('div');
            reflectedTextDiv.classList.add('psi-chi-text');
            if (!reflectedParams || !reflectedParams.isValid || isNaN(reflectedParams.psiDeg) || isNaN(reflectedParams.chiDeg)) { reflectedTextDiv.innerHTML = `ψ: NaN<br>χ: NaN`; }
            else { reflectedTextDiv.innerHTML = `ψ: ${reflectedParams.psiDeg.toFixed(1)}°<br>χ: ${reflectedParams.chiDeg.toFixed(1)}°`; }
            reflectedTextDiv.title = "Polarization ellipse after this element (reflected path)";
            ellipseCell.appendChild(reflectedTextDiv);

            // Output Vector (Reflected)
            const outVecCell = row.insertCell(); let refOutExStr = 'NaN', refOutEyStr = 'NaN', refIntensityMag = NaN;
            if (hasReflectedData && Array.isArray(el.reflectedOutputVector) && el.reflectedOutputVector.length === 2 && el.reflectedOutputVector[0] instanceof Complex && el.reflectedOutputVector[1] instanceof Complex && !isNaN(el.reflectedOutputVector[0].re) && !isNaN(el.reflectedOutputVector[1].re)) { const [refOutEx, refOutEy] = el.reflectedOutputVector; refOutExStr = `${refOutEx.re.toFixed(3)}${refOutEx.im >= 0 ? '+' : ''}${refOutEx.im.toFixed(3)}i`; refOutEyStr = `${refOutEy.re.toFixed(3)}${refOutEy.im >= 0 ? '+' : ''}${refOutEy.im.toFixed(3)}i`; refIntensityMag = Math.sqrt(el.reflectedIntensity || 0); }
            outVecCell.textContent = `[${refOutExStr}, ${refOutEyStr}]`; outVecCell.title = `Reflected Path Output Vector. Magnitude: ${isNaN(refIntensityMag) ? 'NaN' : refIntensityMag.toFixed(3)}`;

            // Intensity (Reflected)
            const intensityCell = row.insertCell(); intensityCell.textContent = (el.reflectedIntensity !== undefined && !isNaN(el.reflectedIntensity)) ? el.reflectedIntensity.toFixed(4) : 'NaN'; intensityCell.style.textAlign = 'right'; intensityCell.title = "Reflected Path Intensity";

            row.insertCell().textContent = '—'; // Actions
        }
    }
    // --- Sync Canvas ---
    if (window.syncCanvasElements) { window.syncCanvasElements(opticalElements); }
    else { console.warn("window.syncCanvasElements not ready for updateTable"); }
}

// Helper to populate element properties cell (Includes Mirror Checkbox)
function populateElementProperties(cell, element) {
    cell.innerHTML = ''; // Clear existing content

    const createPropertyInput = (label, paramName, value, unit = '(deg)', disabled = false, step = 1) => {
        const div = document.createElement('div'); div.classList.add('property-item');
        const labelSpan = document.createElement('span'); labelSpan.textContent = `${label}: `;
        const input = document.createElement('input'); input.type = 'number'; input.value = value?.toFixed(1) ?? '0.0'; input.step = step; input.dataset.param = paramName; input.disabled = disabled; if (disabled) input.style.backgroundColor = '#eee';
        const unitSpan = document.createElement('span'); unitSpan.textContent = ` ${unit}`; unitSpan.classList.add('unit-label');
        div.appendChild(labelSpan); div.appendChild(input); div.appendChild(unitSpan);
        return div;
    };

    if (['Linear Polarizer', 'HWP', 'QWP', 'General Waveplate', 'Arbitrary Birefringent'].includes(element.type)) { cell.appendChild(createPropertyInput('Angle θ', 'angle', element.angle ?? 0, '(deg)', false, 1)); }

    switch (element.type) {
        case 'Linear Polarizer': break;
        case 'Mirror':
            const mirrorDiv = document.createElement('div'); mirrorDiv.classList.add('property-item'); const checkbox = document.createElement('input'); checkbox.type = 'checkbox'; checkbox.checked = element.parameters?.isDoublePass ?? false; checkbox.id = `double-pass-checkbox-${element.id}`; checkbox.dataset.param = 'isDoublePass'; checkbox.addEventListener('change', handleTableCheckboxChange); const label = document.createElement('label'); label.htmlFor = `double-pass-checkbox-${element.id}`; label.textContent = ' 0° (Double Pass)'; label.style.cursor = 'pointer'; label.style.fontWeight = 'normal'; mirrorDiv.appendChild(checkbox); mirrorDiv.appendChild(label); cell.appendChild(mirrorDiv); break;
        case 'HWP': cell.appendChild(createPropertyInput('Retard η', 'eta', element.parameters?.eta ?? 180, '(deg)', true)); break;
        case 'QWP': cell.appendChild(createPropertyInput('Retard η', 'eta', element.parameters?.eta ?? 90, '(deg)', true)); break;
        case 'General Waveplate': cell.appendChild(createPropertyInput('Retard η', 'eta', element.parameters?.eta ?? 0, '(deg)', false, 1)); break;
        case 'Arbitrary Birefringent': cell.appendChild(createPropertyInput('Retard η', 'eta', element.parameters?.eta ?? 0, '(deg)', false, 1)); cell.appendChild(createPropertyInput('Ellipt δ', 'delta', element.parameters?.delta ?? 0, '(deg)', false, 1)); break;
        case 'Faraday Rotator': cell.appendChild(createPropertyInput('Angle θ', 'angle', element.angle ?? 0, '(deg)', true)); cell.appendChild(createPropertyInput('Rot θ_rot', 'theta_rot', element.parameters?.theta_rot ?? -45, '(deg)', false, 1)); break;
        default: cell.textContent = 'N/A';
    }
}


// --- updateIntensityPlot (Handles Mirror Cutoff) ---
function updateIntensityPlot() {
    let zCoords = [0]; let initialIntensity = vecMagSq(initialJonesVector); let intensities = [isNaN(initialIntensity) ? 0 : initialIntensity]; let doublePassMirrorFound = false; let mirrorZ = Infinity; let maxIntensity = isNaN(initialIntensity) ? 0 : initialIntensity;
    const mirrorElement = opticalElements.find(el => el.type === 'Mirror' && el.parameters?.isDoublePass); if (mirrorElement) { doublePassMirrorFound = true; mirrorZ = mirrorElement.position; }
    opticalElements.forEach(el => {
        if (!doublePassMirrorFound || el.position <= mirrorZ) {
            const inputIntensity = vecMagSq(el.inputVector); const lastValidIntensity = intensities.length > 0 ? intensities[intensities.length - 1] : 0; const currentInputIntensity = isNaN(inputIntensity) ? lastValidIntensity : inputIntensity; intensities.push(currentInputIntensity); zCoords.push(el.position - EPSILON); maxIntensity = Math.max(maxIntensity, currentInputIntensity);
            const outputIntensity = el.intensity; const currentOutputIntensity = isNaN(outputIntensity) ? intensities[intensities.length - 1] : outputIntensity; intensities.push(currentOutputIntensity); zCoords.push(el.position); maxIntensity = Math.max(maxIntensity, currentOutputIntensity);
        }
    });
    let lastZ; if (doublePassMirrorFound) { lastZ = mirrorZ; if (zCoords[zCoords.length - 1] < mirrorZ - EPSILON) { zCoords.push(mirrorZ); intensities.push(intensities[intensities.length - 1]); } } else { lastZ = opticalElements.length > 0 ? opticalElements[opticalElements.length - 1].position : 0; const lastIntensity = intensities.length > 0 ? intensities[intensities.length - 1] : 0; zCoords.push(lastZ + Math.max(0.1, Math.abs(lastZ * 0.1))); intensities.push(lastIntensity); }
    const validData = zCoords.map((z, i) => ({ x: z, y: intensities[i] })).filter(p => !isNaN(p.x) && !isNaN(p.y) && isFinite(p.x) && isFinite(p.y));
    if (validData.length < 2) { Plotly.react(intensityPlotDiv, [], { title: 'Intensity vs. Position (z) - No Data' }, {responsive: true}); return; }
    const trace = { x: validData.map(p => p.x), y: validData.map(p => p.y), mode: 'lines', type: 'scatter', line: { shape: 'hv' }, name: doublePassMirrorFound ? 'Intensity (Forward Path)' : 'Intensity' };
    const layout = { title: doublePassMirrorFound ? 'Intensity vs. Position (z) - Forward Path (Reflection Active)' : 'Intensity vs. Position (z)', xaxis: { title: 'Position (z)', autorange: true }, yaxis: { title: 'Intensity', range: [0, Math.max(1.1, maxIntensity * 1.1)], autorange: false }, margin: { l: 50, r: 20, t: 40, b: 40 }, hovermode: 'closest', shapes: [], showlegend: doublePassMirrorFound };
    const showElements = showElementsCheckbox.checked; if (showElements) { opticalElements.forEach(el => { if ((!doublePassMirrorFound || el.position <= mirrorZ) && !isNaN(el.position) && isFinite(el.position)) { layout.shapes.push({ type: 'line', x0: el.position, y0: 0, x1: el.position, y1: 1, yref: 'paper', line: { color: 'rgba(255, 0, 0, 0.7)', width: 1.5, dash: 'dot' } }); } }); }
    Plotly.react(intensityPlotDiv, [trace], layout, {responsive: true});
}

// --- Event Handlers ---

// Add Element Button
function handleAddElement() {
    const type = addElementSelect.value; let lastElementPos = 0; if (opticalElements.length > 0) { const validPositions = opticalElements.map(el => el.position).filter(p => !isNaN(p) && isFinite(p)); if (validPositions.length > 0) { lastElementPos = Math.max(0, ...validPositions); } } const newPosition = lastElementPos + 1.0;
    const newElement = { id: nextElementId++, type: type, position: newPosition, angle: 0, parameters: {}, inputVector: null, outputVector: null, jonesMatrix: null, intensity: null, reflectedInputVector: null, reflectedOutputVector: null, reflectedIntensity: null, isReflectedStage: false };
    if (type === 'QWP') newElement.parameters.eta = 90; if (type === 'HWP') newElement.parameters.eta = 180; if (type === 'General Waveplate') newElement.parameters.eta = 0; if (type === 'Arbitrary Birefringent') { newElement.parameters.eta = 0; newElement.parameters.delta = 0; } if (type === 'Faraday Rotator') newElement.parameters.theta_rot = -45; if (type === 'Mirror') newElement.parameters.isDoublePass = false;
    opticalElements.push(newElement); recalculateSystem();
}

// Remove Element Button
function removeElement(id) {
    const elementToRemove = opticalElements.find(el => el.id === id); if (!elementToRemove) return;
    const wasDoublePassMirror = elementToRemove.type === 'Mirror' && elementToRemove.parameters?.isDoublePass;
    opticalElements = opticalElements.filter(el => el.id !== id); window.opticalElements = opticalElements;
    if(window.removeCanvasElement) { window.removeCanvasElement(id); } else { console.warn("window.removeCanvasElement not available"); }
    recalculateSystem();
}

// Initial Polarization Dropdown Change
function handleInitialPolarizationChange(isInline = false) {
    const selectElement = document.getElementById('initial-polarization-inline'); if (!selectElement) return;
    const selection = selectElement.value;
    initialStateType = selection; // Update state type first

    const customDiv = document.getElementById('custom-initial-polarization-inline');
    const psiChiInputsDiv = document.querySelector('.input-ellipse-cell .psi-chi-inputs'); // Target the inputs

    if (selection === 'Custom') {
        if(customDiv) customDiv.style.display = 'block';
        if(psiChiInputsDiv) psiChiInputsDiv.style.display = 'block'; // Show Psi/Chi inputs for custom
        // Don't immediately recalculate if "Custom..." is selected, wait for user input
        // Update the table to reflect the current custom state
        updateTable(); // Update table to show current custom values
    } else {
        if(customDiv) customDiv.style.display = 'none';
        if(psiChiInputsDiv) psiChiInputsDiv.style.display = 'block'; // Keep Psi/Chi visible for standard states too
        setInitialPolarization(selection); // This sets initialJonesVector and initialStateType
        recalculateSystem(); // Recalculate based on the new standard state
    }
}

// Set Custom Initial Polarization from Jones Vector Inputs
function handleSetCustomInitial(isInline = false) {
    const suffix = '-inline';
    try {
        const ex_re_str = document.getElementById(`initial-ex-real${suffix}`).value;
        const ex_im_str = document.getElementById(`initial-ex-imag${suffix}`).value;
        const ey_re_str = document.getElementById(`initial-ey-real${suffix}`).value;
        const ey_im_str = document.getElementById(`initial-ey-imag${suffix}`).value;

        const ex_re = parseFloat(ex_re_str);
        const ex_im = parseFloat(ex_im_str);
        const ey_re = parseFloat(ey_re_str);
        const ey_im = parseFloat(ey_im_str);

        if ([ex_re, ex_im, ey_re, ey_im].some(isNaN)) {
            throw new Error("Invalid number input. Please enter valid numbers for all parts.");
        }

        const vec = [complex(ex_re, ex_im), complex(ey_re, ey_im)];
        const normVec = normalizeVector(vec); // Always normalize
        initialJonesVector = normVec;
        window.initialJonesVector = initialJonesVector;
        initialStateType = 'Custom'; // Set state to custom
        recalculateSystem(); // Update everything

    } catch (error) {
        alert("Invalid custom polarization input: " + error.message);
        // Optionally revert input fields in table if error occurs
         updateTable();
    }
}

// --- Handle changes to Psi/Chi Inputs ---
function handleInitialPsiChiChange(event, isEnterKey = false) {
    const psiInput = document.getElementById('initial-psi-input');
    const chiInput = document.getElementById('initial-chi-input');
    if (!psiInput || !chiInput) return;

    const psiDegStr = psiInput.value;
    const chiDegStr = chiInput.value;

    const psiDeg = parseFloat(psiDegStr);
    const chiDeg = parseFloat(chiDegStr);

    // --- Basic Validation ---
    if (isNaN(psiDeg) || isNaN(chiDeg)) {
        // If called by Enter key, alert. If by blur, maybe silently revert?
        // For simplicity, let's revert and log warning for now.
        console.warn("Invalid Psi/Chi input (non-numeric). Reverting.");
        // Revert inputs to reflect the *current* state before the invalid change
        const currentParams = calculateEllipseParameters(initialJonesVector);
        psiInput.value = isNaN(currentParams.psiDeg) ? '' : currentParams.psiDeg.toFixed(1);
        chiInput.value = isNaN(currentParams.chiDeg) ? '' : currentParams.chiDeg.toFixed(1);
        if (isEnterKey) {
            alert("Invalid input: Psi and Chi must be numbers.");
        }
        return;
    }

    // --- Optional: Angle Range Clamping/Warning (more user-friendly) ---
    let clampedPsi = psiDeg;
    let clampedChi = chiDeg;
    // Psi is typically 0-180
    if (psiDeg < 0 || psiDeg > 180) {
        // console.warn(`Psi value ${psiDeg}° outside typical range [0, 180]. Using value as entered.`);
        // Or clamp: clampedPsi = Math.max(0, Math.min(180, psiDeg));
    }
     // Chi is typically -45 to 45
    if (chiDeg < -45 || chiDeg > 45) {
       // console.warn(`Chi value ${chiDeg}° outside typical range [-45, 45]. Using value as entered.`);
       // Or clamp: clampedChi = Math.max(-45, Math.min(45, chiDeg));
    }
    // Update input fields if clamped (optional)
    // psiInput.value = clampedPsi.toFixed(1);
    // chiInput.value = clampedChi.toFixed(1);


    // --- Check if values actually changed ---
    const currentParams = calculateEllipseParameters(initialJonesVector);
    const psiDiff = Math.abs(currentParams.psiDeg - clampedPsi);
    const chiDiff = Math.abs(currentParams.chiDeg - clampedChi);

    // Use a tolerance for floating point comparison
    if (psiDiff < 0.05 && chiDiff < 0.05) {
        // console.log("Psi/Chi values didn't change significantly.");
        // Ensure inputs show the correctly formatted current value
        psiInput.value = isNaN(currentParams.psiDeg) ? '' : currentParams.psiDeg.toFixed(1);
        chiInput.value = isNaN(currentParams.chiDeg) ? '' : currentParams.chiDeg.toFixed(1);
        if (isEnterKey) event.target.blur(); // Defocus on Enter if no change
        return; // No change needed
    }

    // --- Values are valid and changed: Calculate new Jones Vector ---
    console.log(`Updating initial state from Psi=${clampedPsi.toFixed(1)}°, Chi=${clampedChi.toFixed(1)}°`);
    const psiRad = clampedPsi * DEG_TO_RAD;
    const chiRad = clampedChi * DEG_TO_RAD;

    const newVector = calculateJonesVectorFromPsiChi(psiRad, chiRad);

    // Update Global State
    initialJonesVector = newVector; // Already normalized from calculation
    window.initialJonesVector = initialJonesVector;
    initialStateType = 'Custom'; // Force state to Custom

    // Recalculate and Update UI
    recalculateSystem();

    // Defocus if change was triggered by Enter key
    if (isEnterKey) {
        event.target.blur();
    }
}

// Handle Enter Key specifically for Psi/Chi inputs
function handleInitialPsiChiEnter(event) {
    if (event.key === 'Enter') {
        event.preventDefault(); // Prevent potential form submission/other default actions
        handleInitialPsiChiChange(event, true); // Pass flag indicating Enter key
    }
}

// Helper to set Jones vector from predefined types
function setInitialPolarization(type) {
    let vec; switch(type) { case 'H': vec = [complex(1), complex(0)]; break; case 'V': vec = [complex(0), complex(1)]; break; case 'D': vec = [complex(sqrt2Inv), complex(sqrt2Inv)]; break; case 'A': vec = [complex(sqrt2Inv), complex(-sqrt2Inv)]; break; case 'R': vec = [complex(sqrt2Inv), complex(0, sqrt2Inv)]; break; case 'L': vec = [complex(sqrt2Inv), complex(0, -sqrt2Inv)]; break; default: console.warn("Unknown polarization type:", type, "defaulting to H."); vec = [complex(1), complex(0)]; type = 'H'; }
    initialJonesVector = normalizeVector(vec); // Ensure normalized
    window.initialJonesVector = initialJonesVector;
    initialStateType = type; // Set the state type name
}


// --- START: Export/Import Handlers ---

// Export current setup to XML
function handleExport() {
    try {
        let xmlString = '<?xml version="1.0" encoding="UTF-8"?>\n';
        xmlString += '<polarizationSetup>\n';

        // Initial State
        xmlString += '  <initialState type="' + initialStateType + '"';
        if (initialStateType === 'Custom') {
            xmlString += ' ex_re="' + (initialJonesVector[0]?.re ?? 0).toFixed(6) + '"';
            xmlString += ' ex_im="' + (initialJonesVector[0]?.im ?? 0).toFixed(6) + '"';
            xmlString += ' ey_re="' + (initialJonesVector[1]?.re ?? 0).toFixed(6) + '"';
            xmlString += ' ey_im="' + (initialJonesVector[1]?.im ?? 0).toFixed(6) + '"';
        }
        xmlString += ' />\n';

        // Elements
        xmlString += '  <elements>\n';
        opticalElements.forEach(el => {
            xmlString += `    <element type="${el.type}" position="${el.position.toFixed(4)}"`;
            if (typeof el.angle === 'number') {
                xmlString += ` angle="${el.angle.toFixed(4)}"`;
            }
            if (el.parameters) {
                if (typeof el.parameters.eta === 'number') {
                    xmlString += ` eta="${el.parameters.eta.toFixed(4)}"`;
                }
                if (typeof el.parameters.delta === 'number') {
                    xmlString += ` delta="${el.parameters.delta.toFixed(4)}"`;
                }
                if (typeof el.parameters.theta_rot === 'number') {
                    xmlString += ` theta_rot="${el.parameters.theta_rot.toFixed(4)}"`;
                }
                if (typeof el.parameters.isDoublePass === 'boolean') {
                    xmlString += ` isDoublePass="${el.parameters.isDoublePass}"`; // Save as "true" or "false"
                }
            }
            xmlString += ' />\n';
        });
        xmlString += '  </elements>\n';
        xmlString += '</polarizationSetup>\n';

        // Create Blob and Download Link
        const blob = new Blob([xmlString], { type: 'application/xml' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'polarization_setup.xml';
        document.body.appendChild(a); // Required for Firefox
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
        console.log("Setup exported successfully.");

    } catch (error) {
        console.error("Error during export:", error);
        alert("Failed to export setup. See console for details.");
    }
}

// Trigger hidden file input for import
function handleImport() {
    if (importFileInput) {
        importFileInput.click(); // Open file dialog
    } else {
        console.error("Import file input element not found.");
    }
}

// Process the selected file
function handleFileSelect(event) {
    const file = event.target.files[0];
    if (!file) {
        console.log("No file selected for import.");
        return;
    }

    // Optional: Check file type (basic check)
    if (!file.name.toLowerCase().endsWith('.xml') && file.type !== 'application/xml') {
         console.warn(`Selected file "${file.name}" might not be an XML file. Proceeding anyway.`);
         // alert(`Warning: Selected file "${file.name}" might not be an XML file.`);
    }

    const reader = new FileReader();

    reader.onload = (e) => {
        const xmlText = e.target.result;
        try {
            const parser = new DOMParser();
            const xmlDoc = parser.parseFromString(xmlText, "application/xml");

            // Check for parsing errors
            const parserError = xmlDoc.querySelector("parsererror");
            if (parserError) {
                console.error("XML Parsing Error:", parserError.textContent);
                throw new Error("Failed to parse XML file. Content may be invalid.");
            }

            // --- Process Initial State ---
            const initialStateNode = xmlDoc.querySelector("polarizationSetup > initialState");
            if (!initialStateNode) {
                throw new Error("Invalid format: <initialState> tag not found.");
            }
            const stateTypeAttr = initialStateNode.getAttribute("type");
            if (!stateTypeAttr) {
                throw new Error("Invalid format: 'type' attribute missing in <initialState>.");
            }

            if (stateTypeAttr === 'Custom') {
                const ex_re = parseFloat(initialStateNode.getAttribute("ex_re") || '0');
                const ex_im = parseFloat(initialStateNode.getAttribute("ex_im") || '0');
                const ey_re = parseFloat(initialStateNode.getAttribute("ey_re") || '0');
                const ey_im = parseFloat(initialStateNode.getAttribute("ey_im") || '0');
                 if ([ex_re, ex_im, ey_re, ey_im].some(isNaN)) {
                     throw new Error("Invalid custom state values in XML.");
                 }
                 const vec = [complex(ex_re, ex_im), complex(ey_re, ey_im)];
                 initialJonesVector = normalizeVector(vec); // Ensure normalized
                 window.initialJonesVector = initialJonesVector;
                 initialStateType = 'Custom';
            } else {
                 // Use existing function for standard types
                 setInitialPolarization(stateTypeAttr);
                 // Ensure initialStateType reflects the loaded type
                 initialStateType = stateTypeAttr;
            }
            console.log("Initial state loaded:", initialStateType);


            // --- Process Elements ---
            const elementsNode = xmlDoc.querySelector("polarizationSetup > elements");
            if (!elementsNode) {
                throw new Error("Invalid format: <elements> tag not found.");
            }

            // Clear existing elements
            opticalElements.length = 0; // Clear array in place
            window.opticalElements = opticalElements; // Update global reference
            nextElementId = 0; // Reset ID counter

            const elementNodes = elementsNode.querySelectorAll("element");
            let activeDoublePassMirrorFound = false; // Track if a double-pass mirror is loaded

            elementNodes.forEach(node => {
                const type = node.getAttribute("type");
                const positionStr = node.getAttribute("position");
                const angleStr = node.getAttribute("angle"); // Might be null

                if (!type || !positionStr) {
                    console.warn("Skipping element due to missing type or position attribute:", node.outerHTML);
                    return; // Skip this element
                }

                const position = parseFloat(positionStr);
                if (isNaN(position)) {
                     console.warn("Skipping element due to invalid position attribute:", node.outerHTML);
                     return;
                }

                const angle = (angleStr !== null) ? parseFloat(angleStr) : 0; // Default angle to 0 if missing
                 if (isNaN(angle)) {
                     console.warn(`Invalid angle "${angleStr}" for element, defaulting to 0.`, node.outerHTML);
                     angle = 0;
                 }

                const newElement = {
                    id: nextElementId++,
                    type: type,
                    position: position,
                    angle: angle,
                    parameters: {},
                    // Reset calculated fields
                    inputVector: null, outputVector: null, jonesMatrix: null, intensity: null,
                    reflectedInputVector: null, reflectedOutputVector: null, reflectedIntensity: null,
                    isReflectedStage: false, reflectedJonesMatrix: null
                };

                // Load parameters based on type
                const etaStr = node.getAttribute("eta");
                const deltaStr = node.getAttribute("delta");
                const thetaRotStr = node.getAttribute("theta_rot");
                const isDoublePassStr = node.getAttribute("isDoublePass");

                if (etaStr !== null && !['Linear Polarizer', 'Mirror', 'Faraday Rotator'].includes(type)) {
                    newElement.parameters.eta = parseFloat(etaStr);
                    if (isNaN(newElement.parameters.eta)) {
                        console.warn(`Invalid eta "${etaStr}" for ${type}, parameter ignored.`, node.outerHTML);
                         delete newElement.parameters.eta;
                    }
                }
                if (deltaStr !== null && type === 'Arbitrary Birefringent') {
                    newElement.parameters.delta = parseFloat(deltaStr);
                     if (isNaN(newElement.parameters.delta)) {
                         console.warn(`Invalid delta "${deltaStr}" for ${type}, parameter ignored.`, node.outerHTML);
                         delete newElement.parameters.delta;
                     }
                }
                if (thetaRotStr !== null && type === 'Faraday Rotator') {
                    newElement.parameters.theta_rot = parseFloat(thetaRotStr);
                     if (isNaN(newElement.parameters.theta_rot)) {
                         console.warn(`Invalid theta_rot "${thetaRotStr}" for ${type}, parameter ignored.`, node.outerHTML);
                         delete newElement.parameters.theta_rot;
                     }
                }
                if (isDoublePassStr !== null && type === 'Mirror') {
                    const isDoublePass = isDoublePassStr.toLowerCase() === 'true';
                    if (isDoublePass && activeDoublePassMirrorFound) {
                        console.warn(`Multiple double-pass mirrors found in imported file. Only the first one encountered will be active. Disabling double-pass for mirror at z=${position}.`);
                        newElement.parameters.isDoublePass = false;
                    } else {
                        newElement.parameters.isDoublePass = isDoublePass;
                        if (isDoublePass) {
                            activeDoublePassMirrorFound = true;
                        }
                    }
                }

                // Add default parameters if needed after loading
                if (type === 'QWP' && typeof newElement.parameters.eta !== 'number') newElement.parameters.eta = 90;
                if (type === 'HWP' && typeof newElement.parameters.eta !== 'number') newElement.parameters.eta = 180;
                // Add others if necessary...

                opticalElements.push(newElement);
            });

            console.log(`Successfully imported ${opticalElements.length} elements.`);

            // Recalculate the system with the new setup
            recalculateSystem();

            // Reset file input value so the change event fires again if the same file is selected
            event.target.value = null;

        } catch (error) {
            console.error("Error during import processing:", error);
            alert("Failed to import setup: " + error.message);
             // Optionally reset to a default state if import fails badly
             // opticalElements.length = 0;
             // setInitialPolarization('H');
             // recalculateSystem();
        }
    };

    reader.onerror = (e) => {
        console.error("Error reading file:", e.target.error);
        alert("Error reading the selected file.");
        event.target.value = null; // Reset input
    };

    reader.readAsText(file); // Start reading the file
}

// --- END: Export/Import Handlers ---


// --- Table Input Handlers ---

// Handle Blur or Enter Key on Number Inputs in Table
function handleTableInputBlurOrEnter(event, isEnterKey = false) {
    const input = event.target; if (input.disabled || input.type !== 'number') return;
    const row = input.closest('tr'); if (!row || !row.dataset.id || row.classList.contains('reflected-row')) return;
    const id = parseInt(row.dataset.id); const param = input.dataset.param; let valueStr = input.value; let value = parseFloat(valueStr);
    const element = opticalElements.find(el => el.id === id); if (!element) { console.error(`Element with ID ${id} not found for input update.`); return; }
    let currentValue; if (param === 'position' || param === 'angle') { currentValue = element[param] ?? 0; } else if (['eta', 'delta', 'theta_rot'].includes(param)) { currentValue = element.parameters?.[param] ?? 0; } else { console.warn(`Unknown parameter '${param}' for input update.`); return; }
    if (isNaN(value) || !isFinite(value)) { input.value = currentValue.toFixed(param === 'position' ? 2 : 1); console.warn(`Invalid input '${valueStr}' for param '${param}', reverting.`); return; }
    if (param === 'position' && value < 0) { value = 0; input.value = value.toFixed(2); }
    const tolerance = 1e-5; if (Math.abs(value - currentValue) < tolerance) { input.value = currentValue.toFixed(param === 'position' ? 2 : 1); if (isEnterKey) input.blur(); return; }
    console.log(`Updating element ${id}, param ${param} from ${currentValue.toFixed(3)} to ${value.toFixed(3)}`);
    if (param === 'position' || param === 'angle') { element[param] = value; } else if (['eta', 'delta', 'theta_rot'].includes(param)) { element.parameters = element.parameters || {}; element.parameters[param] = value; }
    input.value = value.toFixed(param === 'position' ? 2 : 1); recalculateSystem(); if (isEnterKey) { input.blur(); }
}

// Handle Enter Key specifically for number inputs
function handleTableInputEnter(event) { if (event.key === 'Enter' && event.target.type === 'number') { event.preventDefault(); handleTableInputBlurOrEnter(event, true); } }

// Handle Mirror Checkbox Change
function handleTableCheckboxChange(event) {
    const checkbox = event.target; if (checkbox.disabled || checkbox.type !== 'checkbox') return;
    const row = checkbox.closest('tr'); if (!row || !row.dataset.id || row.classList.contains('reflected-row')) return;
    const id = parseInt(row.dataset.id); const param = checkbox.dataset.param; const value = checkbox.checked;
    const element = opticalElements.find(el => el.id === id); if (!element || element.type !== 'Mirror') { console.warn(`Checkbox change ignored: Element ${id} not found or not a Mirror.`); return; }
    console.log(`Updating element ${id}, param ${param} to ${value}`);
    if (value === true) { opticalElements.forEach(el => { if (el.type === 'Mirror' && el.id !== id && el.parameters?.isDoublePass) { console.log(`Disabling double-pass for mirror ${el.id} because ${id} was enabled.`); el.parameters.isDoublePass = false; } }); }
    element.parameters = element.parameters || {}; element.parameters[param] = value; recalculateSystem();
}


// --- Animation Control ---
function toggleAnimation() { if (animationState.running) { stopAnimation(); } else { startAnimation(); } }
function startAnimation() { if (animationState.running) return; animationState.running = true; playPauseBtn.innerHTML = '<i class="fa-solid fa-pause"></i> Pause'; animationState.lastTimestamp = performance.now(); if (typeof animationState.time !== 'number' || isNaN(animationState.time)) { animationState.time = 0; } animationState.requestId = requestAnimationFrame(animateFrame); }
function stopAnimation() { if (!animationState.running) return; animationState.running = false; playPauseBtn.innerHTML = '<i class="fa-solid fa-play"></i> Play'; if (animationState.requestId) { cancelAnimationFrame(animationState.requestId); animationState.requestId = null; } }
function animateFrame(timestamp) { if (!animationState.running) return; let deltaTime = 0; if (animationState.lastTimestamp > 0) { deltaTime = (timestamp - animationState.lastTimestamp) / 1000.0; } animationState.lastTimestamp = timestamp; animationState.time += deltaTime; if (window.updateCanvasVisualization) { const visVector = (!isNaN(initialJonesVector[0].re)) ? initialJonesVector : [complex(1),complex(0)]; window.updateCanvasVisualization(opticalElements, visVector, animationState); } if (animationState.running) { animationState.requestId = requestAnimationFrame(animateFrame); } }

// --- Initialization ---
function init() {
    console.log("Script.js: Initializing...");
    Plotly.newPlot(intensityPlotDiv, [], {}, {responsive: true});
    addElementBtn.addEventListener('click', handleAddElement);
    playPauseBtn.addEventListener('click', toggleAnimation);
    showElementsCheckbox.addEventListener('change', updateIntensityPlot);

    // --- ADDED Event Listeners for Export/Import ---
    if (exportBtn) exportBtn.addEventListener('click', handleExport);
    if (importBtn) importBtn.addEventListener('click', handleImport);
    if (importFileInput) importFileInput.addEventListener('change', handleFileSelect);
    // --- END Added Listeners ---

    setInitialPolarization(initialStateType);
    if (window.initCanvas) { const canvasElement = document.getElementById('visualizationCanvas'); if (canvasElement) { window.initCanvas(canvasElement, opticalElements); console.log("Script.js: Canvas initialization called."); } else { console.error("Canvas element #visualizationCanvas not found!"); } } else { console.error("initCanvas function not found on window object. Check script load order."); }
    recalculateSystem();
    console.log("Script.js: Initialization complete.");
}

// --- Global Access Function (Canvas -> Script) ---
window.updateElementParameterFromCanvas = (id, param, value) => {
    const element = opticalElements.find(el => el.id === id); if (element) { let needsRecalculate = false; const tolerance = 1e-4; if (param === 'position') { const numericValue = parseFloat(value); if (!isNaN(numericValue)) { const clampedValue = Math.max(0, numericValue); if (Math.abs((element.position ?? 0) - clampedValue) > tolerance) { element.position = clampedValue; needsRecalculate = true; } } } else if (param === 'angle') { const numericValue = parseFloat(value); if (!isNaN(numericValue)) { const currentAngle = element.angle ?? 0; if (Math.abs(currentAngle - numericValue) > tolerance) { element.angle = numericValue; needsRecalculate = true; } } } if (needsRecalculate) { recalculateSystem(); } } else { console.warn(`Script: Element with id ${id} not found for update from canvas.`); }
};

// --- Run Initialization ---
document.addEventListener('DOMContentLoaded', init);
// interactiveCanvas.js (Corrected Label Clearing & Animation Direction)
// Manages the Three.js 3D visualization: elements, light field, labels, interaction.

// Ensure Three.js, OrbitControls, CSS2DRenderer are loaded before this script

// =============================================================================
// Global Variables and Constants
// =============================================================================

// --- Three.js Core Components ---
let scene, camera, renderer, controls;
let css2DRenderer; // For HTML labels

// --- Interaction Variables ---
let raycaster, mouse;
let draggableObjects = [];     // Holds meshes and handles for raycasting
let selectedObject = null;     // Currently interacted object (mesh or handle)
let plane;                   // Interaction plane for dragging/rotation
let isDraggingZ = false;       // True if dragging element along Z-axis
let isRotating = false;        // True if rotating element handle
let rotatingElementMesh = null; // Reference to the main mesh being rotated
let dragOffsetZ = 0;           // Offset for smoother Z dragging

// --- Scene Object Groups ---
let lightFieldGroup;          // Contains E-field arrows and envelope lines
let elementGroup;             // Contains element meshes, axes, labels, handles
let reflectedCoordHelper = null; // Group for the reflected coordinate system axes

// --- Visibility State (UI Controlled) ---
let showLabels = true;
let showEField = true;
let showEnvelope = true;
let isDarkMode = false;        // State for canvas theme (false = light, true = dark)

// --- Constants ---
// Visuals
const ELEMENT_RADIUS = 0.5;
const ELEMENT_THICKNESS = 0.05;
const AXIS_EXTENSION = 0.3;     // How far axis line extends beyond radius
const AXIS_COLOR = 0xeeeeee;    // Color of element's orientation axis
const AXIS_LINEWIDTH = 2;
const HANDLE_RADIUS = 0.08;    // Size of the rotation handle sphere
const HANDLE_COLOR = 0xff8c00; // Orange for handle
const HANDLE_HOVER_COLOR = 0xffcc00; // Lighter orange on hover/drag
const LIGHT_COLOR = 0x00ffff;    // Cyan for forward E-field vectors
const ENVELOPE_COLOR = 0xff00ff; // Magenta for forward envelope line
const REFLECTED_LIGHT_COLOR = 0xffa500; // Orange for reflected light
const REFLECTED_ENVELOPE_COLOR = 0xffd700; // Gold for reflected envelope
const CANVAS_BG_LIGHT = 0xf0f0f0; // White background
const CANVAS_BG_DARK = 0x252525;  // Black background
// Math & Physics
const AXIS_LENGTH = 10;        // Visual length of the main Z-axis guide
const WORLD_Z_AXIS = new THREE.Vector3(0, 0, 1); // World Z-axis vector
const EPSILON_INTERACTION = 1e-4; // Tolerance for detecting changes in drag/rotation
// Element Types
const ROTATABLE_ELEMENT_TYPES = [
    'Linear Polarizer',
    'HWP',
    'QWP',
    'General Waveplate',
    'Arbitrary Birefringent'
];

// --- DOM Element References ---
let showLabelsCheckbox, showEFieldCheckbox, showEnvelopeCheckbox;
let canvasContainer;


// =============================================================================
// Initialization (initCanvas)
// =============================================================================
window.initCanvas = function(canvasElement, initialElements) {
    console.log("InteractiveCanvas: Initializing...");
    if (!canvasElement) { console.error("Canvas element not provided!"); return; }
    canvasContainer = canvasElement.parentElement;
    if (!canvasContainer) { console.error("Canvas container not found!"); return; }

    scene = new THREE.Scene();
    // Set initial background color (default white)
    scene.background = new THREE.Color(CANVAS_BG_DARK);
    isDarkMode = true; // Ensure initial state matches background

    const aspect = canvasContainer.clientWidth / canvasContainer.clientHeight;
    camera = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000);
    camera.position.set(-3, 4, 6);
    // Look towards the middle of the typical element area
    const lookAtTarget = new THREE.Vector3(AXIS_LENGTH * 0.25, 0, AXIS_LENGTH * 0.25); // Adjusted lookAt
    camera.lookAt(lookAtTarget); // Initial look at

    renderer = new THREE.WebGLRenderer({ canvas: canvasElement, antialias: true });
    renderer.setSize(canvasContainer.clientWidth, canvasContainer.clientHeight);
    renderer.setPixelRatio(window.devicePixelRatio);

    const css2dContainer = document.getElementById('css2d-renderer');
    if (!css2dContainer) { console.error("CSS2D Renderer container div ('#css2d-renderer') not found!"); return; }
    css2DRenderer = new THREE.CSS2DRenderer();
    css2DRenderer.setSize(canvasContainer.clientWidth, canvasContainer.clientHeight);
    css2dContainer.appendChild(css2DRenderer.domElement);
    const css2dElement = css2DRenderer.domElement;
    css2dElement.style.position = 'absolute';
    css2dElement.style.top = '0px';
    css2dElement.style.left = '0px';
    css2dElement.style.pointerEvents = 'none';

    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.1;
    controls.screenSpacePanning = false;
    controls.target.copy(lookAtTarget); // Set OrbitControls target
    controls.update();

    const ambientLight = new THREE.AmbientLight(0xcccccc, 0.7);
    scene.add(ambientLight);
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight.position.set(3, 5, 4);
    scene.add(directionalLight);

    raycaster = new THREE.Raycaster();
    mouse = new THREE.Vector2();
    plane = new THREE.Plane();

    lightFieldGroup = new THREE.Group();
    scene.add(lightFieldGroup);
    elementGroup = new THREE.Group();
    scene.add(elementGroup);

    // --- World Axes Helper ---
    const axesSize = AXIS_LENGTH / 2; // Use the variable defined earlier
    const worldAxesHelper = new THREE.AxesHelper(axesSize);
    scene.add(worldAxesHelper);

    // --- Add World Axes Labels --- START ---
    const labelOffset = 0.2; // Distance beyond the axis line tip

    const createWorldAxisLabel = (text, color, position) => {
        const div = document.createElement('div');
        div.textContent = text;
        div.style.color = 'white';
        div.style.backgroundColor = color; // Use axis color for background
        div.style.fontSize = '12px';      // Slightly larger for world axes
        div.style.fontWeight = 'bold';
        div.style.fontFamily = 'sans-serif';
        div.style.padding = '1px 4px';
        div.style.borderRadius = '3px';
        div.style.pointerEvents = 'none'; // Ensure labels don't block interaction
        const label = new THREE.CSS2DObject(div);
        label.position.copy(position);
        label.layers.set(0); // Render on the same layer as other labels
        return label;
    };

    // Position labels slightly beyond the end of each axis line
    const xLabelPos = new THREE.Vector3(axesSize + labelOffset, 0, 0);
    const yLabelPos = new THREE.Vector3(0, axesSize + labelOffset, 0);
    const zLabelPos = new THREE.Vector3(0, 0, axesSize + labelOffset);

    const xLabel = createWorldAxisLabel("X", 'red', xLabelPos);
    const yLabel = createWorldAxisLabel("Y", 'green', yLabelPos);
    const zLabel = createWorldAxisLabel("Z", 'blue', zLabelPos);

    scene.add(xLabel);
    scene.add(yLabel);
    scene.add(zLabel);
    // --- Add World Axes Labels --- END ---


    // --- Main Z-Axis Guide Line ---
    const zAxisMaterial = new THREE.LineBasicMaterial({ color: 0xaaaaaa, transparent: true, opacity: 0.5 });
    const zAxisGeometry = new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(0, 0, -AXIS_LENGTH * 0.1),
        new THREE.Vector3(0, 0, AXIS_LENGTH * 1.5) // Extend further if needed
    ]);
    const zAxisLine = new THREE.Line(zAxisGeometry, zAxisMaterial);
    scene.add(zAxisLine);

    setupVisibilityControls();
    setupThemeToggleButton(); // <-- Add this line

    renderer.domElement.addEventListener('pointerdown', onPointerDown, false);
    renderer.domElement.addEventListener('pointermove', onPointerMove, false);
    renderer.domElement.addEventListener('pointerup', onPointerUp, false);
    renderer.domElement.addEventListener('pointerleave', onPointerUp, false); // Handle leaving canvas
    window.addEventListener('resize', onWindowResize, false);

    draggableObjects = [];
    syncCanvasElements(initialElements || []);

    animateCanvas();
    console.log("InteractiveCanvas: Initialization Complete.");
}

// =============================================================================
// UI Control Setup (Checkboxes)
// =============================================================================
function setupVisibilityControls() {
    showLabelsCheckbox = document.getElementById('show-labels-checkbox');
    showEFieldCheckbox = document.getElementById('show-efield-checkbox');
    showEnvelopeCheckbox = document.getElementById('show-envelope-checkbox');

    if (!showLabelsCheckbox || !showEFieldCheckbox || !showEnvelopeCheckbox) {
        console.error("One or more visualization checkboxes not found!");
        return;
    }

    showLabelsCheckbox.checked = showLabels;
    showEFieldCheckbox.checked = showEField;
    showEnvelopeCheckbox.checked = showEnvelope;

    showLabelsCheckbox.addEventListener('change', () => {
        showLabels = showLabelsCheckbox.checked;
        elementGroup.children.forEach(mesh => {
            if (mesh.isMesh && mesh.userData.labelObject instanceof THREE.CSS2DObject) {
                mesh.userData.labelObject.visible = showLabels;
            }
        });
    });

    showEFieldCheckbox.addEventListener('change', () => {
        showEField = showEFieldCheckbox.checked;
        triggerLightFieldUpdate();
    });

    showEnvelopeCheckbox.addEventListener('change', () => {
        showEnvelope = showEnvelopeCheckbox.checked;
        triggerLightFieldUpdate();
    });
}

// =============================================================================
// NEW: Theme Toggle Button Setup
// =============================================================================
function setupThemeToggleButton() {
    const themeToggleBtn = document.getElementById('theme-toggle-btn');
    if (!themeToggleBtn) {
        console.error("Theme toggle button ('#theme-toggle-btn') not found!");
        return;
    }

    themeToggleBtn.addEventListener('click', () => {
        if (!scene) return; // Safety check
        isDarkMode = !isDarkMode; // Toggle the state
        if (isDarkMode) {
            scene.background.setHex(CANVAS_BG_DARK); // Set to black
            themeToggleBtn.title = "Switch to Light Canvas Theme";
        } else {
            scene.background.setHex(CANVAS_BG_LIGHT); // Set to white
            themeToggleBtn.title = "Switch to Dark Canvas Theme";
        }
        // The animateCanvas loop will handle rendering the change
    });
}


// Helper function to request a light field redraw
function triggerLightFieldUpdate() {
     if (typeof window.updateCanvasVisualization === 'function' &&
        window.opticalElements !== undefined &&
        window.initialJonesVector !== undefined &&
        window.animationState !== undefined)
    {
        window.updateCanvasVisualization(
            window.opticalElements,
            window.initialJonesVector,
            window.animationState
        );
    } else {
        console.warn("Cannot trigger light field update: Crucial data/function missing from window scope.");
         console.log("Window state at failure:", { // Log state
             opticalElementsExists: window.opticalElements !== undefined,
             initialJonesVectorExists: window.initialJonesVector !== undefined,
             animationStateExists: window.animationState !== undefined,
             updateCanvasVisualizationExists: typeof window.updateCanvasVisualization === 'function'
         });
    }
}


// =============================================================================
// Light Field Visualization (Called by script.js)
// =============================================================================
window.updateCanvasVisualization = (elements, initialVector, animState) => {
    if (!scene || !lightFieldGroup) {
        console.warn("updateCanvasVisualization called before scene/group ready.");
        return;
    }

    // --- Clear previous light field objects ---
    while (lightFieldGroup.children.length > 0) {
        const object = lightFieldGroup.children[0];
        if (object.geometry) object.geometry.dispose();
        if (object.material) {
            if (Array.isArray(object.material)) object.material.forEach(m => m.dispose());
            else object.material.dispose();
        }
        lightFieldGroup.remove(object);
    }
    // --- Remove reflected coord helper if it exists ---
    removeReflectedCoordSystem(); // Call helper to remove and clean up

    // --- Basic Input Validation ---
    if (!Array.isArray(elements)) { console.error("Invalid 'elements' array", elements); elements = []; }
    if (!Array.isArray(initialVector) || initialVector.length !== 2 || !(initialVector[0] instanceof Complex) || !(initialVector[1] instanceof Complex)) { console.error("Invalid 'initialVector'", initialVector); initialVector = [window.complex(1), window.complex(0)]; }
    if (!animState) { console.error("Invalid 'animState'"); animState = { running: false, time: 0, omega: 2 * Math.PI * 0.5, k: 2 * Math.PI / 1 }; }

    const segmentsPerUnit = 20;

    // --- Find Double-Pass Mirror ---
    let doublePassMirrorFound = false;
    let doublePassMirrorIndex = -1;
    let doublePassMirrorPosition = -1;
    let vectorAtMirrorInput = null;

    elements.forEach((el, index) => {
        if (el.type === 'Mirror' && el.parameters?.isDoublePass) {
            doublePassMirrorFound = true;
            doublePassMirrorIndex = index;
            doublePassMirrorPosition = el.position;
            vectorAtMirrorInput = el.inputVector;
        }
    });

    // --- Draw Forward Path ---
    let currentZ = 0;
    let currentVector = initialVector;
    for (let i = 0; i < elements.length; i++) {
        const el = elements[i];
        if (el.position === undefined || isNaN(el.position)) { console.warn(`Skipping forward segment draw: Invalid position for element ${el.id}`); continue; }
        const startZ = currentZ;
        const endZ = el.position;
        const vectorForSegment = currentVector;

        if (vectorForSegment && startZ < endZ) {
            drawLightSegment(startZ, endZ, vectorForSegment, segmentsPerUnit, animState, false); // isReflected = false
        } else if (startZ > endZ) { console.warn(`Skipping forward segment draw: startZ ${startZ.toFixed(3)} >= endZ ${endZ.toFixed(3)} for element ${el.id}`); }

        currentVector = el.outputVector;
        currentZ = endZ;

        if (i === doublePassMirrorIndex) break;
    }

    // Final forward segment
    if (!doublePassMirrorFound) {
         const finalStartZ = currentZ;
         const finalEndZ = finalStartZ + AXIS_LENGTH * 0.3;
         const finalVector = currentVector;
         if (finalVector && finalStartZ < finalEndZ) {
              drawLightSegment(finalStartZ, finalEndZ, finalVector, segmentsPerUnit, animState, false);
         }
    }

    // --- Draw Reflected Path ---
if (doublePassMirrorFound && vectorAtMirrorInput) {
    if (!Array.isArray(vectorAtMirrorInput) || vectorAtMirrorInput.length !== 2 ||
         !(vectorAtMirrorInput[0] instanceof Complex) || !(vectorAtMirrorInput[1] instanceof Complex) ||
         isNaN(vectorAtMirrorInput[0].re) || isNaN(vectorAtMirrorInput[1].re)) {
        console.warn("Cannot draw reflected path: Invalid vectorAtMirrorInput.", vectorAtMirrorInput);
    } else {
         let currentReflectedVector = window.multiplyMatrixVector(window.REFLECTION_COORD_TRANSFORM, vectorAtMirrorInput);
         // Use the mirror's absolute position for the reflected beam origin.
         const mirrorPos = doublePassMirrorPosition;
         let currentReflectedZ = 0; // Reflected coordinates start at 0 (relative to the mirror)
         drawReflectedCoordSystem(mirrorPos); // Optionally draw helper axes at the mirror

         // Compute phaseOffset so that at z=0 the phase becomes (k*mirrorPos + π - ω*t).
         const phaseOffset = animState.k * mirrorPos;

         if (!isNaN(currentReflectedVector[0]?.re) && !isNaN(currentReflectedVector[1]?.re)) {
             // Iterate backwards through elements before the mirror
             for (let i = doublePassMirrorIndex - 1; i >= 0; i--) {
                 const el = elements[i];
                 if (el.position === undefined || isNaN(el.position)) {
                     console.warn(`Skipping reflected segment draw: Invalid position for element ${el.id}`);
                     continue;
                 }
                 const startZ_reflected = currentReflectedZ;
                 // Convert absolute z to relative z (relative to mirror)
                 const endZ_reflected = el.position - mirrorPos;
                 const vectorForSegment = currentReflectedVector;

                 if (vectorForSegment && startZ_reflected > endZ_reflected) {
                     // Pass mirrorPos as the extra offset so that the reflected beam is drawn starting at the mirror's position
                     drawLightSegment(startZ_reflected, endZ_reflected, vectorForSegment, segmentsPerUnit, animState, true, phaseOffset, mirrorPos);
                 } else if (startZ_reflected < endZ_reflected) {
                     console.warn(`Skipping reflected segment draw: startZ ${startZ_reflected.toFixed(3)} <= endZ ${endZ_reflected.toFixed(3)} for element ${el.id}`);
                 }
                 currentReflectedVector = el.reflectedOutputVector;
                 currentReflectedZ = endZ_reflected;

                 if (!currentReflectedVector || isNaN(currentReflectedVector[0]?.re) || isNaN(currentReflectedVector[1]?.re)) {
                     console.warn(`Stopping reflected path drawing due to invalid vector after element ${el.id}`);
                     break;
                 }
             }
             // Final reflected segment
             const finalReflectedStartZ = currentReflectedZ;
             const finalReflectedEndZ = Math.min(-0.1, finalReflectedStartZ - AXIS_LENGTH * 0.3);
             if (currentReflectedVector && finalReflectedStartZ > finalReflectedEndZ &&
                 !isNaN(currentReflectedVector[0]?.re)) {
                 drawLightSegment(finalReflectedStartZ, finalReflectedEndZ, currentReflectedVector, segmentsPerUnit, animState, true, phaseOffset, mirrorPos);
             }
         } else {
             console.warn("Cannot draw reflected path: Result of reflection transform is invalid.", currentReflectedVector);
         }
    }
}


};


// =============================================================================
// Draw Single Light Segment Helper (Corrected Animation Direction)
// =============================================================================
function drawLightSegment(startZ, endZ, jonesVector, pointsPerUnit, animState, isReflected = false, phaseOffset = 0, zOffset = 0) {
    // --- Input Validation ---
    if (jonesVector === undefined || jonesVector === null || !Array.isArray(jonesVector) ||
        jonesVector.length < 2 ||
        !(jonesVector[0] instanceof Complex) || !(jonesVector[1] instanceof Complex)) { return; }
    const [Ex, Ey] = jonesVector;
    if (isNaN(Ex.re) || isNaN(Ex.im) || isNaN(Ey.re) || isNaN(Ey.im)) { return; }
    if (Math.abs(endZ - startZ) < EPSILON_INTERACTION) { return; }
    if ((!isReflected && startZ > endZ) || (isReflected && startZ < endZ)) { return; }
    if (isNaN(startZ) || isNaN(endZ) || isNaN(pointsPerUnit) || !animState) {
       console.error("drawLightSegment: Invalid numeric parameter or animState."); return;
    }

    // --- Setup ---
    const length = Math.abs(endZ - startZ);
    const numPoints = Math.max(2, Math.ceil(length * pointsPerUnit));
    const zStep = (endZ - startZ) / (numPoints - 1);
    const envelopePoints = [];
    const segmentLightColor    = isReflected ? REFLECTED_LIGHT_COLOR    : LIGHT_COLOR;
    const segmentEnvelopeColor = isReflected ? REFLECTED_ENVELOPE_COLOR : ENVELOPE_COLOR;
    const arrowDensityFactor = Math.max(1, Math.floor(pointsPerUnit / 10));

    // --- Generate Points and Arrows ---
    for (let i = 0; i < numPoints; i++) {
        // Compute the current z in the segment (relative coordinate)
        const z = startZ + i * zStep;
        // For reflected beams add the provided offset so that z=0 becomes mirrorPos
        const effectiveZ = isReflected ? (z + zOffset) : z;

        const timePhaseTerm = (animState.omega || 0) * (animState.time || 0);
        // The phase remains computed with the z relative to the mirror
        const spatialPhaseTerm = (isReflected ? -animState.k : animState.k) * z + (isReflected ? phaseOffset : 0);
        const phase = spatialPhaseTerm - timePhaseTerm;
        const expPhase = window.complex(0, phase).exp();
        if (isNaN(expPhase.re) || isNaN(expPhase.im)) continue;

        // Multiply the Jones components (already adjusted by the reflection sign as needed)
        const ex_inst = isReflected ? -Ex.mul(expPhase).re : Ex.mul(expPhase).re;
        const ey_inst = Ey.mul(expPhase).re;
        // IMPORTANT: Use effectiveZ here so that the drawn point is shifted by zOffset when reflected
        const fieldVectorTip = new THREE.Vector3(ex_inst, ey_inst, effectiveZ);
        if (isNaN(fieldVectorTip.x) || isNaN(fieldVectorTip.y)) continue;

        // Add point to envelope line if enabled
        if (showEnvelope) { envelopePoints.push(fieldVectorTip); }

        // Conditionally draw E-field arrows on selected points
        const drawArrowCondition = (i % arrowDensityFactor === 0) || (i === numPoints - 1);
        if (showEField && drawArrowCondition) {
            const originPoint = new THREE.Vector3(0, 0, effectiveZ);
            const directionVector = fieldVectorTip.clone().sub(originPoint);
            const vecLength = directionVector.length();
            if (vecLength > EPSILON_INTERACTION && !isNaN(vecLength)) {
                const normalizedDirection = directionVector.normalize();
                const arrowLength = vecLength;
                const headLength = Math.max(0.02, Math.min(arrowLength * 0.2, 0.15));
                const headWidth  = Math.max(0.01, Math.min(arrowLength * 0.1, 0.1));
                const arrowHelper = new THREE.ArrowHelper(normalizedDirection, originPoint, arrowLength, segmentLightColor, headLength, headWidth);
                lightFieldGroup.add(arrowHelper);
            }
        }
    }
    // --- Draw Envelope Line ---
    if (showEnvelope && envelopePoints.length >= 2) {
        if (envelopePoints.some(p => isNaN(p.x) || isNaN(p.y) || isNaN(p.z))) {
            console.error("NaN detected in envelopePoints, skipping envelope draw.");
        } else {
            const envelopeMaterial = new THREE.LineBasicMaterial({ color: segmentEnvelopeColor });
            const envelopeGeometry = new THREE.BufferGeometry().setFromPoints(envelopePoints);
            const envelopeLine = new THREE.Line(envelopeGeometry, envelopeMaterial);
            lightFieldGroup.add(envelopeLine);
        }
    }
}



// =============================================================================
// Element Visual Creation and Management
// =============================================================================
function createElementVisual(element) {
    if (!element || element.id === undefined || element.type === undefined || element.position === undefined || isNaN(element.position)) { console.error("Invalid element data for visual creation:", element); return null; }

    const material = new THREE.MeshPhongMaterial({ color: getColorForElement(element.type), opacity: 0.75, transparent: true, side: THREE.DoubleSide });
    const geometry = new THREE.CylinderGeometry(ELEMENT_RADIUS, ELEMENT_RADIUS, ELEMENT_THICKNESS, 32);
    const mesh = new THREE.Mesh(geometry, material);
    mesh.rotation.order = 'XYZ';
    mesh.rotation.x = Math.PI / 2;
    mesh.position.z = element.position;
    mesh.userData = { id: element.id, type: element.type, isElementBody: true };

    const axisMaterial = new THREE.LineBasicMaterial({ color: AXIS_COLOR, linewidth: AXIS_LINEWIDTH });
    const axisLength = ELEMENT_RADIUS + AXIS_EXTENSION;
    const axisPoints = [ new THREE.Vector3(-axisLength, 0, ELEMENT_THICKNESS / 2 + 0.01), new THREE.Vector3(axisLength, 0, ELEMENT_THICKNESS / 2 + 0.01) ];
    const axisGeometry = new THREE.BufferGeometry().setFromPoints(axisPoints);
    const axisLine = new THREE.Line(axisGeometry, axisMaterial);
    axisLine.userData.isAxisLine = true;
    mesh.add(axisLine);

    const labelDiv = document.createElement('div');
    labelDiv.className = 'element-label';
    labelDiv.textContent = element.type;
    const label = new THREE.CSS2DObject(labelDiv);
    label.position.set(0, ELEMENT_RADIUS * 0.7, ELEMENT_THICKNESS / 2 + 0.05);
    label.visible = showLabels;
    label.layers.set(0);
    mesh.add(label);
    mesh.userData.labelObject = label;

    if (ROTATABLE_ELEMENT_TYPES.includes(element.type)) {
        const handleGeometry = new THREE.SphereGeometry(HANDLE_RADIUS, 16, 8);
        const handleMaterial = new THREE.MeshBasicMaterial({ color: HANDLE_COLOR, depthTest: false, transparent: true });
        const handle = new THREE.Mesh(handleGeometry, handleMaterial);
        handle.position.set(axisLength, 0, ELEMENT_THICKNESS / 2 + 0.01);
        handle.renderOrder = 1;
        handle.userData = { parentElementId: element.id, isRotationHandle: true, baseColor: HANDLE_COLOR, hoverColor: HANDLE_HOVER_COLOR };
        mesh.add(handle);
        mesh.userData.rotationHandle = handle;
        draggableObjects.push(handle);
    }

    elementGroup.add(mesh);
    draggableObjects.push(mesh);
    return mesh;
}

function getColorForElement(type) {
    if (type.includes('Polarizer')) return 0xAAAAAA;
    if (type === 'Mirror') return 0xC0C0C0;
    if (type === 'HWP') return 0xff8888;
    if (type === 'QWP') return 0x88ff88;
    if (type === 'General Waveplate') return 0x8888ff;
    if (type === 'Arbitrary Birefringent') return 0xff88ff;
    if (type === 'Faraday Rotator') return 0x88ffff;
    return 0xcccccc;
}

window.syncCanvasElements = (elements) => {
    if (!elementGroup || !draggableObjects) { console.warn("syncCanvasElements called before canvas fully initialized."); return; }
    if (!Array.isArray(elements)) { console.error("Invalid 'elements' array passed to syncCanvasElements:", elements); elements = []; }

    const currentIds = elements.map(el => el.id);
    const existingMeshes = elementGroup.children.filter(c => c.isMesh && c.userData.isElementBody);
    const existingMeshIds = existingMeshes.map(mesh => mesh.userData.id);

    // Remove visuals
    const idsToRemove = existingMeshIds.filter(id => !currentIds.includes(id));
    idsToRemove.forEach(id => { window.removeCanvasElement(id); });

    // Add/Update visuals
    elements.forEach(el => {
        if (el.id === undefined || el.type === undefined || el.position === undefined || isNaN(el.position)) { console.warn("Skipping sync for element with missing/invalid data:", el); return; }

        let visualMesh = existingMeshes.find(mesh => mesh.userData.id === el.id);

        if (!visualMesh) {
            visualMesh = createElementVisual(el);
            if (!visualMesh) return;
        }
        if (!visualMesh) { console.error(`Sync: Could not find or create visual for element ${el.id}`); return; }

        // Update Position
        if (Math.abs(visualMesh.position.z - el.position) > EPSILON_INTERACTION) { visualMesh.position.z = el.position; }

        // Update Rotation
        const targetAngleRad = (el.angle || 0) * DEG_TO_RAD;
        if (Math.abs(visualMesh.rotation.y - targetAngleRad) > EPSILON_INTERACTION) { visualMesh.rotation.y = targetAngleRad; }
        visualMesh.rotation.x = Math.PI / 2; visualMesh.rotation.z = 0;

        // Update Color
        const targetColorHex = getColorForElement(el.type);
        if (visualMesh.material.color.getHex() !== targetColorHex) {
            visualMesh.material.color.setHex(targetColorHex);
            visualMesh.userData.type = el.type;
        }

        // Update Label
        if (visualMesh.userData.labelObject) {
            const labelElement = visualMesh.userData.labelObject.element;
            if (labelElement.textContent !== el.type) { labelElement.textContent = el.type; }
            visualMesh.userData.labelObject.visible = showLabels;
        }

        // Update Handle Presence (Recreate if mismatch)
        const shouldHaveHandle = ROTATABLE_ELEMENT_TYPES.includes(el.type);
        const hasHandle = !!visualMesh.userData.rotationHandle;
        if (shouldHaveHandle !== hasHandle) {
             console.warn(`Handle presence mismatch for element ${el.id}. Recreating visual.`);
             window.removeCanvasElement(el.id);
             createElementVisual(el);
        }
    });
};

window.removeCanvasElement = (id) => {
    if (!elementGroup || !draggableObjects) {
        console.warn("removeCanvasElement called before canvas initialized or elementGroup missing.");
        return;
    }
    const meshToRemove = elementGroup.children.find(child => child.isMesh && child.userData.isElementBody && child.userData.id === id);

    if (meshToRemove) {
        // 1. Remove from draggable objects list
        draggableObjects = draggableObjects.filter(d => d !== meshToRemove);
        if (meshToRemove.userData.rotationHandle) {
            // Also remove the handle if it exists
            draggableObjects = draggableObjects.filter(d => d !== meshToRemove.userData.rotationHandle);
        }
        const labelObject = meshToRemove.userData.labelObject;
        if (labelObject && labelObject instanceof THREE.CSS2DObject) {
            const labelDiv = labelObject.element; // Get the underlying HTML div
            if (labelDiv && labelDiv.parentNode) {
                // Remove the div from its parent in the HTML DOM
                labelDiv.parentNode.removeChild(labelDiv);
                // console.log("Removed label div:", labelDiv.textContent); // Optional debug log
            } else {
                // console.warn(`Could not remove label div for element ${id}. ParentNode missing?`, labelDiv); // Optional warning
            }
            // No need to dispose geometry/material for CSS2DObject itself
        }

        // 3. Dispose WebGL resources of other children (axis lines, handles, etc.)
        meshToRemove.children.forEach(child => {
            // Skip the label object we just handled
            if (child !== labelObject) {
                if (child.geometry) {
                    child.geometry.dispose();
                }
                if (child.material) {
                    // Handle both single material and arrays of materials
                    if (Array.isArray(child.material)) {
                        child.material.forEach(m => m.dispose());
                    } else {
                        child.material.dispose();
                    }
                }
            }
        });

        // 4. Remove the main mesh (and its children, including the label *object*) from the Three.js scene graph
        elementGroup.remove(meshToRemove);

        // 5. Dispose the main mesh's own geometry and material
        if (meshToRemove.geometry) {
            meshToRemove.geometry.dispose();
        }
        if (meshToRemove.material) {
            if (Array.isArray(meshToRemove.material)) {
                meshToRemove.material.forEach(m => m.dispose());
            } else {
                meshToRemove.material.dispose();
            }
        }
        // console.log(`Removed element ${id} visual from canvas.`); // Optional log
    } else {
         // console.warn(`removeCanvasElement: Mesh with ID ${id} not found.`); // Optional warning
    }
};


// =============================================================================
// Reflected Coordinate System Helper (Corrected Label Removal)
// =============================================================================
function drawReflectedCoordSystem(zPos) {
     if (!scene) return;
     removeReflectedCoordSystem(); // Ensure any old one is gone first

     reflectedCoordHelper = new THREE.Group();
     reflectedCoordHelper.name = "reflectedCoordHelper";
     const arrowLength = 0.6;
     const headLength = 0.12;
     const headWidth = 0.06;
     const origin = new THREE.Vector3(0, 0, zPos);

     reflectedCoordHelper.add(new THREE.ArrowHelper(new THREE.Vector3(-1, 0, 0), origin, arrowLength, 0xff0000, headLength, headWidth)); // X'
     reflectedCoordHelper.add(new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), origin, arrowLength, 0x00ff00, headLength, headWidth)); // Y'
     reflectedCoordHelper.add(new THREE.ArrowHelper(new THREE.Vector3(0, 0, -1), origin, arrowLength, 0x0000ff, headLength, headWidth)); // Z'

     // Add CSS2D labels
     if (typeof THREE.CSS2DObject !== 'undefined' && css2DRenderer) {
          const createAxisLabel = (text, color, positionOffset) => {
              const div = document.createElement('div');
              div.textContent = text;
              div.style.color = 'white'; div.style.backgroundColor = color;
              div.style.fontSize = '10px'; div.style.fontFamily = 'sans-serif';
              div.style.padding = '0px 3px'; div.style.borderRadius = '2px';
              const label = new THREE.CSS2DObject(div);
              label.position.copy(origin).add(positionOffset);
              label.layers.set(0);
              // Store reference to the div for later removal
              label.userData.isAxisLabelDiv = true;
              label.userData.labelDiv = div;
              return label;
          };
          const labelDist = arrowLength * 1.1;
          reflectedCoordHelper.add(createAxisLabel("x'", 'red', new THREE.Vector3(-labelDist, 0, 0)));
          reflectedCoordHelper.add(createAxisLabel("y'", 'green', new THREE.Vector3(0, labelDist, 0)));
          reflectedCoordHelper.add(createAxisLabel("z'", 'blue', new THREE.Vector3(0, 0, -labelDist)));
     }
     scene.add(reflectedCoordHelper);
}

function removeReflectedCoordSystem() {
    if (reflectedCoordHelper) {
         // --- Explicitly remove CSS2D label divs from DOM ---
         reflectedCoordHelper.children.forEach(child => {
             if (child instanceof THREE.CSS2DObject && child.userData.isAxisLabelDiv) {
                 const labelDiv = child.userData.labelDiv;
                 if (labelDiv && labelDiv.parentNode) {
                     labelDiv.parentNode.removeChild(labelDiv);
                      // console.log("Removed CSS2D label div:", labelDiv.textContent); // Debug log
                 }
             }
         });
         // ----------------------------------------------------

         // Dispose geometries/materials of THREE objects (ArrowHelpers)
         reflectedCoordHelper.children.forEach(child => {
             if (child instanceof THREE.ArrowHelper) {
                 if (child.line && child.line.geometry) child.line.geometry.dispose();
                 if (child.line && child.line.material) child.line.material.dispose();
                 if (child.cone && child.cone.geometry) child.cone.geometry.dispose();
                 if (child.cone && child.cone.material) child.cone.material.dispose();
             }
         });
        scene.remove(reflectedCoordHelper);
        reflectedCoordHelper = null;
    }
}


// =============================================================================
// User Interaction Handlers (Pointer Events)
// =============================================================================
function onPointerDown(event) {
    if (isDraggingZ || isRotating) return;
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(draggableObjects, false);

    if (intersects.length > 0) {
        const intersectedObject = intersects[0].object;
        const intersectionPoint = intersects[0].point;
        controls.enabled = false;

        if (intersectedObject.userData.isRotationHandle) {
            selectedObject = intersectedObject;
            rotatingElementMesh = selectedObject.parent;
            if (!rotatingElementMesh || !rotatingElementMesh.userData.isElementBody) { console.error("Could not find parent mesh for rotation handle!"); controls.enabled = true; return; }
            isRotating = true; isDraggingZ = false;
            plane.setFromNormalAndCoplanarPoint(WORLD_Z_AXIS, rotatingElementMesh.position);
            renderer.domElement.style.cursor = 'grabbing';
            if (selectedObject.material) selectedObject.material.color.set(selectedObject.userData.hoverColor);
        } else if (intersectedObject.userData.isElementBody) {
            selectedObject = intersectedObject;
            rotatingElementMesh = null; isRotating = false; isDraggingZ = true;
            plane.setFromNormalAndCoplanarPoint(camera.getWorldDirection(plane.normal).negate(), intersectionPoint);
            dragOffsetZ = selectedObject.position.z - intersectionPoint.z;
            renderer.domElement.style.cursor = 'move';
        } else { controls.enabled = true; }
    } else {
        selectedObject = null; isDraggingZ = false; isRotating = false; rotatingElementMesh = null;
        renderer.domElement.style.cursor = 'auto'; controls.enabled = true;
    }
}

function onPointerMove(event) {
    if (!selectedObject) return;
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);
    const intersectionPoint = new THREE.Vector3();

    if (isRotating && rotatingElementMesh) {
        if (raycaster.ray.intersectPlane(plane, intersectionPoint)) {
            const elementCenterWorld = rotatingElementMesh.position.clone();
            const directionVector = intersectionPoint.clone().sub(elementCenterWorld);
            let angleRad = Math.atan2(directionVector.y, directionVector.x);
            let angleDeg = angleRad * RAD_TO_DEG;
            angleDeg = (angleDeg + 360) % 360;
            const roundedAngleDeg = Math.round(angleDeg);
            const finalAngleDeg = (roundedAngleDeg + 360) % 360;

            if (window.updateElementParameterFromCanvas) {
                 const elementData = window.opticalElements?.find(el => el.id === rotatingElementMesh.userData.id);
                 const currentIntAngleDeg = Math.round(elementData?.angle ?? 0);
                 const normalizedCurrentIntAngle = (currentIntAngleDeg + 360) % 360;
                 if (normalizedCurrentIntAngle !== finalAngleDeg) {
                     window.updateElementParameterFromCanvas(rotatingElementMesh.userData.id, 'angle', finalAngleDeg);
                 }
            }
        }
        renderer.domElement.style.cursor = 'grabbing';
    }
    else if (isDraggingZ) {
        if (raycaster.ray.intersectPlane(plane, intersectionPoint)) {
            const newZ = intersectionPoint.z + dragOffsetZ;
            const clampedZ = Math.max(0, newZ);
            if (Math.abs(selectedObject.position.z - clampedZ) > EPSILON_INTERACTION) {
                selectedObject.position.z = clampedZ;
                if (window.updateElementParameterFromCanvas) {
                    const newPosition = parseFloat(clampedZ.toFixed(3));
                    window.updateElementParameterFromCanvas(selectedObject.userData.id, 'position', newPosition);
                }
            }
        }
        renderer.domElement.style.cursor = 'move';
    }
}

function onPointerUp(event) {
    if (isRotating && selectedObject && selectedObject.userData.isRotationHandle) {
        if (selectedObject.material) selectedObject.material.color.set(selectedObject.userData.baseColor);
    }
    isDraggingZ = false; isRotating = false; selectedObject = null; rotatingElementMesh = null; dragOffsetZ = 0;
    controls.enabled = true; renderer.domElement.style.cursor = 'auto';
}


// =============================================================================
// Window Resize Handler
// =============================================================================
function onWindowResize() {
    if (!renderer || !camera || !canvasContainer || !css2DRenderer) return;
    const width = canvasContainer.clientWidth; const height = canvasContainer.clientHeight;
    camera.aspect = width / height; camera.updateProjectionMatrix();
    renderer.setSize(width, height); css2DRenderer.setSize(width, height);
}

// =============================================================================
// Animation Loop
// =============================================================================
function animateCanvas() {
    requestAnimationFrame(animateCanvas);
    controls.update();
    renderer.render(scene, camera);
    if (css2DRenderer) { css2DRenderer.render(scene, camera); }
}

// =============================================================================
// Expose Key Functions
// =============================================================================
// window.initCanvas, updateCanvasVisualization, syncCanvasElements, removeCanvasElement assigned above

console.log("InteractiveCanvas.js loaded.");
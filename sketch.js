let pDynamics, pHiC, pPsPlot; // p5 instances for each canvas

// Simulation parameters with default values
let params = {
    numChromosomes: 3,
    numMonomers: 50,
    simulationSpeed: 30,
    bondLength: 8,
    endPullStrength: 0.02,
    endRandomStep: 1,
    contactRadius: 12,
    hicDecayRate: 0.99,
    monomerSelfRepulsionStrength: 0.1,
    obstacleRepulsionStrength: 0.5,
    extruderForceStrength: 0.5,
    relaxationIterations: 5,
    isPlaying: true,
};

let chromosomes = [];
let hicData; 
let maxHiCValue = 1;
let psData = []; // For P(s) plot: [{s: genomic_dist, P: probability}]

let monomerAnchors = []; 
let circularObstacles = []; 
let nextObstacleId = 0;
let pendingObstacleCenter = null; 

let loopExtruders = [];
let nextExtruderId = 0;

let lastUpdateTime = 0;

// Cached CSS colors
let cachedColors = {};

// DOM Elements
let numChromosomesInput, numMonomersInput, simulationSpeedInput;
let bondLengthInput, endPullStrengthInput, endRandomStepInput, contactRadiusInput, hicDecayRateInput;
let monomerSelfRepulsionInput, obstacleRepulsionInput, extruderForceStrengthInput;
let playPauseButton, resetButton;
let anchorChrIdxInput, anchorMonomerIdxInput, anchorTypeSelect, addAnchorButton, clearAnchorsButton, anchorListUl;
let obstacleRadiusInput, addObstacleButton, clearObstaclesButton, obstacleListUl;
let extruderChrIdxInput, extruderTypeSelect, extruderSpeedInput, extruderMaxSizeInput;
let extruderStaticMonomer1Input, extruderStaticMonomer2Input, extruderStaticIndicesContainer;
let addExtruderButton, clearExtrudersButton, extruderListUl;
let endPullStrengthValueSpan, endRandomStepValueSpan, hicDecayRateValueSpan, monomerSelfRepulsionValueSpan, obstacleRepulsionValueSpan, extruderForceStrengthValueSpan;


const dynamicsSketch = (p) => {
    pDynamics = p;
    p.setup = function() {
        const container = document.getElementById('chromosome-canvas-container');
        const canvas = p.createCanvas(container.offsetWidth, container.offsetHeight);
        canvas.parent('chromosome-canvas-container');
        canvas.mousePressed(handleDynamicsCanvasClick);
        p.pixelDensity(1);
        if (Object.keys(cachedColors).length > 0) p.redraw(); // Redraw if colors are ready
    };

    p.draw = function() {
        if (!cachedColors.panelBgColor) { // Wait for colors to be cached
            p.background(50); return;
        }
        p.background(cachedColors.panelBgColor);
        
        chromosomes.forEach(chr => chr.display(p));
        
        monomerAnchors.forEach(anchor => {
            if (anchor.position) {
                p.fill(255, 255, 0, 200); // Yellow for explicit user anchors
                p.noStroke();
                p.ellipse(anchor.position.x, anchor.position.y, 10, 10);
            }
        });

        circularObstacles.forEach(obs => {
            p.stroke(cachedColors.secondaryAccentColor);
            p.strokeWeight(2);
            p.noFill();
            p.ellipse(obs.center.x, obs.center.y, obs.radius * 2, obs.radius * 2);
        });

        if (pendingObstacleCenter) {
            let tertiaryColorWithAlpha = p.color(cachedColors.tertiaryAccentColor);
            tertiaryColorWithAlpha.setAlpha(100);
            p.fill(tertiaryColorWithAlpha);
            p.noStroke();
            p.ellipse(pendingObstacleCenter.x, pendingObstacleCenter.y, 8, 8);
            p.textSize(12);
            p.fill(cachedColors.tertiaryAccentColor);
            p.text("Obstacle center selected", pendingObstacleCenter.x + 10, pendingObstacleCenter.y);
        }
    };
    
    p.windowResized = function() {
        const container = document.getElementById('chromosome-canvas-container');
        p.resizeCanvas(container.offsetWidth, container.offsetHeight);
    };

    function handleDynamicsCanvasClick() {
        if (!pDynamics || p.mouseX < 0 || p.mouseX > p.width || p.mouseY < 0 || p.mouseY > p.height) return;
        let clickPos = p.createVector(p.mouseX, p.mouseY);

        if (anchorTypeSelect.value === 'fix_to_click') {
            addMonomerAnchor(parseInt(anchorChrIdxInput.value), parseInt(anchorMonomerIdxInput.value), clickPos, 'click_fixed');
            return; 
        }
        
        pendingObstacleCenter = clickPos;
        p.redraw();
    }
};

const hiCSketch = (p) => {
    pHiC = p;
    p.setup = function() {
        const container = document.getElementById('hic-canvas-container');
        const canvas = p.createCanvas(container.offsetWidth, container.offsetHeight);
        canvas.parent('hic-canvas-container');
        p.pixelDensity(1);
        p.colorMode(p.RGB, 255);
        if (Object.keys(cachedColors).length > 0) p.redraw();
        p.noLoop(); 
    };

    p.draw = function() { 
        if (!cachedColors.panelBgColor || !cachedColors.secondaryAccentColor) {
             p.background(50); return;
        }
        p.background(cachedColors.panelBgColor);
        if (!hicData || params.numMonomers === 0) return;

        let M = params.numMonomers;
        let cellSize = Math.min(p.width / M, p.height / M);
        if (M === 0 || cellSize <= 0) return;

        let mapSizeW = cellSize * M;
        let mapSizeH = cellSize * M; 
        
        let offsetX = (p.width - mapSizeW) / 2;
        let offsetY = (p.height - mapSizeH) / 2;

        p.push();
        p.translate(offsetX, offsetY);
        p.noStroke();
        
        const lowColor = p.color(255); // White for low
        const highColor = cachedColors.secondaryAccentColor;

        for (let i = 0; i < M; i++) {
            for (let j = 0; j < M; j++) {
                let val = hicData[i][j];
                let intensity = p.map(val, 0, maxHiCValue > 0 ? maxHiCValue : 1, 0, 1);
                intensity = p.constrain(intensity, 0, 1);

                let c = p.lerpColor(lowColor, highColor, intensity * intensity); 
                p.fill(c);
                p.rect(i * cellSize, j * cellSize, cellSize, cellSize);
            }
        }
        p.pop();
    };
    
    p.windowResized = function() {
        const container = document.getElementById('hic-canvas-container');
        p.resizeCanvas(container.offsetWidth, container.offsetHeight);
        if (pHiC) pHiC.redraw();
    };
};

const psPlotSketch = (p) => {
    pPsPlot = p;
    p.setup = function() {
        const container = document.getElementById('ps-plot-container');
        const canvas = p.createCanvas(container.offsetWidth, container.offsetHeight);
        canvas.parent('ps-plot-container');
        p.pixelDensity(1);
        p.colorMode(p.RGB, 255);
        if (Object.keys(cachedColors).length > 0) p.redraw();
        p.noLoop(); 
    };

    p.draw = function() { 
        if (!cachedColors.panelBgColor || !cachedColors.textColor || !cachedColors.primaryAccentColor) {
            p.background(50); return;
        }
        p.background(cachedColors.panelBgColor);
        if (!psData || psData.length === 0) return;

        const margin = { top: 20, right: 20, bottom: 40, left: 50 };
        const plotWidth = p.width - margin.left - margin.right;
        const plotHeight = p.height - margin.top - margin.bottom;

        if (plotWidth <= 0 || plotHeight <= 0) return;

        let minS = Math.log10(psData[0].s > 0 ? psData[0].s : 1e-9);
        let maxS = Math.log10(psData[psData.length-1].s > 0 ? psData[psData.length-1].s : 1);
        if (minS === maxS) { minS -=1; maxS +=1; }
        
        let validP = psData.filter(d => d.P > 1e-9).map(d => d.P); // Use 1e-9 consistent with storage
        let minP = validP.length > 0 ? Math.log10(Math.min(...validP)) : -9; 
        let maxP = validP.length > 0 ? Math.log10(Math.max(...validP)) : 0;   
        if (minP === maxP) { minP -=1; maxP +=1; }

        p.push();
        p.translate(margin.left, margin.top);

        p.stroke(cachedColors.textColor);
        p.strokeWeight(1);
        p.line(0, plotHeight, plotWidth, plotHeight); 
        p.line(0, 0, 0, plotHeight); 

        p.fill(cachedColors.textColor);
        p.noStroke();
        p.textAlign(p.CENTER, p.TOP);
        p.text("log10(s)", plotWidth / 2, plotHeight + 10);
        p.textAlign(p.CENTER, p.CENTER); 
        p.push();
        p.translate(-margin.left + 15, plotHeight / 2);
        p.rotate(-p.HALF_PI);
        p.text("log10(P(s))", 0, 0);
        p.pop();
        
        for(let i = Math.floor(minS); i <= Math.ceil(maxS); i++) {
            let x = p.map(i, minS, maxS, 0, plotWidth);
            if (x >=0 && x <= plotWidth) {
                p.stroke(cachedColors.textColor);
                p.line(x, plotHeight, x, plotHeight + 5);
                p.noStroke(); p.fill(cachedColors.textColor);
                p.textAlign(p.CENTER);
                p.text(Math.pow(10,i).toExponential(0), x, plotHeight + 15);
            }
        }
         for(let i = Math.floor(minP); i <= Math.ceil(maxP); i++) {
            let y = p.map(i, minP, maxP, plotHeight, 0); 
             if (y >=0 && y <= plotHeight) {
                p.stroke(cachedColors.textColor);
                p.line(0, y, -5, y);
                p.noStroke(); p.fill(cachedColors.textColor);
                p.textAlign(p.RIGHT, p.CENTER);
                p.text(Math.pow(10,i).toExponential(0), -8, y);
            }
        }

        p.stroke(cachedColors.primaryAccentColor);
        p.strokeWeight(2);
        p.noFill();
        p.beginShape();
        psData.forEach(d => {
            if (d.P > 1e-9 && d.s > 0) { 
                let x = p.map(Math.log10(d.s), minS, maxS, 0, plotWidth);
                let y = p.map(Math.log10(d.P), minP, maxP, plotHeight, 0); 
                p.vertex(x, y);
            }
        });
        p.endShape();
        p.pop();
    };
    
    p.windowResized = function() {
        const container = document.getElementById('ps-plot-container');
        p.resizeCanvas(container.offsetWidth, container.offsetHeight);
        if (pPsPlot) pPsPlot.redraw();
    };
};


function globalUpdate() {
    const currentTime = Date.now();
    const timeStepDuration = 1000 / params.simulationSpeed;

    if (params.isPlaying && (currentTime - lastUpdateTime > timeStepDuration)) {
        let stepsToRun = Math.floor((currentTime - lastUpdateTime) / timeStepDuration);
        stepsToRun = Math.min(stepsToRun, 5); 

        const deltaTimePerStep = 1.0 / params.simulationSpeed; // for physics updates

        for (let i = 0; i < stepsToRun; i++) {
            runSimulationStep(deltaTimePerStep);
        }
        lastUpdateTime = currentTime;

        if (pDynamics) pDynamics.redraw();
        if (pHiC) pHiC.redraw();
        if (pPsPlot) pPsPlot.redraw();
    }
}
setInterval(globalUpdate, 1000 / 60); 

class Chromosome {
    constructor(id, numMonomers, p_instance) {
        this.p = p_instance; 
        this.id = id;
        this.numMonomers = numMonomers;
        this.monomers = []; 
        // Each element: { fixed: boolean, position: p5.Vector | null }
        this.isFixed = new Array(numMonomers).fill(null).map(() => ({ fixed: false, position: null }));
        
        this.p.push();
        this.p.colorMode(this.p.HSB, 360, 100, 100, 100);
        this.color = this.p.color((id * 60 + Math.random()*30) % 360, 80, 90, 85); 
        this.p.pop();

        this.initializeMonomers();
    }

    initializeMonomers() {
        this.monomers = [];
        this.isFixed = new Array(this.numMonomers).fill(null).map(() => ({ fixed: false, position: null }));
        const startX = this.p.random(this.p.width * 0.2, this.p.width * 0.8);
        const startY = this.p.random(this.p.height * 0.2, this.p.height * 0.8);
        
        for (let i = 0; i < this.numMonomers; i++) {
            let pos;
            if (i === 0) {
                pos = this.p.createVector(startX, startY);
            } else {
                pos = p5.Vector.add(this.monomers[i-1].pos, p5.Vector.random2D().mult(params.bondLength * 0.5));
            }
            this.monomers.push({pos: pos.copy(), prev_pos: pos.copy()});
        }
        for(let k=0; k<10; k++) this.applyInternalConstraints(); // Initial relaxation
    }

    setFixed(monomerIndex, position) {
        if (monomerIndex >= 0 && monomerIndex < this.numMonomers) {
            this.isFixed[monomerIndex] = { fixed: true, position: position.copy() };
            this.monomers[monomerIndex].pos.set(position);
            this.monomers[monomerIndex].prev_pos.set(position);
        }
    }

    clearFixed(monomerIndex) {
        if (monomerIndex >= 0 && monomerIndex < this.numMonomers) {
            this.isFixed[monomerIndex] = { fixed: false, position: null };
        }
    }
    
    update() { // No deltaTime needed here as it's handled by Verlet integration's implicit timestep
        for (let i = 0; i < this.numMonomers; i++) {
            if (this.isFixed[i].fixed) {
                this.monomers[i].pos.set(this.isFixed[i].position);
                this.monomers[i].prev_pos.set(this.isFixed[i].position);
                continue;
            }

            let m = this.monomers[i];
            let velocity = p5.Vector.sub(m.pos, m.prev_pos);
            let nextPos = p5.Vector.add(m.pos, velocity);
            nextPos.add(p5.Vector.random2D().mult(0.1)); 

            if (i === 0 || i === this.numMonomers - 1) {
                nextPos.add(p5.Vector.random2D().mult(params.endRandomStep));
                if (this.numMonomers > 1) {
                    let targetIndex = (i === 0) ? 
                        Math.min(this.numMonomers - 1, Math.floor(this.numMonomers / 3)) : 
                        Math.max(0, Math.floor(2 * this.numMonomers / 3));
                    if (this.monomers[targetIndex]) { // Ensure target exists
                        let pullDir = p5.Vector.sub(this.monomers[targetIndex].pos, m.pos);
                        nextPos.add(pullDir.normalize().mult(params.endPullStrength * params.bondLength));
                    }
                }
            }
            
            m.prev_pos.set(m.pos);
            m.pos.set(nextPos);
        }

        for (let k = 0; k < params.relaxationIterations; k++) {
            this.applyInternalConstraints();
            this.applyExtruderForces(); // Apply extruder forces
            this.applySelfRepulsion();
            this.applyObstacleRepulsion();
            this.confineToCanvas(); 
        }
    }

    applyInternalConstraints() { 
        for (let i = 0; i < this.numMonomers - 1; i++) {
            let m1 = this.monomers[i];
            let m2 = this.monomers[i+1];
            let delta = p5.Vector.sub(m2.pos, m1.pos);
            let dist = delta.mag();
            if (dist === 0) { delta = p5.Vector.random2D().mult(0.01); dist = 0.01; }
            let diff = (params.bondLength - dist) / dist;

            let m1Fixed = this.isFixed[i].fixed;
            let m2Fixed = this.isFixed[i+1].fixed;

            if (!m1Fixed && !m2Fixed) {
                m1.pos.sub(delta.copy().mult(diff * 0.5));
                m2.pos.add(delta.copy().mult(diff * 0.5));
            } else if (!m1Fixed) {
                m1.pos.sub(delta.copy().mult(diff));
            } else if (!m2Fixed) {
                m2.pos.add(delta.copy().mult(diff));
            }
        }
    }

    applyExtruderForces() {
        if (params.extruderForceStrength <= 0) return;

        loopExtruders.forEach(extruder => {
            if (extruder.chrId !== this.id || !extruder.isActive || extruder.currentMonomerIdx1 === -1 || extruder.currentMonomerIdx2 === -1 || extruder.currentMonomerIdx1 === extruder.currentMonomerIdx2) {
                return;
            }
    
            const idx1 = extruder.currentMonomerIdx1;
            const idx2 = extruder.currentMonomerIdx2;
    
            if (idx1 < 0 || idx1 >= this.numMonomers || idx2 < 0 || idx2 >= this.numMonomers) return;

            let m1 = this.monomers[idx1];
            let m2 = this.monomers[idx2];
    
            let delta = p5.Vector.sub(m2.pos, m1.pos);
            let dist = delta.mag();
            let targetDist = params.bondLength * 0.2; // Loop anchors pulled very close
    
            if (dist > targetDist) { 
                let diffRatio = (targetDist - dist) / dist; 
                
                let m1Fixed = this.isFixed[idx1].fixed;
                let m2Fixed = this.isFixed[idx2].fixed;

                let move1 = delta.copy().mult(diffRatio * params.extruderForceStrength * 0.5);
                let move2 = delta.copy().mult(diffRatio * params.extruderForceStrength * 0.5);

                if (!m1Fixed && !m2Fixed) {
                    m1.pos.sub(move1);
                    m2.pos.add(move2);
                } else if (!m1Fixed) {
                    m1.pos.sub(move1.mult(2));
                } else if (!m2Fixed) {
                    m2.pos.add(move2.mult(2));
                }
            }
        });
    }

    applySelfRepulsion() { 
        if (params.monomerSelfRepulsionStrength <= 0) return;
        const monomerRadius = params.bondLength * 0.3; 
        for (let i = 0; i < this.numMonomers; i++) {
            if (this.isFixed[i].fixed) continue;
            for (let j = i + 2; j < this.numMonomers; j++) { 
                if (this.isFixed[j].fixed) continue;

                let m1 = this.monomers[i];
                let m2 = this.monomers[j];
                let delta = p5.Vector.sub(m2.pos, m1.pos);
                let dist = delta.mag();
                let minDist = monomerRadius * 2;

                if (dist < minDist && dist > 0) {
                    let overlap = minDist - dist;
                    let repulsion = delta.copy().normalize().mult(-overlap * params.monomerSelfRepulsionStrength * 0.5);
                    m1.pos.add(repulsion);
                    m2.pos.sub(repulsion);
                }
            }
        }
    }

    applyObstacleRepulsion() {
        if (params.obstacleRepulsionStrength <= 0 || circularObstacles.length === 0) return;
        const monomerEffectiveRadius = params.bondLength * 0.2; 

        for (let i = 0; i < this.numMonomers; i++) {
            if (this.isFixed[i].fixed) continue;
            let m = this.monomers[i];

            circularObstacles.forEach(obs => {
                let delta = p5.Vector.sub(m.pos, obs.center);
                let dist = delta.mag();
                let minDist = obs.radius + monomerEffectiveRadius;

                if (dist < minDist && dist > 0) {
                    let overlap = minDist - dist;
                    let repulsionForce = delta.normalize().mult(overlap * params.obstacleRepulsionStrength);
                    m.pos.add(repulsionForce);
                }
            });
        }
    }
    
    confineToCanvas() {
        const margin = params.bondLength * 0.5;
        this.monomers.forEach((m, idx) => {
            if (this.isFixed[idx].fixed) return;
            m.pos.x = this.p.constrain(m.pos.x, margin, this.p.width - margin);
            m.pos.y = this.p.constrain(m.pos.y, margin, this.p.height - margin);
        });
    }

    display(p_instance) {
        p_instance.stroke(this.color);
        p_instance.strokeWeight(p_instance.map(params.bondLength, 1, 20, 1.5, 4)); 
        p_instance.noFill();
        p_instance.beginShape();
        this.monomers.forEach(m => p_instance.vertex(m.pos.x, m.pos.y));
        p_instance.endShape();

        this.monomers.forEach((m, idx) => {
            if (this.isFixed[idx].fixed) { // User-set spatial anchors
                 p_instance.fill(255, 255, 0); 
                 p_instance.noStroke();
                 p_instance.ellipse(m.pos.x, m.pos.y, 6, 6);
            }
        });
        // Visualize loop extruder anchors
        loopExtruders.forEach(extruder => {
            if (extruder.chrId === this.id && extruder.isActive) {
                [extruder.currentMonomerIdx1, extruder.currentMonomerIdx2].forEach(monIdx => {
                    if (monIdx !== -1 && monIdx < this.numMonomers && !this.isFixed[monIdx].fixed) { // Don't overdraw user anchors
                        const monomerPos = this.monomers[monIdx].pos;
                        p_instance.fill(255, 0, 255, 150); // Magenta for extruder anchors
                        p_instance.noStroke();
                        p_instance.ellipse(monomerPos.x, monomerPos.y, 5, 5);
                    }
                });
            }
        });
    }
}

function runSimulationStep(deltaTime) { // deltaTime in seconds
    loopExtruders.forEach(extruder => {
        const chr = chromosomes.find(c => c.id === extruder.chrId);
        const numMonomersOnChr = chr ? chr.numMonomers : 0;
        extruder.update(deltaTime, numMonomersOnChr);
    });
    chromosomes.forEach(chr => chr.update());
    updateHiCData();
    calculatePsData();
}

function updateHiCData() {
    if (!hicData || params.numMonomers <= 0) return; // Guard

    for (let i = 0; i < params.numMonomers; i++) {
        for (let j = 0; j < params.numMonomers; j++) {
             if (hicData[i] && typeof hicData[i][j] === 'number') { // Check if hicData[i] and hicData[i][j] are valid
                hicData[i][j] *= params.hicDecayRate;
            } else {
                // This case should ideally not be reached if initialization is correct
                if (!hicData[i]) hicData[i] = Array(params.numMonomers).fill(0);
                else hicData[i][j] = 0;
            }
        }
    }
    maxHiCValue *= params.hicDecayRate;

    chromosomes.forEach(chr => {
        if (chr.numMonomers !== params.numMonomers) return; // Safety check if a chromosome has a different monomer count
        for (let i = 0; i < chr.numMonomers; i++) {
            for (let j = i; j < chr.numMonomers; j++) { 
                if (!chr.monomers[i] || !chr.monomers[j]) continue; // Monomers might not be initialized yet
                let m1 = chr.monomers[i].pos;
                let m2 = chr.monomers[j].pos;
                let dist = p5.Vector.dist(m1, m2);
                
                if (dist < params.contactRadius) {
                    let increment = 1.0; 
                    hicData[i][j] += increment;
                    if (i !== j) hicData[j][i] += increment; 
                    if (hicData[i][j] > maxHiCValue) maxHiCValue = hicData[i][j];
                }
            }
        }
    });
    if (maxHiCValue < 1e-6) maxHiCValue = 1e-6; 
}

function calculatePsData() {
    psData = [];
    if (!hicData || params.numMonomers <= 0) return; // Guard

    for (let s = 1; s < params.numMonomers; s++) { 
        let contactSum = 0;
        let numPairs = 0;
        for (let i = 0; i < params.numMonomers - s; i++) {
            if (hicData[i] && typeof hicData[i][i+s] === 'number') { // Check validity
                 contactSum += hicData[i][i+s];
            }
            numPairs++;
        }
        if (numPairs > 0) {
            let P_s = contactSum / numPairs; 
            psData.push({ s: s, P: P_s > 0 ? P_s : 1e-9 }); 
        } else {
            psData.push({ s: s, P: 1e-9 });
        }
    }
}

function resetSimulation() {
    loadParamsFromUI(); 
    
    chromosomes = [];
    if (!pDynamics || !pDynamics.width || !pDynamics.height) { 
        console.warn("Dynamics canvas not ready for reset. Using default sizes for initialization if necessary.");
    }

    for (let i = 0; i < params.numChromosomes; i++) {
        chromosomes.push(new Chromosome(i, params.numMonomers, pDynamics)); 
    }

    if (params.numMonomers > 0) {
        hicData = Array(params.numMonomers).fill(null).map(() => Array(params.numMonomers).fill(0));
    } else {
        hicData = []; // Handle numMonomers = 0 case
    }
    maxHiCValue = 1e-6;
    psData = [];
    
    // Clear and reapply anchors (they are persistent across resets if still valid)
    reapplyMonomerAnchors(); 
    updateMonomerAnchorListUI();
    
    // Obstacles are also persistent, re-check positions if canvas size changes
    // For now, absolute positions are kept.
    updateObstacleListUI();

    // Loop extruders are typically reset
    loopExtruders = [];
    updateLoopExtruderListUI();

    lastUpdateTime = Date.now();
    if (!params.isPlaying) { 
        params.isPlaying = true; 
        updatePlayPauseButton();
    }
    if (pDynamics) pDynamics.loop();
    if (pHiC) pHiC.redraw();
    if (pPsPlot) pPsPlot.redraw();
}

function setupApplicationColors(p) {
    if (!p || typeof p.color !== 'function') {
        console.error("Valid p5 instance with color function needed for color setup.");
        // Fallback colors if p5 instance isn't ready (should not happen with current flow)
        const style = document.documentElement.style;
        cachedColors.bgColor = style.getPropertyValue('--bg-color').trim() || '#1c1e22';
        cachedColors.panelBgColor = style.getPropertyValue('--panel-bg-color').trim() || '#282c34';
        cachedColors.textColor = style.getPropertyValue('--text-color').trim() || '#abb2bf';
        cachedColors.primaryAccentColor = style.getPropertyValue('--primary-accent-color').trim() || '#61afef';
        cachedColors.secondaryAccentColor = style.getPropertyValue('--secondary-accent-color').trim() || '#e06c75';
        cachedColors.tertiaryAccentColor = style.getPropertyValue('--tertiary-accent-color').trim() || '#98c379';
        return;
    }
    const rootStyle = getComputedStyle(document.documentElement);
    cachedColors.bgColor = p.color(rootStyle.getPropertyValue('--bg-color').trim());
    cachedColors.panelBgColor = p.color(rootStyle.getPropertyValue('--panel-bg-color').trim());
    cachedColors.textColor = p.color(rootStyle.getPropertyValue('--text-color').trim());
    cachedColors.primaryAccentColor = p.color(rootStyle.getPropertyValue('--primary-accent-color').trim());
    cachedColors.secondaryAccentColor = p.color(rootStyle.getPropertyValue('--secondary-accent-color').trim());
    cachedColors.tertiaryAccentColor = p.color(rootStyle.getPropertyValue('--tertiary-accent-color').trim());

    // Ensure CSS var --tertiary-accent-color-rgb is set if needed by .instruction background
    // Example: --tertiary-accent-color-rgb: 152, 195, 121;
    // For now, the CSS has a fallback.
}

function setupControls() {
    numChromosomesInput = document.getElementById('numChromosomes');
    numMonomersInput = document.getElementById('numMonomers');
    simulationSpeedInput = document.getElementById('simulationSpeed');
    bondLengthInput = document.getElementById('bondLength');
    endPullStrengthInput = document.getElementById('endPullStrength');
    endRandomStepInput = document.getElementById('endRandomStep');
    contactRadiusInput = document.getElementById('contactRadius');
    hicDecayRateInput = document.getElementById('hicDecayRate');
    monomerSelfRepulsionInput = document.getElementById('monomerSelfRepulsion');
    obstacleRepulsionInput = document.getElementById('obstacleRepulsion');
    extruderForceStrengthInput = document.getElementById('extruderForceStrength');
    
    endPullStrengthValueSpan = document.getElementById('endPullStrengthValue');
    endRandomStepValueSpan = document.getElementById('endRandomStepValue');
    hicDecayRateValueSpan = document.getElementById('hicDecayRateValue');
    monomerSelfRepulsionValueSpan = document.getElementById('monomerSelfRepulsionValue');
    obstacleRepulsionValueSpan = document.getElementById('obstacleRepulsionValue');
    extruderForceStrengthValueSpan = document.getElementById('extruderForceStrengthValue');

    playPauseButton = document.getElementById('playPauseButton');
    resetButton = document.getElementById('resetButton');
    
    anchorChrIdxInput = document.getElementById('anchorChrIdx');
    anchorMonomerIdxInput = document.getElementById('anchorMonomerIdx');
    anchorTypeSelect = document.getElementById('anchorType');
    addAnchorButton = document.getElementById('addAnchorButton');
    clearAnchorsButton = document.getElementById('clearAnchorsButton');
    anchorListUl = document.querySelector('#anchorList ul');

    obstacleRadiusInput = document.getElementById('obstacleRadius');
    addObstacleButton = document.getElementById('addObstacleButton');
    clearObstaclesButton = document.getElementById('clearObstaclesButton');
    obstacleListUl = document.querySelector('#obstacleList ul');

    extruderChrIdxInput = document.getElementById('extruderChrIdx');
    extruderTypeSelect = document.getElementById('extruderType');
    extruderStaticIndicesContainer = document.getElementById('extruderStaticIndicesContainer');
    extruderStaticMonomer1Input = document.getElementById('extruderStaticMonomer1');
    extruderStaticMonomer2Input = document.getElementById('extruderStaticMonomer2');
    extruderSpeedInput = document.getElementById('extruderSpeed');
    extruderMaxSizeInput = document.getElementById('extruderMaxSize');
    addExtruderButton = document.getElementById('addExtruderButton');
    clearExtrudersButton = document.getElementById('clearExtrudersButton');
    extruderListUl = document.querySelector('#extruderList ul');

    // Event listeners
    numChromosomesInput.addEventListener('change', () => { loadParamsFromUI(); resetSimulation(); });
    numMonomersInput.addEventListener('change', () => { loadParamsFromUI(); resetSimulation(); });
    
    simulationSpeedInput.addEventListener('change', loadParamsFromUI);
    bondLengthInput.addEventListener('change', loadParamsFromUI); 
    contactRadiusInput.addEventListener('change', loadParamsFromUI);

    const addSliderListener = (inputEl, spanEl, paramName, fixedDigits) => {
        inputEl.addEventListener('input', () => { 
            loadParamsFromUI(); 
            spanEl.textContent = params[paramName].toFixed(fixedDigits); 
        });
    };

    addSliderListener(endPullStrengthInput, endPullStrengthValueSpan, 'endPullStrength', 3);
    addSliderListener(endRandomStepInput, endRandomStepValueSpan, 'endRandomStep', 1);
    addSliderListener(hicDecayRateInput, hicDecayRateValueSpan, 'hicDecayRate', 3);
    addSliderListener(monomerSelfRepulsionInput, monomerSelfRepulsionValueSpan, 'monomerSelfRepulsionStrength', 2);
    addSliderListener(obstacleRepulsionInput, obstacleRepulsionValueSpan, 'obstacleRepulsionStrength', 2);
    addSliderListener(extruderForceStrengthInput, extruderForceStrengthValueSpan, 'extruderForceStrength', 2);

    playPauseButton.addEventListener('click', togglePlayPause);
    resetButton.addEventListener('click', resetSimulation);

    addAnchorButton.addEventListener('click', () => {
        let chrIdx = parseInt(anchorChrIdxInput.value);
        let monomerIdx = parseInt(anchorMonomerIdxInput.value);
        let type = anchorTypeSelect.value;
        
        if (type === 'fix_to_center') {
            if (!pDynamics || !pDynamics.width) { alert("Dynamics canvas not ready."); return; }
            let position = pDynamics.createVector(pDynamics.width / 2, pDynamics.height / 2);
            addMonomerAnchor(chrIdx, monomerIdx, position, type);
        } else if (type === 'fix_to_click') {
            alert("Click on the Chromosome Dynamics panel to set the anchor position.");
            // Position will be set by handleDynamicsCanvasClick
        }
    });
    clearAnchorsButton.addEventListener('click', clearMonomerAnchors);
    
    addObstacleButton.addEventListener('click', () => {
        if (!pendingObstacleCenter) {
            alert("Please click on the Dynamics Canvas first to set the obstacle's center.");
            return;
        }
        let radius = parseFloat(obstacleRadiusInput.value);
        if (radius <= 0) { alert("Obstacle radius must be positive."); return; }
        addCircularObstacle(pendingObstacleCenter, radius);
        pendingObstacleCenter = null; 
        if(pDynamics) pDynamics.redraw();
    });
    clearObstaclesButton.addEventListener('click', clearCircularObstacles);

    extruderTypeSelect.addEventListener('change', (e) => {
        if (e.target.value === 'static') {
            extruderStaticIndicesContainer.classList.remove('hidden');
            extruderSpeedInput.disabled = true;
            extruderSpeedInput.value = 0; // Static loops have 0 speed
        } else {
            extruderStaticIndicesContainer.classList.add('hidden');
            extruderSpeedInput.disabled = false;
        }
        loadParamsFromUI(); // Update speed param if changed
    });
    addExtruderButton.addEventListener('click', addLoopExtruder);
    clearExtrudersButton.addEventListener('click', clearLoopExtruders);
    
    // Initialize span values and parameter validation
    loadParamsFromUI(); 
    endPullStrengthValueSpan.textContent = params.endPullStrength.toFixed(3);
    endRandomStepValueSpan.textContent = params.endRandomStep.toFixed(1);
    hicDecayRateValueSpan.textContent = params.hicDecayRate.toFixed(3);
    monomerSelfRepulsionValueSpan.textContent = params.monomerSelfRepulsionStrength.toFixed(2);
    obstacleRepulsionValueSpan.textContent = params.obstacleRepulsionStrength.toFixed(2);
    extruderForceStrengthValueSpan.textContent = params.extruderForceStrength.toFixed(2);
    
    extruderTypeSelect.dispatchEvent(new Event('change')); // Trigger initial state for static indices visibility
}

function loadParamsFromUI() {
    params.numChromosomes = parseInt(numChromosomesInput.value) || 1;
    params.numChromosomes = Math.max(1, Math.min(params.numChromosomes, parseInt(numChromosomesInput.max) || 20));
    numChromosomesInput.value = params.numChromosomes;

    params.numMonomers = parseInt(numMonomersInput.value) || 10;
    params.numMonomers = Math.max(10, Math.min(params.numMonomers, parseInt(numMonomersInput.max) || 200));
    numMonomersInput.value = params.numMonomers;

    params.simulationSpeed = parseInt(simulationSpeedInput.value) || 30;
    params.simulationSpeed = Math.max(1, params.simulationSpeed);
    simulationSpeedInput.value = params.simulationSpeed;

    params.bondLength = parseFloat(bondLengthInput.value) || 8;
    params.endPullStrength = parseFloat(endPullStrengthInput.value);
    params.endRandomStep = parseFloat(endRandomStepInput.value);
    params.contactRadius = parseFloat(contactRadiusInput.value) || 10;
    params.hicDecayRate = parseFloat(hicDecayRateInput.value);
    params.monomerSelfRepulsionStrength = parseFloat(monomerSelfRepulsionInput.value);
    params.obstacleRepulsionStrength = parseFloat(obstacleRepulsionInput.value);
    params.extruderForceStrength = parseFloat(extruderForceStrengthInput.value);

    // Update max values for index inputs
    anchorChrIdxInput.max = Math.max(0, params.numChromosomes - 1);
    anchorMonomerIdxInput.max = Math.max(0, params.numMonomers - 1);
    extruderChrIdxInput.max = Math.max(0, params.numChromosomes - 1);
    extruderStaticMonomer1Input.max = Math.max(0, params.numMonomers - 1);
    extruderStaticMonomer2Input.max = Math.max(0, params.numMonomers - 1);
    extruderMaxSizeInput.max = Math.max(1, params.numMonomers -1);
}

function togglePlayPause() {
    params.isPlaying = !params.isPlaying;
    if (params.isPlaying) {
        lastUpdateTime = Date.now(); 
        if(pDynamics) pDynamics.loop(); 
    } else {
        if(pDynamics) pDynamics.noLoop(); 
    }
    updatePlayPauseButton();
}

function updatePlayPauseButton() {
     playPauseButton.innerHTML = params.isPlaying ? '<span class="icon">❚❚</span> Pause' : '<span class="icon">▶</span> Play';
}

// Monomer Anchor Functions
function addMonomerAnchor(chrIndex, monomerIndex, position, type) {
    if (chrIndex >= params.numChromosomes || chrIndex < 0 || monomerIndex >= params.numMonomers || monomerIndex < 0) {
        alert("Invalid chromosome or monomer index for current settings.");
        return;
    }
    if (!position && type === 'fix_to_click') {
        alert("Click position for anchor not yet defined for fix_to_click.");
        return;
    }
    if (monomerAnchors.some(a => a.chrIndex === chrIndex && a.monomerIndex === monomerIndex)) {
        alert(`Monomer ${monomerIndex} of Chromosome ${chrIndex} is already anchored. Clear existing anchor first.`);
        return;
    }

    const newAnchor = { chrIndex, monomerIndex, position, type, id: `${chrIndex}-${monomerIndex}` };
    monomerAnchors.push(newAnchor);
    
    if (chromosomes[chrIndex]) {
      chromosomes[chrIndex].setFixed(monomerIndex, position);
    }
    updateMonomerAnchorListUI();
}

function clearMonomerAnchors() {
    monomerAnchors.forEach(anchor => {
        if (chromosomes[anchor.chrIndex]) {
            chromosomes[anchor.chrIndex].clearFixed(anchor.monomerIndex);
        }
    });
    monomerAnchors = [];
    updateMonomerAnchorListUI();
}

function reapplyMonomerAnchors() {
    // Existing anchors are stored in monomerAnchors array.
    // Filter out invalid ones and re-apply valid ones.
    const validAnchors = [];
    monomerAnchors.forEach(anchor => {
        if (anchor.chrIndex < params.numChromosomes && anchor.monomerIndex < params.numMonomers) {
            if(anchor.type === 'fix_to_center' && pDynamics && pDynamics.width) { 
                anchor.position = pDynamics.createVector(pDynamics.width / 2, pDynamics.height / 2);
            }
            // For 'fix_to_click', position should already be set. If not, it might be problematic or needs re-click.
            // For simplicity, assume 'click_fixed' positions persist.

            if (chromosomes[anchor.chrIndex] && anchor.position) { 
                chromosomes[anchor.chrIndex].setFixed(anchor.monomerIndex, anchor.position);
                validAnchors.push(anchor);
            }
        }
    });
    monomerAnchors = validAnchors;
}

function updateMonomerAnchorListUI() {
    anchorListUl.innerHTML = ''; 
    monomerAnchors.forEach(anchor => {
        const li = document.createElement('li');
        let posStr = anchor.position ? `(${anchor.position.x.toFixed(0)}, ${anchor.position.y.toFixed(0)})` : '(pos pending)';
        li.textContent = `Chr ${anchor.chrIndex}, Mon ${anchor.monomerIndex} @ ${posStr} (${anchor.type})`;
        
        const removeBtn = document.createElement('button');
        removeBtn.innerHTML = '<span class="icon">🗑️</span>';
        removeBtn.onclick = () => {
            if (chromosomes[anchor.chrIndex]) {
                chromosomes[anchor.chrIndex].clearFixed(anchor.monomerIndex);
            }
            monomerAnchors = monomerAnchors.filter(a => a.id !== anchor.id);
            updateMonomerAnchorListUI();
        };
        li.appendChild(removeBtn);
        anchorListUl.appendChild(li);
    });
}

// Circular Obstacle Functions
function addCircularObstacle(center, radius) {
    if (!center || radius <=0) {
        alert("Invalid obstacle parameters.");
        return;
    }
    circularObstacles.push({ center: center.copy(), radius: radius, id: nextObstacleId++ });
    updateObstacleListUI();
}

function clearCircularObstacles() {
    circularObstacles = [];
    updateObstacleListUI();
}

function updateObstacleListUI() {
    obstacleListUl.innerHTML = '';
    circularObstacles.forEach(obs => {
        const li = document.createElement('li');
        li.textContent = `ID ${obs.id}: Center (${obs.center.x.toFixed(0)}, ${obs.center.y.toFixed(0)}), R: ${obs.radius.toFixed(0)}`;
        
        const removeBtn = document.createElement('button');
        removeBtn.innerHTML = '<span class="icon">🗑️</span>';
        removeBtn.onclick = () => {
            circularObstacles = circularObstacles.filter(o => o.id !== obs.id);
            updateObstacleListUI();
        };
        li.appendChild(removeBtn);
        obstacleListUl.appendChild(li);
    });
     if(pDynamics) pDynamics.redraw(); 
}

// Loop Extruder Class and Functions
class LoopExtruder {
    constructor(config) {
        this.id = config.id;
        this.chrId = config.chrId;
        this.type = config.type; // 'static', 'one_sided_start', 'one_sided_end', 'both_sided_center'
        this.extrusionSpeed = config.extrusionSpeed; // monomers per second
        this.maxLoopSize = config.maxLoopSize; // genomic distance
        
        this.currentMonomerIdx1 = -1;
        this.currentMonomerIdx2 = -1;
        this.currentGenomicSize = 0;
        this.isActive = true;
        this.timeSinceLastExtrusionStep = 0; // accumulates deltaTime

        this.initializeIndices(config.m1_initial, config.m2_initial, config.numMonomersOnChr);
    }

    initializeIndices(m1_initial, m2_initial, numMonomersOnChr) {
        if (numMonomersOnChr <= 0) { this.isActive = false; return; }
        const M = numMonomersOnChr;

        switch (this.type) {
            case 'static':
                this.currentMonomerIdx1 = Math.max(0, Math.min(m1_initial, M - 1));
                this.currentMonomerIdx2 = Math.max(0, Math.min(m2_initial, M - 1));
                if (this.currentMonomerIdx1 === this.currentMonomerIdx2) this.isActive = false; // Invalid static loop
                break;
            case 'one_sided_start':
                this.currentMonomerIdx1 = 0;
                this.currentMonomerIdx2 = 0;
                break;
            case 'one_sided_end':
                this.currentMonomerIdx1 = M - 1;
                this.currentMonomerIdx2 = M - 1;
                break;
            case 'both_sided_center':
                const center = Math.floor(M / 2);
                this.currentMonomerIdx1 = center;
                this.currentMonomerIdx2 = center;
                break;
            default:
                this.isActive = false;
        }
        this.currentGenomicSize = Math.abs(this.currentMonomerIdx2 - this.currentMonomerIdx1);
    }

    update(deltaTimeSeconds, numMonomersOnChr) {
        if (!this.isActive || this.type === 'static' || this.extrusionSpeed <= 0 || numMonomersOnChr <= 0) {
            return;
        }
        const M = numMonomersOnChr;
        this.timeSinceLastExtrusionStep += deltaTimeSeconds;
        const timePerMonomerExtrusion = 1.0 / this.extrusionSpeed;

        while (this.timeSinceLastExtrusionStep >= timePerMonomerExtrusion) {
            if (this.currentGenomicSize >= this.maxLoopSize) break;

            let extrudedOneStep = false;
            switch (this.type) {
                case 'one_sided_start':
                    if (this.currentMonomerIdx2 < M - 1) {
                        this.currentMonomerIdx2++;
                        extrudedOneStep = true;
                    }
                    break;
                case 'one_sided_end':
                    if (this.currentMonomerIdx1 > 0) {
                        this.currentMonomerIdx1--;
                        extrudedOneStep = true;
                    }
                    break;
                case 'both_sided_center':
                    // Prioritize symmetric extrusion if possible
                    let canExtrudeLeft = this.currentMonomerIdx1 > 0;
                    let canExtrudeRight = this.currentMonomerIdx2 < M - 1;
                    if (this.currentGenomicSize + 2 <= this.maxLoopSize && canExtrudeLeft && canExtrudeRight) {
                        this.currentMonomerIdx1--;
                        this.currentMonomerIdx2++;
                        extrudedOneStep = true;
                    } else if (this.currentGenomicSize + 1 <= this.maxLoopSize) { // Try asymmetric if symmetric not possible or exceeds max size
                        if (canExtrudeRight && (!canExtrudeLeft || Math.random() < 0.5)) { // Randomly pick if both options are singly available
                             this.currentMonomerIdx2++; extrudedOneStep = true;
                        } else if (canExtrudeLeft) {
                             this.currentMonomerIdx1--; extrudedOneStep = true;
                        }
                    }
                    break;
            }

            if (extrudedOneStep) {
                this.currentGenomicSize = Math.abs(this.currentMonomerIdx2 - this.currentMonomerIdx1);
                this.timeSinceLastExtrusionStep -= timePerMonomerExtrusion;
            } else {
                break; // Cannot extrude further
            }
        }
    }
}

function addLoopExtruder() {
    const chrId = parseInt(extruderChrIdxInput.value);
    const type = extruderTypeSelect.value;
    const speed = parseFloat(extruderSpeedInput.value);
    const maxSize = parseInt(extruderMaxSizeInput.value);
    
    if (chrId < 0 || chrId >= params.numChromosomes) { alert("Invalid chromosome index for extruder."); return; }
    if (speed < 0) { alert("Extrusion speed cannot be negative."); return; }
    if (maxSize <= 0 && type !== 'static') { alert("Max loop size must be positive for dynamic extruders."); return; }

    const M = chromosomes[chrId] ? chromosomes[chrId].numMonomers : params.numMonomers;
    let m1_initial = 0, m2_initial = 0;

    if (type === 'static') {
        m1_initial = parseInt(extruderStaticMonomer1Input.value);
        m2_initial = parseInt(extruderStaticMonomer2Input.value);
        if (m1_initial < 0 || m1_initial >= M || m2_initial < 0 || m2_initial >= M || m1_initial === m2_initial) {
            alert("Invalid monomer indices for static loop. Ensure they are different and within bounds [0, M-1].");
            return;
        }
        if (speed !== 0) {
             console.warn("Static extruder type selected, speed will be treated as 0.");
        }
    }

    const config = {
        id: nextExtruderId++,
        chrId: chrId,
        type: type,
        extrusionSpeed: (type === 'static') ? 0 : speed,
        maxLoopSize: maxSize,
        m1_initial: m1_initial,
        m2_initial: m2_initial,
        numMonomersOnChr: M
    };
    
    const newExtruder = new LoopExtruder(config);
    if (!newExtruder.isActive) {
        alert("Failed to initialize extruder, possibly due to invalid configuration or chromosome state.");
        return;
    }
    loopExtruders.push(newExtruder);
    updateLoopExtruderListUI();
}

function clearLoopExtruders() {
    loopExtruders = [];
    updateLoopExtruderListUI();
}

function updateLoopExtruderListUI() {
    extruderListUl.innerHTML = '';
    loopExtruders.forEach(extruder => {
        const li = document.createElement('li');
        let typeStr = extruder.type.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
        let details = `Chr ${extruder.chrId}, Type: ${typeStr}`;
        if (extruder.type === 'static') {
            details += `, Anchors: ${extruder.currentMonomerIdx1}-${extruder.currentMonomerIdx2}`;
        } else {
            details += `, Speed: ${extruder.extrusionSpeed.toFixed(1)} M/s, Max Size: ${extruder.maxLoopSize}`;
            details += `, Current Anchors: ${extruder.currentMonomerIdx1}-${extruder.currentMonomerIdx2} (Size: ${extruder.currentGenomicSize})`;
        }

        li.textContent = details;
        
        const removeBtn = document.createElement('button');
        removeBtn.innerHTML = '<span class="icon">🗑️</span>';
        removeBtn.onclick = () => {
            loopExtruders = loopExtruders.filter(e => e.id !== extruder.id);
            updateLoopExtruderListUI();
        };
        li.appendChild(removeBtn);
        extruderListUl.appendChild(li);
    });
}


// Initialize p5 instances and then global setup
document.addEventListener('DOMContentLoaded', () => {
    new p5(dynamicsSketch);
    new p5(hiCSketch);
    new p5(psPlotSketch);

    setTimeout(() => {
        if (pDynamics) { // pDynamics is primary for color parsing
            setupApplicationColors(pDynamics);
        } else { // Fallback if pDynamics isn't ready (shouldn't happen)
            console.warn("pDynamics not ready for initial color setup, using fallback.");
            setupApplicationColors(null); 
        }
        
        setupControls(); // Depends on DOM elements
        loadParamsFromUI(); // Sets initial params from UI, validates them
        resetSimulation(); // Initializes simulation state based on params

        // Redraw all canvases once colors and initial data are ready
        if (pDynamics) {
            if (!params.isPlaying) pDynamics.noLoop(); else pDynamics.loop();
            pDynamics.redraw();
        }
        if (pHiC) pHiC.redraw();
        if (pPsPlot) pPsPlot.redraw();

    }, 150); // Increased delay slightly for safety
});
PartSys.prototype.ForceField = function () {
    this.partCount = partsPerAxis * ((partsPerAxis - 1) * 2);
    this.maxDensity = 1;
    this.minDensity = 1000000;
    var doit = 1;
    this.s0 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s0dot = new Float32Array(this.partCount * PART_MAXVAR);
    this.s1 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s1dot = new Float32Array(this.partCount * PART_MAXVAR);
    this.buffer = new Float32Array(this.partCount * PART_MAXVAR);
    this.sM = new Float32Array(this.partCount * PART_MAXVAR);
    this.sMdot = new Float32Array(this.partCount * PART_MAXVAR);
    this.s2dot = new Float32Array(this.partCount * PART_MAXVAR);
    this.s2 = new Float32Array(this.partCount * PART_MAXVAR);
    this.first = 0;

    for (var i = 0; i < this.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        var xcyc = roundRand3D();
        if (doit == 1) {
            console.log('xc,yc= ' + xcyc[0] + ', ' + xcyc[1]);
            doit = 0;
        }
        xcyc = roundRand3D();
        this.s0[pOff + PART_XVEL] = 0;//INIT_VEL * (this.s0[pOff + PART_YPOS]/2);//(0.4 + 0.2 * xcyc[0]);
        this.s0[pOff + PART_YVEL] = 0;//INIT_VEL * (this.s0[pOff + PART_XPOS]/2);//(0.4 + 0.2 * xcyc[1]);
        this.s0[pOff + PART_ZVEL] = 0;//INIT_VEL * (0.4 + 0.2 * xcyc[2]);
        this.s0[pOff + PART_X_FTOT] = 0;//15 * (.9 - this.s0[pOff + PART_YPOS]);
        this.s0[pOff + PART_Y_FTOT] = 0;//-15 * (.9 - this.s0[pOff + PART_XPOS]);
        this.s0[pOff + PART_Z_FTOT] = 0.0;
        this.s0[pOff + PART_R] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_G] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_B] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_MASS] = 1;//+ 0.2 * Math.random();
        this.s0[pOff + PART_DIAM] = 15;// + 10.0 * Math.random();
        this.s0[pOff + PART_RENDMODE] = Math.floor(4.0 * Math.random()); // 0,1,2 or 3.
        this.s0[pOff + PART_AGE] = 0;  // # of frame-times since creation/initialization
        this.s0[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        this.s0[pOff + PART_MASS_VEL] = 0;  // time-rate-of-change of mass.
        this.s0[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        this.s0[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
        this.s0[pOff + PART_G_VEL] = 0;  // time-rate-of-change of color:grn
        this.s0[pOff + PART_B_VEL] = 0;  // time-rate-of-change of color:blu
        this.s0[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        this.s0[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        this.s0[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
        this.s0[pOff + PART_DENSITY] = 1;
        this.s0[pOff + PART_VOX_X] = Math.max(Math.floor(this.s0[pOff + PART_XPOS] / (h * 2)), 0);
        this.s0[pOff + PART_VOX_Y] = Math.max(Math.floor(this.s0[pOff + PART_YPOS] / (h * 2)), 0);
        this.s0[pOff + PART_VOX_Z] = Math.max(Math.floor(this.s0[pOff + PART_ZPOS] / (h * 2)), 0);
    }
    let index = 0;
    let odd = 0;
    for (var i = 0; i < 1; i++) {
        for (var k = 0; k < (partsPerAxis * 2); k++) {
            for (var j = 0; j < partsPerAxis; j++) {
                let pOff = index * PART_MAXVAR;
                this.s0[pOff + PART_XPOS] = (odd === 0) ? interval * j : interval * j + half_interval;// + (epsilon * (1 * ((partsPerAxis - 1) / 2) - i));
                this.s0[pOff + PART_YPOS] = half_interval * k;// + (epsilon * (1 * ((partsPerAxis - 1) / 2) - j));
                this.s0[pOff + PART_ZPOS] = 0;//(axisSize / (partsPerAxis - 1)) * k + (epsilon * (1 * ((partsPerAxis - 1) / 2) - k));
                index++;
            }
            odd = 1 - odd;
        }
    }
    this.s1 = this.s0.slice(0);
    this.sm = this.s0.slice(0);
    this.s2 = this.s0.slice(0);

    this.cforcer = [
        function (count, s) {
        },
        function (count, s) {
        },
        function (count, s) {
            for (var i = 0; i < count; i++) {
                var pOff = i * PART_MAXVAR;
                s[pOff + PART_X_FTOT] = 0;
                s[pOff + PART_Y_FTOT] = 0;
                s[pOff + PART_Z_FTOT] = -9.832 * s[pOff + PART_MASS];
            }
        },
        function (count, s) {
        },
        function (count, s) {
        },
        function (count, s) {
        },
        function (count, s) {
            for (var i = 0; i < count; i++) {
                var pOff = i * PART_MAXVAR;
                s[pOff + PART_X_FTOT] -= (s[pOff + PART_DIAM] * s[pOff + PART_XVEL] / 2);
                s[pOff + PART_Y_FTOT] -= (s[pOff + PART_DIAM] * s[pOff + PART_YVEL] / 2);
                s[pOff + PART_Z_FTOT] -= (s[pOff + PART_DIAM] * s[pOff + PART_ZVEL] / 2);
            }
        },
        function (count, s) {
        },
        function (count, s) {
            let densities = findDensity(s, count);
            findPressureForce(s, count);
            //      findViscosity(s, count);
            let maxDensity = densities[1];
            let minDensity = densities[0];
            return [minDensity, maxDensity];
        },
        function (count, s) {
        },
        function (count, s) {
        }
    ];


    this.c0 = [
        function (i, s1, s0) {
            if (g_bounce == 0) { //--------------------------------------------------------
                if (s1[i * PART_MAXVAR + PART_XPOS] < 0.0 && s1[i * PART_MAXVAR + PART_XVEL] < 0.0			// simple velocity-reversal
                ) {		// bounce on left wall.
                    s1[i * PART_MAXVAR + PART_XVEL] = -g_resti * s1[i * PART_MAXVAR + PART_XVEL];
                } else if (s1[i * PART_MAXVAR + PART_XPOS] > axisSize && s1[i * PART_MAXVAR + PART_XVEL] > 0.0
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_XVEL] = -g_resti * s1[i * PART_MAXVAR + PART_XVEL];
                }
                if (s1[i * PART_MAXVAR + PART_YPOS] < 0.0 && s1[i * PART_MAXVAR + PART_YVEL] < 0.0
                ) {		// bounce on floor
                    s1[i * PART_MAXVAR + PART_YVEL] = 0;//-g_resti * s1[i * PART_MAXVAR + PART_YVEL];
                } else if (s1[i * PART_MAXVAR + PART_YPOS] > axisSize && s1[i * PART_MAXVAR + PART_YVEL] > 0.0
                ) {		// bounce on ceiling
                    s1[i * PART_MAXVAR + PART_YVEL] = -g_resti * s1[i * PART_MAXVAR + PART_YVEL];
                }
                if (s1[i * PART_MAXVAR + PART_ZPOS] < 0 && s1[i * PART_MAXVAR + PART_ZVEL] < 0.0			// simple velocity-reversal
                ) {		// bounce on left wall.
                    s1[i * PART_MAXVAR + PART_ZVEL] = -g_resti * s1[i * PART_MAXVAR + PART_ZVEL];
                } else if (s1[i * PART_MAXVAR + PART_ZPOS] > axisSize && s1[i * PART_MAXVAR + PART_ZVEL] > 0.0
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_ZVEL] = -g_resti * s1[i * PART_MAXVAR + PART_ZVEL];
                }
                //  -- hard limit on 'floor' keeps y position >= 0;
                if (s1[i * PART_MAXVAR + PART_ZPOS] < 0.0) s1[i * PART_MAXVAR + PART_ZPOS] = 0.0;
            }
            if (g_bounce == 1) { //---------------------------------------------------------------------------
                if (s1[i * PART_MAXVAR + PART_XPOS] < 0.0 && s1[i * PART_MAXVAR + PART_XVEL] < 0.0 // collision!  left wall...
                ) {		// bounce on left wall.
                    s1[i * PART_MAXVAR + PART_XPOS] = 0.0;
                    s1[i * PART_MAXVAR + PART_XVEL] = 0.0;
                } else if (s1[i * PART_MAXVAR + PART_XPOS] > axisSize && s1[i * PART_MAXVAR + PART_XVEL] > 0.0		// collision! right wall...
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_XPOS] = axisSize;					// 1) resolve contact: put particle at wall.
                    s1[i * PART_MAXVAR + PART_XVEL] = 0.0;
                }
                if (s1[i * PART_MAXVAR + PART_YPOS] < 0.0 && s1[i * PART_MAXVAR + PART_YVEL] < 0.0		// collision! floor
                ) {		// bounce on floor

                    s1[i * PART_MAXVAR + PART_YPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    s1[i * PART_MAXVAR + PART_YVEL] = 0;
                } else if (s1[i * PART_MAXVAR + PART_YPOS] > axisSize && s1[i * PART_MAXVAR + PART_YVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling
                    s1[i * PART_MAXVAR + PART_YPOS] = axisSize;					// 1) resolve contact: put particle at wall.
                    s1[i * PART_MAXVAR + PART_YVEL] = 0;

                }
                //	console.log('z = ' + s0[PART_ZPOS] + '  zVel = ' +s0[PART_ZVEL]);
                if (s1[i * PART_MAXVAR + PART_ZPOS] < 0.0 && s1[i * PART_MAXVAR + PART_ZVEL] < 0.0		// collision! left wall...
                ) {		// bounce on floor

                    s1[i * PART_MAXVAR + PART_ZPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    s1[i * PART_MAXVAR + PART_ZVEL] = 0.0;
                } else if (s1[i * PART_MAXVAR + PART_ZPOS] > axisSize && s1[i * PART_MAXVAR + PART_ZVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling

                    s1[i * PART_MAXVAR + PART_ZPOS] = axisSize;					// 1) resolve contact: put particle at wall.
                    s1[i * PART_MAXVAR + PART_ZVEL] = 0.0;
                }
            } else {
                console.log('?!?! unknown constraint: g_bounce==' + g_bounce);
                return;
            }
        },
        function (pOff, s1) {
        },
        function (pOff, s1) {
        }

    ]
}
PartSys.prototype.BouncyS0 = function () {
    this.partCount = 1000;
    var doit = 1;
    this.s0 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s1 = new Float32Array(this.partCount * PART_MAXVAR);
    this.sdot = new Float32Array(this.partCount * PART_MAXVAR);
    this.buffer = new Float32Array(this.partCount * PART_MAXVAR);
    this.sM = new Float32Array(this.partCount * PART_MAXVAR);
    this.s2dot = new Float32Array(this.partCount * PART_MAXVAR);
    this.s2 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s1dot2 = new Float32Array(this.partCount * PART_MAXVAR);
    this.first = 0;

    for (var i = 0; i < this.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        var xcyc = roundRand3D();
        if (doit == 1) {
            console.log('xc,yc= ' + xcyc[0] + ', ' + xcyc[1]);
            doit = 0;
        }
        this.s0[pOff + PART_XPOS] = 0.2 + 0.2 * xcyc[0];		// 0.0 <= randomRound() < 1.0
        this.s0[pOff + PART_YPOS] = 0.2 + 0.2 * xcyc[1];
        this.s0[pOff + PART_ZPOS] = 0.2 + 0.2 * xcyc[2];
        xcyc = roundRand3D();
        this.s0[pOff + PART_XVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[0]);
        this.s0[pOff + PART_YVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[1]);
        this.s0[pOff + PART_ZVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[2]);
        this.s0[pOff + PART_X_FTOT] = 0.0;
        this.s0[pOff + PART_Y_FTOT] = 0.0;
        this.s0[pOff + PART_Z_FTOT] = 0.0;
        this.s0[pOff + PART_R] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_G] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_B] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_MASS] = 0.9 + 0.2 * Math.random();
        this.s0[pOff + PART_DIAM] = 1.0 + 10.0 * Math.random();
        this.s0[pOff + PART_RENDMODE] = Math.floor(4.0 * Math.random()); // 0,1,2 or 3.
        this.s0[pOff + PART_AGE] = 0;  // # of frame-times since creation/initialization
        this.s0[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        this.s0[pOff + PART_MASS_VEL] = 0;  // time-rate-of-change of mass.
        this.s0[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        this.s0[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
        this.s0[pOff + PART_G_VEL] = 0;  // time-rate-of-change of color:grn
        this.s0[pOff + PART_B_VEL] = 0;  // time-rate-of-change of color:blu
        this.s0[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        this.s0[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        this.s0[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

        this.s1 = this.s0.slice(0);
        this.sm = this.s0.slice(0);
        this.s2 = this.s0.slice(0);

    }
    this.cforcer = [
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            s0[pOff + PART_Z_FTOT] = -9.832 * s0[pOff + PART_MASS];
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            s0[pOff + PART_X_FTOT] -= (s0[pOff + PART_DIAM] * s0[pOff + PART_XVEL] / 2);
            s0[pOff + PART_Y_FTOT] -= (s0[pOff + PART_DIAM] * s0[pOff + PART_YVEL] / 2);
            s0[pOff + PART_Z_FTOT] -= (s0[pOff + PART_DIAM] * s0[pOff + PART_ZVEL] / 2);
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        }
    ];

    this.c0 = [
        function (i, s1, s0) {
            if (g_bounce == 0) { //--------------------------------------------------------
                if (s1[i * PART_MAXVAR + PART_XPOS] < 0.0 && s1[i * PART_MAXVAR + PART_XVEL] < 0.0			// simple velocity-reversal
                ) {		// bounce on left wall.
                    s1[i * PART_MAXVAR + PART_XVEL] = -g_resti * s1[i * PART_MAXVAR + PART_XVEL];
                } else if (s1[i * PART_MAXVAR + PART_XPOS] > axisSize && s1[i * PART_MAXVAR + PART_XVEL] > 0.0
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_XVEL] = -g_resti * s1[i * PART_MAXVAR + PART_XVEL];
                }
                if (s1[i * PART_MAXVAR + PART_YPOS] < 0.0 && s1[i * PART_MAXVAR + PART_YVEL] < 0.0
                ) {		// bounce on floor
                    s1[i * PART_MAXVAR + PART_YVEL] = -g_resti * s1[i * PART_MAXVAR + PART_YVEL];
                } else if (s1[i * PART_MAXVAR + PART_YPOS] > axisSize && s1[i * PART_MAXVAR + PART_YVEL] > 0.0
                ) {		// bounce on ceiling
                    s1[i * PART_MAXVAR + PART_YVEL] = -g_resti * s1[i * PART_MAXVAR + PART_YVEL];
                }
                if (s1[i * PART_MAXVAR + PART_ZPOS] < 0 && s1[i * PART_MAXVAR + PART_ZVEL] < 0.0			// simple velocity-reversal
                ) {		// bounce on left wall.
                    s1[i * PART_MAXVAR + PART_ZVEL] = -g_resti * s1[i * PART_MAXVAR + PART_ZVEL];
                } else if (s1[i * PART_MAXVAR + PART_ZPOS] > axisSize && s1[i * PART_MAXVAR + PART_ZVEL] > 0.0
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_ZVEL] = -g_resti * s1[i * PART_MAXVAR + PART_ZVEL];
                }
                //  -- hard limit on 'floor' keeps y position >= 0;
                if (s1[i * PART_MAXVAR + PART_ZPOS] < 0.0) s1[i * PART_MAXVAR + PART_ZPOS] = 0.0;
            } else if (g_bounce == 1) { //---------------------------------------------------------------------------
                if (s1[i * PART_MAXVAR + PART_XPOS] < 0.0 && s1[i * PART_MAXVAR + PART_XVEL] < 0.0 // collision!  left wall...
                ) {		// bounce on left wall.
                    s1[i * PART_MAXVAR + PART_XPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    s1[i * PART_MAXVAR + PART_XVEL] = s0[i * PART_MAXVAR + PART_XVEL];			// we had a the START of the timestep.
                    s1[i * PART_MAXVAR + PART_XVEL] *= g_drag;			// **BUT** velocity during our timestep is STILL
                    if (s1[i * PART_MAXVAR + PART_XVEL] < 0.0) s1[i * PART_MAXVAR + PART_XVEL] = -g_resti * s1[i * PART_MAXVAR + PART_XVEL]; // no sign change--bounce!
                    else s1[i * PART_MAXVAR + PART_XVEL] = g_resti * s1[i * PART_MAXVAR + PART_XVEL];			// sign changed-- don't need another.
                } else if (s1[i * PART_MAXVAR + PART_XPOS] > axisSize && s1[i * PART_MAXVAR + PART_XVEL] > 0.0		// collision! right wall...
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_XPOS] = axisSize;					// 1) resolve contact: put particle at wall.
                    // 2) remove all x velocity gained from forces as
                    // ball moved thru wall in this timestep. HOW?
                    // Assume ball reached wall at START of
                    // the timestep, thus: return to the orig.
                    s1[i * PART_MAXVAR + PART_XVEL] = s0[i * PART_MAXVAR + PART_XVEL];			// velocity we had at the start of timestep;
                    s1[i * PART_MAXVAR + PART_XVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    if (s1[i * PART_MAXVAR + PART_XVEL] > 0.0) s1[i * PART_MAXVAR + PART_XVEL] = -g_resti * s1[i * PART_MAXVAR + PART_XVEL]; // no sign change--bounce!
                    else s1[i * PART_MAXVAR + PART_XVEL] = g_resti * s1[i * PART_MAXVAR + PART_XVEL];			// sign changed-- don't need another.
                }
                if (s1[i * PART_MAXVAR + PART_YPOS] < 0.0 && s1[i * PART_MAXVAR + PART_YVEL] < 0.0		// collision! left wall...
                ) {		// bounce on floor

                    s1[i * PART_MAXVAR + PART_YPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru floor in this timestep. HOW?
                    // Assume ball reached floor at START of
                    // the timestep, thus: return to the orig.
                    s1[i * PART_MAXVAR + PART_YVEL] = s0[i * PART_MAXVAR + PART_YVEL];			// velocity we had at the start of timestep;
                    s1[i * PART_MAXVAR + PART_YVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    if (s1[i * PART_MAXVAR + PART_YVEL] < 0.0) s1[i * PART_MAXVAR + PART_YVEL] = -g_resti * s1[i * PART_MAXVAR + PART_YVEL]; // no sign change--bounce!
                    else s1[i * PART_MAXVAR + PART_YVEL] = g_resti * s1[i * PART_MAXVAR + PART_YVEL];			// sign changed-- don't need another.
                } else if (s1[i * PART_MAXVAR + PART_YPOS] > axisSize && s1[i * PART_MAXVAR + PART_YVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling
                    s1[i * PART_MAXVAR + PART_YPOS] = axisSize;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru ceiling in this timestep. HOW?
                    // Assume ball reached ceiling at START of
                    // the timestep, thus: return to the orig.
                    s1[i * PART_MAXVAR + PART_YVEL] = s0[i * PART_MAXVAR + PART_YVEL];			// velocity we had at the start of timestep;
                    s1[i * PART_MAXVAR + PART_YVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    if (s1[i * PART_MAXVAR + PART_YVEL] > 0.0) s1[i * PART_MAXVAR + PART_YVEL] = -g_resti * s1[i * PART_MAXVAR + PART_YVEL]; // no sign change--bounce!
                    else s1[i * PART_MAXVAR + PART_YVEL] = g_resti * s1[i * PART_MAXVAR + PART_YVEL];			// sign changed-- don't need another.

                }
                //	console.log('z = ' + s0[PART_ZPOS] + '  zVel = ' +s0[PART_ZVEL]);
                if (s1[i * PART_MAXVAR + PART_ZPOS] < 0.0 && s1[i * PART_MAXVAR + PART_ZVEL] < 0.0		// collision! left wall...
                ) {		// bounce on floor

                    s1[i * PART_MAXVAR + PART_ZPOS] = 0.0;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru floor in this timestep. HOW?
                    // Assume ball reached floor at START of
                    // the timestep, thus: return to the orig.
                    s1[i * PART_MAXVAR + PART_ZVEL] = s0[i * PART_MAXVAR + PART_ZVEL];			// velocity we had at the start of timestep;
                    s1[i * PART_MAXVAR + PART_ZVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    if (s1[i * PART_MAXVAR + PART_ZVEL] < 0.0) s1[i * PART_MAXVAR + PART_ZVEL] = -g_resti * s1[i * PART_MAXVAR + PART_ZVEL]; // no sign change--bounce!
                    else s1[i * PART_MAXVAR + PART_ZVEL] = g_resti * s1[i * PART_MAXVAR + PART_ZVEL];			// sign changed-- don't need another.
                } else if (s1[i * PART_MAXVAR + PART_ZPOS] > axisSize && s1[i * PART_MAXVAR + PART_ZVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling

                    s1[i * PART_MAXVAR + PART_ZPOS] = axisSize;					// 1) resolve contact: put particle at wall.
                    // 2) remove all y velocity gained from forces as
                    // ball moved thru ceiling in this timestep. HOW?
                    // Assume ball reached ceiling at START of
                    // the timestep, thus: return to the orig.
                    s1[i * PART_MAXVAR + PART_ZVEL] = s0[i * PART_MAXVAR + PART_ZVEL];			// velocity we had at the start of timestep;
                    s1[i * PART_MAXVAR + PART_ZVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces
                    if (s1[i * PART_MAXVAR + PART_ZVEL] > 0.0) s1[i * PART_MAXVAR + PART_ZVEL] = -g_resti * s1[i * PART_MAXVAR + PART_ZVEL]; // no sign change--bounce!
                    else s1[i * PART_MAXVAR + PART_ZVEL] = g_resti * s1[i * PART_MAXVAR + PART_ZVEL];			// sign changed-- don't need another.
                }
            } else {
                console.log('?!?! unknown constraint: g_bounce==' + g_bounce);
                return;
            }
        },
        function (pOff, s1) {
        },
        function (pOff, s1) {
        }

    ]
}
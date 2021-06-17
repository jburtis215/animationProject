PartSys.prototype.BoidS0 = function () {
    this.partCount = 90;
    var doit = 1;
    this.s0 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s1 = new Float32Array(this.partCount * PART_MAXVAR);
    this.sdot = new Float32Array(this.partCount * PART_MAXVAR);
    this.buffer = new Float32Array(this.partCount * PART_MAXVAR);
    this.sM = new Float32Array(this.partCount * PART_MAXVAR);
    this.s2dot = new Float32Array(this.partCount * PART_MAXVAR);
    this.s2 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s1dot2 = new Float32Array(this.partCount * PART_MAXVAR);
//this.cforcer = new Float32Array(PART_MAXVAR);

    for (var i = 0; i < this.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        var xcyc = roundRand3D();
        if (doit == 1) {
            console.log('xc,yc= ' + xcyc[0] + ', ' + xcyc[1]);
            doit = 0;
        }
        this.s0[pOff + PART_XPOS] = (Math.cos(360 / i) / .9) + .9;		// 0.0 <= randomRound() < 1.0
        this.s0[pOff + PART_YPOS] = (Math.sin(360 / i) / .9) + .9;
        this.s0[pOff + PART_ZPOS] = 1.2;
        xcyc = roundRand3D();
        this.s0[pOff + PART_XVEL] = 0;//INIT_VEL*(-0.2 + 0.4*xcyc[0]);
        this.s0[pOff + PART_YVEL] = 0;//INIT_VEL*(-0.2 + 0.4*xcyc[1]);
        this.s0[pOff + PART_ZVEL] = 0;//INIT_VEL*(-0.2 + 0.4*xcyc[2]);
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
    /*
        this.cforcer[F_NONE] = 0.0;
        this.cforcer[F_MOUSE] = 0.0;
        this.cforcer[F_GRAV_E] = function grav(pOff) {
            this.s0[pOff + PART_Z_FTOT] = -2.832 / partSys1.s0[pOff + PART_MASS];
            };
        this.cforcer[F_GRAV_P] = 0.0;
        this.cforcer[F_WIND] = 0.0;
        this.cforcer[F_BUBBLE] = 0.0;
        this.cforcer[F_DRAG] = 0.995;
        this.cforcer[F_SPRING] = 0.0;
        this.cforcer[F_SPH] = 0.0;
        this.cforcer[F_CHARGE] = 0.0;*/
    this.cforcer = [
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            //s0[pOff + PART_Z_FTOT] = -0.832 /s0[pOff + PART_MASS];
        },
        function (pOff, s0) {

            s0[pOff + PART_X_FTOT] = 7 * (s0[PART_MAXVAR + PART_XPOS] - s0[pOff + PART_XPOS]);
            s0[pOff + PART_Y_FTOT] = 7 * (s0[PART_MAXVAR + PART_YPOS] - s0[pOff + PART_YPOS]);
            s0[pOff + PART_Z_FTOT] = 7 * (s0[PART_MAXVAR + PART_ZPOS] - s0[pOff + PART_ZPOS]);

        },
        function (pOff, s0) {
            //	if(s0[pOff + PART_XPOS] > 1) s0[pOff + PART_X_FTOT] = -.1;
            //	if(s0[pOff + PART_XPOS] < .8) s0[pOff + PART_X_FTOT] = .1;
            //	if(s0[pOff + PART_YPOS] > 1) s0[pOff + PART_Y_FTOT] = -.1;
            //	if(s0[pOff + PART_YPOS] < .8) s0[pOff + PART_Y_FTOT] = .1;
            if (pOff == PART_MAXVAR) {
                s0[pOff + PART_X_FTOT] += .5 * (.9 - s0[pOff + PART_YPOS]);
                s0[pOff + PART_Y_FTOT] += -.5 * (.9 - s0[pOff + PART_XPOS]);
                s0[pOff + PART_X_FTOT] += (axisSize - s0[pOff + PART_ZPOS]) * (.9 - s0[pOff + PART_XPOS]);
                s0[pOff + PART_Y_FTOT] += (axisSize - s0[pOff + PART_ZPOS]) * (.9 - s0[pOff + PART_YPOS]);
                if (s0[pOff + PART_ZPOS] > 2.2)
                    s0[pOff + PART_Z_FTOT] -= 1.5;
                if (s0[pOff + PART_ZPOS] < 1.5) s0[pOff + PART_Z_FTOT] += 1.5;
                //	if(s0[pOff + PART_AGE]%2 == 0) s0[pOff + PART_Z_FTOT] -= 1.5;
                //	else(s0[pOff + PART_AGE]%2 == 1)
            } else {
                for (i = PART_MAXVAR; i < pOff; i += PART_MAXVAR) {
                    var xdist = s0[i + PART_XPOS] - s0[pOff + PART_XPOS];
                    var ydist = s0[i + PART_YPOS] - s0[pOff + PART_YPOS];
                    var zdist = s0[i + PART_ZPOS] - s0[pOff + PART_ZPOS];
                    var dist = Math.sqrt(xdist * xdist + ydist * ydist + zdist * zdist);
                    s0[pOff + PART_X_FTOT] -= (.75 / dist) * (xdist / dist);
                    s0[pOff + PART_Y_FTOT] -= (.75 / dist) * (ydist / dist);
                    s0[pOff + PART_Z_FTOT] -= (.75 / dist) * (zdist / dist);

                    s0[pOff + PART_X_FTOT] += .75 * (s0[i + PART_XVEL] - s0[pOff + PART_XVEL]);
                    s0[pOff + PART_Y_FTOT] += .75 * (s0[i + PART_YVEL] - s0[pOff + PART_YVEL]);
                    s0[pOff + PART_Z_FTOT] += .75 * (s0[i + PART_ZVEL] - s0[pOff + PART_ZVEL]);

                    s0[pOff + PART_X_FTOT] += .85 * (xdist);
                    s0[pOff + PART_Y_FTOT] += .85 * (ydist);
                    s0[pOff + PART_Z_FTOT] += .85 * (zdist);
                }

            }
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
//	var xdif = s0[pOff + PART_XPOS] - .9;
//	var ydif = s0[pOff + PART_YPOS] - .9;

            s0[pOff + PART_X_FTOT] -= (0.3 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_XVEL] / 2);
            s0[pOff + PART_Y_FTOT] -= (0.3 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_YVEL] / 2);
//	s0[pOff + PART_Z_FTOT] -= (0.3 *  0.15*3.14 * s0[pOff + PART_DIAM]* s0[pOff + PART_ZVEL] / 2);

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
                } else if (s1[i * PART_MAXVAR + PART_XPOS] > 4 && s1[i * PART_MAXVAR + PART_XVEL] > 0.0		// collision! right wall...
                ) {		// bounce on right wall
                    s1[i * PART_MAXVAR + PART_XPOS] = 4;					// 1) resolve contact: put particle at wall.
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
                } else if (s1[i * PART_MAXVAR + PART_YPOS] > 4 && s1[i * PART_MAXVAR + PART_YVEL] > 0.0 		// collision! front wall...
                ) {		// bounce on ceiling
                    s1[i * PART_MAXVAR + PART_YPOS] = 4;					// 1) resolve contact: put particle at wall.
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
                }/*
			else if( s1[i* PART_MAXVAR + PART_ZPOS] > 4 && s1[i* PART_MAXVAR + PART_ZVEL] > 0.0 		// collision! front wall...
			) {		// bounce on ceiling

				s1[i* PART_MAXVAR + PART_ZPOS] = 4;					// 1) resolve contact: put particle at wall.
																// 2) remove all y velocity gained from forces as
																// ball moved thru ceiling in this timestep. HOW?
																// Assume ball reached ceiling at START of
																// the timestep, thus: return to the orig.
				s1[i* PART_MAXVAR + PART_ZVEL] = s0[i* PART_MAXVAR + PART_ZVEL];			// velocity we had at the start of timestep;
				s1[i* PART_MAXVAR + PART_ZVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces 
			if(s1[i* PART_MAXVAR + PART_ZVEL] > 0.0) s1[i* PART_MAXVAR + PART_ZVEL] = -g_resti*s1[i* PART_MAXVAR + PART_ZVEL]; // no sign change--bounce!
				else s1[i* PART_MAXVAR + PART_ZVEL] = g_resti*s1[i* PART_MAXVAR + PART_ZVEL];			// sign changed-- don't need another.
			}
			if(s1[i* PART_MAXVAR + PART_ZPOS] > 18) s1[i* PART_MAXVAR + PART_ZPOS] = 0;*/
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
PartSys.prototype.FireS0 = function () {
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
//this.cforcer = new Float32Array(PART_MAXVAR);

    for (var i = 0; i < this.partCount; i++) {
        var pOff = i * PART_MAXVAR;			// starting index of each particle
        var xcyc = roundRand3D();
        if (doit == 1) {
            console.log('xc,yc= ' + xcyc[0] + ', ' + xcyc[1]);
            doit = 0;
        }
        this.s0[pOff + PART_XPOS] = 0.2 + 0.2 / i;		// 0.0 <= randomRound() < 1.0
        this.s0[pOff + PART_YPOS] = 0.2 + 0.2 / i;
        this.s0[pOff + PART_ZPOS] = 0.2 + 0.2 / i;
        xcyc = roundRand3D();
        this.s0[pOff + PART_XVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[0]);
        this.s0[pOff + PART_YVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[1]);
        this.s0[pOff + PART_ZVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[2]);
        this.s0[pOff + PART_X_FTOT] = 0.0;
        this.s0[pOff + PART_Y_FTOT] = 0.0;
        this.s0[pOff + PART_Z_FTOT] = 0.0;
        this.s0[pOff + PART_R] = (i / this.partCount);
        this.s0[pOff + PART_G] = (i / this.partCount);
        this.s0[pOff + PART_B] = (i / this.partCount);
        this.s0[pOff + PART_MASS] = 0.9 + 0.2 * Math.random();
        this.s0[pOff + PART_DIAM] = 1.0 + 10.0 * Math.random();
        this.s0[pOff + PART_RENDMODE] = Math.floor(4.0 * Math.random()); // 0,1,2 or 3.
        this.s0[pOff + PART_AGE] = 0;  // # of frame-times since creation/initialization
        this.s0[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        this.s0[pOff + PART_MASS_VEL] = .1;  // time-rate-of-change of mass.
        this.s0[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        this.s0[pOff + PART_R_VEL] = .05;  // time-rate-of-change of color:red
        this.s0[pOff + PART_G_VEL] = -.05;  // time-rate-of-change of color:grn
        this.s0[pOff + PART_B_VEL] = -.05;  // time-rate-of-change of color:blu
        this.s0[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        this.s0[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        this.s0[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
        this.s0[pOff + PART_EXP] = Math.random() * 15;


        this.s1[pOff + PART_XPOS] = this.s0[pOff + PART_XPOS];		// 0.0 <= randomRound() < 1.0
        this.s1[pOff + PART_YPOS] = this.s0[pOff + PART_YPOS];
        this.s1[pOff + PART_ZPOS] = this.s0[pOff + PART_ZPOS];
        this.s1[pOff + PART_XVEL] = this.s0[pOff + PART_XVEL];
        this.s1[pOff + PART_YVEL] = this.s0[pOff + PART_YVEL];
        this.s1[pOff + PART_ZVEL] = this.s0[pOff + PART_ZVEL];
        this.s1[pOff + PART_X_FTOT] = 0.0;
        this.s1[pOff + PART_Y_FTOT] = 0.0;
        this.s1[pOff + PART_Z_FTOT] = 0.0;
        this.s1[pOff + PART_R] = this.s0[pOff + PART_R];
        this.s1[pOff + PART_G] = this.s0[pOff + PART_G];
        this.s1[pOff + PART_B] = this.s0[pOff + PART_B];
        this.s1[pOff + PART_MASS] = this.s0[pOff + PART_MASS];
        this.s1[pOff + PART_DIAM] = this.s0[pOff + PART_DIAM];
        this.s1[pOff + PART_RENDMODE] = this.s0[pOff + PART_RENDMODE]; // 0,1,2 or 3.
        this.s1[pOff + PART_AGE] = 0;  // # of frame-times since creation/initialization
        this.s1[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        this.s1[pOff + PART_MASS_VEL] = 0;  // time-rate-of-change of mass.
        this.s1[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        this.s1[pOff + PART_R_VEL] = this.s0[pOff + PART_R_VEL];  // time-rate-of-change of color:red
        this.s1[pOff + PART_G_VEL] = this.s0[pOff + PART_G_VEL];  // time-rate-of-change of color:grn
        this.s1[pOff + PART_B_VEL] = this.s0[pOff + PART_B_VEL];  // time-rate-of-change of color:blu
        this.s1[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        this.s1[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        this.s1[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
        this.s1[pOff + PART_EXP] = this.s0[pOff + PART_EXP];

        this.sM[pOff + PART_XPOS] = this.s0[pOff + PART_XPOS];		// 0.0 <= randomRound() < 1.0
        this.sM[pOff + PART_YPOS] = this.s0[pOff + PART_YPOS];
        this.sM[pOff + PART_ZPOS] = this.s0[pOff + PART_ZPOS];
        this.sM[pOff + PART_XVEL] = this.s0[pOff + PART_XVEL];
        this.sM[pOff + PART_YVEL] = this.s0[pOff + PART_YVEL];
        this.sM[pOff + PART_ZVEL] = this.s0[pOff + PART_ZVEL];
        this.sM[pOff + PART_X_FTOT] = 0.0;
        this.sM[pOff + PART_Y_FTOT] = 0.0;
        this.sM[pOff + PART_Z_FTOT] = 0.0;
        this.sM[pOff + PART_R] = this.s0[pOff + PART_R];
        this.sM[pOff + PART_G] = this.s0[pOff + PART_G];
        this.sM[pOff + PART_B] = this.s0[pOff + PART_B];
        this.sM[pOff + PART_MASS] = this.s0[pOff + PART_MASS];
        this.sM[pOff + PART_DIAM] = this.s0[pOff + PART_DIAM];
        this.sM[pOff + PART_RENDMODE] = this.s0[pOff + PART_RENDMODE]; // 0,1,2 or 3.
        this.sM[pOff + PART_AGE] = 0;  // # of frame-times since creation/initialization
        this.sM[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        this.sM[pOff + PART_MASS_VEL] = 0;  // time-rate-of-change of mass.
        this.sM[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        this.sM[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
        this.sM[pOff + PART_G_VEL] = 0;  // time-rate-of-change of color:grn
        this.sM[pOff + PART_B_VEL] = 0;  // time-rate-of-change of color:blu
        this.sM[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        this.sM[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        this.sM[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu

        this.s2[pOff + PART_XPOS] = this.s0[pOff + PART_XPOS];		// 0.0 <= randomRound() < 1.0
        this.s2[pOff + PART_YPOS] = this.s0[pOff + PART_YPOS];
        this.s2[pOff + PART_ZPOS] = this.s0[pOff + PART_ZPOS];
        this.s2[pOff + PART_XVEL] = this.s0[pOff + PART_XVEL];
        this.s2[pOff + PART_YVEL] = this.s0[pOff + PART_YVEL];
        this.s2[pOff + PART_ZVEL] = this.s0[pOff + PART_ZVEL];
        this.s2[pOff + PART_X_FTOT] = 0.0;
        this.s2[pOff + PART_Y_FTOT] = 0.0;
        this.s2[pOff + PART_Z_FTOT] = 0.0;
        this.s2[pOff + PART_R] = this.s0[pOff + PART_R];
        this.s2[pOff + PART_G] = this.s0[pOff + PART_G];
        this.s2[pOff + PART_B] = this.s0[pOff + PART_B];
        this.s2[pOff + PART_MASS] = this.s0[pOff + PART_MASS];
        this.s2[pOff + PART_DIAM] = this.s0[pOff + PART_DIAM];
        this.s2[pOff + PART_RENDMODE] = this.s0[pOff + PART_RENDMODE]; // 0,1,2 or 3.
        this.s2[pOff + PART_AGE] = 0;  // # of frame-times since creation/initialization
        this.s2[pOff + PART_CHARGE] = 0;  // for electrostatic repulsion/attraction
        this.s2[pOff + PART_MASS_VEL] = 0;  // time-rate-of-change of mass.
        this.s2[pOff + PART_MASS_FTOT] = 0;  // force-accumulator for mass-change
        this.s2[pOff + PART_R_VEL] = 0;  // time-rate-of-change of color:red
        this.s2[pOff + PART_G_VEL] = 0;  // time-rate-of-change of color:grn
        this.s2[pOff + PART_B_VEL] = 0;  // time-rate-of-change of color:blu
        this.s2[pOff + PART_R_FTOT] = 0;  // force-accumulator for color-change: red
        this.s2[pOff + PART_G_FTOT] = 0;  // force-accumulator for color-change: grn
        this.s2[pOff + PART_B_FTOT] = 0;  // force-accumulator for color-change: blu
    }
    this.cforcer = [
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            s0[pOff + PART_Z_FTOT] -= 2.832 / s0[pOff + PART_MASS];
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            s0[pOff + PART_X_FTOT] -= (0.45 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_XVEL] / 2);
            s0[pOff + PART_Y_FTOT] -= (0.45 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_YVEL] / 2);
            s0[pOff + PART_Z_FTOT] -= (0.45 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_ZVEL] / 2);
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (particleOffset, s0) {
            s0[particleOffset + PART_R_FTOT] = -.1;
            s0[particleOffset + PART_G_FTOT] = -.1;
            s0[particleOffset + PART_B_FTOT] = -.1;
            s0[particleOffset + PART_MASS_FTOT] = .1;
        }
    ];
    this.c0 = [
        function (i, s1, s0) {
        },
        function (particleOffset, s1) {
            var xcyc = roundRand3D();
            if (s1[particleOffset + PART_AGE] > s1[particleOffset + PART_EXP]) {
                s1[particleOffset + PART_XPOS] = 0.2;
                s1[particleOffset + PART_YPOS] = 0.2;
                s1[particleOffset + PART_ZPOS] = 0.2;
                s1[particleOffset + PART_XVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[0]);
                s1[particleOffset + PART_YVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[1]);
                s1[particleOffset + PART_ZVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[2]);
                s1[particleOffset + PART_R] = 1;
                s1[particleOffset + PART_G] = 1;
                s1[particleOffset + PART_B] = 1;
                s1[particleOffset + PART_R_VEL] = 0;
                s1[particleOffset + PART_G_VEL] = 0;
                s1[particleOffset + PART_B_VEL] = 0;
                s1[particleOffset + PART_MASS] = 0.9 + 0.2 * Math.random();
                s1[particleOffset + PART_AGE] = 0;
            }
        },
        function (pOff, s1) {
        }
    ];
}
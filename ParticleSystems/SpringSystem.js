PartSys.prototype.SpringS0 = function () {
    this.partCount = 10;
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
        this.s0[pOff + PART_XPOS] = 2		// 0.0 <= randomRound() < 1.0
        this.s0[pOff + PART_YPOS] = 0.2;
        this.s0[pOff + PART_ZPOS] = 5.2 - (.5 * i);
        ;
        xcyc = roundRand3D();
        this.s0[pOff + PART_XVEL] = 0;//INIT_VEL*(0.4 + 0.2*xcyc[0]);
        this.s0[pOff + PART_YVEL] = 0;//INIT_VEL*(0.4 + 0.2*xcyc[1]);
        this.s0[pOff + PART_ZVEL] = 0;//INIT_VEL*(0.4 + 0.2*xcyc[2]);
        this.s0[pOff + PART_X_FTOT] = 0.0;
        this.s0[pOff + PART_Y_FTOT] = 0.0;
        this.s0[pOff + PART_Z_FTOT] = 0.0;
        this.s0[pOff + PART_R] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_G] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_B] = 0.2 + 0.8 * Math.random();
        this.s0[pOff + PART_MASS] = 0.9;
        this.s0[pOff + PART_DIAM] = 10;
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
        },//if(pOff != 0)s0[pOff + PART_Z_FTOT] = -2.832 /s0[pOff + PART_MASS];},
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
//	s0[pOff + PART_X_FTOT] -= (0.45 *  0.15*3.14 * s0[pOff + PART_DIAM]* s0[pOff + PART_XVEL] / 2);
//	s0[pOff + PART_Y_FTOT] -= (0.45 *  0.15*3.14 * s0[pOff + PART_DIAM]* s0[pOff + PART_YVEL] / 2);
//	s0[pOff + PART_Z_FTOT] -= (0.45 *  0.15*3.14 * s0[pOff + PART_DIAM]* s0[pOff + PART_ZVEL] / 2);

        },
        function (pOff, s0) {
            if (pOff != 0) {
                //Find the force from the previous particle

                var zcenter = s0[pOff - PART_MAXVAR + PART_ZPOS] - (.25 * (s0[(pOff - PART_MAXVAR) + PART_ZPOS] - s0[pOff + PART_ZPOS]));
                //if(s0[pOff + PART_XPOS] < center){
                s0[pOff + PART_Z_FTOT] += -4 * (s0[pOff + PART_ZPOS] - zcenter + .12);
                var xcenter = s0[pOff - PART_MAXVAR + PART_XPOS] + (.5 * (s0[pOff + PART_XPOS] - s0[(pOff - PART_MAXVAR) + PART_XPOS]));
                s0[pOff + PART_X_FTOT] += -4 * (s0[pOff + PART_XPOS] - xcenter);
                var ycenter = s0[pOff - PART_MAXVAR + PART_YPOS] + (.5 * (s0[pOff + PART_YPOS] - s0[(pOff - PART_MAXVAR) + PART_YPOS]));
                s0[pOff + PART_Y_FTOT] += -4 * (s0[pOff + PART_YPOS] - ycenter);

                s0[pOff + PART_Z_FTOT] += -1 * s0[pOff + PART_ZVEL]

            }
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
        },
        function (pOff, s1) {
        },
        function (pOff, s1) {
        }

    ]
}
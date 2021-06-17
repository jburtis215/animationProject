function PartSys() {
    this.s1 = new Float32Array(this.partCount * PART_MAXVAR);
    this.s0dot = new Float32Array(this.partCount * PART_MAXVAR);
}

findDensity = function (s, partCount) {
    let maxDensity = 0;
    let minDensity = 100000;
    for (let j = 0; j < partCount; j++) {
      //  console.log("new particle");
        var totDensity = 1;
        let p = j * PART_MAXVAR;
        for (let i = 0; i < partCount; i++) {
            let pOff = i * PART_MAXVAR;
            if (sameVoxel(s, p, pOff)) {
                if (pOff !== p) {
                    let xDist = s[p + PART_XPOS] - s[pOff + PART_XPOS];
                    let yDist = s[p + PART_YPOS] - s[pOff + PART_YPOS];
                    let zDist = s[p + PART_ZPOS] - s[pOff + PART_ZPOS];
                    let realDistance = Math.sqrt(Math.pow(xDist, 2) +
                        Math.pow(yDist, 2) +
                        Math.pow(zDist, 2))
                    realDistance = realDistance === 0 ? epsilon : realDistance;
                    if (realDistance < 2 * h) {
                        let kernel = kernelCoefficient * Math.pow((2 - (realDistance / h)), 3);
                        totDensity += kernel * s[pOff + PART_MASS];
            //            console.log("distance:" + realDistance);
              //          console.log(kernel * s[pOff + PART_MASS]);
                    }
                }
            }
        }
        s[p + PART_DENSITY] = totDensity;
      //  console.log(totDensity);
        if (maxDensity < totDensity) {
            maxDensity = totDensity;
        }
        if (minDensity > totDensity) {
            minDensity = totDensity;
        }
    }
    return [minDensity, maxDensity];
}

findPressureForce = function (s, partCount) {

    for (let j = 0; j < partCount; j++) {
        let p = j * PART_MAXVAR;
        var firstComponent = [0, 0, 0];
        var secondComponent = [0, 0, 0];
        for (let i = 0; i < partCount; i++) {
            let pOff = i * PART_MAXVAR;
            if (sameVoxel(s, p, pOff)) {
                if (pOff !== p) {
                    let xDist = s[p + PART_XPOS] - s[pOff + PART_XPOS];
                    //   xDist = xDist === 0? epsilon: xDist;
                    let yDist = s[p + PART_YPOS] - s[pOff + PART_YPOS];
                    // yDist = yDist === 0? epsilon: yDist;
                    let zDist = s[p + PART_ZPOS] - s[pOff + PART_ZPOS];
                    // zDist = zDist === 0? epsilon: zDist;
                    let realDistance = Math.sqrt(Math.pow(xDist, 2) +
                        Math.pow(yDist, 2) +
                        Math.pow(zDist, 2))
                    let totalDistance = Math.abs(xDist) + Math.abs(yDist) + Math.abs(zDist);
                    if (realDistance < 2 * h) {
                        let kernelDer = -45 * (2 - Math.pow((realDistance / h), 2)) / (64 * h * Math.PI);

                        for (var pos = 0; pos < 3; pos++) {
                            let distance = s[p + pos] - s[pOff + pos];
                            let distanceComponent = distance / realDistance;
                         //   if (distance > 0) {
                                firstComponent[pos] += s[pOff + PART_MASS] * kernelDer * (-1 * distanceComponent);
                                secondComponent[pos] += s[pOff + PART_MASS] * ((s[pOff + PART_DENSITY] - DENSITY_CONST) / ((s[pOff + PART_DENSITY] * (s[pOff + PART_DENSITY])))) * kernelDer * (-1 * distanceComponent);
                           //      }
                         /*   else{
                                 firstComponent[pos] -= s[pOff + PART_MASS] * kernelDer * distanceComponent;
                               secondComponent[pos] -= s[pOff + PART_MASS] * ((s[pOff + PART_DENSITY] - DENSITY_CONST) / (s[pOff + PART_DENSITY])) * kernelDer * distanceComponent;
                            }*/

                        }
                    }
                }
            }
        }
        let coefficient = k * s[p + PART_MASS];
        let firstDensCo = ((s[p + PART_DENSITY] - DENSITY_CONST) / (s[p + PART_DENSITY]* s[p + PART_DENSITY]));
        s[p + PART_X_FTOT] += coefficient * ((firstDensCo * firstComponent[0]) + secondComponent[0]);
        s[p + PART_Y_FTOT] += coefficient * ((firstDensCo * firstComponent[1]) + secondComponent[1]);
        s[p + PART_Z_FTOT] += coefficient * ((firstDensCo * firstComponent[2]) + secondComponent[2]);
  //            console.log('x ' + s[p + PART_X_FTOT] + '\ny ' + s[p + PART_Y_FTOT] + '\nz ' + s[p + PART_Z_FTOT]);
    }
}

findViscosity = function (s, partCount) {
    for (let j = 0; j < partCount; j++) {
        var p = j * PART_MAXVAR;
        var totViscosity = [0, 0, 0];
        for (let i = 0; i < partCount; i++) {
            var pOff = i * PART_MAXVAR;
            if (sameVoxel(s, p, pOff)) {
                if (pOff !== p) {
                    let xDist = s[p + PART_XPOS] - s[pOff + PART_XPOS];
                    xDist = xDist === 0 ? epsilon : xDist;
                    let yDist = s[p + PART_YPOS] - s[pOff + PART_YPOS];
                    yDist = yDist === 0 ? epsilon : yDist;
                    let zDist = s[p + PART_ZPOS] - s[pOff + PART_ZPOS];
                    zDist = zDist === 0 ? epsilon : zDist;
                    let realDistance = Math.sqrt(Math.pow(xDist, 2) +
                        Math.pow(yDist, 2) +
                        Math.pow(zDist, 2))
                    let u = h * (s[p + PART_XVEL] - s[pOff + PART_XVEL]) * realDistance / ((Math.pow(realDistance, 2)) + h * h / 100)
                    if (u < 0) {
                        let totalOfInverseDistances = (xDist === 0 ? 1 : 1 / Math.abs(xDist))
                            + (yDist === 0 ? 1 : 1 / Math.abs(yDist))
                            + (zDist === 0 ? 1 : 1 / Math.abs(zDist));

                        let totalDistance = Math.abs(xDist) + Math.abs(yDist) + Math.abs(zDist);
                        let kernelDer = -45 * (2 - Math.pow((realDistance / h), 2)) / (64 * h * Math.PI);

                        for (var pos = 0; pos < 3; pos++) {
                            let distance = s[p + pos] - s[pOff + pos];
                            let distanceComponent = distance / totalDistance;
                            let pForce = (s[p + PART_DENSITY] + s[pOff + PART_DENSITY]) / 2
                            let force = (-speedOfSound * u + 2 * u * u) / pForce
                            totViscosity[pos] += s[pOff + PART_MASS] * force * kernelDer * distanceComponent;
                        }
                        /*
                                                for (let pos = 0; pos < 3; pos++) {
                                                    let distance = s[p + pos] - s[pOff + pos];
                                                    let inverseDistance = (distance === 0 ? 1 : 1 / Math.abs(distance));
                                                    let distanceComponent = inverseDistance / totalOfInverseDistances;
                                                    let pForce = (s[p + PART_DENSITY] + s[pOff + PART_DENSITY]) / 2
                                                    let force = (-speedOfSound * u + 2 * u * u) / pForce
                                                    if (distance < 0) {
                                                        totViscosity[pos] -= s[pOff + PART_MASS] * force * kernelDer * distanceComponent;
                                                    }else{
                                                        totViscosity[pos] += s[pOff + PART_MASS] * force * kernelDer * distanceComponent;
                                                    }
                                                }*/

                    }
                }
            }
        }
        s[p + PART_X_FTOT] -= s[p + PART_MASS] * totViscosity[0];
        s[p + PART_Y_FTOT] -= s[p + PART_MASS] * totViscosity[1];
        s[p + PART_Z_FTOT] -= s[p + PART_MASS] * totViscosity[2];
        console.log('x ' + s[p + PART_MASS] * totViscosity[0] + 'y ' + s[p + PART_MASS] * totViscosity[1] + 'z ' + s[p + PART_MASS] * totViscosity[2])
    }
}

sameVoxel = function (s0, ind1, ind2) {
    return (Math.abs(s0[ind1 + PART_VOX_X] - s0[ind2 + PART_VOX_X]) <= 1) &&
        (Math.abs(s0[ind1 + PART_VOX_Y] - s0[ind2 + PART_VOX_Y]) <= 1) &&
        (Math.abs(s0[ind1 + PART_VOX_Z] - s0[ind2 + PART_VOX_Z]) <= 1)
}


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
PartSys.prototype.WindyS0 = function () {
    this.partCount = 500;
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
        this.s0[pOff + PART_ZPOS] = 0;
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

    this.cforcer = [
        function (pOff, s0) {
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            s0[pOff + PART_Z_FTOT] = -2.832 / s0[pOff + PART_MASS];
        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
            //	if(s0[pOff + PART_XPOS] > 1) s0[pOff + PART_X_FTOT] = -.1;
            //	if(s0[pOff + PART_XPOS] < .8) s0[pOff + PART_X_FTOT] = .1;
            //	if(s0[pOff + PART_YPOS] > 1) s0[pOff + PART_Y_FTOT] = -.1;
            //	if(s0[pOff + PART_YPOS] < .8) s0[pOff + PART_Y_FTOT] = .1;

            s0[pOff + PART_X_FTOT] = 1.5 * (.9 - s0[pOff + PART_YPOS]);
            s0[pOff + PART_Y_FTOT] = -1.5 * (.9 - s0[pOff + PART_XPOS]);
            s0[pOff + PART_X_FTOT] += 5 * (axisSize - s0[pOff + PART_ZPOS]) * (.9 - s0[pOff + PART_XPOS]);
            s0[pOff + PART_Y_FTOT] += 5 * (axisSize - s0[pOff + PART_ZPOS]) * (.9 - s0[pOff + PART_YPOS]);
            s0[pOff + PART_Z_FTOT] += 10 * (Math.abs(.9 - s0[pOff + PART_XPOS])) + 10 * (.2 - s0[pOff + PART_ZPOS]) +
                10 * (Math.abs(.9 - s0[pOff + PART_YPOS]));

        },
        function (pOff, s0) {
        },
        function (pOff, s0) {
//	var xdif = s0[pOff + PART_XPOS] - .9;
//	var ydif = s0[pOff + PART_YPOS] - .9;

            s0[pOff + PART_X_FTOT] -= (0.3 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_XVEL] / 2);
            s0[pOff + PART_Y_FTOT] -= (0.3 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_YVEL] / 2);
            s0[pOff + PART_Z_FTOT] -= (0.3 * 0.15 * 3.14 * s0[pOff + PART_DIAM] * s0[pOff + PART_ZVEL] / 2);

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
                }/*
			else if( s1[i* PART_MAXVAR + PART_ZPOS] > axisSize && s1[i* PART_MAXVAR + PART_ZVEL] > 0.0 		// collision! front wall...
			) {		// bounce on ceiling

				s1[i* PART_MAXVAR + PART_ZPOS] = axisSize;					// 1) resolve contact: put particle at wall.
																// 2) remove all y velocity gained from forces as
																// ball moved thru ceiling in this timestep. HOW?
																// Assume ball reached ceiling at START of
																// the timestep, thus: return to the orig.
				s1[i* PART_MAXVAR + PART_ZVEL] = s0[i* PART_MAXVAR + PART_ZVEL];			// velocity we had at the start of timestep;
				s1[i* PART_MAXVAR + PART_ZVEL] *= g_drag;			// **BUT** reduced by drag (and any other forces 
			if(s1[i* PART_MAXVAR + PART_ZVEL] > 0.0) s1[i* PART_MAXVAR + PART_ZVEL] = -g_resti*s1[i* PART_MAXVAR + PART_ZVEL]; // no sign change--bounce!
				else s1[i* PART_MAXVAR + PART_ZVEL] = g_resti*s1[i* PART_MAXVAR + PART_ZVEL];			// sign changed-- don't need another.
			}*/
                if (s1[i * PART_MAXVAR + PART_ZPOS] > 18) s1[i * PART_MAXVAR + PART_ZPOS] = 0;
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
        function (pOff, s0) {
            s0[pOff + PART_R_FTOT] = -.1;
            s0[pOff + PART_G_FTOT] = -.1;
            s0[pOff + PART_B_FTOT] = -.1;
            s0[pOff + PART_MASS_FTOT] = .1;
        }
    ];
    this.c0 = [
        function (i, s1, s0) {
        },
        function (pOff, s1) {
            var xcyc = roundRand3D();
            if (s1[pOff + PART_AGE] > s1[pOff + PART_EXP]) {
                s1[pOff + PART_XPOS] = 0.2;
                s1[pOff + PART_YPOS] = 0.2;
                s1[pOff + PART_ZPOS] = 0.2;
                s1[pOff + PART_XVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[0]);
                s1[pOff + PART_YVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[1]);
                s1[pOff + PART_ZVEL] = INIT_VEL * (0.4 + 0.2 * xcyc[2]);
                s1[pOff + PART_R] = 1;
                s1[pOff + PART_G] = 1;
                s1[pOff + PART_B] = 1;
                s1[pOff + PART_R_VEL] = 0;
                s1[pOff + PART_G_VEL] = 0;
                s1[pOff + PART_B_VEL] = 0;
                s1[pOff + PART_MASS] = 0.9 + 0.2 * Math.random();
                s1[pOff + PART_AGE] = 0;
            }
        },
        function (pOff, s1) {
        }
    ];
}


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
/*const F_NONE      0       // Non-existent force: ignore this CForcer object
const F_MOUSE     1       // Spring-like connection to the mouse cursor; lets
                            // you 'grab' and 'wiggle' one particle(or several).
const F_GRAV_E    2       // Earth-gravity: pulls all particles 'downward'.
const F_GRAV_P    3       // Planetary-gravity; particle-pair (e0,e1) attract
                            // each other with force== grav* mass0*mass1/ dist^2
const F_WIND      4       // Blowing-wind-like force-field;fcn of 3D position
const F_BUBBLE    5       // Constant inward force towards centerpoint if
                            // particle is > max_radius away from centerpoint.
const F_DRAG      6       // Viscous drag -- proportional to neg. velocity.
const F_SPRING    7       // ties together 2 particles; distance sets force
const F_SPH 8       // a big collection of identical springs; lets you
                            // make cloth & rubbery shapes as one force-making
                            // object, instead of many many F_SPRING objects.
const F_CHARGE    9       // attract/repel by charge and inverse distance;
                            // applies to all charged particles.
const F_FIRE     10
const F_MAXKINDS  11      // 'max' is always the LAST name in our list;
*/

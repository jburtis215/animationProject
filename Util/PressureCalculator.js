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
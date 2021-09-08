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
        s[p + PART_PRESSURE] = Math.pow(speedOfSound, 2) * (totDensity - DENSITY_CONST);
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
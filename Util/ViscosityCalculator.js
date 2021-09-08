function calculateDistances(s, p, pOff) {
    let randomDistance = getRandomInt(1, 2) % 2 === 0 ? epsilon : epsilon * -1;
    let xDist = s[p + PART_XPOS] - s[pOff + PART_XPOS];
    xDist = xDist === 0 ? randomDistance : xDist;
    let yDist = s[p + PART_YPOS] - s[pOff + PART_YPOS];
    yDist = yDist === 0 ? randomDistance : yDist;
    let zDist = s[p + PART_ZPOS] - s[pOff + PART_ZPOS];
    zDist = zDist === 0 ? randomDistance : zDist;
    let realDistance = Math.sqrt(Math.pow(xDist, 2) +
        Math.pow(yDist, 2) +
        Math.pow(zDist, 2))
    return {xDist, yDist, zDist, realDistance};
}

findViscosity = function (s, partCount) {
    for (let j = 0; j < partCount; j++) {
        var p = j * PART_MAXVAR;
        var totViscosity = [0, 0, 0];
        for (let i = 0; i < partCount; i++) {
            var pOff = i * PART_MAXVAR;
            if (sameVoxel(s, p, pOff)) {
                if (pOff !== p) {
                    let {xDist, yDist, zDist, realDistance} = calculateDistances(s, p, pOff);
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

                        let xdistanceComponent = xDist / totalDistance;
                        let ydistanceComponent = yDist / totalDistance;
                        let zdistanceComponent = zDist / totalDistance;
                        let pForce = (s[p + PART_DENSITY] + s[pOff + PART_DENSITY]) / 2
                        let force = (-speedOfSound * u + 2 * u * u) / pForce
                        totViscosity[0] += s[pOff + PART_MASS] * force * kernelDer * xdistanceComponent;
                        totViscosity[1] += s[pOff + PART_MASS] * force * kernelDer * ydistanceComponent;
                        totViscosity[2] += s[pOff + PART_MASS] * force * kernelDer * zdistanceComponent;
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